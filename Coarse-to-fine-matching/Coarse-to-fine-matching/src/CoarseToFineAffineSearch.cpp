#include "CoarseToFineAffineSearch.h"

/*intialize the affine search and downsample the image and the template*/
int CoarseToFineAffineSearch::initialize(Mat _img, Mat _region, Params _params, double _T1){

	params = _params;
	T1 = _T1;
#ifdef RESIZE

	
	bool doResize = true;
	if(_region.rows <= params.patchRows || _region.cols <= params.patchCols){
		doResize = false;
	}

	if(doResize){
		double dR, dC;
		dR = _region.rows*1.0 / params.patchRows;
		dC = _region.cols*1.0 / params.patchCols;
#ifdef USE_PYRAMID
		dR = log(dR)/log(2.0);
		dC = log(dC)/log(2.0);
		widthResizeFactor = 1<<(int)dC;
		heightResizeFactor = 1<<(int)dR;
		region = _region;
		img = _img;
		for(int i = 0; i < (int)dR; i++){
			pyrDown(region, region);
			pyrDown(img, img);
		}
#else
		widthResizeFactor = dC;
		heightResizeFactor = dR;


		int inter_method = 1;
		region.create(Size(_region.rows*1.0/dR, _region.cols*1.0/dC ), _region.type());

		
		resize(_region, region, Size(params.patchCols, params.patchRows),0,0,inter_method);
		Utils::doAntiAliasing(region, region);

		img.create(Size(_img.cols/dC, _img.rows/dR), _img.type());
		resize(_img, img, Size(_img.cols/dC, _img.rows/dR),0,0,inter_method);
		Utils::doAntiAliasing(img, img);
#endif


	}else{
		widthResizeFactor = heightResizeFactor = 1;
		img = _img;
		region = _region;
	}

#else
	widthResizeFactor = heightResizeFactor = 1;
	img = _img;
	region = _region;
#endif
	maxScore = secondMaxScore = 0;	
	return true;
}

/*
computeBestAffineTransform:
	this function goes through the affine hierarchy to reach the leaf node with the best affine transform and localizes the location of the corresponding region
*/
int CoarseToFineAffineSearch::computeBestAffineTransform(double& cS, double& cLambda, double& cTheta, double& cH, int& numNCC, bool matchUsingPrevTrans, AffineTransform affTrans){
	int N = params.N/2;
	int L = params.affineLevels;
	int l = 1;

	int n = 1;

	#ifndef SEARCH_ALL_IMG
	n = params.numLocalMaxima;
	#endif
	
	vector<LocalMaximum> subsets;
	int numSubsets;
	if(!matchUsingPrevTrans){
		numSubsets= 1 << (N-L+l); // number of subsets at the first level;
		numSubsets *= numSubsets;
		for(int i = 0; i < numSubsets; i++){
			LocalMaximum lm;
			lm.loc = Point(0,0);
			lm.regionAnchor = Point(0,0);
			lm.p = i;
			lm.q = i;
			lm.region = img;
			subsets.push_back(lm);
		}
	}else{
		numSubsets= 1; 
		LocalMaximum lm;
		lm.loc = Point(0,0);
		lm.regionAnchor = Point(0,0);
		lm.p = affTrans.pIdx;
		lm.q = affTrans.pIdx;
		lm.region = img;
		subsets.push_back(lm);
		
		l = L;
	}
	timer total_timer;
	total_timer.tic();
	double topNCC, secondTopNCC;
	vector<vector<PQ> > pqVec;
	vector<vector<PQ> > tmpPQVec;
	timer sampling_timer;
	for(; l <= L; l++){		
		pqVec = params.samplesInLevels[l-1];
		tmpPQVec = pqVec;

		Mat img_disp;
		img.copyTo(img_disp);

		// apply the set of transformations in subsets array (subset of the affine transfomration) at level l 
		vector<LocalMaximum> topLocalMaxima = computeBestMatches(subsets, l, pqVec, numNCC);
		
		vector<LocalMaximum>::iterator localMaximaIt = topLocalMaxima.begin();	

		int k = 0; 
		subsets.clear();
		set<int> tI, tJ, locX, locY;
		bool topLM = true;

		while(k < n && localMaximaIt != topLocalMaxima.end()){ // store the top local maxima

			Scalar color;
			int thickness = 1;
			if(k == 0){
				color = Scalar(0,0,255,0);
				thickness = 2;
			}
			else
				color = Scalar::all(0);
		
			LocalMaximum lm = *localMaximaIt;
			if(tI.find(lm.p) == tI.end() || tJ.find(lm.q) == tJ.end()|| locX.find(lm.loc.x) == locX.end() || locY.find(lm.loc.y) == locY.end()){
				Point regionAnchor;
				#ifdef RESIZE
					Rect roiRect = Utils::getROI(lm, img, params.patchCols, params.patchRows, regionAnchor); // get ROI around the local maxima to be used as a search region in the next levels in the affine heirarchy
				#else
					Rect roiRect = Utils::getROI(lm, img, region.cols, region.rows, regionAnchor);
				#endif
				Mat roi = img(roiRect);
				int mul = 4;
				for(int i = 0; i < mul; i++){ 
					LocalMaximum localMax;
					int cp = (*localMaximaIt).p*mul + i;
					int cq = (*localMaximaIt).q*mul + i;
					localMax.p = cp;
					localMax.q = cq;
					localMax.loc = (*localMaximaIt).loc;
					#ifdef SEARCH_ALL_IMG
					localMax.region = img;
					#else
						localMax.region = roi; // set search region to roi
					#endif

					localMax.regionAnchor = regionAnchor;
					subsets.push_back(localMax);
				}
				k++;
				tI.insert(lm.p);
				tJ.insert(lm.q);
				locX.insert(lm.loc.x);
				locY.insert(lm.loc.y);
				topLM = false;
			}
			localMaximaIt++;
		}
	}

	/*the affine matrix elements */
	cS = tmpPQVec[topLocalMaximum.p][0].p.s; 
	cLambda = tmpPQVec[topLocalMaximum.p][0].p.lambda;
	cTheta = tmpPQVec[topLocalMaximum.p][0].q.theta;
	cH = tmpPQVec[topLocalMaximum.p][0].q.h;

	return true;
}
/*
	compute the average of the affine transformations in each node in the hierarchy at level l
	returns a vector of the top local maxima that were localized with their corresponding node in the hierarhcy
*/
vector<LocalMaximum> CoarseToFineAffineSearch::computeBestMatches(vector<LocalMaximum> subsets, int l, vector<vector<PQ> > pqV, int& numNCC){

	int N = params.N;
	int L = params.affineLevels;
	
	int numElems = (1 << L-l);
	int numSamples = ceil(sqrt(numElems*1.0f));

	int *pSet, *qSet;
	int pSetSize = numElems;
	pSet = new int[pSetSize];
	qSet = new int[pSetSize];
	int incr = numElems/((numSamples>1)?numSamples:1);
	int kk = 0;
	for(int k = 0; k < numElems; k++){
		pSet[kk] = k;
		qSet[kk++] = k;
	}

	int cType = CV_32FC3;

	Mat initialMask(region.size(),cType);
	for(int i =0; i < initialMask.rows; i++){
		for(int j =0; j < initialMask.cols; j++){
			for( int cn = 0; cn < initialMask.channels(); cn++ )
				initialMask.at<Vec3f>(i,j)[cn] = 1;
		}
	}

	vector<LocalMaximum> topLocalMaxima;
	
	timer affine_timer;
	double tTime =0;	
	double minScore = 1000000;
	maxScore = secondMaxScore =0;
	long affine_time = 0;
	int itrIndex = -1;

	int i = 0;
	timer subset_timer;
	subset_timer.tic();
	
	int numThreads = omp_get_max_threads();
	omp_set_dynamic(1);
	int* localNCCs = new int[numThreads];
	vector<LocalMaximum>* topLocalMaximaArr = new vector<LocalMaximum>[numThreads];
	for(int i =0 ;i < numThreads; i++){
		localNCCs[i] = 0;
		topLocalMaximaArr[i] = vector<LocalMaximum>();
	}
#ifdef USE_OMP
	#pragma omp parallel shared(topLocalMaximaArr, subsets, pSet, qSet, initialMask, localNCCs) 
	{
	#pragma  omp for 
#endif	
		
	for(i = 0; i < subsets.size(); i++){ // loop over the affine parameters subsets

		itrIndex++;
		int pIndex = subsets[i].p;
		int qIndex = subsets[i].q;

		int stP = (1 << L-l)*(pIndex-1) + 1;
		int enP = (1 << L-l)*(pIndex) + 1;
		int stQ= (1 << L-l)*(qIndex-1) + 1;
		int enQ= (1 << L-l)*(qIndex) + 1;

		Mat *transformsB, *transformsMask;
		transformsB = new Mat[pSetSize];
		transformsMask = new Mat[pSetSize];

		double s, lambda, theta, h;
		int minWidth, minHeight;
		minWidth = minHeight = 1000000;

		int maxWidth, maxHeight;
		maxWidth = maxHeight = 0;

		for(int k = 0; k < pSetSize; k++){ /*apply each affine transformation on both the region and a mask*/
			s = pqV[pIndex][pSet[k]].p.s;
			lambda = pqV[pIndex][pSet[k]].p.lambda;
			theta = pqV[pIndex][pSet[k]].q.theta;
			h = pqV[pIndex][pSet[k]].q.h;

			
			Mat bTmp = AffineTransformGenerator::generateTransform(region, s, lambda, theta, h); 
			Mat maskTemp = AffineTransformGenerator::generateTransform(initialMask, s, lambda, theta, h);

			transformsB[k] = bTmp;
			transformsMask[k] = maskTemp;
			if(bTmp.rows < minHeight)
				minHeight = bTmp.rows;
			if(bTmp.cols < minWidth)
				minWidth = bTmp.cols;

		}	

		Mat Bi = Mat::zeros(minHeight, minWidth, cType);
		Mat maskT = Mat::zeros(minHeight, minWidth, cType);		

		// compute the average of all transformed patches intersected together
		for(int k = 0 ; k < pSetSize; k++){
			Mat currR = transformsB[k];
			Mat currM = transformsMask[k];
			Mat imRoi =  currR(Rect(currR.cols/2-minWidth/2, currR.rows/2-minHeight/2, minWidth,minHeight));
			Mat maskRoi = currM(Rect(currM.cols/2-minWidth/2, currM.rows/2-minHeight/2, minWidth,minHeight));
			
			transformsB[k] = imRoi.clone();
			transformsMask[k] = maskRoi.clone();
			Bi = Bi + transformsB[k];
			maskT = maskT + transformsMask[k];

		}
		Bi /= (pSetSize*1.0);
		maskT /= (pSetSize*1.0);

		Mat tmp (maskT.size(), CV_8UC3);
		threshold(maskT, maskT, 0.999,1.,THRESH_BINARY);
		maskT.convertTo(tmp,CV_8UC3);
		Mat BiU;
		Bi.copyTo(BiU, tmp);			
	
		Mat mask32FC3(tmp.size(), CV_32FC3);
		tmp.convertTo(mask32FC3, CV_32FC3);
		
		Mat imgSqr; // compute the square of the region (needed for the NCC computation)
		subsets[i].region.convertTo(imgSqr,CV_32FC3);
		imgSqr = imgSqr.mul(imgSqr);

		bool debug = false;


#ifdef DEBUG
		debug = true;
#else
		debug = false;
#endif

		int lNCC = 0;
		Mat response;
		affine_timer.tic();


		response = TemplateMatcher::computeResponseMap(subsets[i].region, BiU, mask32FC3, mask32FC3, imgSqr, debug, lNCC); // perform NCC
		double affineInSub = affine_timer.toc();

		tTime += affineInSub;

		LocalMaximum* localMaxima = new LocalMaximum[params.numPatchLocalMaxima];
		if(response.cols != 0){			 
			// find top n local maxima
			Mat tmp = img.clone();
			
			int returned = 0;
			long elapsedTime = 0;
			localMaxima = TemplateMatcher::findTopNLocalMaxima(response,params.numPatchLocalMaxima,1,returned,elapsedTime );
			
			int j = 0;
			Point tmpLoc = localMaxima[0].loc;
			for(j = 0; j < returned; j++){
				localMaxima[j].p = pIndex;
				localMaxima[j].q = qIndex;
				Point currTmpLoc = localMaxima[j].loc;
				
				localMaxima[j].regionAnchor = localMaxima[j].loc;
				#ifndef SEARCH_ALL_IMG
				localMaxima[j].regionAnchor.x += subsets[i].regionAnchor.x;
				localMaxima[j].regionAnchor.y += subsets[i].regionAnchor.y;
					
				localMaxima[j].loc.x += (BiU.cols/2 + subsets[i].regionAnchor.x);
				localMaxima[j].loc.y += (BiU.rows/2 + subsets[i].regionAnchor.y);
				localMaxima[j].patch = BiU.clone();
			
				#endif	
				topLocalMaximaArr[omp_get_thread_num()].push_back(localMaxima[j]);
			}

			
			#pragma omp critical (REG2)
			{
				localNCCs[omp_get_thread_num()]++;
				if(localMaxima[0].val > maxScore){
					maxScore = localMaxima[0].val;
					topLocalMaximum = localMaxima[0];
					if(returned > 1){
						secondMaxScore = localMaxima[1].val;
						secondMaxLoc = localMaxima[1].regionAnchor;
					}else{
						secondMaxScore = -1;
						secondMaxLoc = Point(-1,-1);
					}
				}
			}
		}

			delete[] localMaxima; localMaxima=0;
			delete[] transformsB; transformsB=0;
			delete[] transformsMask; transformsMask=0;
		
	}
#ifdef USE_OMP
	}
#endif



	for(int i = 0; i < numThreads; i++){
		numNCC += localNCCs[i];
		for(int j = 0; j < topLocalMaximaArr[i].size();j++)
			topLocalMaxima.push_back(topLocalMaximaArr[i][j]);
	}

	if(topLocalMaxima.size() > 0){
		sort(topLocalMaxima.begin(), topLocalMaxima.end(), LocalMaximumComparator());
		topLocalMaximum = topLocalMaxima[0];
	}else{
		topLocalMaximum.loc = Point(-1,-1);
		topLocalMaximum.val = -1;
	}

	return topLocalMaxima;
}

void CoarseToFineAffineSearch::deallocate(){
	img.release();
	region.release();
}
