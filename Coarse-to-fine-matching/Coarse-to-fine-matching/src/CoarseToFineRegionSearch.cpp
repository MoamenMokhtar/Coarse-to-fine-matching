#include "CoarseToFineRegionSearch.h"
#define NUM_SPLITS 4
void CoarseToFineRegionSearch::initialize(Mat _I1, Mat _I2, Params _params, int _i1,int _i2, int _L, int _hasMask, string _seq, double _thr1, double _thr2, Mat _mask){
	I1 = _I1;
	I2 = _I2;
	i1 = _i1;
	i2 = _i2;
	merged = Mat::zeros(I1.rows,I1.cols+I2.cols,CV_8UC3);

	Rect I1ROI(Point(0,0),Point(I1.cols, I1.rows));
	Mat I1RoiM = merged(I1ROI);
	Mat I1C, I2C;
	I1.convertTo(I1C, CV_8UC3);
	I2.convertTo(I2C, CV_8UC3);
	addWeighted(I1RoiM,0,I1C,1,0.0,I1RoiM);
	
	Rect I2ROI(Point(I1.cols,0),Point(I1.cols+I2.cols, I2.rows));
	Mat I2RoiM = merged(I2ROI);
	
	addWeighted(I2RoiM,0,I2C,1,0.0,I2RoiM);
	L = _L;
	params = _params;	
	hasMask = _hasMask;
	seq = _seq;
	t1 = _thr1;
	t2 = _thr2;
	mask = _mask;
}

void CoarseToFineRegionSearch::matchImages(){
	

	char dirName[256];
	char systemCall[256];
	
	sprintf(dirName,"resultFiles\\%s",seq.c_str());
	sprintf(systemCall,"mkdir %s",dirName);
	system(systemCall);


		ofstream r1; // regions in the 
		ofstream r2;
		ofstream pF;
		
		char paramsName[512];
		sprintf(paramsName, "%s\\params_%.2f_%.2f_%.2f_%.2f_%d_lms_%d.txt", dirName, t1, t2, params.sMin,params.thetaMin, params.affineLevels, params.numLocalMaxima);
		pF.open(paramsName);
		pF<<"N="<<params.N<<endl;
		pF<<"L="<<params.affineLevels<<endl;
		pF<<"num Lms ="<<params.numLocalMaxima<<endl;
		pF<< "sMin="<<params.sMin<<endl;
		pF<< "sMax="<<params.sMax<<endl;
		pF<< "lambdaMin="<<params.lambdaMin<<endl;
		pF<< "lambdaMax="<<params.lambdaMax<<endl;
		pF<< "thetaMin="<<params.thetaMin*180.0/PI<<endl;
		pF<< "thetaMax="<<params.thetaMax*180.0/PI<<endl;
		pF<< "hMin="<<params.hMin<<endl;
		pF<< "hMax="<<params.hMax<<endl;
		pF<< "patchCols="<<params.patchCols<<endl;
		pF<< "patchRows="<<params.patchRows<<endl;

		char f1Name[512];
		char f2Name[512];
		sprintf(f1Name,"%s\\img%d_%.2f_%.2f.txt",dirName,1,t1,t2);
		sprintf(f2Name,"%s\\img%d_%.2f_%.2f.txt",dirName,2,t1,t2);
		r1.open(f1Name);
		r2.open(f2Name);
		

		int numNCC = 0;
		merged = Utils::mergeTwoImages(I1,I2); 

		timer algo_timer;
		algo_timer.tic();
		vector<RegionSplit> regions, region2;
		RegionSplit s1;
		s1.anchor = Point(0,0);
		s1.parentAnchor = Point(0,0);
		s1.img = I1;
		s1.parentImg= I1;
		RegionSplit splits[NUM_SPLITS];
		RegionSplit splits2[NUM_SPLITS];

		Utils::splitPatch(s1, splits,false);
		for(int i = 0; i < NUM_SPLITS; i++){
			splits[i].parentImg = splits[i].img;
			splits[i].parentAnchor = splits[i].anchor;
			regions.push_back(splits[i]);
		}

		for(int i = 0; i < NUM_SPLITS; i++){
			RegionSplit r = regions[regions.size()-1];
			regions.pop_back();
			Utils::splitPatch(r, splits2,false);
			for(int j = 0; j < NUM_SPLITS; j++){
				splits2[j].parentImg = splits2[j].img;
				splits2[j].parentAnchor = splits2[j].anchor;
				splits2[j].regionLevel = 1;
				region2.push_back(splits2[j]);
			}
		}
		regions = region2;

		int inum = 1;
		int matchingId = -1;
		int attemptedPatchID = -1;

		int numPassed = 0;
			// start region matching
		while(!regions.empty()){
			RegionSplit r = regions[regions.size()-1];
			regions.pop_back();
			if(r.regionLevel > params.regionLevels)
				continue;
			Mat split = r.img;		
			CoarseToFineAffineSearch affineSearch;
			int canUse;
			if(r.match){ // the patch is not under rematching
				canUse = affineSearch.initialize(I2, split, params, t1); // initialize the affine heirarchy search				
			}else{
				// transform the current patch with its parent affine transformation.
				Point ps2[5];
				ps2[0] = Point(r.anchor.x-r.parentAnchor.x, r.anchor.y-r.parentAnchor.y);							
				ps2[1] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols, r.anchor.y-r.parentAnchor.y);	
				ps2[2] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols, r.anchor.y-r.parentAnchor.y+r.img.rows);
				ps2[3] = Point(r.anchor.x-r.parentAnchor.x, r.anchor.y-r.parentAnchor.y+r.img.rows);
				ps2[4] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols/2, r.anchor.y-r.parentAnchor.y+r.img.rows/2);
				Mat transformed = AffineTransformGenerator::generateTransform2(r.parentImg,r.affTrans.s,r.affTrans.lambda,r.affTrans.theta,r.affTrans.h, ps2,5);
				int minX, minY, maxX, maxY;
				minX = minY = 10000;
				maxX = maxY = -10000;

				for(int i = 0; i < 5; i++){
					ps2[i].x += r.affTrans.tx;
					ps2[i].y += r.affTrans.ty;	
					if(ps2[i].x < minX)
						minX = ps2[i].x;
					if(ps2[i].y < minY)
						minY = ps2[i].y;
					if(ps2[i].x >= maxX)
						maxX = ps2[i].x;
					if(ps2[i].y >= maxY)
						maxY = ps2[i].y;
				}	
				int offset = 30;
				int sX = (minX-offset >= 0) ? minX-offset : 0;
				int sY = (minY-offset >= 0)?minY-offset : 0;
				int width = (maxX-sX)+offset*2;
				int height = (maxY-sY)+offset*2;
				
				width = (sX+width > I2.cols)? I2.cols-sX : width;
				height = (sY+height> I2.rows)? I2.rows-sY : height;
				Rect roi(sX, sY, width, height);
				Mat searchRegion = I2(roi).clone();
				canUse = affineSearch.initialize(searchRegion, r.img.clone(), params, t1);
					
			}
			if(!canUse){
				continue;
			}
					

			attemptedPatchID++;
			/*store patch in source image info*/
			StoredRegion sr;
			sr.p1 = r.anchor;
			sr.p2 = Point(r.anchor.x + r.img.cols, r.anchor.y);
			sr.p3 = Point(r.anchor.x, r.anchor.y + r.img.rows);
			sr.p4 = Point(r.anchor.x + r.img.cols, r.anchor.y + r.img.rows);
			sr.p5 = Point(r.anchor.x + r.img.cols/2, r.anchor.y + r.img.rows/2);
			sr.size = r.img.rows*r.img.cols;
			AffineTransform affTrans = AffineTransform();					
			sr.matchingID = -1;
			sr.transform = affTrans;

			if(hasMask && Utils::isBackground(mask, r)){
				regions.push_back(r);
				continue;
			}
			double s, lambda, theta, h;
			int canFind;
			if(r.match){
				canFind = affineSearch.computeBestAffineTransform(s, lambda, theta, h, numNCC);
				/*retain the transformation parameters*/
				affTrans.s = s;
				affTrans.lambda = lambda;
				affTrans.theta = theta;
				affTrans.h = h;							
				affTrans.tx = affineSearch.topLocalMaximum.regionAnchor.x*affineSearch.widthResizeFactor; 
				affTrans.ty = affineSearch.topLocalMaximum.regionAnchor.y*affineSearch.heightResizeFactor;
				affTrans.pIdx = affineSearch.topLocalMaximum.p;
				r.affTrans = affTrans;
			}
			else{
				canFind = affineSearch.computeBestAffineTransform(s, lambda, theta, h, numNCC, 1, r.affTrans);
			}
						
			cout << "Max score " <<affineSearch.maxScore << " Second Max score "  << affineSearch.secondMaxScore << " second/top " << affineSearch.secondMaxScore/affineSearch.maxScore << endl;
			double topScore = affineSearch.maxScore ;

#ifdef USE_SQDIFF_NORMED
				double ratio = affineSearch.maxScore/affineSearch.secondMaxScore;
				if(!canFind || topScore > T1){
#else
				double ratio = affineSearch.secondMaxScore/affineSearch.maxScore;

				if(topScore < t1){ // didn't pass the first threshold test
#endif								

					if(r.match){ // split it and try to match its children with the same parameters
						r1<< sr.print() << endl;
						Utils::splitPatch(r,splits,1);
						for(int i = 0; i < NUM_SPLITS;i++){
							splits[i].affTrans = r.affTrans;
							splits[i].parentImg = r.img.clone();
							splits[i].parentAnchor = r.anchor;
							splits[i].regionLevel = r.regionLevel+1;
							regions.push_back(splits[i]);
						}
					}
					else{ // rematch the patch
						r.match = true;
						r.parentImg = r.img;
						r.parentAnchor = r.anchor;
						regions.push_back(r);
					}
				}else{ // passed the first threshold test
					matchingId++;
					bool accepted = false;
					int red, green, blue, thickness;
					if(ratio < t2){
						red = rand()%255;
						green = rand()%255;
						blue = rand()%255;
						accepted = true;
						thickness = 3;
					}else{
						red = 0;
						green = 0;
						blue = 0;
						thickness = 1;
						accepted = false;
					}
					LocalMaximum selected = affineSearch.topLocalMaximum;
					sr.matchingID = matchingId;
					sr.transform = r.affTrans;

					r1<< sr.print() << endl;	

					/*For the display of the corresponding region*/
					Point ps[5];
					ps[0] = Point(r.anchor.x-r.parentAnchor.x, r.anchor.y-r.parentAnchor.y);							
					ps[1] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols, r.anchor.y-r.parentAnchor.y);	
					ps[2] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols, r.anchor.y-r.parentAnchor.y+r.img.rows);
					ps[3] = Point(r.anchor.x-r.parentAnchor.x, r.anchor.y-r.parentAnchor.y+r.img.rows);
					ps[4] = Point(r.anchor.x-r.parentAnchor.x+r.img.cols/2, r.anchor.y-r.parentAnchor.y+r.img.rows/2);

					Mat tmp = AffineTransformGenerator::generateTransform2(r.parentImg,s,lambda,theta,h, ps,5);

					for(int i = 0; i < 5; i++){
						ps[i].x += r.affTrans.tx;
						ps[i].y += r.affTrans.ty;	
						if(accepted)
							r2 <<ps[i].x <<" " << ps[i].y <<" ";
					}	

					if(accepted){ // store the patch if it passed the ratio test
						r2<< topScore <<" " << ratio<<" " << attemptedPatchID;
						r2<<endl;
					}
					for(int i = 0; i < 4; i++){
						ps[i].x += I1.cols;
					}	

					
					if(accepted){
						for(int i = 0; i < 4; i++)
							line(merged,ps[i], ps[(i+1)%4],Scalar(red, green, blue), thickness,CV_AA);								
						rectangle(merged, r.anchor, Point(r.anchor.x + r.img.cols, r.anchor.y + r.img.rows), Scalar(red,green,blue,0), thickness, 8, 0);
					}
					imwrite("chosen\\split__"+Utils::IntToStr(selected.val)+".png", r.img);
					char mergedName[512];
					sprintf(mergedName,"%s\\merged%d_%.2f_%.2f_%d_Lms_%d.png",dirName,i1,t1,t2,L, params.numLocalMaxima);
					imwrite(mergedName, merged);
					
				}
		}
		
		r1.close();
		r2.close();
		pF<<"MatchingTime="<<algo_timer.toc()/1000 <<endl;
		pF<<"Num NCC ="<<numNCC <<endl;
		pF.close();
		

		// deallocate resources
		for(int i = 0; i < regions.size(); i++)
			regions[i].img.release();
		regions.clear();

}


void CoarseToFineRegionSearch::deallocate(){
	I1.release();
	I2.release();
}