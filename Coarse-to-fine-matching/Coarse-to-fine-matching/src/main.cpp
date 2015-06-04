#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>
#include <time.h>
#include "timer.h"
#include <set>
#include "Params.h"
#include "CoarseToFineAffineSearch.h"
#include "CoarseToFineRegionSearch.h"
#include "Utils.h"
#include "timer.h"

#define PI 3.14159265359
#define EPS 0.0001
using namespace std;
using namespace cv;

int main(){

	ifstream paramFile("..\\..\\Params.txt", ifstream::in);
	char inputLine[500];
	char baseDir[500];
	char seqName[500];
	char outputDir[500];
	char srcImgName[500];
	char targetImgName[500];
	char maskName[500];
	int numSamples, affineLevels, regionLevels, numLocalMaxima, downSampleFactor, numPatchLocalMaxima;
	float sMin,sMax,lambdaMin,lambdaMax,hMin,hMax,thetaMin,thetaMax;
	int hasMask;
	float thr1, thr2;
	int samplingCriteria;
	int initialTiling;
	samplingCriteria=-1;
	numSamples = -1;
	affineLevels = -1;
	regionLevels = -1;
	numLocalMaxima = -1;
	sMin = -1.0;
	sMax = -1.0;
	lambdaMin = -1.0;
	lambdaMax = -1.0;
	hMin = -1;
	hMax = -1;
	thetaMin = -100000;
	thetaMax = -100000;
	thr1 = -1;
	thr2 = -1;
	downSampleFactor = -1;
	numPatchLocalMaxima = -1;
	initialTiling = -1;

	if(paramFile.good()){
		while(!paramFile.eof()){
			paramFile.getline(inputLine, 300);	
			sscanf_s (inputLine,"base_folder=%[^\n]s", baseDir,300);
			sscanf_s (inputLine,"output_base_folder=%[^\n]s", outputDir,300);
			sscanf_s (inputLine,"seq_name=%s", seqName,300);
			sscanf_s (inputLine,"reference_img_name=%s", srcImgName, 300);
			sscanf_s (inputLine,"target_img_name=%s", targetImgName, 300);
 			sscanf_s (inputLine,"n=%d", &numSamples);
			sscanf_s (inputLine,"s_min=%f", &sMin);
			sscanf_s (inputLine,"s_max=%f", &sMax);
 			sscanf_s (inputLine,"lambda_min=%f", &lambdaMin);
			sscanf_s (inputLine,"lambda_max=%f", &lambdaMax);
 			sscanf_s (inputLine,"h_min=%f", &hMin);
 			sscanf_s (inputLine,"h_max=%f", &hMax);
 			sscanf_s (inputLine,"theta_min=%f", &thetaMin);
 			sscanf_s (inputLine,"theta_max=%f", &thetaMax);
			sscanf_s (inputLine,"a_l=%d", &affineLevels);
			sscanf_s (inputLine,"r_l=%d", &regionLevels);
			sscanf_s (inputLine,"sampling_criteria=%d", &samplingCriteria);
			sscanf_s (inputLine,"down_sample=%d", &downSampleFactor);
			sscanf_s (inputLine,"initial_tiling=%d", &initialTiling);
			sscanf_s (inputLine,"thrF=%f", &thr1);
			sscanf_s (inputLine,"thrR=%f", &thr2);
			sscanf_s (inputLine,"num_local_maxima=%d", &numLocalMaxima);
			sscanf_s (inputLine,"num_patch_local_maxima=%d", &numPatchLocalMaxima);
			sscanf_s (inputLine,"has_mask=%d", &hasMask);			
			sscanf_s (inputLine,"mask_name=%s", maskName, 300);
		}
		paramFile.close();
	}else{
		cerr << "ERROR: Can't open params file" << endl;
		return -2;
	}

	if(numSamples == -1 || affineLevels ==-1 || numLocalMaxima ==-1 || sMin == -1 || sMax == -1 || lambdaMin == -1 || lambdaMax == -1 || hMin == -1 || hMax == -1 || thetaMin == -1 || thetaMax == -1 || thr1 == -1  || thr2 == -1 || downSampleFactor == -1 || samplingCriteria == -1 || initialTiling == -1){
		cerr << "error in reading the parameters" << endl;
		exit(1);
	}


	thetaMin = thetaMin*PI/180.0; // convert angle from deg to rad
	thetaMax = thetaMax*PI/180.0; // convert angle from deg to rad
	
	char srcImgDir[256];
	sprintf(srcImgDir, "%s\\%s", baseDir, srcImgName);

	Mat baseImg = imread(srcImgDir,1);

	Mat mask;
	if(hasMask){
		char maskDir[256];
		sprintf(maskDir, "%s\\%s", baseDir, maskName);
		mask = imread(maskName,0);
	}else{
		 mask = Mat::ones(baseImg.size(),CV_8UC1);
	}
	
	char secondImgDir[256];
	sprintf(secondImgDir, "%s\\%s", baseDir, targetImgName);
	cout << secondImgDir << endl;
	Mat img2 = imread(secondImgDir, 1);

	Params params = Utils::computeParamsSteps(sMin,sMax,lambdaMin,lambdaMax,hMin,hMax,thetaMin, thetaMax, numSamples); // compute the step size for each parameter
	params.N = numSamples;
	params.affineLevels = affineLevels;
	params.regionLevels = regionLevels;
	params.patchCols = baseImg.cols/(1<<(downSampleFactor));
	params.patchRows = baseImg.rows/(1<<(downSampleFactor));
	params.numLocalMaxima = numLocalMaxima;
	params.numPatchLocalMaxima = numPatchLocalMaxima;
	params.samplingCriteria = samplingCriteria;
	params.initialTiling = initialTiling;
	int l = 1;

	/*build the affine heirarchy*/
	vector<vector<PQ> > pqVec;
	params.samplesInLevels.push_back(Utils::getSubsets(params, 1 << (numSamples/2-affineLevels+l), pqVec,l,affineLevels));
	for(l = 2; l <= affineLevels; l++){
		params.samplesInLevels.push_back(Utils::getSubsets(params, 1 << (numSamples/2-affineLevels+l), params.samplesInLevels[params.samplesInLevels.size()-1],l,affineLevels));
	}

	baseImg.convertTo(baseImg,CV_32FC3);
	img2.convertTo(img2, CV_32FC3);
	CoarseToFineRegionSearch regionSearch;
	string seqStr(seqName);
	int imgNum = 6;
	regionSearch.initialize(baseImg, img2, params, 1000, imgNum*1000+120, affineLevels, hasMask,seqStr, thr1, thr2, mask);
	regionSearch.matchImages();

	

	return 0;
}

