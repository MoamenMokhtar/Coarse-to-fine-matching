#ifndef PARAMS_H_
#define PARAMS_H_


#include <opencv2\core\core.hpp>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "PQ.h"
using namespace std;
using namespace cv;
/*
This struct holds the parameters used in the affine heirarchy
*/

struct Params{
	int N;
	int numLocalMaxima, numPatchLocalMaxima;
	int numSamples, affineLevels, regionLevels;
	double sMin,sStep, sMax;
	double lambdaMin, lambdaStep, lambdaMax;
	double hMin, hStep, hMax;
	double thetaMin, thetaStep, thetaMax;
	int samplingCriteria;
	int sAxis, thetaAxis, hAxis, lambdaAxis;
	int patchRows;
	int patchCols;
	int initialTiling;

	vector<vector<vector<PQ> > > samplesInLevels;

};


#endif