#ifndef COARSE_TO_FINE_AFFINE_SEARCH_H_
#define COARSE_TO_FINE_AFFINE_SEARCH_H_
//#define DEBUG
#ifndef DEBUG
#define USE_OMP
#endif
#define RESIZE
#define SEARCH_ALL_IMG
#undef SEARCH_ALL_IMG
#define USE_PYRAMID
//#undef USE_PYRAMID

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include <time.h>
#include "Params.h"
#include "LocalMaximum.h"
#include "AffineTransformGenerator.h"
#include "timer.h"
#include "TemplateMatcher.h"
#include <set>
#include "Utils.h"
#include <omp.h>
#include <math.h>
#include "PQ.h"
#include "P.h"
#include "Q.h"
#include "AffineTransform.h"
#include "AffineSolution.h"
#include <unordered_set>
using namespace std;
using namespace cv;


class CoarseToFineAffineSearch{
public:
	Params params;
	Mat img, region;
	double T1;
	double maxScore, secondMaxScore;
	Point secondMaxLoc;
	LocalMaximum topLocalMaximum;
	unordered_set<int> settingSet;
	double widthResizeFactor, heightResizeFactor;
	CoarseToFineAffineSearch(){
	}
	~CoarseToFineAffineSearch(){
	}

	int initialize(Mat _img, Mat _region, Params _params, double _T1);

	int computeBestAffineTransform(double& cS, double& cLambda, double& cTheta, double& cH, int& numNCC, int cL, bool usePrevAffTrans=0, AffineTransform affTrans = AffineTransform());
	vector<LocalMaximum> computeBestMatches(vector<LocalMaximum> subsets, int l, vector<vector<PQ> > pqV, int& numNCC);
	int backtrackTheAffineHierarchy(double& cS, double& cLambda, double& cTheta, double& cH, int& numNCC, int sL, int eL, AffineSolution& prevSolution);
	void deallocate();
};
#endif