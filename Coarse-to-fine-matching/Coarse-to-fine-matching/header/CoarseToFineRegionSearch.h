#ifndef COARSE_TO_FINE_REGION_SEARCH_H_
#define COARSE_TO_FINE_REGION_SEARCH_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
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
#include "CoarseToFineAffineSearch.h"
#include "RegionSplit.h"

#include "StoredRegion.h"
class CoarseToFineRegionSearch{
public:
	Params params;
	Mat I1, I2, merged;
	Mat mask;
	int i1,i2;
	int L;
	int hasMask;
	string seq;
	double t1, t2;
	CoarseToFineRegionSearch(){}
	~CoarseToFineRegionSearch(){}
	void deallocate();

	void initialize(Mat, Mat, Params, int, int, int, int, string, double, double, Mat);
	void matchImages();

};


#endif