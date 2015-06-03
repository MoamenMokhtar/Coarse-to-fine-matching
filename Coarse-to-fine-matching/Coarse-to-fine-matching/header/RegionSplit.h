#ifndef REGION_SPLIT_H_
#define REGION_SPLIT_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include "AffineTransform.h"

using namespace std;
using namespace cv;
struct RegionSplit{
public:
	Point anchor;
	Mat img, parentImg, parentRegion;
	AffineTransform affTrans;
	Mat affineMat;
	int regionLevel;
	Point parentAnchor;
	bool match;
};

#endif