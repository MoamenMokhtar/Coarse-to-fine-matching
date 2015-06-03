#ifndef LOCAL_MAXIMUM_H_
#define LOCAL_MAXIMUM_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>

using namespace std;
using namespace cv;
struct LocalMaximum{
	int p;
	int q;
	float val;
	Point loc;
	Mat region;
	Mat patch;
	Point regionAnchor;
};

#endif