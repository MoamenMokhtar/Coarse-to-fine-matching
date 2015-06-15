#ifndef AFFINE_SOLUTION_H_
#define AFFINE_SOLUTION_H_

#include "AffineTransform.h"

struct AffineSolution{
public:
	Rect roi;
	Point loc;
	AffineTransform affineTrans;
	double matchScore;
	double ratio;
	Mat region;
	float dist;
};


#endif