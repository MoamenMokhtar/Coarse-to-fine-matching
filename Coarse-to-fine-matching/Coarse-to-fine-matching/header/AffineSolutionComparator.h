#ifndef AFFINE_SOLUTION_COMPARATOR_H_
#define AFFINE_SOLUTION_COMPARATOR_H_
#include "AffineSolution.h"

using namespace std;
using namespace cv;


struct AffineSolutionComparator{
	bool operator()(const AffineSolution& a, const AffineSolution& b){
        return a.dist<b.dist;
    }
};
#endif