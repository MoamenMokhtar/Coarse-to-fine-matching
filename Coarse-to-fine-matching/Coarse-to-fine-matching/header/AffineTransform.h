#ifndef AFFINE_TRANSFORM_H_
#define AFFINE_TRANSFORM_H_
#define INIT_VAL -1000000

using namespace std;

class AffineTransform{
public:
	double s, lambda, theta, h, tx, ty;
	int pIdx;
	AffineTransform(){
		s = INIT_VAL;
		lambda = INIT_VAL;
		theta = INIT_VAL;
		h = INIT_VAL;
		tx = INIT_VAL;
		ty = INIT_VAL;
		pIdx = -1;
	}
	~AffineTransform(){}
};

#endif