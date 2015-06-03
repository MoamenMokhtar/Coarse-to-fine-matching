#ifndef STORED_REGION_H_
#define STORED_REGION_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include "AffineTransform.h"
using namespace std;

struct StoredRegion{
	Point p1,p2,p3,p4,p5;
	AffineTransform transform;
	int size;
	int matchingID;

	char* print(){
		char text[3000];
		sprintf(text, "%d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %d",  p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y, p5.x, p5.y, size, transform.s, transform.lambda, transform.theta, transform.h, transform.tx, transform.ty, matchingID);
		//sprintf(text, "%d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %d",  p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y, p5.x, p5.y, size, transform.s, transform.lambda, transform.theta, transform.h, transform.tx, transform.ty, matchingID);
		return text;
	}
};

#endif