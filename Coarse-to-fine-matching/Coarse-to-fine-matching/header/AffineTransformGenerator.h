#ifndef AFFINE_TRANSFORM_GENERATOR_H_
#define AFFINE_TRANSFORM_GENERATOR_H_
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include "Utils.h"
#ifndef PI
#define PI 3.14159265359
#endif



using namespace std;
using namespace cv;


class AffineTransformGenerator{
public:
	static Mat generateTransform(Mat src, double s, double lambda, double theta, double h){
		int numPoints = 4;
		Point ps[4];
		ps[0] = Point(0,0);
		ps[1] = Point(src.cols,0);
		ps[2] = Point(src.cols,src.rows);
		ps[3] = Point(0,src.rows);
		Mat* hps = new Mat[numPoints];

		for(int i = 0; i < numPoints; i++){
			hps[i].create(3, 1, CV_32FC1);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
			hps[i].at<float>(2,0) = 1;
		}

		Mat affMat(2,3,CV_64FC1);   
		affMat.at<double>(0,0) = s;
		affMat.at<double>(0,1) = 0;
		affMat.at<double>(0,2) = /*(-src.cols*s+src.cols)/2*/0;
		affMat.at<double>(1,0) = 0;
		affMat.at<double>(1,1) =s*lambda;
		affMat.at<double>(1,2) = /*(-src.rows*s*lambda +src.rows)/2*/0 ;


		for(int i = 0; i < numPoints; i++){
			Mat tmp(2,1,CV_32FC1);
			
			for(int r = 0; r < affMat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < affMat.cols; c++){
					tmp.at<float>(r,0) += affMat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		Mat res;
		warpAffine( src, res, affMat, Size(src.cols*s,src.rows*s*lambda) );

		affMat.at<double>(0,0) = 1;
		affMat.at<double>(0,1) = h;
		affMat.at<double>(0,2) = h<0?res.rows*abs(h):0;
		affMat.at<double>(1,0) = 0;
		affMat.at<double>(1,1) = 1;
		affMat.at<double>(1,2) = 0;


		for(int i = 0; i < numPoints; i++){
			Mat tmp(2,1,CV_32FC1);
			for(int r = 0; r < affMat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < affMat.cols; c++){
					tmp.at<float>(r,0) += affMat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		Mat res2 = Mat::zeros(res.rows, res.cols+(res.rows*abs(h)),res.type());
		warpAffine( res, res2, affMat, res2.size() );


		double angle = theta*180.0/PI;
		int width = max(res2.rows*abs(sin(angle*PI/180.0)) + res2.cols*cos(angle*PI/180.0), res2.cols);
		int height = max(res2.rows*cos(angle*PI/180.0) + res2.cols*abs(sin(angle*PI/180.0)), res2.rows);
		Mat res3 = Mat::zeros(height,width,res2.type());
		int shiftX = (res3.cols - res2.cols)/2;
		int shiftY = (res3.rows - res2.rows)/2;
		shiftY = shiftY < 0 ? 0 : shiftY;
		Rect rct(Point(shiftX , shiftY),Point(res2.cols + shiftX, res2.rows + shiftY));
		Mat rr1 = res3(rct);
		addWeighted(rr1,0,res2,1,0,rr1);

		Point center = Point( res3.cols/2, res3.rows/2);
		Mat rot_mat = getRotationMatrix2D( center, angle, 1 );
		
		for(int i = 0; i < numPoints; i++){
			hps[i].at<float>(0,0) += shiftX;
			hps[i].at<float>(1,0) += shiftY;
			Mat tmp(2,1,CV_32FC1);
			for(int r = 0; r < rot_mat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < rot_mat.cols; c++){
					tmp.at<float>(r,0) += rot_mat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		int minX, minY, maxX, maxY;
		minX = minY = 10000;
		maxX = maxY = -10000;

		for(int i = 0; i < 4; i++){
			if(ps[i].x < minX)
				minX = ps[i].x;
			if(ps[i].y < minY)
				minY = ps[i].y;
			if(ps[i].x >= maxX)
				maxX = ps[i].x;
			if(ps[i].y >= maxY)
				maxY = ps[i].y;
		}	
		minX = (minX < 0) ? 0 : minX;
		minY = (minY < 0) ? 0 : minY;
		warpAffine( res3, res3, rot_mat, res3.size());
		Rect roi(Point(minX, minY), Point(maxX, maxY));
		res3 = res3(roi);
		delete[] hps;
		return res3;
	}


static Mat generateTransform2(Mat src, double s, double lambda, double theta, double h, Point* ps, int numPoints){
		
		Mat* hps = new Mat[numPoints];

		for(int i = 0; i < numPoints; i++){
			hps[i].create(3, 1, CV_32FC1);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
			hps[i].at<float>(2,0) = 1;
		}

		Mat affMat(2,3,CV_64FC1);   
		affMat.at<double>(0,0) = s;
		affMat.at<double>(0,1) = 0;
		affMat.at<double>(0,2) = /*(-src.cols*s+src.cols)/2*/0;
		affMat.at<double>(1,0) = 0;
		affMat.at<double>(1,1) =s*lambda;
		affMat.at<double>(1,2) = /*(-src.rows*s*lambda +src.rows)/2*/0 ;


		for(int i = 0; i < numPoints; i++){
			Mat tmp(2,1,CV_32FC1);
			
			for(int r = 0; r < affMat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < affMat.cols; c++){
					tmp.at<float>(r,0) += affMat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		Mat res;
		warpAffine( src, res, affMat, Size(src.cols*s,src.rows*s*lambda) );

		affMat.at<double>(0,0) = 1;
		affMat.at<double>(0,1) = h;
		affMat.at<double>(0,2) = h<0?res.rows*abs(h):0;
		affMat.at<double>(1,0) = 0;
		affMat.at<double>(1,1) = 1;
		affMat.at<double>(1,2) = 0;


		for(int i = 0; i < numPoints; i++){
			Mat tmp(2,1,CV_32FC1);
			for(int r = 0; r < affMat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < affMat.cols; c++){
					tmp.at<float>(r,0) += affMat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		Mat res2 = Mat::zeros(res.rows, res.cols+(res.rows*abs(h)),res.type());
		warpAffine( res, res2, affMat, res2.size() );


		double angle = theta*180.0/PI;
		int width = max(res2.rows*abs(sin(angle*PI/180.0)) + res2.cols*cos(angle*PI/180.0), res2.cols);
		int height = max(res2.rows*cos(angle*PI/180.0) + res2.cols*abs(sin(angle*PI/180.0)), res2.rows);
		Mat res3 = Mat::zeros(height,width,res2.type());
		int shiftX = (res3.cols - res2.cols)/2;
		int shiftY = (res3.rows - res2.rows)/2;
		shiftY = shiftY < 0 ? 0 : shiftY;
		Rect rct(Point(shiftX , shiftY),Point(res2.cols + shiftX, res2.rows + shiftY));
		Mat rr1 = res3(rct);
		addWeighted(rr1,0,res2,1,0,rr1);

		Point center = Point( res3.cols/2, res3.rows/2);
		Mat rot_mat = getRotationMatrix2D( center, angle, 1 );
		
		for(int i = 0; i < numPoints; i++){
			hps[i].at<float>(0,0) += shiftX;
			hps[i].at<float>(1,0) += shiftY;
			Mat tmp(2,1,CV_32FC1);
			for(int r = 0; r < rot_mat.rows; r++){
				tmp.at<float>(r,0) = 0;
				for(int c = 0; c < rot_mat.cols; c++){
					tmp.at<float>(r,0) += rot_mat.at<double>(r,c)*hps[i].at<float>(c,0);
				}
			}
			ps[i].x = tmp.at<float>(0,0);
			ps[i].y = tmp.at<float>(1,0);
			hps[i].at<float>(0,0) = ps[i].x;
			hps[i].at<float>(1,0) = ps[i].y;
		}

		int minX, minY, maxX, maxY;
		minX = minY = 10000;
		maxX = maxY = -10000;

		for(int i = 0; i < 4; i++){
			if(ps[i].x < minX)
				minX = ps[i].x;
			if(ps[i].y < minY)
				minY = ps[i].y;
			if(ps[i].x >= maxX)
				maxX = ps[i].x;
			if(ps[i].y >= maxY)
				maxY = ps[i].y;
		}	
		warpAffine( res3, res3, rot_mat, res3.size());
		minX = minX >=0 ? minX:0;
		minY = minY >=0 ? minY:0;
		maxX = maxX < res3.cols ? maxX : res3.cols-1;
		maxY = maxY < res3.rows ? maxY : res3.rows-1;
		Rect roi(Point(minX, minY), Point(maxX, maxY));
		res3 = res3(roi);
		delete[] hps;
		return res3;
	}
};


#endif