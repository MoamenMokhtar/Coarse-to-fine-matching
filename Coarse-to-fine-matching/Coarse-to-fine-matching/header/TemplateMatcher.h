#ifndef TEMPLATE_MATCHER_H_
#define TEMPLATE_MATCHER_H_

//#define USE_CCORR_NORMED
//#define USE_COEFF_NORMED
//#define USE_SQDIFF_NORMED
//#define USE_OMP2
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <time.h>
#include "Params.h"
#include "LocalMaximum.h"
#include "timer.h"
#include "LocalMaximumComparator.h"
#include "Utils.h"

using namespace std;
using namespace cv;

class TemplateMatcher{
public:
	static Mat doNCC(Mat img, Mat templ, Mat mask, Mat mask2, Mat imgSqr){


		int matchMethod = CV_TM_CCORR;
		int resultRows = img.rows - templ.rows + 1;
		int resultCols = img.cols - templ.cols + 1;
		if(resultRows <= 0 || resultCols <= 0){
			cout <<" >>> " << "resultRows " << resultRows <<" resultCols " << resultCols <<endl;
			Mat res(0,0,0);
			return res;
		}

		Mat CC1, CC2, imgDenom1;
#ifdef USE_OMP2
		#pragma omp parallel
		{
			#pragma omp sections nowait
			{
				#pragma omp section
					matchTemplate(img, templ, CC1, CV_TM_CCORR);
				#pragma omp section				
					matchTemplate(img, mask,CC2,CV_TM_CCORR);
				#pragma omp section					
					matchTemplate(imgSqr,mask2,imgDenom1,CV_TM_CCORR);
				}

		}
#else

		matchTemplate(img, templ, CC1, matchMethod );
		matchTemplate(img, mask,CC2, matchMethod);
		matchTemplate(imgSqr,mask2,imgDenom1, matchMethod);
		
#endif


		double templMean, templStd;
		double sum = 0;
		templMean = 0;
		int numElems = 0;
		for(int i = 0; i < templ.rows; i++){
			for(int j = 0; j < templ.cols; j++){
				if(mask.at<Vec3f>(i,j)[0] == 1){
					for(int c = 0; c < templ.channels(); c++){
						sum += img.at<Vec3f>(i,j)[c]* templ.at<Vec3f>(i,j)[c];
						templMean += templ.at<Vec3f>(i,j)[c];
						numElems++;
					}
				}
			}
		}
		templMean /= (numElems);
		templStd = 0;
		for(int i = 0; i < templ.rows; i++){
			for(int j = 0; j < templ.cols; j++){
				if(mask.at<Vec3f>(i,j)[0] == 1){
					for(int c = 0; c < templ.channels(); c++){
						templStd += (((templ.at<Vec3f>(i,j)[c]) - templMean)*((templ.at<Vec3f>(i,j)[c]) - templMean));
					}
				}
			}
		}
		templStd /= (numElems - 1);
		templStd = sqrt(templStd);

		Mat imgDenom2;
		imgDenom2 = CC2.mul(CC2);

		CC2 *= templMean;

		Mat numerator;
		numerator = CC1 - CC2;
		imgDenom2 /= numElems;

		Mat imgDenom;
		imgDenom = imgDenom1 - imgDenom2;

		double templDenom = sqrt((numElems-1)*1.0f)*templStd;

		Mat result(resultRows, resultCols, CV_32FC1);
		Mat denom;
		for(int i = 0; i < result.rows; i++){
			for(int j = 0; j < result.cols; j++){
				float imgD = sqrt(max(imgDenom.at<float>(i,j),0));	
				result.at<float>(i,j) = imgD==0?0:(numerator.at<float>(i,j))/(imgD * templDenom);
			}
		}
		return result;
	}

	static Mat computeResponseMap(Mat img, Mat templ,Mat mask,  Mat mask2, Mat imgSqr, bool debug, int& numNCC){
		int matchMethod = CV_TM_CCORR;
		numNCC++;
		int resultRows = img.rows - templ.rows + 1;
		int resultCols = img.cols - templ.cols + 1;
		if(resultRows <= 0 || resultCols <= 0){
			Mat res(0,0,0);
			return res;
		}
		return doNCC(img, templ, mask, mask2, imgSqr);
	}

	static LocalMaximum* findTopNLocalMaxima(Mat res, int n, bool blackout, int& numReturnedLocalMaxima, long& timeE){
		
		int numMaxima = 0;
		LocalMaximum* localMaxima = new LocalMaximum[res.rows*res.cols];
		float* p, *pU, *pD;
		timer e_timer;
		
		for(int i = 0; i < res.rows; i++){
			e_timer.tic();		
			p = res.ptr<float>(i);
			if(i > 0)
				pU = res.ptr<float>(i-1);
			if(i < res.rows-1)
				pD = res.ptr<float>(i+1);
			
			for(int j = 0; j < res.cols; j++){
				
				int numNeighbors = 0;
				bool isLocalMaximum = true;
				float val = p[j];
				if(val == 0)
					continue;
				if(i > 0 && pU[j] > val){
					continue;
				}
				if(i < res.rows-1 && pD[j] > val){
					continue;
				}
				if(j > 0 && p[j-1] > val){
					continue;
				}
				if(j < res.cols-1 && p[j + 1] > val){
					continue;
				}				
				LocalMaximum lm;
				lm.val = val;
				lm.loc = Point(j,i);
				localMaxima[numMaxima++] = lm;
			}
			timeE += e_timer.toc();
		}

		LocalMaximum* topNLocalMaxima = new LocalMaximum[n];
		numReturnedLocalMaxima = 0;

		sort(localMaxima, localMaxima+numMaxima,LocalMaximumComparator());

		int k = 0;
		for(int i = 0; k < n && i < numMaxima; i++){
			LocalMaximum lm = localMaxima[i];
			if(res.at<float>(lm.loc.y, lm.loc.x) != -1){ //  not blacked out
				topNLocalMaxima[k] = lm;
				if(blackout)
					Utils::blackoutAroundMaximum(res, lm);
				k++;
			}
		}
		numReturnedLocalMaxima = k;
		
		delete[] localMaxima; localMaxima = 0;
		return topNLocalMaxima;
	}
};

#endif
