#ifndef UTILS_H_
#define UTILS_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\core\core.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "LocalMaximum.h"
#include <math.h>
#include <time.h>
#include "timer.h"
#include <set>
#include "Params.h"
#include "PQ.h"
#include "RegionSplit.h"
using namespace std;
class Utils{
public:
	/*in the case of segmented foreground object, this function checks whether the patch in hand belongs to the background or not*/
	static bool isBackground(Mat mask, RegionSplit split){
		
		int offset = 5;
		int shiftX[] = {0, offset, -offset, 0, 0};
		int shiftY[] = {0, 0, 0, offset, -offset};
		for(int i = 0; i < 5; i++){
			int xLoc = split.anchor.x + split.img.cols/2+shiftX[i];
			int yLoc = split.anchor.y + split.img.rows/2+shiftY[i];
			if(mask.at<uchar>(yLoc, xLoc) == 0)
				return false;
		}
		
		return true;
	}
	static void doAntiAliasing(Mat& src, Mat& dst){

		GaussianBlur(src, dst,Size(3,3),0.2,0.2);

	}

	static Rect getROI(Point p, Mat img, int templCols, int templRows, Point& regionAnchor, float factor){

		int width = templCols*factor;
		int height = templRows*factor;
		int x = p.x - width/2;
		int y = p.y - height/2;
		regionAnchor.x = x;
		regionAnchor.y = y;
		if(x < 0){
			x = 0;
			regionAnchor.x = 0;
		}
		if(y < 0){
			y = 0;
			regionAnchor.y = 0;
		}
		if(x + width > img.cols)
			width = img.cols - x;
		if(y + height > img.rows)
			height = img.rows - y;
		return Rect(x, y, width, height);

	}
	/*given an image and a local maximum, extracts an ROI from the image around the location of the local maximum*/
	static Rect getROI(LocalMaximum lm, Mat img, int templCols, int templRows, Point& regionAnchor){

		int width = templCols*2.5;
		int height = templRows*2.5;
		int x = lm.loc.x - width/2;
		int y = lm.loc.y - height/2;
		regionAnchor.x = x;
		regionAnchor.y = y;
		if(x < 0){
			x = 0;
			regionAnchor.x = 0;
		}
		if(y < 0){
			y = 0;
			regionAnchor.y = 0;
		}
		if(x + width > img.cols)
			width = img.cols - x;
		if(y + height > img.rows)
			height = img.rows - y;
		return Rect(x, y, width, height);

	}
	/*given a local maximum and an image, it draws a cross on the image at the location of the local maximum*/
	static void drawCrossAtMaximum(Mat& result1, LocalMaximum lm, Scalar color, int thickness){
		int width = 10;
		int height = 10;
		int x = lm.loc.x - width/2;
		int y = lm.loc.y - height/2;
		if(x < 0)
			x = 0;
		if(y < 0)
			y = 0;
		if(x + width > result1.cols)
			width = result1.cols - x;
		if(y + height > result1.rows)
			height = result1.rows - y;

		line(result1, Point(x,lm.loc.y), Point(x+width,lm.loc.y),color,thickness);
		line(result1, Point(lm.loc.x,y), Point(lm.loc.x,y+height),color,thickness);
	}

	static void blackoutAroundMaximum(Mat& result1, LocalMaximum lm){
		double ratio = 4;
		int width = (result1.cols < ratio )? result1.cols : (result1.cols/ratio);
		int height =(result1.rows < ratio )? result1.rows : (result1.rows/ratio);
		int x = lm.loc.x - width/2;
		int y = lm.loc.y - height/2;
		if(x < 0)
			x = 0;
		if(y < 0)
			y = 0;
		if(x + width > result1.cols)
			width = result1.cols - x;
		if(y + height > result1.rows)
			height = result1.rows - y;

		rectangle(result1,Rect(x,y,width,height),Scalar::all(-1),-1);
	}

	template <class T>
	static string IntToStr(T number){
		stringstream ss;//create a stringstream
		ss << number;//add number to the stream
		return ss.str();//return a string with the contents of the stream
	}

	/*computes the step size for each parameter*/
	static Params computeParamsSteps(double minS, double maxS, double minLambda, double maxLambda, double minH, double maxH, double minTheta, double maxTheta, int numSamples){

		bool twoParams = false;
		if(minLambda == maxLambda && minH == maxH)
				twoParams = true;

		int pExpo = numSamples/2 ;
		int qExpo = numSamples/2 + numSamples%2;
		double pAxis = 1 << (pExpo);
		double qAxis = 1 << (qExpo);
		double sAxis;
		double lambdaAxis;
		double thetaAxis;
		double hAxis;
		if(!twoParams){
			sAxis = pAxis/2;
			lambdaAxis = pAxis/2;
			thetaAxis = qAxis/2;
			hAxis = qAxis/2;
		}else{
			sAxis = pAxis;
			lambdaAxis = 0;
			thetaAxis = qAxis;
			hAxis = 0;
		}


		double sStep = (maxS - minS)*1.0/(sAxis - 1);
		double lambdaStep = (maxLambda - minLambda)*1.0/(lambdaAxis - 1);
		double hStep = (maxH - minH)*1.0/(hAxis - 1);
		double thetaStep = (maxTheta - minTheta)*1.0/(thetaAxis - 1);	

		Params params;
		params.sMin = minS;
		params.sStep = sStep;
		params.sMax = maxS;
		params.hMin = minH;
		params.hStep = hStep;
		params.hMax = maxH;
		params.lambdaMin = minLambda;
		params.lambdaStep = lambdaStep;
		params.lambdaMax = maxLambda;
		params.thetaMin = minTheta;
		params.thetaStep = thetaStep;
		params.thetaMax = maxTheta;
		params.numSamples = numSamples;
		params.sAxis = sAxis ;
		params.thetaAxis = thetaAxis ;
		params.hAxis = hAxis;
		params.lambdaAxis = lambdaAxis;

		return params;

	}

	/*this function samples the parameters space and creates the subsets that form the affine heirarchy*/
	static vector<vector<PQ> > getSubsets(Params params, int numSubsets,vector<vector<PQ> > pqVec, int l, int L) {
	
		int kp = 0;
		int kq = 0; 
		bool twoParams = false;
		if( l == 1){ // for the top level

			int movesPSize;
			int movesQSize;

			if(params.lambdaMin == params.lambdaMax && params.hMin == params.hMax)
				twoParams = true;

			if(!twoParams){
				movesPSize = params.sAxis + params.lambdaAxis;
				movesQSize = params.thetaAxis+params.hAxis;  
			}else{
				movesPSize = params.sAxis*2;
				movesQSize = params.thetaAxis*2;
			}

			int subsetSize = 1 << (L-l);

			pair<int,int>* movesP = new pair<int,int>[movesPSize];
			pair<int,int>* movesQ = new pair<int,int>[movesQSize];

			int j = 0;

			


			if(params.samplingCriteria){ 
				/*                                            ^
				 sample the space in the shape of 4X	      |XX
															  |XX
															  ----> */

				int xStart = 0;
				int yStart = 0;

				for(int k = 0; k < 4; k++){

					for(int i = 0; i < params.thetaAxis/2; i+=2){
						pair<int, int> thQ(yStart + i, xStart + i);
						pair<int, int> thP(yStart + i, xStart + i);
						movesP[j] = thP;
						movesQ[j] = thQ;
						j++;
					}

					for(int i = params.thetaAxis/2-1; i > 0 ; i-=2){
						pair<int, int> thQ(yStart + i, xStart + params.thetaAxis/2 - i);
						pair<int, int> thP(yStart + i, xStart + params.thetaAxis/2 - i);
						movesP[j] = thP;
						movesQ[j] = thQ;
						j++;
					}
					xStart += params.thetaAxis/2;
					if(xStart >= params.thetaAxis){
						xStart = 0;
						yStart += params.thetaAxis/2;
					}
				}
			}else{ /* or in the shape of an X |X
											  --- */				
				for(int i = 0; i < movesPSize/2; i++){
					pair<int, int> thQ(i,i);
					pair<int, int> thP(i,i);
					movesP[j] = thP;
					movesQ[j] = thQ;
					j++;
				}

				for(int i = movesPSize/2; i >=0; i--){
					pair<int, int> thQ(movesQSize/2-i,i);
					pair<int, int> thP(movesPSize/2-i,i);
					movesP[j] = thP;
					movesQ[j] = thQ;
					j++;				
				}
			}

			// build the P and Q subsets
			vector<PQ> pqv;
			vector<vector< P> > fP;
			vector<vector< Q> > fQ;
			
			int mIdx = 0;
			for(int i =0 ;i < numSubsets; i++){
				vector<P> pV;
				vector<Q> qV;
				for(int j = 0; j < subsetSize; j++){
					P p;
					p.s = (params.sMin) + movesP[mIdx].first*params.sStep;
					p.lambda = params.lambdaMin + movesP[mIdx].second*params.lambdaStep;
					pV.push_back(p);

					Q q;
					q.theta = params.thetaMin + movesQ[mIdx].first*params.thetaStep;
					q.h = params.hMin + movesQ[mIdx].second*params.hStep;
					qV.push_back(q);
					mIdx++;
				}
				fP.push_back(pV);
				fQ.push_back(qV);
				
			}

			for(int i =0; i < fP.size(); i++){
				Vector<Q> tmpQ = fQ[i];
				for(int j = 0; j < fP.size(); j++){
					Vector<P> tmpP = fP[j];					
					vector<PQ> tmpPQ;
					for(int k = 0; k < tmpQ.size(); k++){
						PQ pq;
						pq.p = tmpP[k];
						pq.q = tmpQ[k];
						tmpPQ.push_back(pq);
					}
					pqVec.push_back(tmpPQ);
				}
			}
			return pqVec;
		} // for the next levels in the hierarchy, just subdivide each subset to its children
		else{
			vector<vector<PQ> > tmpPQ;
			for(int k = 0; k < pqVec.size(); k++){
				vector<PQ> cPQ = pqVec[k];
				int halfSize = cPQ.size()/2;
				vector<PQ> tPQ;
				for(int i = 0; i < halfSize; i++){
					PQ pq;
					pq.p = cPQ[i].p;
					pq.q = cPQ[i].q;
					tPQ.push_back(pq);
				}
				tmpPQ.push_back(tPQ);
				tPQ.clear();

				for(int i = 0; i < halfSize; i++){
					PQ pq;
					pq.p = cPQ[i].p;
					pq.q = cPQ[i+halfSize].q;
					tPQ.push_back(pq);
				}
				tmpPQ.push_back(tPQ);
				tPQ.clear();

				for(int i = 0; i < halfSize; i++){
					PQ pq;
					pq.p = cPQ[i+halfSize].p;
					pq.q = cPQ[i].q;
					tPQ.push_back(pq);
				}
				tmpPQ.push_back(tPQ);
				tPQ.clear();

				for(int i = 0; i < halfSize; i++){
					PQ pq;
					pq.p = cPQ[i+halfSize].p;
					pq.q = cPQ[i+halfSize].q;
					tPQ.push_back(pq);
				}
				tmpPQ.push_back(tPQ);
				tPQ.clear();
			}

			return tmpPQ;
		}

	}
	/*merges two images together in a larger one*/
	static Mat mergeTwoImages(Mat I1, Mat I2){
		Mat merged = Mat::zeros(I1.rows,I1.cols+I2.cols,I1.type()); // merged image for display

		Rect I1ROI(Point(0,0),Point(I1.cols, I1.rows));
		Mat I1RoiM = merged(I1ROI);
		addWeighted(I1RoiM,0,I1,1,0.0,I1RoiM);
	
		Rect I2ROI(Point(I1.cols,0),Point(I1.cols+I2.cols, I2.rows));
		Mat I2RoiM = merged(I2ROI);
		addWeighted(I2RoiM,0,I2,1,0.0,I2RoiM);
		return merged;
	}
	/*given a patch, this function subdivides it into 4 equal rectangluar patches*/
	static void splitPatch(RegionSplit s, RegionSplit* regions, bool useParent){
		Mat r = s.img;
		Mat ul(r,Rect(0,0,r.cols/2,r.rows/2));
		Mat ur(r,Rect(r.cols/2,0,r.cols/2,r.rows/2));
		Mat dl(r,Rect(0,r.rows/2,r.cols/2,r.rows/2));
		Mat dr(r,Rect(r.cols/2,r.rows/2,r.cols/2,r.rows/2));
	
		RegionSplit s1;
		s1.anchor = s.anchor;
		s1.img = ul;
		s1.match = !useParent;

		RegionSplit s2;
		s2.anchor = s.anchor;
		s2.anchor.x += r.cols/2;	
		s2.img = ur;
		s2.match = !useParent;

		RegionSplit s3;
		s3.anchor = s.anchor;
		s3.anchor.y += r.rows/2;
		s3.img = dl;
		s3.match = !useParent;

		RegionSplit s4;
		s4.anchor = s.anchor;
		s4.anchor.x += r.cols/2;
		s4.anchor.y += r.rows/2;
		s4.img = dr;
		s4.match = !useParent;

		regions[0] = s1;
		regions[1] = s2;
		regions[2] = s3;
		regions[3] = s4;

		return;

	}
};

#endif