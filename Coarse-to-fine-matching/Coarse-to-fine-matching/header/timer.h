
#pragma once

#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <windows.h>
#ifdef __linux__
#include <sys/time.h>
#endif

/*! Timer class for timing purposes */
class timer
{
    LARGE_INTEGER t;	
public:
	/*! Start timer */
	void tic();

	/*! Stop timer and return elapsed time in milliseconds */
	double toc();

private:
#ifdef __linux__
	long getCurrentTime(){
		struct timeval start;
		long mtime, seconds, useconds;
		gettimeofday(&start, NULL);
		seconds  = start.tv_sec; // seconds since epoch
		useconds = start.tv_usec; // microSeconds since epoch
		mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
		return mtime;
	}
#endif
};

