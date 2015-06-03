#include "timer.h"

// Timer class for timing purposes
void timer::tic()	
{
#ifdef __linux__
	t = getCurrentTime();
#elif _WIN32
	QueryPerformanceCounter(&t);
#endif

}

double timer::toc()
{

#ifdef __linux__
	LARGE_INTEGER t2 = getCurrentTime();
	double elapsedTime = (t2 - t)/1000.0; // elapsed time in ms
	return elapsedTime;
#elif _WIN32
	LARGE_INTEGER frequency;
	LARGE_INTEGER t1;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&t1);
	double	elapsedTime = (t1.QuadPart - t.QuadPart) * 1000.0 / frequency.QuadPart;
	return elapsedTime;
#endif
}


