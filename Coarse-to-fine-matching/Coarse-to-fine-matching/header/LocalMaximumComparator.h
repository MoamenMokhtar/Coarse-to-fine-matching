#ifndef LOCAL_MAXIMUM_COMPARATOR_H_
#define LOCAL_MAXIMUM_COMPARATOR_H_
#include "LocalMaximum.h"
#include "TemplateMatcher.h"

struct LocalMaximumComparator{
	inline bool operator()(LocalMaximum l1, LocalMaximum l2){
#ifdef USE_SQDIFF_NORMED
		return l1.val < l2.val;
#else 
		return l1.val > l2.val;
#endif
	}
};
#endif