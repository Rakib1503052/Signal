#ifndef INTERNALFUNC_H
#define INTERNALFUNC_H

#include "signal.h"

//This header is for declaring internal functions.
//These functions should not be used outside of library.

namespace Internal
{
std::vector<double> hanning(size_t);
std::vector<double> hamming(size_t, double a0 = 0.54);
}

#endif // INTERNALFUNC_H
