#define _USE_MATH_DEFINES
#include "internalFunc.h"

using namespace std;

vector<double> Internal::hamming(size_t len, double a0)
{
    int N = len-1;
    double a1 = 1-a0;
    vector<double> ret(len,0);

    for(size_t n = 0; n < len; n++)
    {
        ret[n] = a0 - (a1 * cos(2 * M_PI * n / N));
    }

    return ret;
}

vector<double> Internal::hanning(size_t len)
{
    return Internal::hamming(len, 0.5);
}
