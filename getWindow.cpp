#include "signal.h"
#include "internalFunc.h"

using namespace std;

vector<string> windowNames = {"hamm", "hann"};

vector<double> Signal::Frequency::get_window(string window, size_t len)
{
    bool flag = false;
    for(auto& i:windowNames)
    {
        if(window == i)
        {
            flag = true;
            break;
        }
    }

    if(!flag){throw runtime_error ("Error: Window name is wrong.");}

    vector<double> ret;

    if(window == "hann")
    {
        ret = Internal::hanning(len);
        return ret;
    }
    if(window == "hamm")
    {
        ret = Internal::hamming(len);
        return ret;
    }
}
