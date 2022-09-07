#include "signal.h"

using namespace std;

vector<double> Signal::DSF::padding(vector<double> data, size_t pad_len)
{
    vector<double> leftPad(data.begin()+1, data.begin()+pad_len);
    reverse(leftPad.begin(), leftPad.end());
    vector<double> rightPad(data.rbegin()+1, data.rbegin()+pad_len);
    double leftEnd = data[0], rightEnd = data[data.size()-1];

    for(size_t i = 0; i < pad_len-1; i++)
    {
        leftPad[i] = 2 * leftEnd - leftPad[i];
        rightPad[i] = 2 * rightEnd - rightPad[i];
    }

    data.insert(data.begin(), leftPad.begin(), leftPad.end());
    data.insert(data.end(), rightPad.begin(), rightPad.end());

    return data;
}

vector<double> Signal::DSF::padRemover(vector<double> data, size_t pad_len)
{
    data.erase(data.begin(), data.begin() + pad_len);
    data.erase(data.end() - pad_len, data.end());

    return data;
}
