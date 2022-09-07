#include "signal.h"

using namespace std;

vector<double> Signal::DSF::lfilter(vector<double> B, vector<double> A, const vector<double> &X, vector<double> &Zi)
{
    //This function is magic. That's right, it's MAGIC. Don't touch it if you don't want to get cursed.
    //I was wrong. This IS black magic.

    if (A.empty())
    {
        throw domain_error("The feedback filter coefficients are empty.");
    }
    if (all_of(A.begin(), A.end(), [](double coef){ return coef == 0; }))
    {
        throw domain_error("At least one of the feedback filter coefficients has to be non-zero.");
    }
    if (A[0] == 0)
    {
        throw domain_error("First feedback coefficient has to be non-zero.");
    }

    //Output vector
    vector<double> Y(X.size(), 0);

    // Normalize feedback coefficients if a[0] != 1;
    auto a0 = A[0];
    if (a0 != 1.0)
    {
        transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
        transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
    }

    size_t input_size = X.size();
    size_t filter_order = max(A.size(), B.size());
    B.resize(filter_order, 0);
    A.resize(filter_order, 0);
    Zi.resize(filter_order, 0);
    Y.resize(input_size);

    const double *x = &X[0];
    const double *b = &B[0];
    const double *a = &A[0];
    double *z = &Zi[0];
    double *y = &Y[0];

    for (size_t i = 0; i < input_size; ++i)
    {
        size_t order = filter_order - 1;
        while (order)
        {
            if (i >= order)
            {
                z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
            }
            --order;
        }
        y[i] = b[0] * x[i] + z[0];
    }
    Zi.resize(filter_order - 1);

    return Y;
}


vector<double> Signal::DSF::filtfilt(vector<double> B, vector<double> A, const vector<double> &X)
{
    //Hopefully you've watched Doctor Strange: Multiverse of Madness. If you haven't then watch it ASAP.
    //And then you'll know why you shouldn't mess with black magic.
    //Because this is BLACK MAGIC. Right, lfilter was magic and this filtfilt is black magic.

    //I was wrong again. This is not black magic, this is outright voodoo. Got it from some dude on SO.

    int len = X.size();     // length of input
    int na = A.size();
    int nb = B.size();
    int nfilt = max(nb, na);
    int nfact = 3 * (nfilt); // length of edge transients

    if (len <= nfact)
    {
        throw domain_error("Input data too short! Data must have length more than 3 times filter order.");
    }

    // set up filter's initial conditions to remove DC offset problems at the
    // beginning and end of the sequence
    B.resize(nfilt, 0);
    A.resize(nfilt, 0);

    double y0;
    vector<double> signal1, signal2, zi;
    signal1 = Signal::DSF::padding(X, nfact);

    // Calculate initial conditions
    vector<double> zzi = lfilter_zi(B, A);
    zi.resize(zzi.size());

    // Do the forward and backward filtering
    y0 = signal1[0];
    transform(zzi.begin(), zzi.begin() + zzi.size(), zi.begin(), [y0](double val){ return val*y0; });
    signal2 = lfilter(B, A, signal1, zi);     //Forward filter

    reverse(signal2.begin(), signal2.end());
    y0 = signal2[0];
    transform(zzi.begin(), zzi.begin() + zzi.size(), zi.begin(), [y0](double val){ return val*y0; });
    signal1 = lfilter(B, A, signal2, zi);     //Reverse filter

    reverse(signal1.begin(), signal1.end());
    vector<double> Y = Signal::DSF::padRemover(signal1, nfact-1);        //Output

    return Y;
}

vector<double> Signal::DSF::lfilter_zi(vector<double> b, vector<double> a)
{
    //I don't know what this is. Scipy had it, I wrote it.

    size_t n = a.size() > b.size() ? a.size() : b.size();

    if(a.size() != b.size()) //Zero padding to match size
    {
        a.insert(a.end(), (n - a.size()), 0.0);
        b.insert(b.end(), (n - a.size()), 0.0);
    }

    if (a[0] != 1) //Normalize to a[0]=1
    {
        double temp_coef = a[0];
        for (size_t i = 0; i < n; i++)
        {
            b[i] = b[i] / temp_coef;
            a[i] = a[i] / temp_coef;
        }
    }

    vector<vector<double>> IminusA = matalg::mat_sub(matalg::mat_I(n-1),
                                                     matalg::mat_T(matalg::companion(a)));

    vector<double> B (n-1, 0);

    for(size_t i = 0; i < B.size(); i++) {B[i] = b[i+1] - (a[i+1] * b[0]);}

    vector<double> zi = matalg::solve_linSystem(IminusA, B);

    return zi;
}
