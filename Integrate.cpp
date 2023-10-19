#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double f(double* x);

double leftSq(vector<double> x, double h, int n);
double rightSq(vector<double> x, double h, int n);
double trap(vector<double> x, double h, int n);
double simp(vector<double> x, double h, int n);

int main()
{
    int n = 100, N = n + 1;
    vector<double> x;
    double b = 1, a = - 1, sumR = 0, sumL = 0, sumT = 0, sumS = 0, h;

    for (int i = 0; i < N; i++)
    {
        x.push_back((double(i) * (b - a)) / double(n) + a);
    }
    h = x[2] - x[1];
    sumL = leftSq(x, h, N);
    sumR = rightSq(x, h, N);
    sumT = trap(x, h, N);
    sumS = simp(x, h, N);

    ofstream outf("integr.txt");

    outf << sumL <<  "/" << "It is left" << endl;
    outf << sumR <<  "/" "It is right" << endl;
    outf << sumT <<  "/" << "It is Trap" << endl;
    outf << sumS <<  "/" << "It is Simps" << endl;
}

double f(double* x)
{
    return 1 / (pow(*x, 3) + *x + 10);
}


double leftSq(vector<double> x, double h, int n)
{
    int i = 0;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += h * (f(&x[i]));
    }
    return sum;
}

double rightSq(vector<double> x, double h, int n)
{
    int i = 0;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += h * (f(&x[i + 1]));
    }
    return sum;
}
double trap(vector<double> x, double h, int n)
{
    int i = 0;
    double sum = 0;
    for (i = 0; i < n; i++)
    {
        sum += (h/2) * (f(&x[i + 1]) + f(&x[i]));
    }
    return sum;
}
double simp(vector<double> x, double h, int n)
{
    int i = 0;
    double sum = 0;
    double sred = (x[n - 1] + x[0]) / 2;
    for (i = 0; i < n; i++)
    {
        sum += (h / 6) * (f(&x[i + 1]) + 4 * f(&sred) + f(&x[i]));
    }
    return sum;
}