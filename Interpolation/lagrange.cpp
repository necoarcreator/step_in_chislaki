#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;


double f(double x)
{
    return sin(x / 3 + exp(sin(x / 3) * sin(x / 3)));
}


double e(int i, double x, double knots[], int n)
{
    double sum = 1, sum_z = 1;
    int j;

    for (j = 0; j < n + 1; j++)
    {
        if (j == i)
        {
            continue;
        }
        else
        {
            sum *= (x - knots[j]);
            sum_z *= (knots[i] - knots[j]);
        }
    }
    return (sum / sum_z);
}

double L(int n, double x, double knots[])
{
    double sum = 0;
    int i;

    for (i = 0; i < n + 1; i++)
    {
        sum += f(knots[i]) * e(i, x, knots, n);

    }

    return sum;

}


int main(void)
{

    int n = 50, i, j, N = 100, k;
    double a = 0.0, b = 10.0;
    double knots[5], x[101];
    //std::vector<double> x;

    for (i = 0; i < n + 1; i++)
    {
        knots[i] = (double(i) * (b - a) / double(n));

        //knots[i] = (a + b) / 2 + (b - a) * cos(M_PI * (2 * i + 1) / (2 * n)) / 2;
    }
    for (j = 0; j < N + 1; j++)
    {
        x[j] = (double(j) * (b - a) / double(N));
    }

    ofstream outf("Interp_Lg.txt");

    //ofstream outf("Interp_Lag_Cheb.txt");

    for (k = 0; k < N + 1; k++)
    {
        outf << x[k] << "/" << L(n, x[k], knots) << endl;
    }

    return 0;
}
