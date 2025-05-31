#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

double f(double x)
{
    return sin(x / 3 + exp(sin(x / 3) * sin(x/3)));
}
/*
double c(int n, double knots[])
{
    double sum = 0, niz = 1;
    int i, j;

    for (i = 0; i < n; i++)
    {

        for (j = 0; j < n; j++)
        {
            if (j != i) {
                niz *= (knots[i] - knots[j]);
            }
        }
        sum += f(knots[i]) / (niz);
        niz = 1;
    }

    return sum;
}

double givePr(double x, double knots[], int n)
{
    int i;
    double pr = 1;
    for (i = 1; i < n; i++)
    {
        pr *= (x - knots[i - 1]);
    }

    return pr;
}
double Newt(int n, double x, double knots[])
{
    int i;
    double pr = givePr(x, knots, n);
    if (n == 0)
    {
        return c(0, knots);
    }
    else
    {
        return c(n, knots) * pr + Newt(n - 1, x, knots);
    }
}
*/





double Newt(int n, double x, double knots[])
{
    int i, j;
    double pr = 1, niz = 1, sum = f(knots[0]);

    for (i = 1; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j != i) {
                niz *= (knots[i] - knots[j]);
            }
        }
        pr *= (x - knots[i - 1]);
        sum += pr * f(knots[i]) / niz;
        niz = 1;
    }
    return sum;
}

/*
double newton(double x, vector<double> X, vector<double> Y, int n) {
    vector<vector<double>> f(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        f[i][0] = Y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            f[i][j] = (f[i + 1][j - 1] - f[i][j - 1]) / (X[i + j] - X[i]);
        }
    }
    double y = f[0][0];
    double w = 1.0;
    for (int j = 1; j < n; j++) {
        w *= (x - X[j - 1]);
        y += f[0][j] * w;
    }
    return y;
}
*/
int main()
{
    int n = 20, i, j, N = 100, k;
    double knots[5], x[100];
    double a = 0.0, b = 10.0;

    for (i = 0; i < n; i++)
    {
        knots[i] = (double(i) * (b - a) / double(n - 1));
        //knots[i] = (a + b) / 2 + (b - a) * cos(M_PI * (2 * i + 1) / (2 * n)) / 2;
        cout << knots[i] << endl;
    }
    for (j = 0; j < N; j++)
    {
        x[j] = (double(j) * (b - a) / double(N - 1));
    }

    ofstream outf("Interp_New.txt");
    //ofstream outf("Interp_New_Cheb.txt");

    for (k = 0; k < N; k++)
    {
        outf << x[k] << "/" << Newt(n, x[k], knots) << endl;
    }
    return 0;
}

