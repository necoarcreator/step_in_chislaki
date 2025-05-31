#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double f(float* x);
double F_tr(float x);

double leftSq(vector<double> x, double h, int n);
double rightSq(vector<double> x, double h, int n);
double trap(vector<double> x, double h, int n);
double simp(vector<double> x, double h, int n);
double err(double a, double b, double h);

int main()
{
    int n = 100, N = n + 1;
    vector<double> x;
    double b = 1, a = - 1, sumR = 0, sumL = 0, sumT = 0, sumS = 0, h;

   /* for (int i = 0; i < N; i++)
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
    */

    ofstream outf("error.txt");

    double hh[10], mx[10];
    int J = size(hh);
    outf << "h   " << "error" << endl;
    for (int i = 0; i < J; i++)
    {
        hh[i] = pow(10, -i);
        mx[i] = err(a, b, hh[i]);
        outf << hh[i] << "/" << mx[i] << endl;
    }

}

double f(double* x)
{
    return 1 / (pow(*x, 3) + *x + 10);
}

double F_tr(double x)
{
    return log(x + 2) / 13 - log(pow(x,2) - 2 * x + 5) / 26 + 3 * atan((x - 1) / 2) / 26;
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

/*double rightSq(vector<double> x, double h, int n)
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
}*/

double err(double a, double b, double h)
{
    vector<double> x;
    double J = (b - a) / h;
    int i;

        for (i = 0; i < J; i++)
        {
            x.push_back((double(i) * (b - a)) / double(J) + a);
        }

    double y, y_tr, err;
    y = leftSq(x, h, J);
    y_tr = F_tr(b) - F_tr(a);
    err = y - y_tr;
    return err;
}