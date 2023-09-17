
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

using namespace std;

double foo(double x)
{
    return sin(x/3);
}


double f(double x)
{
    return foo(x/3 + exp(foo(x) * foo(x)));
}


double c(int n, double knots)
{
    double sum = 0, niz = 1;
    int i, j;
    
    for (i = 0; i < n +  1; i++) 
    { 
        
    for (j = 0; j < n + 1; j++)
    {   
        if (j != i) {
        niz *= (knots[i] - knots[j])
        }
    }
    sum += f(knots[i])/(niz);
    niz = 1;
    }
}



double L(int n, double x, double pr, double knots)
{   
    int i;
    if (n == 0)
    {
        return c(0, knots);
    }
    else
    {   
        return c(n, knots) * pr + L(n-1, x, pr/(x - knots[n - 1]));
    }
}


int main()
{
    int n = 100, i, N = 1000, k;
    double knots[n];
    double a = 0.0, b = 10.0;
    
    for (i = 0; i < n+1; i++)
    {
        knots[i] = (double(i))/(b - a);
        //knots[i] = (a + b)/2 + (b - a) * cos(M_PI * (2*i + 1)/(2*n))/2
    }
    
    ofstream outf("Interp_New.txt");

//pr надо затащить внутрь L!!!//
    for (k = 0; k < N + 1; k++)
    {
    outf << x[k] << "/" << L(n, x[k], pr,  knots) << endl;
    }
    return 0;
}

