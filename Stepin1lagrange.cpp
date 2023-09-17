
#include <iostream>
#include <math.h>
#include <vector>
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


double e(int i, double x, double knots, int n)
{
    double sum = 0, sum_z = 0;
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
    return (sum/sum_z);
}
    
 double L(int n, double x, double knots)   
 {
     double sum = 0;
     int i;
     
     for(i = 0; i < n + 1; i++)
     {
         sum += f(knots[i]) * e(i, x, knots, n);
         
     }
     
     return sum;
     
 }
    

main(void)
{
    
    int n = 100, i, N = 1000, k;
    double a = 0.0, b = 10.0;
    double knots[n + 1], x[N + 1];
    //std::vector<double> x;
    
    for (i = 0; i < n + 1; i++)
    {
        knots[i] = (double(i)/(b - a));
        //knots[i] = (a + b)/2 + (b - a) * cos(M_PI * (2*i + 1)/(2*n))/2
    }
    
    ofstream outf("Interp_Lag.txt");
    
    for (k = 0; k < N + 1; k++)
    {
    outf << x[k] << "/" << L(n, x[k], knots) << endl;
    }
    
    return 0;
}
