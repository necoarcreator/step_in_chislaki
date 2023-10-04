#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>

double f(double x);
double computate_right(double x[], int i); //first right 
double computate_left(double x[], int i); //first left 
double computate_centre(double x[], int i); //first central
double computate_right_precisely(double x[]); //first right boundary with precision O(h^2)
double computate_left_precisely(double x[], int N); //first left boundary with precision O(h^2)
double computate_centre_second(double x[], int i); //second central 
double computate_left_second(double x[], int N); //second left boundary with precision O(h^2)
double computate_right_second(double x[]); //second right boundary with precision O(h^2)
int main()
{   
    int N = 100;
    double a = -5.0, b = 5.0;
    double x[N + 1];

    // first derivative with presicion O(h) //

    std::ofstream outf("first_o(h).txt");

    for (int i = 0; i < N + 1; i++)
    {
        x[i] = (double(i) * (b - a) / double(N)) + a;
    }

    for (int i = 0; i < N + 1; i++)
    {
        if (i != N)
        {
            std::cout << computate_right(x, i) << std::endl;
            outf << x[i] << "/" << computate_right(x, i) << std::endl;
        }
        if (i == N)
        {
            std::cout << "RIGHT BOUND IS " << computate_left(x, N) << std::endl;
            outf << x[i] << "/" << computate_left(x, i) << std::endl;
        }
    }

    // now - with precision O(h^2) //

    //std::ofstream outf("first_o(h2).txt");

    for (int i = 0; i < N + 1; i++)
    {
        if (i == 0)
        {
            std::cout << "LEFT BOUND IS " << computate_left_precisely(x, N) << std::endl;
            //outf << x[i] << "/" << computate_left_precisely(x, N) << std::endl;
        }
        if ((i != N) && (i != 0))
        {
            std::cout << computate_centre(x, i) << std::endl;
            //outf << x[i] << "/" << computate_centre(x, i) << std::endl;
        }
        if (i == N)
        {
            std::cout << "RIGHT BOUND IS " << computate_right_precisely(x) << std::endl;
            //outf << x[i] << "/" << computate_right_precisely(x) << std::endl;
        }
    }
    // now - the second derivative //

    //std::ofstream outf("second_o(h2).txt");

    for (int i = 0; i < N + 1; i++)
    {
        if (i == 0)
        {
            std::cout << "LEFT BOUND IS " << computate_right_second(x) << std::endl;
            //outf << x[i] << "/" << computate_right_second(x) << std::endl;
        }
        if ((i != N) && (i != 0))
        {
            std::cout << computate_centre_second(x, i) << std::endl;
            //outf << x[i] << "/" << computate_centre_second(x, i) << std::endl;
        }
        if (i == N)
        {
            std::cout << "RIGHT BOUND IS " << computate_left_second(x, N) << std::endl;
            //outf << x[i] << "/" << computate_left_second(x, N) << std::endl;
        }
    }

}

double f(double x)
{
    return atan(log(x * x + 1) + 1) * atan(log(x * x + 1) + 1);
}

double computate_right(double x[], int i)
{
    return (f(x[i+1]) - f(x[i])) / (x[i + 1] - x[i]);
}
double computate_left(double x[], int i)
{
    return (f(x[i]) - f(x[i - 1])) / (x[i] - x[i - 1]);
}
double computate_centre(double x[], int i)
{
    return (f(x[i + 1]) - f(x[i - 1])) / (x[i + 1] - x[i - 1]);
}
double computate_right_precisely(double x[])
{
    return (-3 * f(x[0]) - f(x[2]) + 4 * f(x[1])) / (x[2] - x[0]);
}

double computate_left_precisely(double x[], int N)
{
    return (3 * f(x[N]) + f(x[N - 2]) - 4 * f(x[N - 1])) / (x[N] - x[N - 2]);
}
double computate_centre_second(double x[], int i)
{
    return (f(x[i + 1]) - 2 * f(x[i]) + f(x[i - 1])) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]));
}

double computate_left_second(double x[], int N)
{
    return (-1 * f(x[N - 3]) + 4 * f(x[N - 2]) - 5 * f(x[N - 1]) + 2 * f(x[N])) / ((x[N] - x[N - 1]) * (x[N] - x[N - 1]));
}
double computate_right_second(double x[])
{
    return (2 * f(x[0]) - 5 * f(x[1]) + 4 * f(x[2]) - f(x[3])) / ((x[1] - x[0]) * (x[1] - x[0]));
}