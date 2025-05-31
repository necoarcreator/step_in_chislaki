#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;

double f(double t, double x)
{
	return (x + 2 * pow(t, 2) * tanh(x * t)) / (pow(cosh(x * t), 2));
		//-2 - (1 + tan(x - t) + 2 * pow(tan(x - t), 2)) / (2 * cos(x - t));
		//(x + 2 * pow(t, 2) * tanh(x*t)) / (pow(cosh(x*t),2));
}

double u0(double t, double x)
{
	return x + tanh(x * t);
		//x * x + 1 / (2 * cos(x - t));
		//x + tanh(x * t);
}

double phi(double x)
{
	return x;
		//x * x + 1 / (2 * cos(x));
		//x;
}

double gamma0(double t)
{
	return 1 + t;
		//sqrt(2 + t);
		//1 + t;
}
double gammal(double t)
{
	return 1 + t / pow(cosh(t), 2);
		//(t - 2) / sqrt(1 + t);
		//1 + t / pow(cosh(t),2);
}

vector <double> progon(double a, double b, vector<double> A, vector<double> B, vector<double> C, vector<double> d, int I) {

	vector<double> u(I);
	vector<double> a_t;
	vector<double> b_t;

	a_t.push_back(-C[0] / B[0]);
	b_t.push_back(d[0] / B[0]);

	for (int i = 1; i < I; i++) {
		a_t.push_back(-C[i] / (A[i] * a_t[i - 1] + B[i]));
		b_t.push_back((d[i] - A[i] * b_t[i - 1]) / (A[i] * a_t[i - 1] + B[i]));
	}

	u[I-1] = ((d[I - 1] - A[I - 1] * b_t[I - 2]) / (B[I - 1] + A[I - 1] * a_t[I - 2]));



	for (int i = I - 2; i >= 0; i--) {
		u[i] = a_t[i] * u[i + 1] + b_t[i];
	}

	return u;
}

int main()
{
	double a = 0.0, b = 1.0;
	double a_sq = 1.0;
	double h = 0.05, tau = 0.01, T = 1.0;
	double I{ (b - a) / h + 1};
	double N{ T / tau + 1};
	I = int(I);
	N = int(N);

	double alpha0 = 1.0;
		//3.0;
		//1.0;
	double beta0 = 1.0;
		//-1.0;
		//1.0;
	double alphal = 0.0;
		//1.0;
		//0.0;
	double betal = 1.0;
		//0.0;
		//1.0;

	vector<vector<double>> u;
	vector<double> u_temp;
	
	for (int n = 0; n < N + 1; n++)
	{
		for (int i = 0; i < I + 1; i++)
		{
			u_temp.push_back(0);
		}
		u.push_back(u_temp);
		u_temp.clear();
	}

	for (int i = 0; i < I + 1; i++) {
		u[0][i] = phi(i * h);
	}

	double syg = 0.0;
	double A = a_sq * a_sq * (1. - syg) * tau / h / h;
	double B = -2* a_sq * a_sq * (1. - syg) * tau /h /h - 1.;
	double C = a_sq * a_sq * (1. - syg) * tau / h / h;

	double A0 = alpha0 - beta0 / h;
	double B0 = beta0 / h;
	double BI = betal / h;
	double CI = alphal - betal / h;
	

	vector<double> d(I);
	vector<double> d_temp;

	vector<double> a_t(I);
	vector<double> b_t(I);
	vector<double> c_t(I);

	//single presition
	a_t[0] = 0.;
	b_t[0] = alpha0 * h - beta0;
	c_t[0] = beta0;
	
	for (int i = 1; i < I - 1; i++)
	{
		a_t[i] = A;
		b_t[i] = B;
		c_t[i] = C;
	}

	a_t[I - 1] = -betal;
	b_t[I - 1] = alphal * h + betal;
	c_t[I - 1] = 0.;

	bool isSuccess = true;
	try {
		if (a != 0.0)
		{
			throw "The program won't work correctly with x_left != 0.";
		}
		for (int n = 1; n < N; n++) {
			d[0] = gamma0(tau * n) * h;
			for (int i = 1; i < I - 1; i++) {
				d[i] = -f(tau * (n - 1), h * i) * tau - u[n - 1][i] - a_sq * a_sq * syg * tau / h / h * (u[n - 1][i + 1] - 2. * u[n - 1][i] + u[n - 1][i - 1]);
			}
			d[I - 1] = gammal(tau * n) * h;
			u[n] = progon(a, b, a_t, b_t, c_t, d, I);
		}
			
		
		
	}
	catch (const string exeption)
	{	
		isSuccess = false;
		cerr << "Error:" << exeption << endl;
	}
	

	if (isSuccess)
	{
		ofstream out("output_tepl.txt");
		if (out.is_open())
		{
			for (int n = 0; n < N + 1; n++)
			{
				for (int i = 0; i < I; i++) {
					out << i * h << " " << n * tau << " " << u[n][i] << " " << u0(n * tau, i * h) << " " << log(abs(u[n][i] - u0(n * tau, i * h))) << "\n";
				}
				out << "-\n";

			}
			cout << "Program ended successefully!\n";
		}
		else
		{
			cout << "Error while reading file!\n";
		}
		out.close();
	}

	//second presition
	a_t[0] = 0.;
	b_t[0] = 2. * alpha0 * h - 2. * beta0 - beta0 * h * h / tau / a_sq / a_sq;
	c_t[0] = 2. * beta0;

	for (int i = 1; i < I - 1; i++)
	{
		a_t[i] = A;
		b_t[i] = B;
		c_t[i] = C;
	}

	a_t[I - 1] = -2. * betal;
	b_t[I - 1] = 2. * alphal * h + 2. * betal + betal * h * h / tau / a_sq / a_sq;
	c_t[I - 1] = 0.;

	bool isSuccess2 = true;
	try {
		if (a != 0.0)
		{
			throw "The program won't work correctly with x_left != 0.";
		}

		for (int n = 1; n < N; n++) {
			d[0] = gamma0(tau * n) * h * 2. - beta0 * f(tau * n, 0) * h * h / a_sq / a_sq - u[n - 1][0] * h * h * beta0 / tau / a_sq / a_sq;
			for (int i = 1; i < I - 1; i++) {
				d[i] = -f(tau * (n - 1), h * i) * tau - u[n - 1][i] - a_sq * a_sq * syg * tau / h / h * (u[n - 1][i + 1] - 2. * u[n - 1][i] + u[n - 1][i - 1]);
			}
			d[I - 1] = 2. * gammal(tau * n) * h + h * h * betal * u[n - 1][I - 1] / a_sq / a_sq / tau + betal * h * h * f(tau * n, (I - 1) * h) / a_sq / a_sq;
			u[n] = progon(a, b, a_t, b_t, c_t, d, I);
		}

	}
	catch (const string exeption)
	{
		isSuccess2 = false;
		cerr << "Error:" << exeption << endl;
	}


	if (isSuccess2)
	{
		ofstream out("output_tepl2.txt");
		if (out.is_open())
		{
			for (int n = 0; n < N + 1; n++)
			{
				for (int i = 0; i < I; i++) {
					out << i * h << " " << n * tau << " " << u[n][i] << " " << u0(n * tau, i * h) << " " << log(abs(u[n][i] - u0(n * tau, i * h))) << "\n";
				}
				out << "-\n";

			}
			cout << "Program ended successefully!\n";
		}
		else
		{
			cout << "Error while reading file!\n";
		}
		out.close();
	}
	return 0;
}