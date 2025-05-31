#include<stdio.h>
#include<fstream>
#include<math.h>
#include<omp.h>
#include <vector>
#include<iostream>


using namespace std;


double u0(double x, double t) { //analit
	return  3 / 2 * pow(x, 2) + exp(-t) * acos(x / 2);
		//2 * sin(x + t);
}

double f(double x, double t) { //neodnor
	return -3 + x * exp(-t) / pow((4 - pow(x, 2)), 3 / 2) + 2 * exp(-t) * acos(x / 2);
		//-2 * sin(x + t);
}

double phi(double x) { //nach smesh
	return 3 / 2 * pow(x, 2) + acos(x / 2);
		//2 * sin(x);
}

double psi(double x) { //nach skor
	return -acos(x / 2);
		//2 * cos(x);
}

double gamma0(double t) { //neodn v gr
	return (atan(1) * 4 - 1/2) * exp(-t);
		//2 * (sin(t) - cos(t));
}

double gammal(double t) { //neodn v gr sprava
	return 3 / 2 + atan(1) * 4 / 3 * exp(-t);
		//2 * sin(1 + t);
}


int main() {

	// initializing


	double T = 1.0;
	double x = 1.0;
	double h = 0.1;
	double tau = 0.0005;
	double a_sq = 0.5; //skorost volni
	int N = T / tau +1; 
	int I = x / h +1; 

	double alpha0 = 2.0;
		//1.0;
	double beta0 = 1.0;
		//-1.0;
	double alphal = 1.0;
		//1.0;
	double betal = 0.0;
		//0.0;

	vector<vector<double>> u;
	vector<double> u_temp;

	// filling zeros
	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < I; i++) {
			u_temp.push_back(0);
		}
		u.push_back(u_temp);
	}
	// initial condition

	for (int i = 0; i < I; i++) {
		u[0][i] = phi(i * h);
		u[1][i] = u[0][i] + tau * psi(i * h);
	}

	//single precision sol

	for (int n = 2; n < N; n++) {

		u[n][0] = h / (alpha0 * h - beta0) * gamma0(n * tau) - beta0 / (alpha0 * h - beta0) * u[n][1];

		for (int i = 1; i < I; i++) {
			u[n][i] = 2 * u[n-1][i] - u[n-2][i] + a_sq * a_sq * tau * tau / h / h * (u[n-1][i+1] - 2 * u[n-1][i] + u[n - 1][i-1]) + tau * tau * f(i * h, (n - 1) * tau);
		}

		u[n][I - 1] = h / (alphal * h - betal) * gammal(n * tau) - betal / (alphal * h - betal) * u[n][I - 2];

		for (int i = 0; i < I; i++) {
			u[n-2][i] = u[n-1][i];
			u[n-1][i] = u[n][i];
		}


	}

	// output results

	ofstream out("output_v1.txt");
	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < I; i++) {
			out << i * h << " " << n * tau << " " << u[n][i] << " " << u0(i * h, n * tau) << " " << log(abs(u[n][i] - u0(i * h, n * tau))) << "\n";
		}
		out << "-\n";
	}
	out.close();
	cout << "Programm ended successefully!\n";

	for (int i = 0; i < I; i++) {
		u[0][i] = phi(i * h);
		u[1][i] = u[0][i] + tau * psi(i * h);
	}

	//double precision sol

	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < I; i++)
		{
			u[n][i] = 0;
		}
	}

	// initial condition

	for (int i = 0; i < I; i++) {
		u[0][i] = phi(i * h);
		u[1][i] = u[0][i] + tau * psi(i * h);
	}

	for (int n = 2; n < N; n++) {

		u[n][0] = h * gamma0(n * tau) / (alpha0 * h - 3 * beta0 / 2) + beta0 / 2 * (4 * u[n][1] - u[n][2]);

		for (int i = 1; i < I; i++) {
			u[n][i] = 2 * u[n - 1][i] - u[n - 2][i] + a_sq * a_sq * tau * tau / h / h * (u[n - 1][i + 1] - 2 * u[n - 1][i] + u[n - 1][i - 1]) + tau * tau * f(i * h, (n - 1) * tau);
		}

		u[n][I - 1] = h * gammal(n * tau) / (alphal * h - 3 * betal / 2) + betal / 2 * (4 * u[n][I - 2] - u[n][I - 3]);

		for (int i = 0; i < I; i++) {
			u[n - 2][i] = u[n - 1][i];
			u[n - 1][i] = u[n][i];
		}


	}

	// output results

	ofstream out2("output_v2.txt");
	for (int n = 0; n < N; n++)
	{
		for (int i = 0; i < I; i++) {
			out2 << i * h << " " << n * tau << " " << u[n][i] << " " << u0(i * h, n * tau) << " " << log(abs(u[n][i] - u0(i * h, n * tau))) << "\n";
		}
		out2 << "-\n";
	}
	out2.close();
	cout << "Programm ended successefully!\n";


	return 0;
}
