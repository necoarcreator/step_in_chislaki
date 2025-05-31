#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

double trueFunc(double x)
{
	return x * atan(x) - exp(x);
}

double trueDeriv(double x)
{
	return atan(x) + x / (1 + pow(x, 2)) - exp(x);
}
double f(double x, double u, double v)
{
	return ((-2) * x * v + 2 * u + 2 + (1 - 2 * x - pow(x, 2)) * exp(x)) / (1 + pow(x, 2));
}

void euler(double a, double b, vector<vector<double>>& res, double h, bool is_error_mode = false)
{
	vector <double> x{ a };

	while (round(x.back() * 1000) / 1000 != b)
	{
		x.push_back(x.back() + h);
	}


	double N = x.size();

	for (int i = 1; i < N; i++)
	{
		//res[0][i - 1] === u[i - 1]
		//res[1][i - 1] === v[i - 1]
		res[0].push_back(res[0][i - 1] + h * res[1][i - 1]);
		res[1].push_back(res[1][i - 1] + h * f(x[i - 1], res[0][i - 1], res[1][i - 1]));
	}
	if (is_error_mode)
	{
		for (int k = 0; k < N; k++)
		{
			res[0][k] = abs(trueFunc(x[k]) - res[0][k]);
			res[1][k] = abs(trueDeriv(x[k]) - res[1][k]);
		}
	}
	
	return;
}

void RK4(double a, double b, vector<vector<double>>& res, double h, bool is_error_mode = false)
{
	vector <double> x{ a };

	while (round(x.back() * 1000) / 1000 != b)
	{
		x.push_back(x.back() + h);
	}


	double N = x.size();

	double ku0, ku1, ku2, ku3;
	double kv0, kv1, kv2, kv3;

	for (int i = 1; i < N; i++)
	{
		//res[0][i - 1] === u[i - 1]
		//res[1][i - 1] === v[i - 1]

		ku0 = res[1][i - 1]; 
		kv0 = f(x[i - 1], res[0][i - 1], res[1][i - 1]);

		ku1 = res[1][i - 1] + h * kv0 / 2;
		kv1 = f(x[i - 1] + h / 2, res[0][i - 1] + h * ku0 / 2, res[1][i - 1] + h * kv0 / 2);

		ku2 = res[1][i - 1] + h * kv1 / 2;
		kv2 = f(x[i - 1] + h / 2, res[0][i - 1] + h * ku1 / 2, res[1][i - 1] + h * kv1 / 2);

		ku3 = res[1][i - 1] + h * kv2;
		kv3 = f(x[i - 1] + h, res[0][i - 1] + h * ku2, res[1][i - 1] + h * kv2);


		res[0].push_back(res[0][i - 1] + h * (ku0 + 2 * ku1 + 2 * ku2 + ku3)/ 6);
		res[1].push_back(res[1][i - 1] + h * (kv0 + 2 * kv1 + 2 * kv2 + kv3) / 6);
	}
	if (is_error_mode)
	{
		for (int k = 0; k < N; k++)
		{
			res[0][k] = abs(trueFunc(x[k]) - res[0][k]);
			res[1][k] = abs(trueDeriv(x[k]) - res[1][k]);
		}
	}
	return;


 }

void Adams3(double a, double b, vector<vector<double>>& res, double h, bool is_error_mode = false)
{
	vector <double> x{ a };

	while (round(x.back() * 1000) / 1000 != b)
	{
		x.push_back(x.back() + h);
	}


	double N = x.size();
	double ku0, ku1, ku2, ku3;
	double kv0, kv1, kv2, kv3;

	for (int i = 1; i < 4; i++)
	{
		ku0 = res[1][i - 1];
		kv0 = f(x[i - 1], res[0][i - 1], res[1][i - 1]);

		ku1 = res[1][i - 1] + h * kv0 / 2;
		kv1 = f(x[i - 1] + h / 2, res[0][i - 1] + h * ku0 / 2, res[1][i - 1] + h * kv0 / 2);

		ku2 = res[1][i - 1] + h * kv1 / 2;
		kv2 = f(x[i - 1] + h / 2, res[0][i - 1] + h * ku1 / 2, res[1][i - 1] + h * kv1 / 2);

		ku3 = res[1][i - 1] + h * kv2;
		kv3 = f(x[i - 1] + h, res[0][i - 1] + h * ku2, res[1][i - 1] + h * kv2);


		res[0].push_back(res[0][i - 1] + h * (ku0 + 2 * ku1 + 2 * ku2 + ku3) / 6);
		res[1].push_back(res[1][i - 1] + h * (kv0 + 2 * kv1 + 2 * kv2 + kv3) / 6);
	}

	for (int k = 4; k < N; k++)
	{
		res[1].push_back(res[1][k - 1] + (h / 24) * (55 * f(x[k - 1], res[0][k - 1], res[1][k - 1])
			- 59 * f(x[k - 2], res[0][k - 2], res[1][k - 2]) 
			+ 37 * f(x[k - 3], res[0][k - 3], res[1][k - 3])
			- 9 * f(x[k - 4], res[0][k - 4], res[1][k - 4])));
		
		res[0].push_back(res[0][k - 1] + (h / 24) * (55 * res[1][k - 1]
			- 59 * res[1][k - 2] + 37 * res[1][k - 3] - 9 * res[1][k - 4]));	
	}
	if (is_error_mode)
	{
		for (int k = 0; k < N; k++)
		{
			res[0][k] = abs(trueFunc(x[k]) - res[0][k]);
			res[1][k] = abs(trueDeriv(x[k]) - res[1][k]);
		}
	}
	return;

}

int main()
{
	double a = 0.0, b = 1.0;

	double h = 0.05;

	double u0 = -1.0, u0_der = -1.0;

	vector <double> x{ a };

	while (round(x.back() * 1000)/1000 != b)
	{
		x.push_back(x.back() + h);
	}


	double N = x.size();

	vector<vector<double>> u_eul{ {u0}, {u0_der} };

	euler(a,b, u_eul, h);

	ofstream file_euler("euler.txt");
	if (file_euler.is_open())
	{
		for (int j = 0; j < N; j++)
		{
			file_euler << x[j] << " " << u_eul[0][j] << " " << u_eul[1][j] << endl;
		}
	}
	else
	{
		cout << "programm ended during Euler method with code 1";
	}
	file_euler.close();

	vector<vector<double>> u_rk4{ {u0}, {u0_der} };

	RK4(a, b, u_rk4, h);
	ofstream file_rk4("rk4.txt");
	if (file_rk4.is_open())
	{
		for (int j = 0; j < N; j++)
		{
			file_rk4 << x[j] << " " << u_rk4[0][j] << " " << u_rk4[1][j] << endl;
		}
	}
	else
	{
		cout << "programm ended during RK4 method with code 1";
	}
	file_rk4.close();

	vector<vector<double>> u_eul_err_small{ {u0}, {u0_der} };
	vector<vector<double>> u_eul_err_big{ {u0}, {u0_der} };

	euler(a, b, u_eul_err_big, h, true);
	euler(a, b, u_eul_err_small, h/2, true);

	ofstream file_eul_err("eul_err.txt");
	if (file_eul_err.is_open())
	{
		for (int j = 0; j < N; j++)
		{
			file_eul_err << log(x[j]) << " " << u_eul_err_big[0][j] << " "
				<< u_eul_err_small[0][j] << " " << u_eul_err_big[1][j] <<
				" " << u_eul_err_big[1][j] << endl;
		}
	}
	else
	{
		cout << "programm ended during Euler method error finding with code 1";
	}
	file_eul_err.close();

	vector<vector<double>> u_add{ {u0}, {u0_der} };

	Adams3(a, b, u_add, h);
	ofstream file_add("addams.txt");
	if (file_add.is_open())
	{
		for (int j = 0; j < N; j++)
		{
			file_add << x[j] << " " << u_add[0][j] << " " << u_add[1][j] << endl;
		}
	}
	else
	{
		cout << "programm ended during Addams method with code 1";
	}
	file_add.close();

};



