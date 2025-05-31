#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double P(double x) {  // коэфы перед u'
    return cos(x);
}

double Q(double x) {   // коэфы перед u
    return sin(x);
}

double F(double x) {
    return (1 - sin(x) - cos(x));
}

double U0(double x) {  // точное решение
    return cos(x) + sin(x);
}

vector <double> progon(vector<double> X, vector<double> A, vector<double> B, vector<double> C, vector<double> d, int N) {

    vector<double> u(N, 0);
    vector<double> a;
    vector<double> b;

    a.push_back(-C[0] / B[0]);
    b.push_back(d[0] / B[0]);

    for (int i = 1; i <= N; i++) {
        a.push_back(-C[i] / (A[i] * a[i - 1] + B[i]));
        b.push_back((d[i] - A[i] * b[i - 1]) / (A[i] * a[i - 1] + B[i]));
    }

    u.push_back((d[N] - A[N] * b[N - 1]) / (B[N] + A[N] * a[N - 1]));

    for (int i = N - 1; i >= 0; i--) {
        u[i] = a[i] * u[i + 1] + b[i];
    }

    return u;
}

vector<double> sol1(vector<double> X, double gammaa, double gammab, double nua, double nub, double etaa, double etab, int N, double h) {
    vector<double> p;
    for (int i = 0; i <= N; i++) {
        p.push_back(P(X[i]));
    }

    vector<double> q;
    for (int i = 0; i <= N; i++) {
        q.push_back(Q(X[i]));
    }

    vector<double> f;
    f.push_back(gammaa);
    for (int i = 1; i < N; i++) {
        f.push_back(F(X[i]));
    }
    f.push_back(gammab);

    vector<double> A;
    vector<double> B;
    vector<double> C;

    A.push_back(0);
    B.push_back(nua - etaa / h);
    C.push_back(etaa / h);

    for (int i = 0; i < N - 1; i++) {
        A.push_back(1 / h / h - p[i] / 2 / h);
        B.push_back(q[i] - 2 / h / h);
        C.push_back(1 / h / h + p[i] / 2 / h);
    }

    A.push_back(-etab / h);
    B.push_back(nub + etab / h);
    C.push_back(0);

    vector<double> sol = progon(X, A, B, C, f, N);
    return sol;
}

vector<double> sol2(vector<double> X, double gammaa, double gammab, double nua, double nub, double etaa, double etab, int N, double h) {
    vector<double> p;
    for (int i = 0; i <= N; i++) {
        p.push_back(P(X[i]));
    }

    vector<double> q;
    for (int i = 0; i <= N; i++) {
        q.push_back(Q(X[i]));
    }

    vector<double> f;
    f.push_back(gammaa + F(X[0]) * etaa * h / (2. - p[0] * h));
    for (int i = 1; i < N; i++) {
        f.push_back(F(X[i]));
    }
    f.push_back(gammab - etab * F(X[N]) * h / (2. + p[N] * h));

    vector<double> A;
    vector<double> B;
    vector<double> C;

    A.push_back(0);
    B.push_back(nua - etaa / 2. / h * ((4. - 2 * q[0] * h * h) / (2. - p[0] * h)));
    C.push_back(etaa / 2. / h * (1. + (2. + p[0] * h) / (2. - p[0] * h)));

    for (int i = 0; i < N - 1; i++) {
        A.push_back(1 / h / h - p[i] / 2. / h);
        B.push_back(q[i] - 2. / h / h);
        C.push_back(1 / h / h + p[i] / 2. / h);
    }

    A.push_back(-etab / 2. / h * (1. - (p[N] * h - 2.) / (p[N] * h + 2.)));
    B.push_back(nub + etab / 2. / h * ((4. - 2. * q[N] * h * h) / (2. + p[N] * h)));
    C.push_back(0);

    vector<double> sol = progon(X, A, B, C, f, N);
    return sol;
}

int main()
{
    double start = 0.; //левая граница
    double fin = 1.; //правая граница
    double h = 0.05; // шаг тут куча символов, чтоб не пропустить efdhfpaerughvh alsdjfhweriufhasdljfh aryfhq ifgsedrghaerklufhdafugsoeriawelrkfgarlkjgheraghwaerghsdflghsdlfuighsldfghspdruihypero
    int N = int((fin - start) / h) + 1;

    vector <double> X;
    for (int i = 0; i <= N; i++) {
        X.push_back(start + i * h);
    }

    const double nua = 1.;
    const double etaa = -1.;
    const double nub = 1.;
    const double etab = 0.;
    const double gammaa = 0;
    const double gammab = 1.3818;

    vector<double> U1 = sol1(X, gammaa, gammab, nua, nub, etaa, etab, N, h);
    vector<double> U2 = sol2(X, gammaa, gammab, nua, nub, etaa, etab, N, h);

    ofstream csv_file;
    const char csv_file_name1[64] = "Up1.csv";
    csv_file.open(csv_file_name1);
    csv_file << "X,U1,U2,Y\n";
    for (int i = 0; i < N; i++) {
        csv_file << X[i] << ", " << U1[i] << ", " << U2[i] << ", " << U0(X[i]) << endl;
    }
    csv_file.close();

    vector <long double> H; // массив шагов
    vector <long double> O1; // массив ошибок
    vector <long double> O2; // массив ошибок
    vector <long double> Hl; // массив шагов отлагорифмированных
    vector <long double> Ol1; // массив ошибок отлагорифмированнаых
    vector <long double> Ol2; // массив ошибок отлагорифмированнаых

    for (int i = 1; i <= 5; i++) { //сюда вписать до какого шага мы идём (10^-наше число)
        long double hl = pow(10, -i);
        H.push_back(hl);
        int N1 = int(double(fin - start) / pow(10, -i)) + 1;
        vector <double> X1;
        vector <double> Y1;
        for (int i = 0; i <= N1; i++) {
            X1.push_back(start + i * hl);
            Y1.push_back(U0(X1[i]));
        }
        vector <double> U11 = sol1(X1, gammaa, gammab, nua, nub, etaa, etab, N1, hl);
        vector <double> U22 = sol2(X1, gammaa, gammab, nua, nub, etaa, etab, N1, hl);
        long double maxDiff1 = 0;
        long double maxDiff2 = 0;
        long double diff1 = 0;
        long double diff2 = 0;
        for (int i = 0; i <= N1; i++) {
            diff1 = abs(Y1[i] - U11[i]);
            diff2 = abs(Y1[i] - U22[i]);
            if (diff1 > maxDiff1) maxDiff1 = diff1;
            if (diff2 > maxDiff2) maxDiff2 = diff2;
        }
        O1.push_back(maxDiff1);
        O2.push_back(maxDiff2);
        Hl.push_back(log(pow(10, -i)));
        Ol1.push_back(log(maxDiff1));
        Ol2.push_back(log(maxDiff2));
        cout << i << "-iy progon gotov" << endl;
    }
    const char csv_file_name2[64] = "Up2.csv";
    csv_file.open(csv_file_name2);
    csv_file << "H,O1,O2,Hl,Ol1,Ol2\n";
    for (int i = 0; i < H.size(); i++) {
        csv_file << H[i] << ", " << O1[i] << ", " << O2[i] << ", " << Hl[i] << ", " << Ol1[i] << ", " << Ol2[i] << endl;
        //cout << H[i] << ", " << O1[i] << ", " << O2[i] << ", " << Hl[i] << ", " << Ol1[i] << ", " << Ol2[i] << endl;
    }
    csv_file.close();

    cout << "gotovo\n";

    return 0;
}