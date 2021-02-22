#include <iostream>
#include "Diff_equation.h"
using namespace std;

float Diff_equation::f(float x, float y) { return x-y; }

void Diff_equation :: solve_rungekutta(float x, float y, float h, int m) {
    float k1, k2, k3, k4;
    int i = 1;
    X[0] = x; Y[0] = y;
    
    do
    {
        k1 = h * f(X[i - 1], Y[i - 1]);
        k2 = h * f(X[i - 1] + h / 2, Y[i - 1] + k1 / 2);
        k3 = h * f(X[i - 1] + h / 2, Y[i - 1] + k2 / 2);
        k4 = h * f(X[i - 1] + h, Y[i - 1] + k3);
        Y[i] = Y[i - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        X[i] = X[i - 1] + h;
        i = i + 1;
    } while (i <= m);
}

void Diff_equation::solve_predictor(float x, float y, float h, int m) {
    if (this -> m > 3) {
        solve_rungekutta(x, y, h, 3);
        for (int i = 4; i <= this->m; i++) {
            Y[i] = Y[i - 1] + h * (55 * f(X[i - 1], Y[i - 1])
				- 59 * f(X[i - 2], Y[i - 2]) + 37 * f(X[i - 3], Y[i - 3])
				- 9 * f(X[i - 4], Y[i - 4])) / 24;
        }
    }
    else solve_rungekutta(x, y, h, m);
}

void Diff_equation::solve_corrector(float x, float y, float h, int m) {
    if (this -> m > 3) {
        solve_rungekutta(x, y, h, 3);
        for (int i = 4; i <= this -> m; i++) {
            Y[i] = method_adams(i, h);
        }
    }
    else solve_rungekutta(x, y, h, m);
}

float Diff_equation :: method_adams(int i, float h) {
    float delta;
    float delta1 = f(X[i - 1], Y[i - 1]) - f(X[i - 2], Y[i - 2]);
    float delta2 = f(X[i - 1], Y[i - 1]) - 2 * f(X[i - 2], Y[i - 2])
		+ f(X[i - 3], Y[i - 3]);
    float delta3 = f(X[i - 1], Y[i - 1]) - 3 * f(X[i - 2], Y[i - 2])
		+ 3 * f(X[i - 3], Y[i - 3]) - f(X[i - 4], Y[i - 4]);
    delta = Y[i - 1] + h * Y[i - 1] + pow(h, 2) / 2 * delta1
		+ 5 / 12 * pow(h, 3) * delta2 + 3 / 8 * pow(h, 4) * delta3;
    return delta;
}

float Diff_equation::accuracy(float h) { return pow(h, 5); }
void Diff_equation::print_method() {
    for (int i = 0; i <= this -> m; i++) {
        cout << "x[" << i << "] = " << X[i] << " y[" << i << "] = " << Y[i] << "\n";
    }
}

float Diff_equation :: approximation(int k, float x, float y) {
	float result = 0;
	float* coef = new float[k + 1];
	float* coef2 = new float[k + 1];
	int n = 1;

	for (int i = 0; i <= k + 1; i++) { coef[i] = 0; coef2[i] = 0; };
	coef[0] = 1;

	for (int j = 1; j <= k + 1; j++) {
		for (int i = 0; i < n; i++) {
			coef[i] *= -1;
			if (i == 1) coef[i] += 1;
			coef2[i + 1] = coef[i] / (i + 1);
		}
		coef2[0] = 1;
		n++;
		for (int i = 0; i <= k + 1; i++) { coef[i] = coef2[i]; }
	}
	for (int i = 0; i <= k + 1; i++) {
		result += coef2[i] * pow(x, i);
	}
	return result;
}

void Diff_equation :: get_equation_approximation(int k, float x, float y) {
	float result = 0;
	float* coef = new float[k + 1];
	float* coef2 = new float[k + 1];
	int n = 1;

	for (int i = 0; i <= k + 1; i++) { coef[i] = 0; coef2[i] = 0; };
	coef[0] = 1;

	for (int j = 1; j <= k + 1; j++) {
		for (int i = 0; i < n; i++) {
			coef[i] *= -1;
			if (i == 1) coef[i] += 1;
			coef2[i + 1] = coef[i] / (i + 1);
		}
		coef2[0] = 1;
		n++;
		for (int i = 0; i <= k + 1; i++) { coef[i] = coef2[i]; }
	}
	cout << "y = ";
	for (int i = 0; i <= k + 1; i++) {
		cout << "(" << coef2[i] << "*x^" << i << ")";
		if (i != k + 1) cout << " + ";
	}
}

void Diff_equation::solve_approximation(float x, float y) {
	int i = 1;
	while (abs(approximation(i, x, y) - approximation(i - 1, x, y)) > pow(10, -10))
	{
		i++;
	};
	get_equation_approximation(i, x, y);

}