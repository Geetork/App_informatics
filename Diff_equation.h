#include <iostream>

class Diff_equation {
private:
    int m;
public:
    float* Y;
    float* X;

    Diff_equation(int m) {
        this -> m = m;
        X = new float[this -> m];
        Y = new float[this -> m];
    }

    float f(float, float);
    void solve_rungekutta(float, float, float, int);
    void solve_predictor(float, float, float, int);
    void solve_corrector(float, float, float, int);
    void solve_approximation(float, float);

    float accuracy(float);
    float method_adams(int, float);

    float approximation(int, float, float);
    void get_equation_approximation(int, float, float);

    void print_method();
};




