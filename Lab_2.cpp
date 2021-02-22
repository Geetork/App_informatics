#include <iostream>
#include "Diff_equation.h"
using namespace std;

int main() {
    setlocale(LC_ALL, "Russian");
    float x0, y0, h;
    int m;

    cout << "Метод Рунге-Кутта: \nx0 = "; cin >> x0;
    cout << "y0 = "; cin >> y0;
    cout << "h = "; cin >> h;
    cout << "m = "; cin >> m;

    Diff_equation diff_equation_rungekutta(m);
    diff_equation_rungekutta.solve_rungekutta(x0, y0, h, m);
    diff_equation_rungekutta.print_method();
    cout << "\nТочность метода = " << diff_equation_rungekutta.accuracy(h);

    cout << "\nМетод прогноза: \nx0 = "; cin >> x0;
    cout << "y0 = "; cin >> y0;
    cout << "h = "; cin >> h;
    cout << "m = "; cin >> m;

    Diff_equation diff_equation_predictor(m);
    diff_equation_rungekutta.solve_predictor(x0, y0, h, m);
    diff_equation_rungekutta.print_method();

    cout << "Метод коррекции: \nx0 = "; cin >> x0;
    cout << "y0 = "; cin >> y0;
    cout << "h = "; cin >> h;
    cout << "m = "; cin >> m;

    Diff_equation diff_equation_corrector(m);
    diff_equation_rungekutta.solve_corrector(x0, y0, h, m);
    diff_equation_rungekutta.print_method();

    cout << "Метод приближений: \nx0 = "; cin >> x0;
    cout << "y0 = "; cin >> y0;

    diff_equation_rungekutta.solve_approximation(x0, y0);

    return 0;
}