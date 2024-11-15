/*********************************************
Kod stanowi uzupe?nienie materia?ów do ?wicze?
w ramach przedmiotu metody optymalizacji.
Kod udost?pniony na licencji CC BY-SA 3.0
Autor: dr in?. ?ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

int main() {
    try {
        lab1();
    } catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0() {
    //Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1() {
    double epsilon = 1e-5, gamma = 1e-200;
    int Nmax = 1000;
    double d = 1.0;
    double alpha = 1.5;

    // exp
    cout << "Optymalizacja funkcji testowej metoda ekspansji:" << endl;
    double x0 = 100;
    double *interval = expansion(ff1T, x0, d, alpha, Nmax);
    double a = interval[0];
    double b = interval[1];
    cout << "Przedzial znaleziony metoda ekspansji: [" << a << ", " << b << "]" << endl;
    cout << "ile fcall = " << solution::f_calls << endl;
    delete[] interval;

    // fib
    solution::clear_calls();
    solution opt;
    opt = fib(ff1T, a, b, epsilon);
    cout << "\nMetoda Fibonacciego\n" << opt << endl;

    // lag
    solution::clear_calls();
    opt = lag(ff1T, a, b, epsilon, gamma, Nmax);
    cout << "Metoda Lagrange'a\n" << opt << endl;
}

void lab2() {
}

void lab3() {
}

void lab4() {
}
