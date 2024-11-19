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
#define SEPARATOR ";"

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

bool is_global_min(double x, double y) {
    const double global_x = 62.7482;
    const double global_y = -0.921148;

    const double tolerance = 0.001;

    return (std::abs(x - global_x) < tolerance) && (std::abs(y - global_y) < tolerance);
}

string to_string_with_comma(double value) {
    ostringstream oss;
    oss << value;
    string result = oss.str();
    replace(result.begin(), result.end(), '.', ',');
    return result;
}

void save_single_column_matrix_to_file(matrix m, ofstream& file) {
    for (int i = 0; i < get_len(m); i++) {
        file << (i + 1) << SEPARATOR << to_string_with_comma(m(i, 0)) << "\n";
    }
}

void lab1() {
    double epsilon = 1e-5, gamma = 1e-200;
    int Nmax = 1000;
    double d = 1.0;
    double alpha = 5.0;

    // Z ZAW??ANIEM PRZEDZIA?U (TABELA 1, TABELA 2)
    /*
    srand(time(0));

    ofstream expFile("expansion_results-" + format("{:.2f}", alpha) + ".csv");
    ofstream fibFile("fibonacci_results-" + format("{:.2f}", alpha) + ".csv");
    ofstream lagFile("lagrange_results-" + format("{:.2f}", alpha) + ".csv");

    expFile << "alpha" << SEPARATOR << "x0" << SEPARATOR << "a" << SEPARATOR << "b" << SEPARATOR << "f_calls\n";
    fibFile << "alpha" << SEPARATOR << "x*" << SEPARATOR << "y*" << SEPARATOR << "f_calls" << SEPARATOR << "min(local/global)\n";
    lagFile << "alpha" << SEPARATOR << "x*" << SEPARATOR << "y*" << SEPARATOR << "f_calls" << SEPARATOR << "min(local/global)\n";

    for (int i = 0; i < 100; i++) {
        // EXPANSION
        double x0 = m2d(200 * rand_mat() - 100);
        double *interval = expansion(ff1T, x0, d, alpha, Nmax);
        double a = interval[0];
        double b = interval[1];
        expFile << to_string_with_comma(alpha) << SEPARATOR << to_string_with_comma(x0) << SEPARATOR << to_string_with_comma(a) << SEPARATOR << to_string_with_comma(b) << SEPARATOR << solution::f_calls << "\n";
        delete[] interval;

        // FIBONACCI
        solution::clear_calls();
        solution opt = fib(ff1T, a, b, epsilon);
        std::string minType = is_global_min(m2d(opt.x), m2d(opt.y)) ? "globalne" : "lokalne";
        fibFile << to_string_with_comma(alpha) << SEPARATOR << to_string_with_comma(m2d(opt.x)) << SEPARATOR << to_string_with_comma(m2d(opt.y)) << SEPARATOR << solution::f_calls << SEPARATOR << minType << "\n";

        // LAGRANGE
        solution::clear_calls();
        opt = lag(ff1T, a, b, epsilon, gamma, Nmax);
        minType = is_global_min(m2d(opt.x), m2d(opt.y)) ? "globalne" : "lokalne";
        lagFile << to_string_with_comma(alpha) << SEPARATOR << to_string_with_comma(m2d(opt.x)) << SEPARATOR << to_string_with_comma(m2d(opt.y)) << SEPARATOR << solution::f_calls << SEPARATOR << minType << "\n";
    }

    expFile.close();
    fibFile.close();
    lagFile.close();
    */

    // BEZ ZAW??ANIA PRZEDZIA?U (WYKRES)
    double a = -100, b = 100;

    ofstream fibFile("fibonacci_results-no-expansion.csv");
    ofstream lagFile("lagrange_results-no-expansion.csv");

    fibFile << "i" << SEPARATOR << "(b-a)" << "\n";
    lagFile << "i" << SEPARATOR << "(b-a)" << "\n";

    solution optFib = fib(ff1T, a, b, epsilon);
    save_single_column_matrix_to_file(optFib.ud, fibFile);
    solution::clear_calls();

    solution optLag = lag(ff1T, a, b, epsilon, gamma, Nmax);
    save_single_column_matrix_to_file(optLag.ud, lagFile);
    solution::clear_calls();

    fibFile.close();
    lagFile.close();
}


void lab2() {
}

void lab3() {
}

void lab4() {
}
