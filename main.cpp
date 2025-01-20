/*********************************************
Kod stanowi uzupe?nienie materia?ów do ?wicze?
w ramach przedmiotu metody optymalizacji.
Kod udost?pniony na licencji CC BY-SA 3.0
Autor: dr in?. ?ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include <float.h>

#include"opt_alg.h"
#define SEPARATOR ";"

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

int main() {
    try {
        lab3();
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

void save_single_column_matrix_to_file(matrix m, ofstream &file) {
    for (int i = 0; i < get_len(m); i++) {
        file << (i + 1) << SEPARATOR << to_string_with_comma(m(i, 0)) << "\n";
    }
}

void lab1() {
    double epsilon = 1e-10, gamma = 1e-200;
    int Nmax = 1000;
    double d = 1.0;
    double alpha = 5.0;

    // Z ZAWEZANIEM PRZEDZIALU (TABELA 1, TABELA 2)
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


    // BEZ ZAWEZANIA PRZEDZIALU (WYKRES)
    /*
    double a = -100, b = 100;

    ofstream fibFile("fibonacci_results-no-expansion.csv");
    ofstream lagFile("lagrange_results-no-expansion.csv");

    fibFile << "i" << SEPARATOR << "(b-a)" << "\n";
    lagFile << "i" << SEPARATOR << "(b-a)" << "\n";

    solution optFib = fib(ff1T, a, b, epsilon);
    save_single_column_matrix_to_file(optFib.ud, fibFile);
    cout << optFib << endl;
    solution::clear_calls();


    solution optLag = lag(ff1T, a, b, epsilon, gamma, Nmax);
    save_single_column_matrix_to_file(optLag.ud, lagFile);
    cout << optLag << endl;
    solution::clear_calls();

    fibFile.close();
    lagFile.close();
    */

    double a = 1e-4; // 1 cm ^ 2
    double b = 1e-2; // 100 cm^2
    // double epsilon = 1e-5; // jakas dokladnosc
    // double gamma = 1e-200; // jakas dokladnosc ale wieksza
    int maxIter = 1000;

    ofstream fibFile("simulation-fib.csv");
    fibFile << "i" << SEPARATOR << "vA" << SEPARATOR << "vB" << SEPARATOR << "tB" << endl;
    solution opt = fib(ff1R, a, b, epsilon);
    cout << opt << endl;

    matrix *simDataFib = getSimulationData1R(opt.x);
    fibFile << hcat(simDataFib[0], simDataFib[1]);
    fibFile.close();
    solution::clear_calls();


    ofstream lagFile("simulation-lag.csv");
    lagFile << "i" << SEPARATOR << "vA" << SEPARATOR << "vB" << SEPARATOR << "tB" << endl;
    opt = lag(ff1R, a, b, epsilon, gamma, maxIter);
    cout << opt;
    matrix *simDataLag = getSimulationData1R(opt.x);
    lagFile << hcat(simDataLag[0], simDataLag[1]);
    lagFile.close();
    solution::clear_calls();
}


bool isLegit(matrix x0, double a) {
    return x0(0) >= 1 && x0(1) >= 1 && pow(x0(0), 2) + pow(x0(1), 2) <= pow(a, 2);
}

void lab2() {
    double a = 4.4934;
    double epsilon = 1e-3;
    int nMax = 1e4;
    double cEx = 1, dcEx = 2;
    double cIn = 10, dcIn = 0.5;

    // PROBLEM TESTOWY
    /*
    ofstream fileExternal("external_results-" + format("{:.2f}", a) + ".csv");
    ofstream fileInternal("internal_results-" + format("{:.2f}", a) + ".csv");

    fileExternal <<
            "i" << SEPARATOR <<
            "x1_0" << SEPARATOR <<
            "x2_0" << SEPARATOR <<
            "x1*" << SEPARATOR <<
            "x2*" << SEPARATOR <<
            "r*" << SEPARATOR <<
            "y*" << SEPARATOR <<
            "f_calls" << endl;
    fileInternal <<
            "i" << SEPARATOR <<
            "x1_0" << SEPARATOR <<
            "x2_0" << SEPARATOR <<
            "x1*" << SEPARATOR <<
            "x2*" << SEPARATOR <<
            "r*" << SEPARATOR <<
            "y*" << SEPARATOR <<
            "f_calls" << endl;

    for (int i = 0; i < 100; ++i) {
        matrix x0;
        solution optEx, optIn;

        do {
            x0 = 5 * rand_mat(2, 1) + 1;
        } while (!isLegit(x0, a));

        optEx = pen(ff2Tz, x0, cEx, dcEx, epsilon, nMax, a);
        fileExternal <<
                to_string_with_comma((i + 1)) << SEPARATOR <<
                to_string_with_comma(x0(0)) << SEPARATOR <<
                to_string_with_comma(x0(1)) << SEPARATOR <<
                to_string_with_comma(optEx.x(0)) << SEPARATOR <<
                to_string_with_comma(optEx.x(1)) << SEPARATOR <<
                to_string_with_comma(norm(optEx.x)) << SEPARATOR <<
                to_string_with_comma(m2d(optEx.y)) << SEPARATOR <<
                to_string_with_comma(m2d(solution::f_calls)) << endl;
        solution::clear_calls();

        optIn = pen(ff2Tw, x0, cIn, dcIn, epsilon, nMax, a);
        fileInternal <<
                to_string_with_comma((i + 1)) << SEPARATOR <<
                to_string_with_comma(x0(0)) << SEPARATOR <<
                to_string_with_comma(x0(1)) << SEPARATOR <<
                to_string_with_comma(optIn.x(0)) << SEPARATOR <<
                to_string_with_comma(optIn.x(1)) << SEPARATOR <<
                to_string_with_comma(norm(optIn.x)) << SEPARATOR <<
                to_string_with_comma(m2d(optIn.y)) << SEPARATOR <<
                to_string_with_comma(m2d(solution::f_calls)) << endl;
        solution::clear_calls();
    }

    fileExternal.close();
    fileInternal.close();
    */

    // PROBLEM RZECZYWISTY
    matrix x0 = matrix(2, 1);
    x0(0) = 20 * m2d(rand_mat()) - 10;
    x0(1) = 30 * m2d(rand_mat()) - 15;
    cout << x0 << endl << endl;

    solution optR = pen(ff2R, x0, cEx, dcEx, epsilon, nMax);
    cout << optR << endl;

    ofstream simFile("simulation-nm.csv");
    simFile << "t" << SEPARATOR << "x" << SEPARATOR << "y" << endl;
    matrix *simData = getSimulationData2R(optR.x);
    for (int i = 0; i < get_len(*simData); i++) {
        simFile << to_string_with_comma(simData[0](i, 0)) << SEPARATOR << to_string_with_comma(simData[1](i, 0)) <<
                SEPARATOR << to_string_with_comma(simData[1](i, 2)) << endl;
    }
    // simFile << hcat(simData[0], simData[1]);
    simFile.close();

    solution::clear_calls();
}

void lab3() {
    solution opt;
    int Nmax = 10e3;
    double h = 0.12;
    double epsilon = 1e-3;

    ofstream sd("sg_results.csv");
    ofstream cg("cg_results.csv");
    ofstream n("n_results.csv");

    for (int i = 0; i < 100; ++i) {
        matrix x0 = 20 * rand_mat(2, 1) - 10;

        solution SD_sol = SD(ff3T, gf3, x0, h, epsilon, Nmax);
        sd << x0(0) << SEPARATOR << x0(1) << SEPARATOR << SD_sol.x(0) << SEPARATOR << SD_sol.x(1) << SEPARATOR
                << SD_sol.y[0] << SEPARATOR << solution::f_calls << SEPARATOR << solution::g_calls << endl;

        solution::clear_calls();

        solution CG_sol = CG(ff3T, gf3, x0, h, epsilon, Nmax);
        cg << x0(0) << SEPARATOR << x0(1) << SEPARATOR << CG_sol.x(0) << SEPARATOR << CG_sol.x(1) << SEPARATOR
                << CG_sol.y[0] << SEPARATOR << solution::f_calls << SEPARATOR << solution::g_calls << endl;

        solution::clear_calls();

        solution N_sol = Newton(ff3T, gf3, hf3, x0, h, epsilon, Nmax);
        n << x0(0) << SEPARATOR << x0(1) << SEPARATOR << N_sol.x(0) << SEPARATOR << N_sol.x(1) << SEPARATOR
                << N_sol.y[0] << SEPARATOR << solution::f_calls << SEPARATOR << solution::g_calls << endl;
    }

    sd.close();
    cg.close();
    n.close();
    solution::clear_calls();
}

void lab4() {
}
