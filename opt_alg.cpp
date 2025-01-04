#include"opt_alg.h"

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1,
            matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    } catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1,
                  matrix ud2) {
    try {
        double *p = new double[2]{0, 0};

        int i = 0;
        solution X0(x0);
        solution X1(x0 + d);

        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);

        if (X1.y == X0.y) {
            p[0] = m2d(X0.y);
            p[1] = m2d(X1.y);
            return p;
        }

        if (X1.y > X0.y) {
            d = -d;
            X1.x = x0 + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y >= X0.y) {
                p[0] = m2d(X1.y);
                p[1] = m2d(X0.y - d);
                return p;
            }
        }

        solution X2;
        while (solution::f_calls <= Nmax) {
            ++i;
            X2.x = x0 + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);

            if (X2.y > X1.y) break;
            X0 = X1;
            X1 = X2;
        }

        p[0] = (d > 0) ? X0.x(0, 0) : X2.x(0, 0);
        p[1] = (d > 0) ? X2.x(0, 0) : X0.x(0, 0);
        return p;
    } catch (std::string &ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}


solution fib(
    matrix (*ff)(matrix, matrix, matrix), // funkcja celu
    double a,
    double b,
    double epsilon,
    matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;

        float range = b - a;

        int k = static_cast<int>(ceil(log2(sqrt(5) * range / epsilon) / log2((1 + sqrt(5)) / 2)));
        Xopt.ud = matrix(range);

        vector<double> fibonaccis(k, 1);
        for (int i = 2; i < k; i++) {
            fibonaccis[i] = fibonaccis[i - 2] + fibonaccis[i - 1];
        }

        solution A(a), B(b);
        solution C(B.x - fibonaccis[k - 2] / fibonaccis[k - 1] * (B.x - A.x));
        solution D(A.x + B.x - C.x);

        for (int i = 0; i < k - 2; i++) {
            if (C.fit_fun(ff, ud1, ud2) < D.fit_fun(ff, ud1, ud2)) {
                B = D;
            } else {
                A = C;
            }

            C.x = B.x - fibonaccis[k - i - 2] / fibonaccis[k - i - 1] * (B.x - A.x);
            D.x = A.x + B.x - C.x;
            Xopt.ud.add_row(B.x - A.x);
        }

        Xopt = (C.y < D.y) ? C : D;
        Xopt.flag = 1;
        return Xopt;
    } catch (std::string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax,
             matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        solution A(a), B(b), C((a + b) / 2);
        solution D_0(a);
        solution D_1;

        Xopt.ud = B.x - A.x;

        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        while (true) {
            matrix l1 = A.y * (pow(B.x, 2) - pow(C.x, 2));
            matrix l2 = B.y * (pow(C.x, 2) - pow(A.x, 2));
            matrix l3 = C.y * (pow(A.x, 2) - pow(B.x, 2));
            double l = m2d(l1 + l2 + l3);

            matrix m1 = A.y * (B.x - C.x);
            matrix m2 = B.y * (C.x - A.x);
            matrix m3 = C.y * (A.x - B.x);
            double m = m2d(m1 + m2 + m3);

            if (m <= 0) {
                D_0.fit_fun(ff, ud1, ud2);
                cout << "m <= 0\n";
                Xopt = D_0;
                return Xopt;
            }

            D_1 = solution(0.5 * l / m);
            D_1.fit_fun(ff, ud1, ud2);

            if (A.x <= D_1.x && D_1.x <= C.x) {
                if (D_1.y < C.y) {
                    B = C;
                    C = D_1;
                } else {
                    A = D_1;
                }
            } else if (C.x <= D_1.x && D_1.x <= B.x) {
                if (D_1.y < C.y) {
                    A = C;
                    C = D_1;
                } else {
                    B = D_1;
                }
            } else {
                cout << "D_1 is outside the interval [A.x, C.x] and [C.x, B.x]\n";
                D_0.fit_fun(ff, ud1, ud2);
                Xopt = D_0;
                return Xopt;
            }

            Xopt.ud.add_row(B.x - A.x);

            if (m2d(B.x - A.x) < epsilon || abs(m2d(D_1.x - D_0.x)) < gamma) {
                Xopt = D_1;
                return Xopt;
            }

            if (solution::f_calls > Nmax) {
                Xopt = D_1;
                return Xopt;
            }

            D_0 = D_1;
        }
    } catch
    (string ex_info) {
        cerr << ex_info;
        throw ("solution lag(...):\n" + ex_info);
    }
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax,
            matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        //Tu wpisz kod funkcji

        return XB;
    } catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon,
               int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,matrix ud2) {
    try {
        solution Xopt;
        double alpha = 1, beta = 0.5, gamma = 2, delta = 0.5, s = 1;
        solution X(x0), X1;

        while (true) {
            X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
            if (norm(X1.x - X.x) < epsilon) {
                Xopt = X1;
                Xopt.flag = 1;
                break;
            }

            if (solution::f_calls > Nmax) {
                Xopt = X1;
                Xopt.flag = 0;
                break;
            }

            c = dc * c;
            X = X1;
        }

        return Xopt;
    } catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma,
                double delta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);
        matrix D = ident_mat(n); // Macierz jednostkowa
        int N = n + 1; // Liczba wierzchołków simpleksu

        solution* S = new solution[N]; // Simpleks

        // Inicjalizacja simpleksu
        S[0].x = x0;
        S[0].fit_fun(ff, ud1, ud2);
        for (int i = 1; i < N; ++i) {
            S[i].x = S[0].x + s * D[i - 1];
            S[i].fit_fun(ff, ud1, ud2);
        }

        solution PO, PE, PZ;
        matrix pc;
        int i_min, i_max;

        while (true) {
            // Znalezienie indeksów punktów minimalnego i maksymalnego
            i_min = i_max = 0;
            for (int i = 1; i < N; ++i) {
                if (S[i_min].y > S[i].y)
                    i_min = i;
                if (S[i_max].y < S[i].y)
                    i_max = i;
            }

            // Obliczenie środka ciężkości (p)
            pc = matrix(n, 1);
            for (int i = 0; i < N; ++i) {
                if (i != i_max)
                    pc = pc + S[i].x;
            }
            pc = pc / (N - 1.0);

            // Refleksja
            PO.x = pc + alpha * (pc - S[i_max].x);
            PO.fit_fun(ff, ud1, ud2);

            if (PO.y < S[i_min].y) {
                // Ekspansja
                PE.x = pc + gamma * (PO.x - pc);
                PE.fit_fun(ff, ud1, ud2);
                S[i_max] = (PE.y < PO.y) ? PE : PO;
            } else if (PO.y < S[i_max].y) {
                // Akceptacja refleksji
                S[i_max] = PO;
            } else {
                // Kontrakcja
                PZ.x = pc + beta * (S[i_max].x - pc);
                PZ.fit_fun(ff, ud1, ud2);
                if (PZ.y < S[i_max].y) {
                    S[i_max] = PZ;
                } else {
                    // Redukcja
                    for (int i = 0; i < N; ++i) {
                        if (i != i_min) {
                            S[i].x = S[i_min].x + delta * (S[i].x - S[i_min].x);
                            S[i].fit_fun(ff, ud1, ud2);
                        }
                    }
                }
            }

            // Warunek stopu
            double max_dist = 0.0;
            for (int i = 0; i < N; ++i) {
                double dist = norm(S[i].x - S[i_min].x);
                if (dist > max_dist)
                    max_dist = dist;
            }

            if (max_dist < epsilon) {
                solution Xopt = S[i_min];
                Xopt.flag = 1;
                delete[] S;
                return Xopt;
            }

            if (solution::f_calls > Nmax) {
                solution Xopt = S[i_min];
                Xopt.flag = 0;
                delete[] S;
                return Xopt;
            }
        }
    } catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}


solution SD(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution CG(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0,
            double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix),
                matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    } catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}
