#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y;

    double x_val = x(0);
    y = -cos(0.1 * x_val) * exp(-pow(0.1 * x_val - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x_val, 2);
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    Y[0].~matrix();
    Y[1].~matrix();
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

/// <summary>
/// 
/// </summary>
matrix df1(
    double t, // 
    matrix Y, // poczatkowe wartosci
    matrix ud1,
    matrix ud2) {
    double a = 0.98, b = 0.63, g = 9.81;
    double PA = 0.5, TA = 90, PB = 1, DB = 0.00365665, Fin = 0.01, Tin = 20;

    matrix dY(3, 1);

    double FAout = (Y(0) > 0)
                       ? a * b * m2d(ud2) * sqrt(2 * g * Y(0) / PA)
                       : 0;

    double FBout = (Y(1) > 0)
                       ? a * b * DB * sqrt(2 * g * Y(1) / PB)
                       : 0;

    dY(0) = -FAout;
    dY(1) = FAout + Fin - FBout;
    dY(2) = Fin / Y(1) * (Tin - Y(2)) + FAout / Y(1) * (TA - Y(2));

    return dY;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(4, 1);

    double r = 0.12;
    double C = 0.47;
    double ro = 1.2;
    double S = M_PI * pow(r, 2);
    double m = 0.6;
    double g = 9.81;

    double Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1));
    double Fmx = ro * Y(3) * m2d(ud2) * M_PI * pow(r, 3);
    double Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3));
    double Fmy = ro * Y(1) * m2d(ud2) * M_PI * pow(r, 3);

    dY(0) = Y(1);
    dY(1) = (-Fmx - Dx) / m;
    dY(2) = Y(3);
    dY(3) = (-Fmy - Dy - m * g) / m;

    return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    double x_val = x(0);
    matrix y = -cos(0.1 * x_val) * exp(-pow(0.1 * x_val - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x_val, 2);
    return y;
}

matrix *getSimulationData1R(matrix x, matrix ud1, matrix ud2) {
    double t0 = 0;
    double dt = 1;
    double t_end = 2000;
    double vA = 5; // obj. wody w A
    double vB = 1; // obj. wody w B
    double tB_0 = 20; // temp. pocz. w B
    matrix initialValues = matrix(3, new double[3]{vA, vB, tB_0});

    return solve_ode(df1, t0, dt, t_end, initialValues, ud1, x);
}

/// <summary>
/// Funkcja celu - LAB 1
/// </summary>
matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    matrix *simulationData = getSimulationData1R(x, ud1, ud2);;

    int dataLength = get_len(simulationData[0]);

    double maxValue = simulationData[1](0, 2);
    for (int index = 1; index < dataLength; ++index) {
        maxValue = std::max(maxValue, simulationData[1](index, 2));
    }

    matrix result = abs(maxValue - 50);
    return result;
}

double ff2T(matrix x, matrix ud1, matrix ud2) {
    double pom = M_PI * std::sqrt(std::pow(x(0) / M_PI, 2) + std::pow(x(1) / M_PI, 2));
    return std::sin(pom) / pom;
}

double boundary(int i, matrix x, double a) {
    switch (i) {
        case 0:
            return -x(0) + 1.0;
        case 1:
            return -x(1) + 1.0;
        case 2:
            return norm(x) - a;
        default:
            return 0.0;
    }
}

matrix ff2Tw(matrix x, matrix ud1, matrix ud2) {
    double y = ff2T(x);

    double gValue = -x(0) + 1;
    if (gValue > 0) {
        y = 1e10;
    } else {
        y -= m2d(ud2 / gValue);
    }
    return y;
}

matrix ff2Tz(matrix x, matrix ud1, matrix ud2) {
    double y = ff2T(x);

    for (int i = 0; i < 3; ++i) {
        double gValue = boundary(i, x(0), x(1), ud1(0));
        if (gValue > 0)
            y += m2d(ud2 * pow(gValue, 2));
    }

    return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0(4, new double[4]{0, x(0), 100, 0});
    matrix *Y = solve_ode(df2, 0, 0.01, 7, Y0, ud1, x(1));
    int n = get_len(Y[0]);
    int i50 = 0, i0 = 0;
    for (int i = 0; i < n; ++i) {
        if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
            i50 = i;
        if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
            i0 = i;
    }

    y = -Y[1](i0, 0);
    if (abs(x(0) - 10) > 0)
        y = y + ud2 * pow(abs(x(0) - 10), 2);
    if (abs(x(1) - 20) > 0)
        y = y + ud2 * pow(abs(x(1) - 20), 2);
    if (abs(Y[1](i50, 0) - 5) - 1 > 0)
        y = y + ud2 * pow(abs(Y[1](i50, 0) - 5) - 1, 2);
    return y;
}
