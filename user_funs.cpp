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

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
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

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    double x_val = x(0);
    matrix y = -cos(0.1 * x_val) * exp(-pow(0.1 * x_val - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x_val, 2);
    return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    matrix result;
    double vA = 5, vB = 1, tB_0 = 20;
    matrix initialValues = matrix(3, new double[3]{vA, vB, tB_0});

    matrix *simulationData = solve_ode(df1, 0, 1, 2000, initialValues, ud1, ud2);
    int dataLength = get_len(simulationData[0]);

    double maxValue = simulationData[1](0, 2);
    for (int index = 1; index < dataLength; ++index) {
        maxValue = std::max(maxValue, simulationData[1](index, 2));
    }

    result = abs(maxValue - 50);
    return result;
}
