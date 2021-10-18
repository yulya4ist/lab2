#include <iostream>
#include "matrix_work.h"
#include "vector"
#include "personal_vector.h"
#define not_square 2
#define square 1
//Данные макросы передаются при считывании матрицы из файла
//Параметр 1 означает, что матрица квадратная, параметр 2 - что не квадратная
#define link_data_in "/home/qw/Рабочий стол/laba2_1/"
// Данный макрос используется для указания пути к файлам с исходными данными
#define link_data_out "/home/qw/Рабочий стол/laba2_1/"
// Данный макрос используется для указания пути к файлам с результатами работы программы

using namespace std;
using namespace pv_id;

double dihotomy(double (*f)(matrix<double> const&, double), double a, double b, double eps, matrix<double> const &A);

template <typename float_type>
vector<float_type> residual (matrix<float_type> const &A, vector<float_type> const &B, vector<float_type> const &X)
//Вычисления вектора невязки
{
    vector<float_type> B1 = A*X;
    vector<float_type> v_residual = B - B1;
    return v_residual;
}

vector<double> yacoby_method(matrix<double> const&A, vector<double> const &B,
                             vector<double> const &X_0, double eps, ofstream& out, int& iter)
{
    matrix<double> C(A);
    vector<double> Y(B);
    for (int i = 0; i < A.lin(); ++i) {
        for (int j = 0; j < A.col(); ++j) {
            if (i != j) {
                C[i][j] = -A[i][j]/A[i][i];
            } else {
                C[i][j] = 0;
            }
        }
        Y[i] = B[i]/A[i][i];
    }
    out << endl << "Матрица C для метода Якоби" << endl;
    C.output_f(out);
    out << endl << "Вектор Y для метода Якоби" << endl;
    tVectorOut(Y,out);
    out << endl << "Октаэдрическая норма матрицы C для метода Якоби" << endl;
    out << C.abs_1() << endl;
    out << endl << "Кубическая норма матрицы C для метода Якоби" << endl;
    out << C.abs_inf() << endl;
    vector<double> X_prev = X_0;
    vector<double> X = C*X_0+Y;
    iter = 0;
    while (abs_1(X-X_prev) > eps*(1-C.abs_1())/C.abs_1()) {
        X_prev = X;
        X = C*X_prev+Y;
        iter++;
    }
    return X;
}

vector<double> zeidel_method(matrix<double> const&A, vector<double> const &B,
                             vector<double> const &X_0, double eps, ofstream& out, int& iter)
{
    matrix<double> U(A);
    matrix<double> D(A);
    matrix<double> L(A);
    for (int i = 0; i < A.lin(); ++i) {
        for (int j = 0; j < A.col(); ++j) {
            if (i>j) {
                U[i][j] = 0;
                D[i][j] = 0;
            } else if (i<j) {
                D[i][j] = 0;
                L[i][j] = 0;
            } else {
                U[i][j] = 0;
                L[i][j] = 0;
            }
        }
    }
    matrix<double> C(A.col());
    C = -((D+L).inverse2()*U);
    vector<double> Y(B);
    Y = (D+L).inverse2()*B;
    out << endl << "Матрица C для метода Зейделя" << endl;
    C.output_f(out);
    out << endl << "Вектор Y для метода Зейделя" << endl;
    tVectorOut(Y,out);
    out << endl << "Октаэдрическая норма матрицы C для метода Зейделя" << endl;
    out << C.abs_1() << endl;
    out << endl << "Кубическая норма матрицы C для метода Зейделя" << endl;
    out << C.abs_inf() << endl;
    matrix<double> Cu(C);
    matrix<double> Cd(C);
    matrix<double> Cl(C);
    for (int i = 0; i < C.lin(); ++i) {
        for (int j = 0; j < C.col(); ++j) {
            if (i>j) {
                Cu[i][j] = 0;
                Cd[i][j] = 0;
            } else if (i<j) {
                Cd[i][j] = 0;
                Cl[i][j] = 0;
            } else {
                Cu[i][j] = 0;
                Cl[i][j] = 0;
            }
        }
    }
    out << endl << "Матрицы C_u, C_d, C_l для метода Зейделя" << endl;
    out << endl << "Матрица C_u" << endl;
    Cu.output_f(out);
    out << endl << "Матрица C_d" << endl;
    Cd.output_f(out);
    out << endl << "Матрица C_l" << endl;
    Cl.output_f(out);
    iter = 0;
    vector<double> X_prev = X_0;
    vector<double> X = C*X_0+Y;
    while (abs_1(X-X_prev) > eps*(1-C.abs_1())/C.abs_1()) {
        iter++;
        X_prev = X;
        X = C*X_prev+Y;
    }
    return X;
}

vector<double> relacs_method(matrix<double> const&A, vector<double> const &B,
                             vector<double> const &X_0, double eps, double w, ofstream& out, int& iter)
{
    matrix<double> U(A);
    matrix<double> D(A);
    matrix<double> L(A);
    for (int i = 0; i < A.lin(); ++i) {
        for (int j = 0; j < A.col(); ++j) {
            if (i>j) {
                U[i][j] = 0;
                D[i][j] = 0;
            } else if (i<j) {
                D[i][j] = 0;
                L[i][j] = 0;
            } else {
                U[i][j] = 0;
                L[i][j] = 0;
            }
        }
    }
    matrix<double> C(A.col());
    C = (D+w*L).inverse2()*(D+w*(L-A));
    vector<double> Y(B);
    Y = w*(D+w*L).inverse2()*B;
    out << endl << "Матрица C для метода релаксации с w=" << w << endl;
    C.output_f(out);
    out << endl << "Вектор Y для метода релаксации с w=" << w << endl;
    tVectorOut(Y,out);
    out << endl << "Октаэдрическая норма матрицы C для метода релаксации с w=" << w << endl;
    out << C.abs_1() << endl;
    out << endl << "Кубическая норма матрицы C для метода релаксации с w=" << w << endl;
    out << C.abs_inf() << endl;
    vector<double> X_prev = X_0;
    vector<double> X = C*X_0+Y;
    iter = 0;
    while (abs_1(X-X_prev) > eps*(1-C.abs_1())/C.abs_1()) {
        iter++;
        X_prev = X;
        X = C*X_prev+Y;
    }
    return X;
}

vector<double> three_diagonal_relacs_method(matrix<double> const&A, vector<double> const &B,
                                            vector<double> const &X_0, double eps, double w)
{
    vector<double> X_prev = X_0;
    vector<double> X = X_prev;
    int i = 0;
    X[i] = (1-w)*X_prev[i]-w*(A[i][i+1]/A[i][i])*X_prev[i+1]+w*B[i]/A[i][i];
    for (i = 1; i < X.size()-1; ++i) {
        X[i] = -w*(A[i][i-1]/A[i][i])*X[i-1]+(1-w)*X_prev[i]-w*(A[i][i+1]/A[i][i])*X_prev[i+1]+w*B[i]/A[i][i];
    }
    i = X.size()-1;
    X[i] = -w*(A[i][i-1]/A[i][i])*X[i-1]+(1-w)*X_prev[i]+w*B[i]/A[i][i];
    while (abs_1(X-X_prev) > eps) {
        X_prev = X;
        i = 0;
        X[i] = (1-w)*X_prev[i]-w*(A[i][i+1]/A[i][i])*X_prev[i+1]+w*B[i]/A[i][i];
        for (i = 1; i < X.size()-1; ++i) {
            X[i] = -w*(A[i][i-1]/A[i][i])*X[i-1]+(1-w)*X_prev[i]-w*(A[i][i+1]/A[i][i])*X_prev[i+1]+w*B[i]/A[i][i];
        }
        i = X.size()-1;
        X[i] = -w*(A[i][i-1]/A[i][i])*X[i-1]+(1-w)*X_prev[i]+w*B[i]/A[i][i];
    }
    return X;
}



vector<double> simple_iteraition(matrix<double> const&A, vector<double> const &B,
                                 vector<double> const &X_0, double eps, ofstream& out, int& iter)
{
    matrix<double> A_new(A);
    vector<double> B_new(B);
    matrix<double> E(A);
    E.init_unit();
    double a = -5;
    double b = -4;
    while (((-a*A)+E).abs_1()>((-b*A)+E).abs_1()) {
        a += 1;
        b += 1;
    }
    auto f = [](matrix<double> const &A, double x){
        matrix<double> E(A);
        E.init_unit();
        return ((-x*A)+E).abs_1();
    };
    double tau = dihotomy(f,a-1,b,eps,A);
    if ((tau*A-E).abs_1()>=1) {
        a = -5;
        b = -4;
        while (((-a*A)+E).abs_inf()>((-b*A)+E).abs_inf()) {
            a += 1;
            b += 1;
        }
        auto g = [](matrix<double> const &A, double x){
            matrix<double> E(A);
            E.init_unit();
            return ((-x*A)+E).abs_inf();
        };
        tau = dihotomy(g,a-1,b,eps,A);
        if ((tau*A-E).abs_inf()>=1) {
            for (int i = 1; i < pow(2,A.lin()); ++i) {
                matrix<double> A_changed(A);
                for (int j = 0; j < A.lin(); ++j) {
                    int sign = i/pow(2,j);
                    sign = sign % 2;
                    if (sign == 1) {
                        for (int u = 0; u < A.col(); u++) {
                            A_changed[j][u] = -A[j][u];
                        }
                    }
                }
                a = -5;
                b = -4;
                while (((-a*A_changed)+E).abs_1()>((-b*A_changed)+E).abs_1()) {
                    a += 1;
                    b += 1;
                }
                tau = dihotomy(f,a-1,b,eps,A_changed);
                if (f(A_changed,tau)<1) {
                    A_new = A_changed;
                    for (int j = 0; j < A.lin(); ++j) {
                        int sign = i/pow(2,j);
                        sign = sign % 2;
                        if (sign == 1) {
                            B_new[j] = -B[j];
                        }
                    }
                    break;
                }
                a = -5;
                b = -4;
                while (((-a*A)+E).abs_inf()>((-b*A)+E).abs_inf()) {
                    a += 1;
                    b += 1;
                }
                tau = dihotomy(g,a-1,b,eps,A);
                if (g(A_changed,tau)<1) {
                    A_new = A_changed;
                    for (int j = 0; j < A.lin(); ++j) {
                        int sign = i/pow(2,j);
                        sign = sign % 2;
                        if (sign == 1) {
                            B_new[j] = -B[j];
                        }
                    }
                    break;
                }
            }
        }
    }
    if (f(A_new,tau)>=1 && (((-tau*A_new)+E).abs_inf()>=1)) {
        cout << "Как бы не меняли строки не сходится" << endl;
        return B;
    }
    iter = 0;
    matrix<double> C(A_new);
    C = -(tau*A_new-E);
    out << endl << "Матрица C для метода простой итерации" << endl;
    C.output_f(out);
    out << endl << "Вектор Y для метода простой итерации" << endl;
    tVectorOut(tau*B_new,out);
    out << endl << "Октаэдрическая норма матрицы C для метода простой итерации" << endl;
    out << C.abs_1() << endl;
    out << endl << "Кубическая норма матрицы C для метода простой итерации" << endl;
    out << C.abs_inf() << endl;
    vector<double> X_prev = X_0;
    vector<double> X = (C*X_0)+(tau*B_new);
    while (abs_1(X-X_prev) > eps*(1-C.abs_1())/C.abs_1()) {
        iter++;
        X_prev = X;
        X = (C*X_prev)+(tau*B_new);
    }
    return X;
}

double dihotomy(double (*f)(matrix<double> const&, double), double a, double b, double eps, matrix<double> const &A){
    double result = (a+b)/2;
    double delta = eps/4;
    double x1 = result-delta;
    double x2 = result+delta;
    while (eps < (b-a)) {
        if (f(A,x1)<f(A,x2)) {
            b = x2;
        } else {
            a = x1;
        }
        result = (a+b)/2;
        x1 = result-delta;
        x2 = result+delta;
    }
    return (a+b)/2;
}

void lab2(char* const &input1, char* const &input2, char* const& output_name)
{
    ofstream n_out;
    n_out.open(output_name);
    ofstream crash;
    crash.open(link_data_out"crash.txt");
    try{
        n_out << "Результаты работы программы:" << endl << endl;
        matrix<double> A(square, input1);
        n_out << "Исходная матрица:" << std::endl;
        A.output_f(n_out);
        n_out << std::endl;
        vector<double> B(read_vector<double>(input2));
        n_out << "Исходная правая часть СЛАУ:" << endl;
        tVectorOut(B,n_out);
        n_out << endl;

        n_out << "I. Для точности eps = 1e-4" << endl;

        int iter_s = 0;
        vector<double> Xs = simple_iteraition(A,B,B,1e-4,n_out,iter_s);

        n_out << endl << "Результат метода простой итерации" << endl;
        tVectorOut(Xs,n_out);
        n_out << endl << "Число итераций для метода простой итерации: " << iter_s << endl;

        int iter0 = 0;
        vector<double> X = yacoby_method(A,B,B,1e-4, n_out, iter0);

        n_out << endl << "Результат метода Якоби" << endl;
        tVectorOut(X,n_out);
        n_out << endl << "Число итераций для метода Якоби: " << iter0 << endl;

        int iter1 = 0;
        vector<double> X1 = zeidel_method(A,B,B,1e-4,n_out,iter1);

        n_out << endl << "Результат метода Зейделя" << endl;
        tVectorOut(X1,n_out);
        n_out << endl << "Число итераций для метода Зейделя: " << iter1 << endl;

        int iter2 = 0;
        vector<double> X2 = relacs_method(A,B,B,1e-4,0.5,n_out,iter2);

        n_out << endl << "Результат метода релаксации при w = 0.5" << endl;
        tVectorOut(X2,n_out);
        n_out << endl << "Число итераций для метода релаксации при w = 0.5: " << iter2 << endl;

        matrix<double> A1(A);
        A1[0][3] = 0;
        A1[3][0] = 0;
        A1[0][2] = 0;
        A1[2][0] = 0;
        A1[1][3] = 0;
        A1[3][1] = 0;

        n_out << endl << "Матрица А1 (трёхдиагональная)" << endl;
        A1.output_f(n_out);

        int iter3 = 0;
        vector<double> X3 = yacoby_method(A1,B,B,1e-4,crash, iter3);

        n_out << endl << "Результат метода Якоби при A1 вместо A" << endl;
        tVectorOut(X3,n_out);

        vector<double> X4 = three_diagonal_relacs_method(A1,B,B,1e-4,1);

        n_out << endl << "Результат специального метода релаксации при A1" << endl;
        tVectorOut(X4,n_out);

        n_out << endl << "II. Для точности eps = 1e-7" << endl;

        int iter_s1 = 0;
        vector<double> Xs1 = simple_iteraition(A,B,B,1e-7,n_out,iter_s1);

        n_out << endl << "Результат метода простой итерации" << endl;
        tVectorOut(Xs1,n_out);
        n_out << endl << "Число итераций для метода простой итерации: " << iter_s1 << endl;

        int iter01 = 0;
        vector<double> X01 = yacoby_method(A,B,B,1e-7, n_out, iter01);

        n_out << endl << "Результат метода Якоби" << endl;
        tVectorOut(X01,n_out);
        n_out << endl << "Число итераций для метода Якоби: " << iter01 << endl;

        int iter11 = 0;
        vector<double> X11 = zeidel_method(A,B,B,1e-7,n_out,iter11);

        n_out << endl << "Результат метода Зейделя" << endl;
        tVectorOut(X11,n_out);
        n_out << endl << "Число итераций для метода Зейделя: " << iter11 << endl;

        int iter21 = 0;
        vector<double> X21 = relacs_method(A,B,B,1e-7,0.5,n_out,iter21);

        n_out << endl << "Результат метода релаксации при w = 0.5" << endl;
        tVectorOut(X21,n_out);
        n_out << endl << "Число итераций для метода релаксации при w = 0.5: " << iter21 << endl;


        n_out << endl << "Матрица А1 (трёхдиагональная)" << endl;
        A1.output_f(n_out);

        int iter31 = 0;
        vector<double> X31 = yacoby_method(A1,B,B,1e-7,crash, iter31);

        n_out << endl << "Результат метода Якоби при A1 вместо A" << endl;
        tVectorOut(X31,n_out);

        vector<double> X41 = three_diagonal_relacs_method(A1,B,B,1e-7,1);

        n_out << endl << "Результат специального метода релаксации при A1" << endl;
        tVectorOut(X41,n_out);
    }
    catch(char const* &s)
    //Обработка ошибок выводит в файл название ошибки
    {
        n_out << endl << "Ошибка: " << s;
    }
}

int main() {
    lab2(link_data_in"matrix26_2.dat",link_data_in"vector26_2.dat",link_data_out"Chist.txt");
    lab2(link_data_in"matrix13_2.dat",link_data_in"vector13_2.dat",link_data_out"Kis.txt");
    return 0;
}
