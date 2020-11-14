#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A1 -2.0
#define A2 3.0
#define B1 -1.0
#define B2 4.0

#define N 10
#define M 10

#define eps 0.000001

#define COUNT 20000


// todo: use realloc or optimize memory usage (priority 1)

double h1(int m) {
    return (double) (A2 - A1) / m;
}

double h2(int n) {
    return (double) (B2 - B1) / n;
}

double H1 = (double) (A2 - A1) / M;
double H2 = (double) (B2 - B1) / N;

double weight_1(int i) {
    if ((i >= 1) && (i <= M - 1)) {
        return 1.0;
    } else if ((i == M) || (i == 0)) {
        return 0.5;
    } else {
        printf("ERROR: weight1 %d\n", i);
    }
    return 0;
}

double weight_2(int j) {
    if ((j >= 1) && (j <= N - 1)) {
        return 1.0;
    } else if ((j == N) || (j == 0)) {
        return 0.5;
    } else {
        printf("ERROR: weight2 %d\n", j);
    }
    return 0;
}


double k_function(double x, double y) {
    return 1 + (x + y) * (x + y);
}


double f_function(double x, double y) {
    double tmp = k_function(x, y);
    return tmp + 2 * tmp * tmp + 4 * tmp * tmp * tmp;
}


double psi_top(double x) {
    double tmp = (x*x + 17);
    return -32 / (tmp * tmp * tmp) + 2 / tmp;
}


double psi_bottom(double x) {
    double tmp = (x*x + 2);
    return -8 / (tmp * tmp * tmp) + 2 / tmp;
}


double** f_fill(double **f, int row, int col, double *x, double *y) {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
//           todo: change with conditions;
            f[i][j] = f_function(x[i], y[j]);
        }
    }
    return f;
}


double** a_fill(double **a, double *x, double *y) {
    for (int i = 1; i <= M; ++i) {
        for (int j = 1; j <= N; ++j) {
            a[i][j] = k_function(x[i] - 0.5 * H1, y[j]);
        }
    }
    return a;
}


double** b_fill(double **b, double *x, double *y) {
    for (int i = 1; i <= M; ++i) {
        for (int j = 1; j <= N; ++j) {
            b[i][j] = k_function(x[i], y[j] - 0.5 * H2);
        }
    }
    return b;
}




double deriv_x_beyond(double **w, int i, int j) {
        return (w[i+1][j] - w[i][j]) / H1;
}

double deriv_x_past(double **w, int i, int j) {
    return (w[i][j] - w[i-1][j]) / H1;
}

double deriv_y_beyond(double **w, int i, int j) {
    return (w[i][j+1] - w[i][j]) / H2;
}

double deriv_y_past(double **w, int i, int j) {
    return (w[i][j] - w[i][j-1]) / H2;
}


double u_function(double x, double y) {
    return 2 / (1 + x*x + y*y);
}

double dot_product(double **u, double **v) {
    double dot_value = 0;
    for (int i = 0; i <= M; ++i) {
        double sum = 0;
        for (int j = 0; j <= N; ++j) {
            sum += H2 * weight_1(i) * weight_2(j) * u[i][j] * v[i][j];
        }
        dot_value += H1 * sum;
    }
    return dot_value;
}


double norm(double **u) {
    double norm_squared = dot_product(u, u);
    return sqrt(norm_squared);
}

double** matrix_multiply(double **a, int row1, int col1, double **b, int row2, int col2) {
    printf("LOG: MULTIPLY MATRICES\n");
    if (col1 != row2) {
        printf("ERROR: wrong row and column in matrix multiplying\n");
        exit(-1);
    }
    double** c = (double **) malloc(row1  * sizeof(double*));
    for (int i = 0; i < row1; i++) {
        c[i] = (double*) malloc(col2 * sizeof(double));
        for (int j = 0; j < col2; j++) {
            c[i][j] = 0;
            for (int k = 0; k < col1; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

double ** matrix_difference(double **a, double **b, int row, int col) {
    printf("LOG: MATRIX DIFFERENCE\n");
    double** c = (double **) malloc(row  * sizeof(double*));
    for (int i = 0; i < row; i++) {
        c[i] = (double*) malloc(col * sizeof(double));
        for (int j = 0; j < col; j++) {
            c[i][j] = a[i][j] - b[i][j];

        }
    }
    return c;
}


double** fill_first_right_border_condition(double **w, double **pertrubed_f, const double *y) {
    for (int j = 0; j <= N; ++j) {
        pertrubed_f[M][j] = w[M][j];
//        w[M][j] += 2.0 / (10.0 + y[j] * y[j]);
    }
    return pertrubed_f;
}


double** fill_first_left_border_condition(double **w, double **pertrubed_f, const double *x) {
    for (int j = 0; j <= N; ++j) {
        pertrubed_f[0][j] = w[0][j];
//        w[0][j] += 2.0 / (5.0 + x[j] * x[j]);
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_top_condition(double **w, double **pertrubed_f,
                                  const double *x, const double *y,
                                  const double q) {
    for (int i = 1; i <= M - 1; ++i) {
        pertrubed_f[i][N] = w[i][N] *( 2 / (H2 * H2) * k_function(x[i], y[N]-0.5 * H2) +
                (q + 2 / H2) +
                (k_function(x[i] + 0.5 * H1, y[N]) + k_function(x[i] - 0.5 * H1, y[N])) / (H1 * H1)) +
                w[i][N-1] * (-2 / (H2 * H2) * k_function(x[i], y[N]-0.5 * H2)) +
                w[i+1][N] * (-k_function(x[i] + 0.5 * H1, y[N])/ (H1 * H1)) +
                w[i-1][N] * (-k_function(x[i] - 0.5 * H1, y[N])/ (H1 * H1));
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_bottom_condition(double **w, double **pertrubed_f,
                                     const double *x, const double *y,
                                     const double q) {
    for (int i = 1; i <= M -1; ++i) {
        pertrubed_f[i][0] = w[i][0] * (2 / (H2 * H2) * k_function(x[i], y[1]-0.5 * H2) + (q + 2 / H2) +
                (k_function(x[i] + 0.5 * H1, y[0]) + k_function(x[i] - 0.5 * H1, y[0])) / (H1 * H1)) +
                w[i][1] * ( -2 / (H2 * H2) *k_function(x[i], y[1]-0.5 * H2)) +
                w[i+1][0] * (k_function(x[i] + 0.5 * H1, y[0])/ (H1 * H1)) +
                w[i-1][0] * (k_function(x[i] - 0.5 * H1, y[0])/ (H1 * H1));
    }
    return  pertrubed_f;
}

//  q === 1
double** fill_operator_part(double **w, double **pertrubed_f,
                            const double *x, const double *y,
                            const double q) {
    for (int i = 1; i <= M - 1; ++i) {
        for (int j = 2; j <= N - 1; ++j) {
            pertrubed_f[i][j] = w[i][j] * (q + 1 / (H1 * H1) * (k_function(x[i] + 0.5 * H1, y[j]) + k_function(x[i] - 0.5 * H1, y[j])) +
                    1 / (H2 * H2) * (k_function(x[i], y[j] + 0.5 * H2) + k_function(x[i], y[j] - 0.5 * H2))) +
                    w[i+1][j] * (-1 / (H1 * H1) * (k_function(x[i] + 0.5 * H1, y[j]))) +
                    w[i-1][j] * (-1 / (H1 * H1) * (k_function(x[i] - 0.5 * H1, y[j]))) +
                    w[i][j+1] * (-1 / (H2 * H2) * (k_function(x[i], y[j] + 0.5 * H2))) +
                    w[i][j-1] * (-1 / (H2 * H2) * (k_function(x[i], y[j] - 0.5 * H2)));
        }
    }
    return pertrubed_f;
}


double** fill_f(double **f, const double *x, const double *y) {

    for (int i = 1; i <= M - 1; ++i) {
        f[i][N] = f_function(x[i], y[N]) + 2 / H2 * psi_top(x[i]);
        f[i][1] = f_function(x[i], y[0]) + 2 / H2 * psi_bottom(x[i]);
    }

    for (int j = 0; j <= N; ++j) {
        f[0][j] =  2.0 / (5.0 + y[j] * y[j]);
        f[N][j] =  2.0 / (10.0 + y[j] * y[j]);
    }

    for (int i = 1; i <= M - 1; ++i) {
        for (int j = 2; j <= N - 1; ++j) {
            f[i][j] = f_function(x[i], y[j]);
        }
    }
    return f;
}


int test_multiply_matrices(void) {
    double **a = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; ++i) {
        a[i] = (double *) malloc(2 * sizeof(double));
    }

    double **b = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; ++i) {
        b[i] = (double *) malloc(2 * sizeof(double));
    }
    int c[3][2];
    a[0][0] = 1; a[1][0] = -1; a[2][0] = 5;
    a[0][1] = 2; a[1][1] = 4; a[2][1] = -3;

    b[0][0] = 0; b[1][0] = -1;
    b[0][1] = 2; b[1][1] = 5;

    c[0][0] = -2; c[1][0] = -4; c[2][0] = 3;
    c[0][1] = 12; c[1][1] = 18; c[2][1] = -5;

    double** res = matrix_multiply(a, 3, 2, b, 2, 2);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (c[i][j] != res[i][j]) {
                printf("LOG: TEST desn't pass\n");

                for (int m = 0; m< 3; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        printf("m=%d, n=%d, result=%f expected=%d| ",m, n, res[m][n], c[m][n]);
                    }
                    printf("\n");
                }
                printf("\n\n");

                exit(-1);
            }
         }
    }
    printf("LOG: TEST passed\n");
    return 1;

}

int test_diff_of_matrices(void) {
    double **a = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; ++i) {
        a[i] = (double *) malloc(2 * sizeof(double));
    }

    double **c = (double **) malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; ++i) {
        c[i] = (double *) malloc(2 * sizeof(double));
    }
    a[0][0] = 1; a[1][0] = -1; a[2][0] = 5;
    a[0][1] = 2; a[1][1] = 4; a[2][1] = -3;

    c[0][0] = -2; c[1][0] = -4; c[2][0] = 3;
    c[0][1] = 12; c[1][1] = 18; c[2][1] = -5;

    double diff[3][2];
    diff[0][0] = 3; diff[1][0] = 3; diff[2][0] = 2;
    diff[0][1] = -10; diff[1][1] = -14; diff[2][1] = 2;

    double** res = matrix_difference(a, c, 3, 2);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (fabs(diff[i][j] - res[i][j]) > eps) {
                printf("LOG: TEST desn't pass\n");

                for (int m = 0; m< 3; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        printf("m=%d, n=%d, result=%f expected=%f| ",m, n, res[m][n], diff[m][n]);
                    }
                    printf("\n");
                }
                printf("\n\n");

                exit(-1);
            }
        }
    }
    printf("LOG: TEST passed\n");
    return 1;

}

//double ** fill_r() {
//
//
//
//}

double** apply_operator(double **w, double **pertrubed_f,
                          const double *x, const double *y,
                          const double q) {
    pertrubed_f = fill_first_left_border_condition(w, pertrubed_f, x);
    pertrubed_f = fill_first_right_border_condition(w, pertrubed_f, y);
    pertrubed_f = fill_third_bottom_condition(w, pertrubed_f, x, y, q);
    pertrubed_f = fill_third_top_condition(w, pertrubed_f, x, y, q);
    pertrubed_f = fill_operator_part(w, pertrubed_f, x, y, q);
    return pertrubed_f;
}

int main() {
    test_multiply_matrices();
    test_diff_of_matrices();
    double **a = (double **) malloc((M + 1) * sizeof(double *));
    for (int i = 0; i < M + 1; ++i) {
        a[i] = (double *) malloc((N + 1) * sizeof(double));
    }

    for (int i = 0; i < M+1; ++i){
        for (int j = 0; j < N + 1; ++j) {
            a[i][j] = 0.0;
        }
    }


    double **b = (double **) malloc((M + 1) * sizeof(double *));
    for (int i = 0; i < M + 1; ++i) {
        b[i] = (double *) malloc((N + 1) * sizeof(double));
    }

    double **r = (double **) malloc((M + 1) * sizeof(double *));
    for (int i = 0; i < M + 1; ++i) {
        r[i] = (double *) malloc((N + 1) * sizeof(double));
    }

    double **pertrubed_f = (double **) malloc((M + 1) * sizeof(double *));
    for (int i = 0; i < M + 1; ++i) {
        pertrubed_f[i] = (double *) malloc((N + 1) * sizeof(double));
    }

    double *x = (double *) malloc ((M + 1) * sizeof(double));
    double *y = (double *) malloc ((N + 1) * sizeof(double));


    for (int i = 0; i <= M; ++i) {
        x[i] = A1 + i * H1;
        printf("%f ", x[i]);
    }
    printf("\n");


    for (int i = 0; i <= N; ++i) {
        y[i] = B1 + i * H2;
        printf("%f ", y[i]);
    }
    printf("\n");

    double q = 1.0;

// FILL MATRIX B
    printf("LOG: FILLING B\n");
    b = fill_f(b, x, y);





    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            printf("i=%d, j=%d, result=%f | ",i, j, a[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");


    double ***w = (double ***) malloc(COUNT * sizeof(double **));
    for (int i = 0; i < COUNT; ++i){
        w[i] = (double **) malloc((M + 1) * sizeof(double *));
        for (int j = 0; j < M + 1; ++j) {
            w[i][j] = (double *) malloc((N + 1) * sizeof(double));
        }
    }
    int k = 0;
    for (int out_ind = 0; out_ind < COUNT; ++out_ind) {
        for (int i = 0; i < M + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                w[out_ind][i][j] = 0;
            }
        }
    }

    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            w[0][i][j] = u_function(x[i], y[j]);
        }
    }


    double** ar = (double **) malloc((M + 1) * sizeof(double *));
    for (int i = 0; i < M + 1; ++i) {
        ar[i] = (double *) malloc((N + 1) * sizeof(double));
    }

    double current_norm;
    do {
        printf("LOG: Start iteration #%d\n", k + 1);

        pertrubed_f = fill_first_left_border_condition(w[k], pertrubed_f, x);
        pertrubed_f = fill_first_right_border_condition(w[k], pertrubed_f, y);
        pertrubed_f = fill_third_bottom_condition(w[k], pertrubed_f, x, y, q);
        pertrubed_f = fill_third_top_condition(w[k], pertrubed_f, x, y, q);
        pertrubed_f = fill_operator_part(w[k], pertrubed_f, x, y, q);

        r = matrix_difference(pertrubed_f,
                              b, M + 1, N + 1);


        // check dimensions
        ar = apply_operator(r, ar, x, y, q);

        double tmp = norm(ar);
        double t = dot_product(ar, r) / (tmp * tmp);

        for (int i = 0; i < M + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {

//                printf("%f %f\n", t*r[i][j], w[k][i][j]);
                w[k+1][i][j] = w[k][i][j] - t * r[i][j];
            }
        }



        current_norm = norm(matrix_difference(w[k+1], w[k], M + 1, N + 1));

        ++k;
        printf("%f\n", current_norm);
    } while ( current_norm > eps);

    for (int i = 0; i < M ; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%f ", w[k-1][i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < M ; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%f ", b[i][j]);
        }
        printf("\n");
    }


    return 0;
}
