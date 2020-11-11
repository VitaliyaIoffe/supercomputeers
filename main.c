#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A1 -2.0
#define A2 3.0
#define B1 -1.0
#define B2 4.0

#define N 20
#define M 20

# define eps 0.000001

#define COUNT 2000

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
    } else if (i == M) {
        return 0.5;
    } else {
        printf("ERROR: weight1");
    }
    return 0;
}

double weight_2(int j) {
    if ((j >= 1) && (j <= N - 1)) {
        return 1.0;
    } else if (j == N) {
        return 0.5;
    } else {
        printf("ERROR: weight2");
    }
    return 0;
}


double k_function(double x, double y) {
    return 1 + (x + y) * (x + y);
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


double norm(double **w_next, double **w_prev) {
    double norm = 0;
    for (int i = 0; i <= M; ++i) {

        for (int j = 0; j <= N; ++j) {

        }
//        todo: fix norm
//        norm += ro(i, j) * (w_next  )
    }
    return norm;
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

double** matrix_multiply(double **a, int row1, int col1, double **b, int row2, int col2) {
    if (col1 != row2) {
        printf("ERROR: wrong row and column in matrix multiplying");
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
    double** c = (double **) malloc(row  * sizeof(double*));
    for (int i = 0; i < row; i++) {
        c[i] = (double*) malloc(col * sizeof(double));
        for (int j = 0; j < col; j++) {
            c[i][j] = a[i][j] - b[i][j];

        }
    }
    return c;
}

//double ** fill_r() {
//
//}

int main() {

    double **a = (double **) malloc((N + 1) * sizeof(double *));
    for (int i = 0; i < N + 1; ++i) {
        a[i] = (double *) malloc((M + 1) * sizeof(double));
    }

    double **b = (double **) malloc((N + 1) * sizeof(double *));
    for (int i = 0; i < N + 1; ++i) {
        b[i] = (double *) malloc((M + 1) * sizeof(double));
    }

    double **r = (double **) malloc((N + 1) * sizeof(double *));
    for (int i = 0; i < N + 1; ++i) {
        r[i] = (double *) malloc((M + 1) * sizeof(double));
    }

    double *x = (double *) malloc ((N + 1) * sizeof(double));
    double *y = (double *) malloc ((M + 1) * sizeof(double));

    for (int i = 0; i <= M; ++i) {
        x[i] = A1 + i * H1;
    }

    for (int i = 0; i <= N; ++i) {
        y[i] = B1 + i * H2;
    }

    a = a_fill(a, x, y);
    b = b_fill(b, x, y);

    double ***w = (double ***) malloc( COUNT * sizeof(double **));
    for (int i = 0; i < COUNT; ++i){
        w[i] = (double **) malloc((N + 1) * sizeof(double *));
        for (int j = 0; j < N + 1; ++j) {
            w[i][j] = (double *) malloc((M + 1) * sizeof(double));
        }
    }
    int k = 0;
    for (int i = 0; i < N+1; ++i) {
        for (int j = 0; j < M+1; ++j) {
            w[0][i][j] = 0;
        }
    }
    do {
        r = matrix_difference(matrix_multiply(a, M, N, w[k], M, N), b, M, N);

        // check dimensions
        double** ar = matrix_multiply(a, M, N, r, M, 1);
        double tmp = norm(ar, ar);
        double t = dot_product(ar, r) / (tmp * tmp);
        for (int i = 0; i < N+1; ++i) {
            for (int j = 0; j < M+1; ++j) {
                w[k+1][i][j] = w[k][i][j] - t * r[i][j];
            }
        }
        ++k;

    } while (norm(w[k+1], w[k]) > eps);

    return 0;
}
