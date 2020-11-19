#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include </usr/include/mpi/mpi.h>

#define A1 -2.0
#define A2 3.0
#define B1 -1.0
#define B2 4.0

#define N 20
#define M 20

#define eps 0.00001


double H1 = (double) (A2 - A1) / M;
double H2 = (double) (B2 - B1) / N;

double weight_1(int i, const int m) {
    if ((i >= 1) && (i <= m - 1)) {
        return 1.0;
    } else if ((i == m) || (i == 0)) {
        return 0.5;
    } else {
        printf("ERROR: weight1 %d\n", i);
    }
    return 0;
}

double weight_2(int j, const int n) {
    if ((j >= 1) && (j <= n - 1)) {
        return 1.0;
    } else if ((j == n) || (j == 0)) {
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


double u_function(double x, double y) {
    return 2 / (1 + x*x + y*y);
}


double dot_product(double **u, double **v, const int m, const int n) {
    double dot_value = 0;
    int i, j;
    for (i = 0; i <= m; ++i) {
        double sum = 0;
        for (j = 0; j <= n; ++j) {
            sum += H2 * weight_1(i, m) * weight_2(j, n) * u[i][j] * v[i][j];
        }
        dot_value += H1 * sum;
    }
    return dot_value;
}


double norm(double **u, const int m, const int n) {

    return sqrt(dot_product(u, u, m, n));
}


double ** matrix_difference(double **a, double **b, double **c, int row, int col) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    return c;
}


double** fill_first_right_border_condition(double **w, double **pertrubed_f, const double *y, const int m, const int n) {
    int i, j;
    for (j = 0; j <= n; ++j) {
        pertrubed_f[m][j] = w[m][j];
    }
    return pertrubed_f;
}


double** fill_first_left_border_condition(double **w, double **pertrubed_f, const double *x, const int n) {
    int i, j;
    for (j = 0; j <= n; ++j) {
        pertrubed_f[0][j] = w[0][j];
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_top_condition(double **w, double **pertrubed_f,
                                  const double *x, const double *y,
                                  const double q, const int m, const int n) {
    int i, j;
    for (i = 1; i <= m - 1; ++i) {
        pertrubed_f[i][n] = w[i][n] *( 2 / (H2 * H2) * k_function(x[i], y[n]-0.5 * H2) + (q + 2 / H2) +
                (k_function(x[i] + 0.5 * H1, y[n]) + k_function(x[i] - 0.5 * H1, y[n])) / (H1 * H1)) +
                w[i][n-1] * (-2 / (H2 * H2) * k_function(x[i], y[n]-0.5 * H2)) +
                w[i+1][n] * (-k_function(x[i] + 0.5 * H1, y[n])/ (H1 * H1)) +
                w[i-1][n] * (-k_function(x[i] - 0.5 * H1, y[n])/ (H1 * H1));
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_bottom_condition(double **w, double **pertrubed_f,
                                     const double *x, const double *y,
                                     const double q, const int m, const int n) {
    int i, j;
    for (i = 1; i <= m -1; ++i) {
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
                            const double q, const int m, const int n) {
    int i, j;
    for (i = 1; i <= m - 1; ++i) {
        for (j = 2; j <= n - 1; ++j) {
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


double** fill_f(double **f, const double *x, const double *y, const int m, const int n) {
    int i, j;
    for (i = 1; i <= m - 1; ++i) {
        f[i][n] = f_function(x[i], y[n]) + 2 / H2 * psi_top(x[i]);
        f[i][1] = f_function(x[i], y[0]) + 2 / H2 * psi_bottom(x[i]);
    }

    for (j = 0; j <= n; ++j) {
        f[0][j] =  2.0 / (5.0 + y[j] * y[j]);
        f[n][j] =  2.0 / (10.0 + y[j] * y[j]);
    }

    for (i = 1; i <= m - 1; ++i) {
        for (j = 2; j <= n - 1; ++j) {
            f[i][j] = f_function(x[i], y[j]);
        }
    }
    return f;
}


double** apply_operator(double **w, double **pertrubed_f,
                          const double *x, const double *y,
                          const double q, const int m, const int n) {
    pertrubed_f = fill_first_left_border_condition(w, pertrubed_f, x, n);
    pertrubed_f = fill_first_right_border_condition(w, pertrubed_f, y, m, n);
    pertrubed_f = fill_third_bottom_condition(w, pertrubed_f, x, y, q,m , n);
    pertrubed_f = fill_third_top_condition(w, pertrubed_f, x, y, q, m , n);
    pertrubed_f = fill_operator_part(w, pertrubed_f, x, y, q, m ,n);
    return pertrubed_f;
}


int main(int argc, char *argv[]) {

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm grid_comm;

    int dimensions[2];
    int wrap_around[2];
    int reorder = 1;
    int coordinates[2];
    int my_grid_rank;
    int coords[2], Pi, Pj;

    int M_local, N_local, i0, i1, j0, j1;
    M_local = (M) / size -1;
    N_local = (N ) / size -1;


    dimensions[0] = 2;
    dimensions[1] = 2;

    wrap_around[0] = 1, wrap_around[1] = 1;
    printf("%d %d\n", N_local, M_local);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, reorder, &grid_comm);




    MPI_Comm_rank(grid_comm, &my_grid_rank);
    MPI_Cart_coords(grid_comm, my_grid_rank, 2, coordinates);
    Pi = coords[0];
    Pj = coords[1];

    i0 = (M + 1) / size * Pi;
    j0 = (N + 1) / size * Pj;
    i1 = i0 + M_local;
    j1 = j0 + N_local;

    printf("pi=%d pj=%d size=%d rank=%d\n", Pi, Pj, size, rank);


    double q = 1.0, current_norm;
    int k = 0, i, j, out_ind;

    double **b = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i < M_local + 1; ++i) {
        b[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **c = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i < M_local + 1; ++i) {
        c[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **r = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i < M_local + 1; ++i) {
        r[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **pertrubed_f = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i < M_local + 1; ++i) {
        pertrubed_f[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double** ar = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i < M_local + 1; ++i) {
        ar[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double *x = (double *) malloc ((M_local + 1) * sizeof(double));
    double *y = (double *) malloc ((N_local + 1) * sizeof(double));


    for (i = 0; i <= M_local; ++i) {
        x[i] = A1 + i * H1;
    }


    for (i = 0; i <= N_local; ++i) {
        y[i] = B1 + i * H2;
    }


    double ***w = (double ***) malloc(2 * sizeof(double **));
    for (i = 0; i < 2; ++i){
        w[i] = (double **) malloc((M_local + 1) * sizeof(double *));
        for (j = 0; j < M_local + 1; ++j) {
            w[i][j] = (double *) malloc((N_local + 1) * sizeof(double));
        }
    }

    b = fill_f(b, x, y, M_local, N_local);

    for (out_ind = 0; out_ind < 2; ++out_ind) {
        for (i = 0; i < M_local + 1; ++i) {
            for (j = 0; j < N_local + 1; ++j) {
                w[out_ind][i][j] = 0;
            }
        }
    }

    for (i = 0; i < M_local + 1; ++i) {
        for (j = 0; j < N_local + 1; ++j) {
            w[0][i][j] = u_function(x[i], y[j]);
        }
    }



    do {
        pertrubed_f = apply_operator(w[0], pertrubed_f, x, y, q, M_local, N_local);

        r = matrix_difference(pertrubed_f, b, c, M_local + 1, N_local + 1);

        ar = apply_operator(r, ar, x, y, q, M_local, N_local);

        double tmp = norm(ar, M_local, N_local);
        double t = dot_product(ar, r, M_local, N_local) / (tmp * tmp);

        for (i = 0; i < M_local + 1; ++i) {
            for (j = 0; j < N_local + 1; ++j) {
                w[1][i][j] = w[0][i][j] - t * r[i][j];
            }
        }

        current_norm = norm(matrix_difference(w[1], w[0], c, M_local + 1, N_local + 1), M_local, N_local);

 //      swap values for minimal memory usages
        for (i = 0; i < M_local + 1; ++i) {
            for (j = 0; j < N_local + 1; ++j) {
                w[0][i][j] = w[1][i][j];
            }
        }
        ++k;
        if (k) {
            printf("process= %d iteration = %d: norm = %f tmp=%f  t=%f\n", rank, k, current_norm, tmp, t);
        }
    } while (current_norm > eps);

    MPI_Finalize();
    printf("%f %d\n", current_norm, k);


    return 0;
}
