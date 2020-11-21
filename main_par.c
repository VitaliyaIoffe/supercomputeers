#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include </usr/include/mpi/mpi.h>

#define A1 (-2.0)
#define A2 3.0
#define B1 (-1.0)
#define B2 4.0

#define N 140
#define M 140

#define eps 0.00001


double H1 = (double) (A2 - A1) / M;
double H2 = (double) (B2 - B1) / N;


double weight_1(int i) {
    if ((i >= 1) && (i <= M - 1)) {
        return 1.0;
    } else if ((i == M) || (i == 0)) {
        return 0.5;
    } else {
        printf("ERROR: weight1 %d\n", i);
        exit(-1);
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
        exit(-1);
    }
    return 0;
}


double k_function(double x, double y) {
    return 1 + (x + y) * (x + y);
}


double f_function(double x, double y) {
    double tmp = k_function(x, y);
    return tmp + 2 * tmp * tmp + 4 * x * y * tmp * tmp * tmp;
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


double dot_product(double **u, double **v,
                   const int m, const int n,
                   const int i0, const int i1,
                   const int j0, const int j1,
                   const int size, const int rank) {
    double dot_value = 0, global_dot_value = 0;
    int i, j;

    for (i = i0; i <= i1; ++i) {
        double sum = 0;

//        todo investigate whyyy
//        if (i==10)
//            v[i-i0][0] = 0;
        for (j = j0; j < j0 + n; ++j) {
            sum += H2 * weight_1(i) * weight_2(j) * u[i - i0][j - j0] * v[i - i0][j - j0];


        }
        dot_value += H1 * sum;
    }
    MPI_Allreduce(&dot_value, &global_dot_value, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);


    return global_dot_value;
}


double norm(double **u,
            const int m, const int n,
            const int i0, const int i1,
            const int j0, const int j1,
            const int size, const int rank) {

    double norm = dot_product(u, u, m, n, i0, i1, j0, j1, size, rank);
    return sqrt(norm);
}


double** matrix_difference(double **a, double **b, double **c,
                           int row, int col) {
    int i, j;
    for (i = 0; i <= row; i++) {
        for (j = 0; j <= col; j++) {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    return c;
}


double** fill_first_right_border_condition(double **w, double **pertrubed_f, const double *y,
                                           const int m, const int n,
                                           const int i0, const int j0,const int i1, const int j1) {
    int i, j;
    if (i1 != M) {
        return pertrubed_f;
    }
    for (j = j0; j <= j1; ++j) {
        pertrubed_f[i1-i0][j-j0] = w[i1-i0][j-j0];
    }
    return pertrubed_f;
}


double** fill_first_left_border_condition(double **w, double **pertrubed_f, const double *x,
                                          const int m, const int n,const int i0, const int j0, const int i1, const int j1) {
    int i, j;
    if (i0 != 0) {
        return pertrubed_f;
    }
    for (j = j0; j <= j1; ++j) {
        pertrubed_f[0][j-j0] = w[0][j-j0];
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_top_condition(double **w, double **pertrubed_f,
                                  const double *x, const double *y,
                                  const double q, const int m, const int n,
                                  const int i0, const int j0, const int i1, const int j1,
                                  const double *left_rcv_buf, const  double *right_rcv_buf,
                                  const double *top_rcv_buf, const double *bottom_rcv_buf) {
    int i, j;
    double w_i_next, w_i_prev;
    if (j1 != N) {
        return pertrubed_f;
    }
    double w_next, w_prev;
    for (i = i0; i <= i1; ++i) {
        if (!(i0 > 1 && i0 <= M - 1)) {
            continue;
        }

        if (i + 1 > i1) {
            w_i_next = top_rcv_buf[j1] * (-k_function(x[i-i0] + 0.5 * H1, y[j1])/ (H1 * H1));
        } else {
            w_i_next = w[i+1-i0][j1] * (-k_function(x[i-i0] + 0.5 * H1, y[j1])/ (H1 * H1));
        }

        if (i - 1 < i0) {
            w_i_prev = bottom_rcv_buf[j1] * (-k_function(x[i-i0] - 0.5 * H1, y[j1])/ (H1 * H1));
        } else {
            w_i_prev = w[i-i0][j1] * (-k_function(x[i-i0] - 0.5 * H1, y[j1])/ (H1 * H1));
        }


        pertrubed_f[i-i0][j1] = w[i-i0][j1] *( 2 / (H2 * H2) * k_function(x[i-i0], y[j1]-0.5 * H2) + (q + 2 / H2) +
                (k_function(x[i-i0] + 0.5 * H1, y[j1]) + k_function(x[i-i0] - 0.5 * H1, y[j1])) / (H1 * H1)) +
                w[i-i0][j1-1] * (-2 / (H2 * H2) * k_function(x[i-i0], y[j1]-0.5 * H2)) +
                w_i_next + w_i_prev;
    }
    return pertrubed_f;
}

//  q === 1
double** fill_third_bottom_condition(double **w, double **pertrubed_f,
                                     const double *x, const double *y,
                                     const double q, const int m, const int n,
                                     const int i0, const int j0,
                                     const int i1, const int j1,
                                     double *left_rcv_buf, double *right_rcv_buf,
                                     double *top_rcv_buf, double *bottom_rcv_buf) {
    int i, j;
    double w_i_next, w_i_prev;
    if (j0 != 0) {
        return pertrubed_f;
    }
    for (i = i0; i <= i1; ++i) {
        if (!(i >= 1 && i <= M - 1)) {
            continue;
        }


        if (i + 1 > i1) {
            w_i_next = top_rcv_buf[0] * (-1 / (H1 * H1) * (k_function(x[i-i0] + 0.5 * H1, y[0])));
        } else {
            w_i_next = w[i+1-i0][0] * (k_function(x[i-i0] + 0.5 * H1, y[0])/ (H1 * H1));
        }

        if (i - 1 < i0) {
            w_i_prev = bottom_rcv_buf[0] * (-1 / (H1 * H1) * (k_function(x[i-i0] - 0.5 * H1, y[0])));
        } else {
            w_i_prev = w[i-1-i0][0] * (k_function(x[i-i0] - 0.5 * H1, y[0])/ (H1 * H1));
        }


        pertrubed_f[i-i0][0] = w[i-i0][0] * (2 / (H2 * H2) * k_function(x[i-i0], y[1]-0.5 * H2) + (q + 2 / H2) +
                (k_function(x[i-i0] + 0.5 * H1, y[0]) + k_function(x[i-i0] - 0.5 * H1, y[0])) / (H1 * H1)) +
                w[i-i0][1] * ( -2 / (H2 * H2) *k_function(x[i-i0], y[1]-0.5 * H2)) +
                w_i_next +
                w_i_prev;
    }
    return  pertrubed_f;
}

//  q === 1
double** fill_operator_part(double **w, double **pertrubed_f,
                            const double *x, const double *y,
                            const double q, const int m, const int n,
                            int i0, int j0, int i1, int j1,
                            const double *left_rcv_buf, const double *right_rcv_buf,
                            const double *top_rcv_buf, const double *bottom_rcv_buf) {
    int i, j;

    double w_i_next, w_j_next, w_i_prev, w_j_prev;
    for (i = i0; i <= i1; ++i) {
        if (!(i >= 1 && i <= M - 1)) {
            continue;
        }
        for (j = j0; j < j1; ++j) {
            if (!(j >= 1 && j <= N - 1)) {
                continue;
            }

            if (i + 1 >= i1) {
                w_i_next = top_rcv_buf[i-i0] * (-1 / (H1 * H1) * (k_function(x[i-i0] + 0.5 * H1, y[j-j0])));
            } else {
                w_i_next = w[i+1-i0][j-j0] * (-1 / (H1 * H1) * (k_function(x[i-i0] + 0.5 * H1, y[j-j0])));
            }

            if (i - 1 < i0) {
                w_i_prev = bottom_rcv_buf[i - i0] * (-1 / (H1 * H1) * (k_function(x[i-i0] - 0.5 * H1, y[j-j0])));
            } else {
                w_i_prev = w[i-1-i0][j-j0] * (-1 / (H1 * H1) * (k_function(x[i-i0] - 0.5 * H1, y[j-j0])));
            }

            if (j + 1 >= j1) {
                w_j_next = right_rcv_buf[j - j0] * (-1 / (H2 * H2) * (k_function(x[i-i0], y[j-j0] + 0.5 * H2)));
            } else {
                w_j_next = w[i-i0][j+1-j0] * (-1 / (H2 * H2) * (k_function(x[i-i0], y[j-j0] + 0.5 * H2)));
            }

            if (j - 1 < j0) {
                w_j_prev = left_rcv_buf[j - j0] * (-1 / (H2 * H2) * (k_function(x[i-i0], y[j-j0] - 0.5 * H2)));
            } else {
                w_j_prev = w[i-i0][j-1-j0] * (-1 / (H2 * H2) * (k_function(x[i-i0], y[j-j0] - 0.5 * H2)));
            }

            pertrubed_f[i-i0][j-j0] = (
                    w[i-i0][j-j0] * (q + 1 / (H1 * H1) * (
                            k_function(x[i-i0] + 0.5 * H1, y[j-j0]) + k_function(x[i-i0] - 0.5 * H1, y[j-j0])) +
                            1 / (H2 * H2) * (
                                    k_function(x[i-i0], y[j-j0] + 0.5 * H2) + k_function(x[i-i0], y[j-j0] - 0.5 * H2)))) +
                                            w_i_next + w_i_prev +
                                            w_j_next + w_j_prev;
            if (i == 1 && j == 10 && isnan(w_i_next)) {
                printf("ALARM %d %d  %f, %f\n", i, j, top_rcv_buf[i-i0], w[i+1-i0][j-j0]);
            }
        }
    }

    return pertrubed_f;
}


double** fill_f(double **f, const double *x, const double *y,
                const int m, const int n,
                const int i0, const int i1,
                const int j0, const int j1) {
    int i, j;
    for (i = i0; i <= i1; ++i) {
        if (i + i0 == 0 || i + i0 > M - 1) {
            continue;
        }
        if (j1 == N){
            f[i-i0][j1-j0] = f_function(x[i-i0], y[j1-j0]) + 2 / H2 * psi_top(x[i-i0]);
        }
        if (j0 == 0) {
            f[i-i0][1] = f_function(x[i-i0], y[0]) + 2 / H2 * psi_bottom(x[i-i0]);
        }
    }

    for (j = j0; j <= j1; ++j) {
        if (i0 == 0) {
            f[0][j-j0] =  2.0 / (5.0 + y[j-j0] * y[j-j0]);
        }
        if (i1 == M) {
            f[i1-i0][j-j0] =  2.0 / (10.0 + y[j-j0] * y[j-j0]);
        }
    }

    for (i = i0; i <= i1; ++i) {
       if (i == 0 || i > M - 1) {
            continue;
        }
        for (j = j0; j <= j1; ++j) {
            if (j < 1 || j > N - 1) {
                continue;
            }
            f[i-i0][j-j0] = f_function(x[i-i0], y[j-j0]);
        }
    }

    return f;
}


double** apply_operator(double **w, double **pertrubed_f,
                        const double *x, const double *y,
                        const double q, const int m, const int n,
                        const int i0, const int j0, const int i1, const int j1,
                        const int *rank,  MPI_Comm comm, int rank_number) {


    int left_rank = rank[0], right_rank = rank[1], top_rank = rank[2], bottom_rank = rank[3];

    double right_snd_buf[n + 1], right_rcv_buf[n + 1];
    double left_snd_buf[n + 1], left_rcv_buf[n + 1];
    double top_snd_buf[m + 1], top_rcv_buf[m + 1];
    double bottom_snd_buf[m + 1], bottom_rcv_buf[m + 1];


    MPI_Request bottom_snd_req, bottom_rcv_req;
    MPI_Status bottom_snd_stat, bottom_rcv_stat;

    MPI_Request right_snd_req, right_rcv_req;
    MPI_Status right_snd_stat, right_rcv_stat;

    MPI_Request top_snd_req, top_rcv_req;
    MPI_Status top_snd_stat, top_rcv_stat;

    MPI_Request left_snd_req, left_rcv_req;
    MPI_Status left_snd_stat, left_rcv_stat;

    int tag = 1;
    int j, i;

    for (j = j0; j <= j1; ++j) {
        left_snd_buf[j-j0] = w[0][j-j0];
        right_snd_buf[j-j0] = w[i1-i0][j-j0];
    }

    for (i = i0; i <= i1; ++i) {
        top_snd_buf[i-i0] = w[i-i0][j1-j0];
        bottom_snd_buf[i-i0] = w[i-i0][0];
    }

    if (left_rank != -1) {
        MPI_Isend(&left_snd_buf[0], n+1, MPI_DOUBLE, left_rank, tag, comm, &left_snd_req);
        MPI_Irecv(&left_rcv_buf[0], n+1, MPI_DOUBLE, left_rank, tag, comm, &left_rcv_req);
    }
//
    if (right_rank != -1) {
        MPI_Isend(&right_snd_buf[0], n+1, MPI_DOUBLE, right_rank, tag, comm, &right_snd_req);
        MPI_Irecv(&right_rcv_buf[0], n+1, MPI_DOUBLE, right_rank, tag, comm, &right_rcv_req);
    }

    if (top_rank != -1) {
        MPI_Isend(&top_snd_buf[0], m+1, MPI_DOUBLE, top_rank, tag, comm, &top_snd_req);
        MPI_Irecv(&top_rcv_buf[0], m+1, MPI_DOUBLE, top_rank, tag, comm, &top_rcv_req);
    }

    if (bottom_rank != -1) {
        MPI_Isend(&bottom_snd_buf[0], m+1, MPI_DOUBLE, bottom_rank, tag, comm, &bottom_snd_req);
        MPI_Irecv(&bottom_rcv_buf[0], m+1, MPI_DOUBLE, bottom_rank, tag, comm, &bottom_rcv_req);
    }

    if (left_rank != -1)
    {
        MPI_Wait(&left_rcv_req, &left_rcv_stat);
        MPI_Wait(&left_snd_req, &left_snd_stat);
   }

    if (right_rank != -1) {
        MPI_Wait(&right_rcv_req, &right_rcv_stat);
        MPI_Wait(&right_snd_req, &right_snd_stat);
    }
    if (top_rank != -1) {
        MPI_Wait(&top_rcv_req, &top_rcv_stat);
        MPI_Wait(&top_snd_req, &top_snd_stat);
    }

    if (bottom_rank != -1) {
        MPI_Wait(&bottom_rcv_req, &bottom_rcv_stat);
        MPI_Wait(&bottom_snd_req, &bottom_snd_stat);
    }

    if (i0 == 0) {
        pertrubed_f = fill_first_left_border_condition(w, pertrubed_f, x, m+1, n+1, i0, j0, i1, j1);
    }

    if (i1 == M) {
        pertrubed_f = fill_first_right_border_condition(w, pertrubed_f, y, m+1, n+1, i0, j0, i1, j1);
    }

    if (j0 != 0) {
        pertrubed_f = fill_third_bottom_condition(w, pertrubed_f, x, y, q, m+1, n+1, i0, j0, i1, j1,
                                                  left_rcv_buf, right_rcv_buf,
                                                  right_rcv_buf, bottom_rcv_buf);
    }

    if (j1 != N) {
        pertrubed_f = fill_third_top_condition(w, pertrubed_f, x, y, q, m+1, n+1, i0, j0, i1, j1,
                                               left_rcv_buf, right_rcv_buf,
                                               top_rcv_buf, bottom_rcv_buf);
    }

    pertrubed_f = fill_operator_part(w, pertrubed_f, x, y, q, m+1, n+1, i0, j0, i1, j1,
                                     left_rcv_buf, right_rcv_buf,
                                     top_rcv_buf, bottom_rcv_buf);
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


    int n_procs = 10;
    dimensions[0] = atoi(argv[argc-2]);
    dimensions[1] = atoi(argv[argc-1]);

    if (n_procs != dimensions[0] * dimensions[1]) {
        printf("ERROR: dim[0] = %d, dim[1] = %d, procs = %d\n", dimensions[0], dimensions[1], n_procs);
        exit(-1);
    }


    M_local = (M) / dimensions[0] - 1;
    N_local = (N) / dimensions[1] - 1;

    wrap_around[0] = 1, wrap_around[1] = 1;
    printf("%d %d\n", N_local, M_local);

    MPI_Cart_create(MPI_COMM_WORLD, 2 , dimensions, wrap_around, reorder, &grid_comm);

    MPI_Cart_coords(grid_comm, rank, 2, coords);
    Pi = coords[0];
    Pj = coords[1];


    i0 = (M) / dimensions[0] * Pi;
    j0 = (N) / dimensions[1] * Pj;
    i1 = i0 + M_local + 1;
    j1 = j0 + N_local + 1;

    int ranks[4], left_rank = -1, right_rank = -1, top_rank = -1, bottom_rank = -1;
    if (Pj > 0) {
        coords[0] = Pi;
        coords[1] = Pj - 1;
        MPI_Cart_rank(grid_comm, coords, &bottom_rank);
    }

    if (Pi > 0) {
        coords[0] = Pi - 1;
        coords[1] = Pj;
        MPI_Cart_rank(grid_comm, coords, &left_rank);
    }

    if (Pj + 1 < dimensions[1]) {
        coords[0] = Pi;
        coords[1] = Pj + 1;
        MPI_Cart_rank(grid_comm, coords, &top_rank);
    }

    if (Pi + 1 < dimensions[0]) {
        coords[0] = Pi + 1;
        coords[1] = Pj;
        MPI_Cart_rank(grid_comm, coords, &right_rank);
    }

    ranks[0] = left_rank;
    ranks[1] = right_rank;
    ranks[2] = top_rank;
    ranks[3] = bottom_rank;
    printf("pi=%d pj=%d i0=%d i1=%d j0=%d j1=%d size=%d rank=%d\n", Pi, Pj, i0, i1, j0, j1, size, rank);
//    printf("pi=%d pj=%d left=%d right=%d top=%d bottom=%d size=%d rank=%d\n", Pi, Pj, left_rank, right_rank, top_rank, bottom_rank, size, rank);


    double q = 1.0, current_norm;
    int k = 0, i, j, out_ind;

    double **b = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i <= M_local + 1; ++i) {
        b[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **c = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i <= M_local + 1; ++i) {
        c[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **r = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i <= M_local + 1; ++i) {
        r[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double **pertrubed_f = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i <= M_local + 1; ++i) {
        pertrubed_f[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double** ar = (double **) malloc((M_local + 1) * sizeof(double *));
    for (i = 0; i <= M_local + 1; ++i) {
        ar[i] = (double *) malloc((N_local + 1) * sizeof(double));
    }

    double *x = (double *) malloc ((M_local + 1) * sizeof(double));
    double *y = (double *) malloc ((N_local + 1) * sizeof(double));


    for (i = i0; i <= i1; ++i) {
        x[i-i0] = A1 + i * H1;
    }


    for (j = j0; j <= j1; ++j) {
        y[j-j0] = B1 + j * H2;
    }


    double ***w = (double ***) malloc(2 * sizeof(double **));
    for (i = 0; i < 2; ++i){
        w[i] = (double **) malloc((M_local + 1) * sizeof(double *));
        for (j = 0; j <= M_local + 1; ++j) {
            w[i][j] = (double *) malloc((N_local + 1) * sizeof(double));
        }
    }


    b = fill_f(b, x, y,
               M_local, N_local,
               i0, i1, j0, j1);


    for (out_ind = 0; out_ind < 2; ++out_ind) {
        for (i = i0; i <= i1; ++i) {
            for (j = j0; j <= j1; ++j) {
                w[out_ind][i-i0][j-j0] = 0;
            }
        }
    }

    for (i = i0; i <= i1; ++i) {
        for (j = j0; j <= j1; ++j) {
            w[0][i-i0][j-j0] = u_function(x[i-i0], y[j-j0]);
            c[i-i0][j-j0] = 0;
            ar[i-i0][j-j0] = 0;
        }
    }


    do {
        pertrubed_f = apply_operator(w[0], pertrubed_f, x, y, q, M_local, N_local, i0, j0, i1, j1, ranks, grid_comm, rank);

        r = matrix_difference(pertrubed_f, b, r, M_local+1 , N_local+1);

        ar = apply_operator(r, ar, x, y, q, M_local, N_local, i0, j0, i1, j1, ranks, grid_comm, rank);

        double tmp = norm(ar, M_local+1, N_local+1, i0, i1, j0, j1, size, rank);

        double t = dot_product(ar, r, M_local+1, N_local+1, i0, i1, j0, j1, size, rank) / (tmp * tmp);
//        printf("t = %f %d %d %d %d tmp %f\n", dot_product(ar, r, M_local+1, N_local+1, i0, i1, j0, j1, size, rank), i0, i1, j0, j1, tmp);
        for (i = i0; i <= i1; ++i) {
            for (j = j0; j <= j1; ++j) {
                w[1][i-i0][j-j0] = w[0][i-i0][j-j0] - t * r[i-i0][j-j0];
            }
        }

        c = matrix_difference(w[1], w[0], c, M_local + 1, N_local + 1);

        current_norm = norm(c,M_local+1, N_local+1,
                            i0, i1, j0, j1,
                            size, rank);

//        //      swap values for minimal memory usages
        for (i = i0; i <= i1; ++i) {
            for (j = j0; j <= j1; ++j) {
                w[0][i-i0][j-j0] = w[1][i-i0][j-j0];
            }
        }
        if (k % 20 == 0) {
            printf("process = %d, iteration = %d: norm = %f, tmp = %f,  t = %f\n", rank, k, current_norm, tmp, t);
        }
        ++k;
//    } while (0);
    } while (current_norm > eps);

    free(x);
    free(y);


    for (i = 0; i < M_local + 1; ++i) {
        free(b[i]);
        free(c[i]);
        free(ar[i]);
        free(r[i]);
        free(pertrubed_f[i]);
    }
    free(b);
    free(c);
    free(ar);
    free(r);
    free(pertrubed_f);


    for (i = 0; i < 2; ++i){
        for (j = 0; j < M_local + 1; ++j) {
            free(w[i][j]);
        }
        free(w[i]);
    }
    free(w);
//    fclose(myfile);

    MPI_Finalize();
    printf("%f %d\n", current_norm, k);


    return 0;
}
