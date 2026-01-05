#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

/* ------------------ mpicontrol (style giống slide) ------------------ */
int myrank, nprocs, ierror, length;
MPI_Status status;
char hostname[MPI_MAX_PROCESSOR_NAME];

/* ------------------ mxmdat (n,m,l + matrices) ------------------ */
/* a(n,m) x b(m,l) = c(n,l) */
int n, m, l;
int ilo, ihi;              /* chunk of columns [ilo..ihi] (1-based, as in slide) */
int chunk_cols;            /* = ihi - ilo + 1 */

double *a = NULL;          /* full A, size n*m (redundant on all ranks) */
double *b = NULL;          /* local B chunk, size m*chunk_cols */
double *c = NULL;          /* local C chunk, size n*chunk_cols */

/* --------------- helpers: index macros (row-major) ---------------- */
#define A(i,j) a[(i)*m + (j)]                       /* i:0..n-1, j:0..m-1 */
#define B(j,k) b[(j)*chunk_cols + (k)]              /* j:0..m-1, k:0..chunk_cols-1 */
#define C(i,k) c[(i)*chunk_cols + (k)]              /* i:0..n-1, k:0..chunk_cols-1 */

/* ------------------ ordered printing (slide style) ------------------ */
void print_ordered(const char *string) {
    /* Slide dùng barrier/loop. Cách này bám sát: barrier + loop in theo rank */
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < nprocs; i++) {
        if (i == myrank) {
            printf("%s\n", string);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/* ------------------ initialize_mpi (slide 0.1/0.2) ------------------ */
void initialize_mpi(void) {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Get_processor_name(hostname, &length);

    char string[256];
    snprintf(string, sizeof(string), "rank=%d nprocs=%d host=%s", myrank, nprocs, hostname);
    print_ordered(string);
}

/* ------------------ workshare (slide version 1.0) ------------------ */
/* (1:n) is partitioned into nprocs pieces. task <rank> gets (ilo:ihi) */
void workshare(int rank, int nprocs_, int N, int *ilo_out, int *ihi_out) {
    int nrem   = N % nprocs_;
    int nchunk = (N - nrem) / nprocs_;

    if (rank < nrem) {
        *ilo_out = 1 + rank * (nchunk + 1);
        *ihi_out = *ilo_out + nchunk;
    } else {
        *ilo_out = 1 + rank * nchunk + nrem;
        *ihi_out = *ilo_out + nchunk - 1;
    }
}

/* ------------------ mat_dimensions (master-only IO style) ----------- */
/* Slide đọc từ stdin ở master, rồi broadcast */
void mat_dimensions(void) {
    if (myrank == 0) {
        printf("Multiplication of 2 Matrices: a(n,m) x b(m,l) = c(n,l)\n");
        printf("Input of matrix dimensions n,m,l : ");
        fflush(stdout);
        if (scanf("%d %d %d", &n, &m, &l) != 3) {
            fprintf(stderr, "Input error.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/* ------------------ mat_alloc (slide v2.0) -------------------------- */
/* allocate ( a(n,m), b(m,ilo:ihi), c(n,ilo:ihi) ) */
void mat_alloc(int ilo_, int ihi_) {
    ilo = ilo_;
    ihi = ihi_;
    chunk_cols = (ihi >= ilo) ? (ihi - ilo + 1) : 0;

    a = (double*)malloc((size_t)n * (size_t)m * sizeof(double));
    b = (double*)malloc((size_t)m * (size_t)chunk_cols * sizeof(double));
    c = (double*)malloc((size_t)n * (size_t)chunk_cols * sizeof(double));

    if (!a || !b || !c) {
        fprintf(stderr, "rank %d: allocation failed\n", myrank);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }
}

/* ------------------ mat_input (slide v2.0) -------------------------- */
/* Each MPI process initializes full A redundantly and its chunk of B & C */
void mat_input(void) {
    /* a(i,1:m) = i  (Fortran 1-based), so in C: a[row,*] = row+1 */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A(i, j) = (double)(i + 1);
        }
    }

    /* b(1:m,i) = i for i in [ilo..ihi] */
    for (int k_global = ilo; k_global <= ihi; k_global++) {
        int k_local = k_global - ilo; /* 0..chunk_cols-1 */
        for (int j = 0; j < m; j++) {
            B(j, k_local) = (double)k_global;
        }
    }

    /* c(1:n,ilo:ihi) = 0 */
    memset(c, 0, (size_t)n * (size_t)chunk_cols * sizeof(double));
}

/* ------------------ matmul_chunk (serial triple loop from slide) ---- */
/* slide serial:
   do i = 1, n
     do k = 1, l
       do j = 1, m
         c(i,k) += a(i,j)*b(j,k)
*/
void matmul_chunk(void) {
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < chunk_cols; k++) {
            double sum = 0.0;
            for (int j = 0; j < m; j++) {
                sum += A(i, j) * B(j, k);
            }
            C(i, k) = sum;
        }
    }
}

/* ------------------ mat_output (slide v2.0) ------------------------- */
/* Each MPI prints its chunk of C ordered */
void mat_output(void) {
    MPI_Barrier(MPI_COMM_WORLD);

    for (int r = 0; r < nprocs; r++) {
        if (r == myrank) {
            printf("Matrix C: (chunk %d) columns [%d..%d]\n", myrank, ilo, ihi);
            for (int i = 0; i < n; i++) {
                /* print c(i,ilo:ihi) */
                for (int k = 0; k < chunk_cols; k++) {
                    printf("%10.2f ", C(i, k));
                }
                printf("\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(void) {
    initialize_mpi();

    mat_dimensions();                               /* slide: read n,m,l */
    workshare(myrank, nprocs, l, &ilo, &ihi);       /* slide v2.0: chunk on l */
    {
        char string[256];
        snprintf(string, sizeof(string),
                 "rank=%d work: ilo=%d ihi=%d chunk=%d",
                 myrank, ilo, ihi, (ihi >= ilo) ? (ihi - ilo + 1) : 0);
        print_ordered(string);
    }

    mat_alloc(ilo, ihi);                             /* slide v2.0 alloc */
    mat_input();                                     /* slide v2.0 init */
    matmul_chunk();                                  /* compute local C */
    mat_output();                                    /* print chunk */

    free(a); free(b); free(c);
    MPI_Finalize();
    return 0;
}