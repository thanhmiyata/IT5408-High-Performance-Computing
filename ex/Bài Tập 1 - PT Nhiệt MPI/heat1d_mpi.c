#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Lưu ý: để dễ đúng slide, chọn M chia hết cho số process */
#define M    100
#define L    1.0f
#define D    0.1f
#define TIME 1.0f
#define DT   0.0005f

float dx;
int   Ntime;

/* DHB2 theo slide trang cuối: dùng Tl/Tr và con trỏ */
void DHB2(float *Ts, float Tl, float Tr, float *dTs, int ms) {
    int i;
    float c, l, r;
    for (i = 0; i < ms; i++) {
        c = *(Ts + i);
        l = (i == 0)      ? Tl : *(Ts + i - 1);
        r = (i == ms - 1) ? Tr : *(Ts + i + 1);
        *(dTs + i) = D * (l - 2.0f * c + r) / (dx * dx);
    }
}

int main(int argc, char **argv) {
    int rank, size;
    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    dx = L / (M - 1);
    Ntime = (int)(TIME / DT);

    int mc = M / size; /* slide: mc = N/NP */

    float *T  = NULL;                 /* chỉ root có */
    float *Ts = (float*)malloc(mc * sizeof(float));
    float *dTs= (float*)malloc(mc * sizeof(float));
    if (!Ts || !dTs) MPI_Abort(MPI_COMM_WORLD, 1);

    if (rank == 0) {
        T = (float*)malloc(M * sizeof(float));
        if (!T) MPI_Abort(MPI_COMM_WORLD, 2);
        for (int i = 0; i < M; i++) *(T + i) = 25.0f; /* init tại root */
    }

    /* Bước 3 slide: Scatter */
    MPI_Scatter(T, mc, MPI_FLOAT, Ts, mc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    FILE *f = NULL;
    if (rank == 0) {
        f = fopen("T_mpi.csv", "w");
        fprintf(f, "t");
        for (int i = 0; i < M; i++) fprintf(f, ",T%d", i);
        fprintf(f, "\n");
    }

    for (int t = 0; t <= Ntime; t++) {
        /* Gather để root lưu file theo yêu cầu bài tập */
        MPI_Gather(Ts, mc, MPI_FLOAT, T, mc, MPI_FLOAT, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            fprintf(f, "%f", t * DT);
            for (int i = 0; i < M; i++) fprintf(f, ",%f", *(T + i));
            fprintf(f, "\n");
        }

        /* B4.1-a: truyền thông Tl (đúng block code slide) */
        float Tl, Tr;

        if (rank == 0) {
            Tl = 100.0f;
            if (size > 1) {
                MPI_Send(Ts + mc - 1, 1, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
            }
        } else if (rank == size - 1) {
            MPI_Recv(&Tl, 1, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &stat);
        } else {
            MPI_Send(Ts + mc - 1, 1, MPI_FLOAT, rank + 1, rank, MPI_COMM_WORLD);
            MPI_Recv(&Tl, 1, MPI_FLOAT, rank - 1, rank - 1, MPI_COMM_WORLD, &stat);
        }

        /* B4.1-b: truyền thông Tr (bám code slide; sửa count nhận = 1 cho đúng) */
        if (rank == size - 1) {
            Tr = 25.0f;
            if (size > 1) {
                MPI_Send(Ts, 1, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
            }
        } else if (rank == 0) {
            MPI_Recv(&Tr, 1, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &stat);
        } else {
            MPI_Send(Ts, 1, MPI_FLOAT, rank - 1, rank, MPI_COMM_WORLD);
            MPI_Recv(&Tr, 1, MPI_FLOAT, rank + 1, rank + 1, MPI_COMM_WORLD, &stat);
        }

        /* B4.2: tính toán */
        DHB2(Ts, Tl, Tr, dTs, mc);

        for (int i = 0; i < mc; i++)
            *(Ts + i) = *(Ts + i) + *(dTs + i) * DT;
    }

    if (rank == 0) {
        fclose(f);
        free(T);
    }
    free(Ts);
    free(dTs);

    MPI_Finalize();
    return 0;
}