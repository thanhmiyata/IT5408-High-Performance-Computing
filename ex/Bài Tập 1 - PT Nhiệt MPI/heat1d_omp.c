#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

/* Tham số */
#define M    101
#define L    1.0f
#define D    0.1f
#define TIME 1.0f
#define DT   0.0005f

float dx;
int   Ntime;

/* DHB2 vẫn giữ đúng style con trỏ như thầy, nhưng tính theo miền [start, stop) */
void DHB2_range(float *T, float *dT, int start, int stop) {
    int i;
    float c, l, r;
    for (i = start; i < stop; i++) {
        c = *(T + i);
        l = (i == 0)   ? 100.0f : *(T + i - 1);
        r = (i == M-1) ? 25.0f  : *(T + i + 1);
        *(dT + i) = D * (l - 2.0f * c + r) / (dx * dx);
    }
}

int main() {
    int i, t;
    float *T, *dT, *Tnew;

    dx = L / (M - 1);
    Ntime = (int)(TIME / DT);

    T    = (float*)malloc(M * sizeof(float));
    dT   = (float*)malloc(M * sizeof(float));
    Tnew = (float*)malloc(M * sizeof(float));
    if (!T || !dT || !Tnew) return 1;

    for (i = 0; i < M; i++) *(T + i) = 25.0f;

    FILE *f = fopen("T_omp.csv", "w");
    if (!f) return 2;

    fprintf(f, "t");
    for (i = 0; i < M; i++) fprintf(f, ",T%d", i);
    fprintf(f, "\n");

    for (t = 0; t <= Ntime; t++) {
        fprintf(f, "%f", t * DT);
        for (i = 0; i < M; i++) fprintf(f, ",%f", *(T + i));
        fprintf(f, "\n");

        /* Chia miền theo slide: start=ID*Mc, stop=start+Mc (Mc=M/NT) */
        #pragma omp parallel
        {
            int ID = omp_get_thread_num();
            int NT = omp_get_num_threads();

            int Mc = M / NT;
            int start = ID * Mc;
            int stop  = (ID == NT - 1) ? M : (start + Mc); /* thread cuối ăn phần dư */

            /* Tính dT trong miền */
            DHB2_range(T, dT, start, stop);

            /* Đồng bộ vì phụ thuộc dữ liệu (slide nhấn mạnh) */
            #pragma omp barrier

            /* Update Tnew trong miền */
            for (int j = start; j < stop; j++) {
                *(Tnew + j) = *(T + j) + *(dT + j) * DT;
            }
        }

        /* Copy Tnew -> T */
        for (i = 0; i < M; i++) *(T + i) = *(Tnew + i);
    }

    fclose(f);
    free(T);
    free(dT);
    free(Tnew);
    return 0;
}