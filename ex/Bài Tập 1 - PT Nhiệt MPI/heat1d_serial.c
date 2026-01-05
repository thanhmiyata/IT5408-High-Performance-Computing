#include <stdio.h>
#include <stdlib.h>

/* Tham số bài toán */
#define M    101
#define L    1.0f
#define D    0.1f
#define TIME 1.0f
#define DT   0.0005f

float dx;
int   Ntime;

/* DHB2 theo đúng style trong slide: dùng con trỏ *(T+i) */
void DHB2(float *T, float *dT) {
    int i;
    float c, l, r;

    /* Slide có chỗ i < M-1 nhưng lại check i==M-1.
       Để đúng bài toán và đúng điều kiện biên, chạy i < M */
    for (i = 0; i < M; i++) {
        c = *(T + i);
        l = (i == 0)   ? 100.0f : *(T + i - 1);
        r = (i == M-1) ? 25.0f  : *(T + i + 1);
        *(dT + i) = D * (l - 2.0f * c + r) / (dx * dx);
    }
}

int main() {
    int i, t;
    float *T, *dT;

    dx = L / (M - 1);
    Ntime = (int)(TIME / DT);

    T  = (float*)malloc(M * sizeof(float));
    dT = (float*)malloc(M * sizeof(float));
    if (!T || !dT) return 1;

    /* 1) Khởi tạo: T(i)=25, 0<i<M-1 (slide) */
    for (i = 0; i < M; i++) *(T + i) = 25.0f;

    /* 2) Lưu mảng T từng bước thời gian ra file (slide yêu cầu) */
    FILE *f = fopen("T_serial.csv", "w");
    if (!f) return 2;

    fprintf(f, "t");
    for (i = 0; i < M; i++) fprintf(f, ",T%d", i);
    fprintf(f, "\n");

    /* 3) Tích hợp theo thời gian (flow đúng slide trang tích hợp) */
    for (t = 0; t <= Ntime; t++) {
        fprintf(f, "%f", t * DT);
        for (i = 0; i < M; i++) fprintf(f, ",%f", *(T + i));
        fprintf(f, "\n");

        DHB2(T, dT);
        for (i = 0; i < M; i++)
            *(T + i) = *(T + i) + *(dT + i) * DT;
    }

    fclose(f);
    free(T);
    free(dT);
    return 0;
}