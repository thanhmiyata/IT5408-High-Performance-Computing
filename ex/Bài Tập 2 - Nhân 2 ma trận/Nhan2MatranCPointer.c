#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
//=========================
void DisplayMatrix(float *A, int row,  int col){
  int i,j;
  for(i=0;i<row;i++){
    for(j=0;j<col;j++) printf("  %0.0f",*(A+i*col+j));
    printf("\n");
  }
}
//=========================
int main()
{
  int m, n, p, i, j, k;
  float tmp,sum = 0.0;
// Khai bao Ma tran A, B, C dang con tro
  float *A, *B, *C;
  A = (float *) malloc ((m*n)*sizeof(float));
  B = (float *) malloc ((n*p)*sizeof(float));
  C = (float *) malloc ((m*p)*sizeof(float));
// Nhap Ma tran A 
  printf("Enter the rows and columns of matrix A:\n");
  printf("The rows:");
  scanf("%d", &m);
  printf("The columns:");
  scanf("%d", &n);
  printf("Enter the elements of matrix A\n");
   for (  i = 0 ; i < m ; i++ )
    for ( j = 0 ; j < n ; j++ ){
      scanf("%f", &tmp);
      *(A+i*n+j) = tmp;
    }
  printf("The Matrix A:\n");      
  DisplayMatrix(A, m, n);

// Nhap Ma tran B
  printf("Enter the columns of matrix B:\n");
  printf("The columns:");
  scanf("%d", &p);
  printf("Enter the elements of matrix B\n");
   for (  i = 0 ; i < n ; i++ )
    for ( j = 0 ; j < p ; j++ ){
      scanf("%f", &tmp);
      *(B+i*p+j) = tmp;
    }
  printf("The Matrix B:\n");
  DisplayMatrix(B, n, p);

// Nhan 2 Ma Tran
  for ( i = 0 ; i < m ; i++ )
    {
      for ( j = 0 ; j < p ; j++ )
      {
        for ( k = 0 ; k < n ; k++ )
        {
          sum = sum + (*(A + i*n + k)) * (*(B + k*p + j));
        }
 
        *(C + i*p + j) = sum;
        sum = 0;
      }
    }

// In ket qua 
  printf("The Matrix C:\n");
  DisplayMatrix(C, m, p);
//
  return 0;
}

