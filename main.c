#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double f (double,double);
double g( double, double);
double jacobi ** (double **, double **, int);
int main(void) {
    int n, i, j;
    double h;
    double **b, **x;
    printf("n=?\n");
    scanf("%d", &n);
    h=(double)1/(n+1);
    printf("%le\n",h);
    b= (double **)malloc(n* sizeof(double*));
    x=(double **)calloc(n, sizeof(double*));
    if(b==NULL || x==NULL){
        printf("error en la memòria\n");
        exit (1);
    }

    // comentari
    for(i=0; i<n; i++) {
        for (j = 0; j < n; j++) {
            b[i][j] = h * h * f((i + 1) * h, (j + 1) * h);
            if (i == 0)
                b[i][j] += g(i * h, (j + 1) * h);
            if (i == n - 1)
                b[i][j] += g((n + 1) * h, (j + 1) * h);
            if (j == 0)
                b[i][j] += g((i + 1) * h, j * h);
            if (j == n - 1)
                b[i][j] += g((i + 1) * h, (n + 1) * h);
        }
    }




    for(i=0;i<n;i++){
        free(b[i]);
        free(x[i]);
    }
    free(b);
    free(x);
    return 0;
}
//hola ana
double f (double x, double y){
    return -2*pow(x,3) + 6*x*pow(y,2);
}
double g (double x, double y){
    return pow(x, 3)*pow(y, 2)-x*pow(y,4);
}
int jacobi (double **b, double ** x0, int n){
    int i, j,k, cont=0;
    double **x1;
    x1=(double**)calloc(n, sizeof(double*));
    for(i=0;i<n;i++){
        for(j=0; j<n; j++){
            x1[j][i]=b[j][i];
            for(k=0;k<n;k++) {
                x1[j][i] -=;
            }
            cont ++;
        }
    }
}