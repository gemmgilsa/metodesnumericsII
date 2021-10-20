#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TOL 10e-10
#define TOP 1000

double f (double,double);
double g( double, double);
int jacobi (double **, double **, int, double, int);
double norma(double**, int);
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
//hola ana guapa
double f (double x, double y){
    return -2*pow(x,3) + 6*x*pow(y,2);
}
//chechiiiiiii
double g (double x, double y){
    return pow(x, 3)*pow(y, 2)-x*pow(y,4);
}
int jacobi (double **b, double ** x0, int n, double tol, int top){
    int i, j,k, iter=0;
    double **x1, err=10;
    x1=(double**)calloc(n, sizeof(double*));
    while(tol<=err && iter<top){
        err=0;
        for(i=0;i<n;i++){
            for(j=0; j<n; j++){
                x1[j][i]=b[j][i];
                if(i!=0)
                    x1[j][i]+=x0[j][i-1];
                if(j!=0)
                    x1[j][i]+=x0[j-1][i];
                if(i!=(n-1))
                    x1[j][i]+=x0[j][i+1];
                if((j!=(n-1)))
                    x1[j][i]+=x0[j+1][i];
                x1[j][i]/=4;
            }
        }
        iter++;
        err=abs(norma(x1, n)-norma(x0, n));
    }
    for(i=0;i<n;i++){
        free(x1[i]);
    }
    free(x1);
    if(iter==top){
        return -1;
    }
    return iter;
}
int gaussseidel (double **b, double ** x, int n, double tol, int top){
    int i, j,k, iter=0;
    double **xaux, err=10, aux;
    xaux= (double **)malloc(n* sizeof(double*));
    for(i=0;i<n; i++){
        xaux[i]=(double*)malloc(n*sizeof(double));
    }
    while(tol<=err && iter<top){
        err=0;
        for(i=0;i<n;i++){
            for(j=0; j<n; j++){
                xaux[j][i]=x[j][i];
                aux=x[j][i];
                x[j][i]=b[j][i];
                if(i!=0)
                    x[j][i]+=x[j][i-1];
                if(j!=0)
                    x[j][i]+=x[j-1][i];
                if(i!=(n-1))
                    x[j][i]+=x[j][i+1];
                if((j!=(n-1)))
                    x[j][i]+=x[j+1][i];
                x[j][i]/=4;
            }
        }
        iter++;
        err=abs(norma(x, n)-norma(xaux, n));
    }
    for(i=0;i<n;i++){
        free(xaux[i]);
    }
    free(xaux);
    if(iter==top){
        return -1;
    }
    return iter;
}
double norma(double **x, int n){
    int i, j;
    double max=0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(abs(x[i][j]> max)){
                max=abs(x[i][j]);
            }
        }
    }
    return max;

}