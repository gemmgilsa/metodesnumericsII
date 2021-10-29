/* ANA VICTORIA GALINDO, GEMMA GIL SAURA  MÈTODES NUMÈRICS II */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 10e-8
#define TOP 100000

double f (double,double);
double g( double, double);
int jacobi (double **, double **, int, double, int);
int gaussseidel (double **, double ** , int , double , int);
int SOR (double **, double ** , int , double , int, double);
double norma(double**, double**, int);
void printar( int, int, double **, double, double);
int main(void) {
    int n, i, j, k,t, min=TOP;
    double h, w, wo;
    double **b, **x;
    printf("n=?\n");
    scanf("%d", &n);
    /*Calculem la h=1/(n+1) i escrivim el que dona*/
    h=(double)1/(n+1);
    printf("h=%le\n",h);
    /*Inicialitzem les matrius i els donem memòria*/
    b=(double **)malloc(n* sizeof(double*));
    if(b==NULL){
        printf("Error en la matriu b\n");
        exit(1);
    }
    for(i=0;i<n;i++){
        b[i]=(double*)malloc(n*sizeof(double));
        if(b[i]==NULL){
            printf("Error en la matriu, fila b[%d]\n", i);
            exit(1);
        }
    }
    x=(double **)calloc(n, sizeof(double*));
    if(x==NULL){
        printf("Error en la matriu x\n");
        for(i=0;i<n;i++){
            free(b[i]);
        }
        free(b);
        exit(1);
    }
    for(i=0;i<n;i++){
        x[i]=(double*)calloc(n,sizeof(double));
        if(x[i]==NULL){
            printf("Error en la matriu x[%d]\n", i);
            for(i=0;i<n;i++){
                free(b[i]);
            }
            free(b);
            exit(1);
        }
    }

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
    /*Calculem el terme independent i inicialitzem l'aproximació inicial a 0 */
   /* printf("terme independent i aproximacio inicial\n");
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            printf(" j = %d i= %d bij=%le   uij=%le\n", j+1, i+1, b[i][j], x[i][j]);
        }
    }*/
    printf("Comencem el proces iteratiu\n");
    /*Cridem a jacobi i calculem les iteracions que li costa en convergir*/
    k=jacobi(b, x, n, TOL, TOP);
    if(k<0){
        printf("El metode de jacobi no convergeix \n");
    }else{
        printf("Jacobi convergeix a la iter %d\n", k);
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            x[i][j]=0;
        }
    }
    /*Cridem a gausseidel i calculem les iteracions tambe*/
    k=gaussseidel(b, x, n, TOL, TOP);
    if(k<0){
        printf("El metode de gausseidel no convergeix \n");
    }else{
        printf("gausseidel convergeix a la iter %d\n", k);
    }
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            x[i][j]=0;
        }
    }

    /*Definim w, anem provant per a diferents w, i mirem quina és la óptima, és a dir, amb quina w SOR convergeix amb menos iteracions */
    for(i=1;i<200;i++){
        w=i*0.01;
        k=SOR(b, x, n, TOL, TOP, w);
     //   printf("w=%le, k=%d \n",w, k);
        if(k<min && k>=0){
            min=k;
            wo=w;
        }
        for(j=0;j<n;j++){
            for(t=0;t<n;t++){
                x[j][t]=0;
            }
        }
    }
    if(k<0){
        printf("El metode de SOR no convergeix \n");
    }else{
        printf("SOR convergeix a la iter %d\n", min);
        printf("La w optima es %le\n",wo);
    }
    for(i=0;i<n;i++){
        free(b[i]);
        free(x[i]);
    }
    free(b);
    free(x);
    return 0;
}
/*double f(double x, double y){
    return (x*x-2)*sin(y)+(4-2*y*y)*cos(x);
}
double g(double x, double y){
    return x*x*sin(y)-2*y*y*cos(x);
}*/
double f (double x, double y){
    return -2*pow(x,3) + 6*x*pow(y,2);
}

double g (double x, double y){
    return pow(x, 3)*pow(y, 2)-x*pow(y,4);
}
int jacobi (double **b, double ** x0, int n, double tol, int top){
    int i, j, iter=0;
    double **x1, err=1, errant, quo;
    x1=(double**)calloc(n, sizeof(double*));
    if(x1==NULL){
        printf("Error en la matriu\n");
        exit(1);
    }
    for(i=0;i<n;i++) {
        x1[i] = (double *) calloc(n, sizeof(double));
        if (x1[i] == NULL) {
            printf("Error en la matriu\n");
            exit(1);
        }
    }
    while(err>tol && iter<top){
        errant=err;
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
                x1[j][i]=x1[j][i]/4;
            }
        }
        iter ++;
        err=fabs(norma(x1,x0, n));
        quo=err/errant;
        // printar(iter, n, x1, err, quo);
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                x0[i][j]=x1[i][j];
            }
        }
    }
    printf("Radi espectral de jacobi = %le\n", quo);
    for(i=0;i<n;i++){
        free(x1[i]);
    }
    free(x1);
    if(iter==top){
        iter=-1;
        printar(iter, n, x1, err, quo);
        return -1;
    }
    return iter;
}
int gaussseidel (double **b, double ** x, int n, double tol, int top){
    int i, j, iter=0;
    double **xaux, err=1, errant, quo;
    xaux= (double **)malloc(n* sizeof(double*));
    for(i=0;i<n; i++){
        xaux[i]=(double*)malloc(n*sizeof(double));
    }
    while(tol<=err && iter<=top){
        errant=err;
        for(i=0;i<n;i++){
            for(j=0; j<n; j++){
                xaux[j][i]=x[j][i];
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
        err=fabs(norma(x,xaux, n));
        quo=err/errant;
        //printf("%le\n",err);
    }
    printf("Radi espectral de Gauss-Seidel = %le\n", quo);
    for(i=0;i<n;i++){
        free(xaux[i]);
    }
    free(xaux);
    if(iter==top){
        return -1;
    }
    return iter;
}

int SOR (double **b, double ** x, int n, double tol, int top, double w){
    int i, j, iter=0;
    double **xaux, err=1,errant, quo, aux;
    xaux= (double **)malloc(n* sizeof(double*));
    for(i=0;i<n; i++){
        xaux[i]=(double*)malloc(n*sizeof(double));
    }
    while(tol<=err && iter<top){
        errant=err;
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
                x[j][i]=(x[j][i]*(double)w/4 + (1-w)*aux);
            }
        }
        iter++;
        err=fabs(norma(x,xaux,n));
        //printf("%le\n",err);
        quo=err/errant;
    }
    printf("Radi espectral de SOR = %le\n", quo);
    for(i=0;i<n;i++){
        free(xaux[i]);
    }
    free(xaux);
    if(iter==top){
        return -1;
    }
    return iter;
}

double norma(double **x, double **y, int n){
    int i, j;
    double max=0;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(fabs(x[i][j]-y[i][j])> max){
                max=fabs(x[i][j]-y[i][j]);
            }
        }
    }
    return max;

}

void printar( int iter, int n, double **x, double err, double quo){
    int i, j;
    if(iter<0){
        printf("err_infinit= %le\n", err);
        return;
    }
    for(j=0;j<n;j++){
        for(i=0;i<n;i++){
            printf("k=%d   j=%d   i=%d   uij = %le\n", iter, j, i, x[i][j]);
        }
    }
    printf("k= %d, error= %le ,  quo=%le\n", iter,err, quo);
    return ;
}