#include <stdio.h>
#include <math.h>
#include <conio.h>

#include "Solvers.h"
#include "Derivative.h"


#define n 31
#define m 31

double T=200;
int input = 100;

double Re = 100;
double a_x = 0, b_x = 3;
double a_y = 0, b_y = 1;
double tau = 0.01;
char name_out[100];

double p1=1,p2=0;


double Norm(double **a,int lN,int lM)
{
    double max=fabs(a[0][0]);

    for(int i=0;i<=lN;i++)
    {
        for(int j=0;j<=lM;j++)
        {


                if(fabs(a[i][j])>max)
                {
                    max=fabs(a[i][j]);
                }

        }
    }
    return max;
}
double Sc(double **x, double **y,int ln,int lm)
{
    double s=0;
    for(int i=0;i<=ln;i++)
    {
        for(int j=0;j<=lm;j++)
        {
            s+=x[i][j]*y[i][j];
        }
    }
    return s;
}
void null(double **u,int ln,int lm)
{
    for(int i=0;i<=ln;i++)
        for(int j=0;j<=lm;j++)
            u[i][j]=0;
}
void Print(FILE *f,int t, double *x, double *y, double **u,double **v, double **p)
{
    fprintf(f, "TITLE = \"PROTEKANIE\"\n");
    fprintf(f, "VARIABLES = \"X\",\"Y\",\"U\",\"V\",\"P\"\n");
    fprintf(f, "ZONE T=\"Zone%d\",I=%d,J=%d,F=POINT\n",10000+t, m, n  );

    for (int i = 0; i < n; i++)
        {
        for (int j = 0; j < m ; j++)
        {
            fprintf(f, "%lf,%lf,%lf,%lf,%lf\n", x[i], y[j], u[i][j], v[i][j], p[i][j]);
        }
    }
}

double F_t()
{
    return (p2-p1)/(b_x-a_x);
}
void FindSpeed(double *u,int N, double h,double t)
{
    double *u_temp = new double[N+1];

    double *A=new double[N+1];
    double *B=new double[N+1];
    double *C=new double[N+1];
    double *F=new double[N+1];


    for(int i = 0; i<=N; i++)
        {
            if(i==0)
            {
                A[i]=0;
                B[i]=1.0/2.0;
                C[i]=1.0/2.0;
                F[i]=0;
            }
            else if(i==N)
            {
                A[i]=1.0/2.0;
                B[i]=1.0/2.0;
                C[i]=0;
                F[i]=0;
            }
            else
            {
                A[i] =  - 1.0 / (h*h) * 1.0/Re;
                B[i] =  1.0/tau + (2.0)/(h*h) * 1.0/Re;
                C[i] =  - 1.0 / (h*h) * 1.0/Re;
                F[i] =  - 1.0 * F_t() + u[i]/tau;
            }
        }
        SolveByScalarRun(N,u_temp,A,B,C,F);

        for(int i=0;i<=N;i++)
            u[i]=u_temp[i];


    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp;
}

double Operator_p(double **p, int i,int j, double hx, double hy)
{

    if( i==0 )
    {
        return -1*(p[i+1][j]-p[i][j])/(hx*hx);
    }
    else if( i==n )
    {
        return (p[i][j]-p[i-1][j])/(hx*hx);
    }
    else if( j==0 )
    {
        return -1*(p[i][j+1]-p[i][j])/(hy*hy);
    }
    else if( j== m )
    {
        return (p[i][j]-p[i][j-1])/(hy*hy);
    }
    else
    {
        return -1*( ( p[i+1][j] - 2.0*p[i][j] + p[i-1][j] )/(hx*hx) + ( p[i][j+1] - 2.0*p[i][j] + p[i][j-1] )/(hy*hy) );
    }
}
double RightPart_p(double **u, double **v, int i, int j, double hx, double hy)
{

    if(i==0)
    {
        //return p1;
        return -1 * F_t() / hx;
    }
    else if(i==n)
    {
        //return p2;
        return F_t() / hx;
    }
    else if(j==0)
    {
        //return p1;
        return -1 * 0 / hy;
    }
    else if(j==m)
    {
        //return p2;
        return 0 / hy;
    }
    else
    {
        return -1 * ( ( (u[i][j]-u[i-1][j])/(hx) + (v[i][j]-v[i][j-1])/(hy) ) / tau );
    }
}
void Solve_p(double **p,double **u,double **v,double hx,double hy)
{
    double **Rn = new double *[n+1];
    for(int i=0;i<=n;i++)
        Rn[i] = new double [m+1];
    double **_Rn = new double *[n+1];
    for(int i=0;i<=n;i++)
        _Rn[i] = new double [m+1];
    double **Pn = new double *[n+1];
    for(int i=0;i<=n;i++)
        Pn[i] = new double [m+1];
    double **Pn1 = new double *[n+1];
    for(int i=0;i<=n;i++)
        Pn1[i] = new double [m+1];
    double **Vn = new double *[n+1];
    for(int i=0;i<=n;i++)
        Vn[i] = new double [m+1];
    double **Sn = new double *[n+1];
    for(int i=0;i<=n;i++)
        Sn[i] = new double [m+1];
    double **Tn = new double *[n+1];
    for(int i=0;i<=n;i++)
        Tn[i] = new double [m+1];

    double pn,an,wn,pn1,bn;

    null(Rn,n,m);
    null(_Rn,n,m);
    null(Pn,n,m);
    null(Pn1,n,m);
    null(Vn,n,m);
    null(Sn,n,m);
    null(Tn,n,m);

    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
            if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                Rn[i][j] = RightPart_p(u,v,i,j,hx,hy) - Operator_p(p,i,j,hx,hy);
    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
            if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                _Rn[i][j] = Rn[i][j];

    pn = 1.;
    an = 1.;
    bn = 1.;
    wn = 1.;
    int n1 = 0;
    double eps = 1e-12;
    while(Norm(Rn,n,m) > eps)
    {
        n1++;
        pn1 = Sc(_Rn,Rn,n,m);
        bn = (pn1/pn)*(an/wn);

        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Pn1[i][j] = Rn[i][j] + bn*(Pn[i][j]-wn*Vn[i][j]);
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Vn[i][j] = Operator_p(Pn1,i,j,hx,hy);
        an = pn1/Sc(_Rn,Vn,n,m);
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Sn[i][j] = Rn[i][j] - an*Vn[i][j];
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Tn[i][j] = Operator_p(Sn,i,j,hx,hy);
        wn = Sc(Tn,Sn,n,m)/Sc(Tn,Tn,n,m);
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    p[i][j] += an*Pn1[i][j] + wn*Sn[i][j];
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Rn[i][j] = Sn[i][j] - wn*Tn[i][j];
        pn = pn1;
        for(int i=0;i<=n;i++)
            for(int j=0;j<=m;j++)
                if(     !(   (i==0 && j==0)   ||   (i==0 && j==m)   ||    (i==n && j==0)   ||   (i==n && j==m)  ||  (i==n/2  && j==0)   )     )
                    Pn[i][j] = Pn1[i][j];
    }
    for(int i=0;i<n+1;i++)
    {
        delete Rn[i];
        delete _Rn[i];
        delete Pn[i];
        delete Pn1[i];
        delete Vn[i];
        delete Sn[i];
        delete Tn[i];

    }
    delete []Rn;
    delete []_Rn;
    delete []Pn;
    delete []Pn1;
    delete []Vn;
    delete []Sn;
    delete []Tn;
    Rn=NULL;
    _Rn=NULL;
    Pn=NULL;
    Pn1=NULL;
    Vn=NULL;
    Sn=NULL;
    Tn=NULL;
}

void SolveTransport(double *S_1,double *S_1_2,double **u,double **u_n,double **v_n,int N,int M,double hx,double hy,double t)
{
    double *u_temp_x=new double[N+1];
    double *u_temp_y=new double[M+1];

    double **u_1_2=new double*[N+1];
    for(int i=0;i<=N;i++)
        u_1_2[i]=new double[M+1];


    for(int i=0;i<=N;i++)
        for(int j=0;j<=M;j++)
        {
            u_1_2[i][j]=0;
        }

    for(int j=1;j<M;j++)
    {
        u_1_2[0][j]=S_1_2[j];
        u_1_2[N][j]=S_1_2[j];
    }



    double *A=new double[M+1];
    double *B=new double[M+1];
    double *C=new double[M+1];
    double *F=new double[M+1];

    for(int j=0;j<=M;j++)
    {
        A[j]=0;
        B[j]=0;
        C[j]=0;
        F[j]=0;
    }

    for(int i=1;i<N;i++)
    {
        for(int j = 0; j<=M; j++)
        {
            if(j == 0)
            {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=S_1_2[j];
            }
            else if(j == M)
            {
                A[j]=0;
                B[j]=1;
                C[j]=0;
                F[j]=S_1_2[j];
            }
            else
            {
                A[j] =  - v_n[i][j] / (2*hy) - 1.0 / (hy*hy)*(1.0/Re);
                B[j] =  2.0/tau+(2.0)/(hy*hy)*(1.0/Re);
                C[j] =  v_n[i][j] / (2*hy) - 1.0 / (hy*hy)*(1.0/Re);
                F[j] =  - u_n[i][j] * Derivative2d_x(u_n,i,j,hx,'c') + (1.0/Re)*Laplas2d_x(u_n,i,j,hx) + 2.0*u_n[i][j]/tau;
            }
        }
        SolveByScalarRun(M,u_temp_y,A,B,C,F);

        for(int j=0;j<=M;j++)
            u_1_2[i][j]=u_temp_y[j];

        }

        for(int i=0;i<=N;i++)
        {
            for(int j=0;j<=M;j++)
            {
                u[i][j]=u_1_2[i][j];
            }
        }

        delete []A;
        delete []B;
        delete []C;
        delete []F;


        A=new double[N+1];
        B=new double[N+1];
        C=new double[N+1];
        F=new double[N+1];


        for(int j=1;j<M;j++)
        {
            for(int i = 0; i<=N; i++)
            {
                if(i==0)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S_1[j];
                }
                else if(i==N)
                {
                    A[i]=0;
                    B[i]=1;
                    C[i]=0;
                    F[i]=S_1[j];
                }
                else
                {
                    A[i] =  - u_n[i][j] / (2*hx) - 1.0 / (hx*hx)*(1.0/Re);
                    B[i] =  2.0/tau+(2.0)/(hx*hx)*(1.0/Re);
                    C[i] =  u_n[i][j] / (2*hx) - 1.0 / (hx*hx)*(1.0/Re);
                    F[i] =  - v_n[i][j] * Derivative2d_y(u_1_2,i,j,hy,'c') + (1.0/Re)*Laplas2d_y(u_1_2,i,j,hy)+2.0*u_1_2[i][j]/tau;
                }
            }
        SolveByScalarRun(N,u_temp_x,A,B,C,F);

        for(int i=0;i<=N;i++)
            u[i][j]=u_temp_x[i];
        }

        for(int i=1;i<N;i++)
        {
            u[i][0]=S_1[0];
            u[i][M]=S_1[M];
        }


    for(int i=0;i<=N;i++)
    {
        delete []u_1_2[i];
    }
    delete []u_1_2;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []u_temp_x;
    delete []u_temp_y;
    u_1_2=NULL;
    A=NULL;
    B=NULL;
    C=NULL;
    F=NULL;
    u_temp_x=NULL;
    u_temp_y=NULL;
}




int main()
{
    double hx = (b_x-a_x)/(n-1);
    double hy = (b_y-a_y)/(m-1);

    printf("%lf %lf\n",hx,hy);

    int vN=n;
    int vM=m-1;

    int uN=n-1;
    int uM=m;

    double **u = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        u[i] = new double [uM+1];

    double **u_check_n = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        u_check_n[i] = new double [uM+1];

    double **u_check_n1 = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        u_check_n1[i] = new double [uM+1];

    double **u_stop = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        u_stop[i] = new double [uM+1];



    double **w_u = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        w_u[i] = new double [uM+1];

    double **v = new double *[vN+1];
    for(int i=0;i<vN+1;i++)
        v[i] = new double [vM+1];

    double **p = new double *[n+1];
    for(int i=0;i<=n;i++)
        p[i] = new double [m+1];


    double **u_n = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        u_n[i] = new double [uM+1];

    double **v_n = new double *[vN+1];
    for(int i=0;i<vN+1;i++)
        v_n[i] = new double [vM+1];

    for(int i=0;i<=n;i++)
        for(int j=0;j<=m;j++)
        {
            p[i][j]=0;
        }
    for(int i=0;i<=vN;i++)
        for(int j=0;j<=vM;j++)
        {
            v[i][j]=0;
            v_n[i][j]=0;
        }
    for(int i=0;i<=uN;i++)
        for(int j=0;j<=uM;j++)
        {
            w_u[i][j]=0;
            u[i][j]=0;
            u_check_n[i][j]=0;
            u_check_n1[i][j]=0;
            u_n[i][j]=0;
        }

    double **v_proec = new double *[uN+1];
    for(int i=0;i<uN+1;i++)
        v_proec[i] = new double [uM+1];

    double **u_proec = new double *[vN+1];
    for(int i=0;i<vN+1;i++)
        u_proec[i] = new double [vM+1];

    double **p_out = new double *[n];
    for(int i=0;i<n;i++)
        p_out[i] = new double [m];
    double **u_out = new double *[n];
    for(int i=0;i<n;i++)
        u_out[i] = new double [m];
    double **ut = new double *[n];
    for(int i=0;i<n;i++)
        ut[i] = new double [m];
    double **v_out = new double *[n];
    for(int i=0;i<n;i++)
        v_out[i] = new double [m];

    double *x = new double [n+1];
    x[0]=a_x;
    for(int i=1;i<n+1;i++)
        x[i]=x[i-1]+hx;
    double *y = new double [m+1];
    y[0]=-b_y/2.0;
    for(int i=1;i<m+1;i++)
        y[i]=y[i-1]+hy;

    double delta_p = p2 - p1;
    double R = b_y/2.0;
    double L = (b_x - a_x);

    for(int j=0;j<m;j++)
        ut[0][j]= ( delta_p*Re / ( 2 * L) ) * ( (y[j]*y[j]) - ( R * R )  );


    for(int j=0;j<m;j++)
    {
        for(int i=1;i<n;i++)
        {
            ut[i][j]=ut[i-1][j];
        }
    }

    double *S=new double[uM+1];
    double *S_g=new double[uM+1];
    double *S_g_1_2=new double[uM+1];

    double *S_gv=new double[vM+1];
    double *S_gv_1_2=new double[vM+1];

    for(int i=0;i<=uM;i++)
    {
        S[i]=0;
        S_g[i]=0;
        S_g_1_2[i]=0;
    }
    for(int i=0;i<=vM;i++)
    {
        S_gv[i]=0;
        S_gv_1_2[i]=0;
    }


    int it = 0;
    for(double t=0;t<T;t+=tau)
    {
        it++;
        for(int i=1;i<uN;i++)
        {
            for(int j=1;j<uM;j++)
            {
                v_proec[i][j]=(v_n[i][j]+v_n[i+1][j]+v_n[i][j-1]+v_n[i+1][j-1])/4.0;
            }
        }
        for(int i=1;i<vN;i++)
        {
            for(int j=1;j<vM;j++)
            {
                u_proec[i][j]=(u_n[i-1][j+1]+u_n[i][j+1]+u_n[i][j]+u_n[i-1][j])/4.0;
            }
        }

        FindSpeed(S,uM,hy,t);
        for(int i=0;i<=uM;i++)
        {
            S_g[i] = S[i] + tau * F_t();
        }
        for(int i=0;i<=uM;i++)
        {
            S_g_1_2[i] = ( S_g[i] + S_g_1_2[i] ) / 2.0;
        }

        for(int i=0;i<=uN;i++)
            for(int j=0;j<=uM;j++)
            {
                u_check_n[i][j]=S_g[j];
            }

        SolveTransport(S_g,S_g_1_2,u,u_n,v_proec,uN,uM,hx,hy,t);

         for(int i=1;i<uN;i++)
            for(int j=1;j<uM;j++)
            {
                u_check_n1[i][j]=u_check_n[i][j]-u[i][j];
            }

        for(int i=0;i<=uM;i++)
        {
            S_g_1_2[i] = S_g[i];
        }

        SolveTransport(S_gv,S_gv_1_2,v,v_n,u_proec,vN,vM,hx,hy,t);

        p[n/2][0]=0;
        Solve_p(p,u,v,hx,hy);

        printf("\nNei = %2.20lf\n",(p[n/2][1]-p[n/2][0])/(hy));
        getch();

        for(int i=0;i<=uN;i++)
            for(int j=0;j<=uM;j++)
                if(!(   (i==0 && j==0)   ||   (i==0 && j==uM)   ||    (i==uN && j==0)   ||   (i==uN && j==uM)     ))
                    u[i][j] = u[i][j] - tau*( (p[i+1][j]-p[i][j])/(hx) );

        for(int i=0;i<=vN;i++)
            for(int j=0;j<=vM;j++)
                if(!(   (i==0 && j==0)   ||   (i==0 && j==vM)   ||    (i==vN && j==0)   ||   (i==vN && j==vM)     ))
                    v[i][j] = v[i][j] - tau*( (p[i][j+1]-p[i][j])/(hy) );

        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                if     (   (i==0 && j==0)   ||   (i==0 && j==m-1)   ||    (i==n-1 && j==0)   ||   (i==n-1 && j==m-1)     )
                {
                    p_out[i][j]=0;
                    v_out[i][j]=0;
                    u_out[i][j]=0;
                }
                else
                {
                    u_out[i][j]=(u[i][j+1]+u[i][j])/2.0;
                    v_out[i][j]=(v[i][j]+v[i+1][j])/2.0;
                    p_out[i][j]=(p[i][j]+p[i][j+1]+p[i+1][j+1]+p[i+1][j])/4.0;
                }

            }
        }

        sprintf(name_out,"./Solve/Zone_%d.dat",it);
        if( it%input == 0 )
        {
            FILE *f = fopen(name_out,"w");
            Print(f,it,x,y,u_out,v_out,p_out);
            fclose(f);
        }

        for(int i=0;i<=uN;i++)
            for(int j=0;j<=uM;j++)
                u_stop[i][j]=u[i][j]-u_n[i][j];

        if(Norm(u_stop,uN,uM)<pow(10,-10))
        {
            FILE *f = fopen("stop.txt","w");
            fprintf(f,"stop on %lf",t);
            fclose(f);
            break;
        }

        for(int i=0;i<=vN;i++)
            for(int j=0;j<=vM;j++)
            {
                v_n[i][j]=v[i][j];
            }
        for(int i=0;i<=uN;i++)
            for(int j=0;j<=uM;j++)
            {
                u_n[i][j]=u[i][j];
            }

        for(int i=0;i<n;i++)
        {
            for(int j=0;j<m;j++)
            {
                w_u[i][j]=ut[i][j]-u_out[i][j];
            }
        }

        double divV=0;
        for(int i=1;i<n;i++)
        {
            for(int j=1;j<m;j++)
            {
                divV+=((u_n[i][j]-u_n[i-1][j])/(hx) + (v_n[i][j]-v_n[i][j-1])/(hy));
            }
        }

        printf("\rtime = %2.5lf div(V*) = %2.10lf Norm(w) = %2.20lf",t,fabs(divV),Norm(w_u,n-1,m-1));
    }

    sprintf(name_out,"./Solve/Zone_%d.dat",it);
    FILE *f = fopen(name_out,"w");
    Print(f,it,x,y,u_out,v_out,p_out);
    fclose(f);

    delete []S;
    delete []S_g;
    delete []S_g_1_2;

    for(int i=0;i<=uN;i++)
    {
        delete u[i];
        delete u_check_n[i];
        delete u_check_n1[i];
        delete u_stop[i];
        delete u_n[i];
        delete v_proec[i];
    }
    delete []u;
    delete []u_check_n;
    delete []u_check_n1;
    delete []u_stop;
    delete []u_n;
    delete []v_proec;

    for(int i=0;i<=vN;i++)
    {
        delete v[i];
        delete v_n[i];
        delete u_proec[i];
    }
    delete []v;
    delete []v_n;
    delete []u_proec;

    for(int i=0;i<=n;i++)
    {
        delete p[i];
    }
    delete []p;

    for(int i=0;i<n;i++)
    {
        delete v_out[i];
        delete p_out[i];
        delete u_out[i];
    }
    delete []u_out;
    delete []v_out;
    delete []p_out;

    u=NULL;
    v=NULL;
    p=NULL;
    u_n=NULL;
    v_n=NULL;
    u_proec=NULL;
    u_stop=NULL;
    u_check_n=NULL;
    u_check_n1=NULL;
    v_proec=NULL;
    v_out=NULL;
    u_out=NULL;
    p_out=NULL;
    S=NULL;
    S_g=NULL;
    S_g_1_2=NULL;

    printf("\nFinalisation work\n");
    return 0;
}
