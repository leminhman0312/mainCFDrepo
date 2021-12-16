/* Exact Solutions to the Navier-Stokes Equation in 2D
   This code was originally written by John Burkardt
   and distributed under the GNU LGPL license.
   References:
     Bennett Fox, Algorithm 647: Implementation and Relative Efficiency of
       Quasirandom Sequence Generators, ACM Transactions on Mathematical
       Software, Volume 12, Number 4, December 1986, pages 362-376.
     Geoffrey Taylor, A E Green, Mechanism for the production of small eddies
       from large ones, Proceedings of the Royal Society of London, Series A,
       Volume 158, 1937, pages 499-521.
     Geoffrey Taylor, On the decay of vortices in a viscous fluid,
       Philosophical Magazine, Volume 46, 1923, pages 671-674.
     Maxim Olshanskii, Leo Rebholz, Application of barycenter refined meshes
       in linear elasticity and incompressible fluid dynamics, ETNA: Electronic
       Transactions in Numerical Analysis, Volume 38, pages 258-274, 2011.
     Paul Bratley, Bennett Fox, Linus Schrage, A Guide to Simulation, Second
       Edition, Springer, 1987, ISBN: 0387964673, LC: QA76.9.C65.B73.
     Peter Lewis, Allen Goodman, James Miller, A Pseudo-Random Number Generator
       for the System/360, IBM Systems Journal, Volume 8, Number 2, 1969, pages 136-143.
     Pierre L'Ecuyer, Random Number Generation, in Handbook of Simulation,
       edited by Jerry Banks, Wiley, 1998, ISBN: 0471134031, LC: T57.62.H37.

   streamlining by dudley.benton@gmail.com
   change output format from GNUplot to generic
 */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

char header[FILENAME_MAX];
char filename[FILENAME_MAX];

void grid_2d(int x_num,double x_lo,double x_hi,int y_num,double y_lo,double y_hi,double x[],double y[])
  {
  int i,j;
  double xi,yj;
  if(x_num==1)
    {
    for(j=0;j<y_num;j++)
      for(i=0;i<x_num;i++)
        x[i+j*x_num]=(x_lo+x_hi)/2.;
    }
  else
    {
    for(i=0;i<x_num;i++)
      {
      xi=((double)(x_num-i-1)*x_lo+(double)(i)*x_hi)/(double)(x_num-1);
      for(j=0;j<y_num;j++)
        x[i+j*x_num]=xi;
      }
    }
  if(y_num==1)
    {
    for(j=0;j<y_num;j++)
      for(i=0;i<x_num;i++)
        y[i+j*x_num]=(y_lo+y_hi)/2.;
    }
  else
    {
    for(j=0;j<y_num;j++)
      {
      yj=((double)(y_num-j-1)*y_lo+(double)(j)*y_hi)/(double)(y_num-1);
      for(i=0;i<x_num;i++)
        y[i+j*x_num]=yj;
      }
    }
  }

void ns2de_V2D(char*header,int n,double x[],double y[],double u[],double v[],double p[],double s)
  {
  int i;
  FILE*fp;
  strcpy(filename,header);
  strcat(filename,".v2d");
  printf("writing data: %s\n",filename);
  fp=fopen(filename,"wt");
  for(i=0;i<n;i++)
    fprintf(fp,"%lG %lG %lG %lG %lG %lG %lG\n",x[i],y[i],u[i],v[i],s*u[i],s*v[i],p[i]);
  fclose(fp);
  }

void uvp_lukas(double nu,double rho,int n,double x[],double y[],double t,double u[],double v[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
    u[i]=-cos(M_PI*x[i])/M_PI;
    v[i]=-y[i]*sin(M_PI*x[i]);
    p[i]=-0.;
    }
  }

void lukas()
  {
  int n,x_num=21,y_num=21;
  double nu,rho,s,t,x_hi,x_lo,y_hi,y_lo,*p,*u,*v,*x,*y;
  printf("Lukas Bystricky Flow:\n");
  printf("Generate a velocity field on a regular grid.\n");
  x_lo=y_lo=0.;
  x_hi=y_hi=1.;
  x=(double*)malloc(x_num*y_num*sizeof(double));
  y=(double*)malloc(x_num*y_num*sizeof(double));
  grid_2d(x_num,x_lo,x_hi,y_num,y_lo,y_hi,x,y);
  nu=rho=1.;
  n=x_num*y_num;
  t=0.;
  u=(double*)malloc(x_num*y_num*sizeof(double));
  v=(double*)malloc(x_num*y_num*sizeof(double));
  p=(double*)malloc(x_num*y_num*sizeof(double));
  uvp_lukas(nu,rho,n,x,y,t,u,v,p);
  strcpy(header,"lukas");
  s=0.25;
  ns2de_V2D(header,n,x,y,u,v,p,s);
  free(p);
  free(u);
  free(v);
  free(x);
  free(y);
  }

void uvp_poiseuille(double nu,double rho,int n,double x[],double y[],double t,double u[],double v[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
    u[i]=1.-y[i]*y[i];
    v[i]=0.;
    p[i]=-2.*rho*nu*x[i];
    }
  }

void poiseuille()
  {
  int n,x_num=21,y_num=21;
  double nu,rho,s,t,x_hi,x_lo,y_hi,y_lo,*p,*u,*v,*x,*y;
  printf("Spiral Flow:\n");
  printf("Generate a velocity field on a regular grid.\n");
  x_lo=0.;
  x_hi=6.;
  y_lo=-1.;
  y_hi=1.;
  x=(double*)malloc(x_num*y_num*sizeof(double));
  y=(double*)malloc(x_num*y_num*sizeof(double));
  grid_2d(x_num,x_lo,x_hi,y_num,y_lo,y_hi,x,y);
  nu=rho=1.;
  n=x_num*y_num;
  t=0.;
  u=(double*)malloc(x_num*y_num*sizeof(double));
  v=(double*)malloc(x_num*y_num*sizeof(double));
  p=(double*)malloc(x_num*y_num*sizeof(double));
  uvp_poiseuille(nu,rho,n,x,y,t,u,v,p);
  strcpy(header,"poiseuille");
  s=5.;
  ns2de_V2D(header,n,x,y,u,v,p,s);
  free(p);
  free(u);
  free(v);
  free(x);
  free(y);
  }

void uvp_spiral(double nu,double rho,int n,double x[],double y[],double t,double u[],double v[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
    u[i]=(1.+nu*t)*2.*x[i]*x[i]*(x[i]-1.)*(x[i]-1.)*y[i]*(2.*y[i]-1.)*(y[i]-1.);
    v[i]=-(1.+nu*t)*2.*x[i]*(2.*x[i]-1.)*(x[i]-1.)*y[i]*y[i]*(y[i]-1.)*(y[i]-1.);
    p[i]=y[i];
    }
  }

void spiral()
  {
  int n,x_num=21,y_num=21;
  double nu,rho,s,t,x_hi,x_lo,y_hi,y_lo,*p,*u,*v,*x,*y;
  printf("Spiral Flow:\n");
  printf("Generate a velocity field on a regular grid.\n");
  x_lo=y_lo=0.;
  x_hi=y_hi=1.;
  x=(double*)malloc(x_num*y_num*sizeof(double));
  y=(double*)malloc(x_num*y_num*sizeof(double));
  grid_2d(x_num,x_lo,x_hi,y_num,y_lo,y_hi,x,y);
  nu=rho=1.;
  n=x_num*y_num;
  t=0.;
  u=(double*)malloc(x_num*y_num*sizeof(double));
  v=(double*)malloc(x_num*y_num*sizeof(double));
  p=(double*)malloc(x_num*y_num*sizeof(double));
  uvp_spiral(nu,rho,n,x,y,t,u,v,p);
  strcpy(header,"spiral");
  s=5.;
  ns2de_V2D(header,n,x,y,u,v,p,s);
  free(p);
  free(u);
  free(v);
  free(x);
  free(y);
  }

void uvp_taylor(double nu,double rho,int n,double x[],double y[],double t,double u[],double v[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
    u[i]=-cos(M_PI*x[i])*sin(M_PI*y[i]);
    v[i]=sin(M_PI*x[i])*cos(M_PI*y[i]);
    p[i]=-0.25*rho*(cos(2.*M_PI*x[i])+cos(2.*M_PI*y[i]));
    u[i]=u[i]*exp(-2.*M_PI*M_PI*nu*t);
    v[i]=v[i]*exp(-2.*M_PI*M_PI*nu*t);
    p[i]=p[i]*exp(-4.*M_PI*M_PI*nu*t);
    }
  }

void taylor()
  {
  int n,x_num=21,y_num=21;
  double nu,rho,s,t,x_hi,x_lo,y_hi,y_lo,*p,*u,*v,*x,*y;
  printf("Taylor Flow:\n");
  printf("Generate a velocity field on a regular grid.\n");
  x_lo=0.5;
  x_hi=2.5;
  y_lo=0.5;
  y_hi=2.5;
  x=(double*)malloc(x_num*y_num*sizeof(double));
  y=(double*)malloc(x_num*y_num*sizeof(double));
  grid_2d(x_num,x_lo,x_hi,y_num,y_lo,y_hi,x,y);
  nu=rho=1.;
  n=x_num*y_num;
  t=0.;
  u=(double*)malloc(x_num*y_num*sizeof(double));
  v=(double*)malloc(x_num*y_num*sizeof(double));
  p=(double*)malloc(x_num*y_num*sizeof(double));
  uvp_taylor(nu,rho,n,x,y,t,u,v,p);
  strcpy(header,"taylor");
  s=0.10;
  ns2de_V2D(header,n,x,y,u,v,p,s);
  free(p);
  free(u);
  free(v);
  free(x);
  free(y);
  }

void uvp_vortex(double nu,double rho,int n,double x[],double y[],double t,double u[],double v[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
    u[i]=-cos(M_PI*x[i])*sin(M_PI*y[i]);
    v[i]=sin(M_PI*x[i])*cos(M_PI*y[i]);
    p[i]=-0.25*rho*(cos(2.*M_PI*x[i])+cos(2.*M_PI*y[i]));
    }
  }

void vortex()
  {
  int n,x_num=21,y_num=21;
  double nu,rho,s,t,x_hi,x_lo,y_hi,y_lo,*p,*u,*v,*x,*y;
  printf("Vortex Flow:\n");
  printf("Generate a velocity field on a regular grid.\n");
  x_lo=0.5;
  x_hi=1.5;
  y_lo=0.5;
  y_hi=1.5;
  x=(double*)malloc(x_num*y_num*sizeof(double));
  y=(double*)malloc(x_num*y_num*sizeof(double));
  grid_2d(x_num,x_lo,x_hi,y_num,y_lo,y_hi,x,y);
  nu=rho=1.;
  n=x_num*y_num;
  t=0.;
  u=(double*)malloc(x_num*y_num*sizeof(double));
  v=(double*)malloc(x_num*y_num*sizeof(double));
  p=(double*)malloc(x_num*y_num*sizeof(double));
  uvp_vortex(nu,rho,n,x,y,t,u,v,p);
  strcpy(header,"vortex");
  s=0.10;
  ns2de_V2D(header,n,x,y,u,v,p,s);
  free(p);
  free(u);
  free(v);
  free(x);
  free(y);
  }

int main(int argc,char**argv,char**envp)
  {
  printf("NS2DE library functions\n");
  lukas();
  poiseuille();
  spiral();
  taylor();
  vortex();
  return(0);
  }
