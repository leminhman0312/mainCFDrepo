/* apply Green's Lemma over a polygon
   dudley.benton@gmail.com
 */

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

typedef struct{double x,y;}POLY;

double GreensLemma(POLY*poly,int n,double a,double b)
  {
  static double A4[]={0.339981043584856,0.861136311594053};
  static double W4[]={0.652145154862546,0.347854845137454};
  double dX,dY,S,Xa,Xb,X1,X2,Ya,Yb,Y1,Y2;
  int i,j;
  X2=poly[n-1].x;
  Y2=poly[n-1].y;
  for(S=i=0;i<n;i++)
    {
    X1=X2;
    Y1=Y2;
    X2=poly[i].x;
    Y2=poly[i].y;
    dX=(X2-X1)/2.;
    dY=(Y2-Y1)/2.;
    Xb=(X2+X1)/2.;
    Yb=(Y2+Y1)/2.;
    for(j=0;j<sizeof(A4)/sizeof(A4[0]);j++)
      {
      Xa=pow(Xb-A4[j]*dX,a);
      Ya=pow(Yb-A4[j]*dY,b+1.);
      S+=W4[j]*dX*Xa*Ya;
      Xa=pow(Xb+A4[j]*dX,a);
      Ya=pow(Yb+A4[j]*dY,b+1.);
      S+=W4[j]*dX*Xa*Ya;
      }
    }
  return(S/(b+1.));
  }

POLY poly[]={
  {0.00,0.00},
  {0.53,0.90},
  {1.35,0.90},
  {1.35,0.97},
  {5.19,0.97},
  {3.60,0.00}};
 
double GQ2D(double f(double,double),double g1(double),double g2(double),double a,double b)
  {
  int i,j,k;
  double cx,cy,dx,dy,q,w,x,y1,y2;
  static double A8[]={0.183434642496,0.525532409916,0.796666477414,0.960289856498};
  static double W8[]={0.362683783378,0.313706645878,0.222381034453,0.101228536290};
  cx=(a+b)/2;
  dx=(b-a)/2;
  q=0;
  for(i=0;i<4;i++)
    {
    for(k=-1;k<=1;k+=2)
      {
      x=cx+k*dx*A8[i];
      y1=g1(x);
      y2=g2(x);
      cy=(y1+y2)/2;
      dy=(y2-y1)/2;
      w=dy*W8[i];
      for(j=0;j<4;j++)
        q+=w*W8[j]*(f(x,cy-dy*A8[j])+f(x,cy+dy*A8[j]));
      }
    }
  return(q*dx);
  }

double a,b;

double f(double y,double x)
  {
  return(pow(x,a)*pow(y,b));
  }

double g1(double y)
  {
  if(y<0.9)
    return(0.53*y/0.9);
  return(1.35);
  }

double g2(double y)
  {
  return(3.6+1.59*y/0.97);
  }

int main(int argc,char**argv,char**envp)
  {
  printf("illustrating Green's Lemma\n");
  printf(" a   b   1D   2D\n");
  for(a=0.;a<3.01;a+=0.5)
    for(b=0.;b<3.01;b+=0.5)
      printf("%3.1lf %3.1lf %4.1lf %4.1lf\n",a,b,GreensLemma(poly,6,a,b),GQ2D(f,g1,g2,0.,0.97));
  return(0);
  }
