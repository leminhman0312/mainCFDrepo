/* compare methods of integration over an element,
   including Green's Lemma around the perimeter
   and Gauss Quadrature (1D and 2D)
   dudley.benton@gmail.com
 */
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#define y1 yyy
#include <math.h>
#undef y1

double fmin(double x,double y)
  {
  if(x<y)
    return(x);
  return(y);
  }

double fmax(double x,double y)
  {
  if(x>y)
    return(x);
  return(y);
  }

typedef struct{double x,y;}POLY;

double PolygonArea(POLY*p,int n)
  {
  int i,j;
  double A;
  for(i=n-1,A=j=0;j<n;i=j++)
    A+=(p[j].x-p[i].x)*(p[i].y+p[j].y);
  return(A/2.);
  }

double GreensLemma(POLY*poly,int n,double f(double,double))
  {
  static double A[]={0.339981043584856,0.861136311594053};
  static double W[]={0.652145154862546,0.347854845137454};
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
    for(j=0;j<sizeof(A)/sizeof(double);j++)
      {
      Xa=Xb-A[j]*dX;
      Ya=Yb-A[j]*dY;
      S+=W[j]*dX*f(Xa,Ya);
      Xa=Xb+A[j]*dX;
      Ya=Yb+A[j]*dY;
      S+=W[j]*dX*f(Xa,Ya);
      }
    }
  return(S);
  }

double GQ2D(double f(double,double),double g1(double),double g2(double),double a,double b)
  {
  int i,j,k;
  double cx,cy,dx,dy,q,w,x1,x2,y;
  static double A[]={0.125233408511,0.367831498998,0.587317954287,
                     0.769902674194,0.904117256370,0.981560634247};
  static double W[]={0.249147045813,0.233492536538,0.203167426723,
                     0.160078328543,0.106939325995,0.047175336387};
  cy=(a+b)/2;
  dy=(b-a)/2;
  for(q=i=0;i<sizeof(A)/sizeof(double);i++)
    {
    for(k=-1;k<=1;k+=2)
      {
      y=cy+k*dy*A[i];
      x1=g1(y);
      x2=g2(y);
      cx=(x1+x2)/2;
      dx=(x2-x1)/2;
      w=dx*W[i];
      for(j=0;j<sizeof(A)/sizeof(double);j++)
        q+=w*W[j]*(f(cx-dx*A[j],y)+f(cx+dx*A[j],y));
      }
    }
  return(q*dy);
  }

double mu,rho;

double dpdt(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  double dudx,dvdy,p,u,v;
  p=p1+p2*x+p3*y;
  u=u1+u2*x+u3*y+u4*x*x+u5*x*y+u6*y*y;
  v=v1+v2*x+v3*y+v4*x*x+v5*x*y+v6*y*y;
  dudx=u2+2.*u4*x+u5*y;
  dvdy=v3+v5*x+2.*v6*y;
  return(-u*p2-v*p3-p*(dudx+dvdy));
  }

double dpdtx(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-p2*(u1*x+u2*x*x/2.+u3*y*x+u4*x*x*x/3.+u5*x*x*y/2.+u6*y*y*x)
         -p3*(v1*x+v2*x*x/2.+v3*y*x+v4*x*x*x/3.+v5*x*x*y/2.+v6*y*y*x)
         -p2*(v5+2.*u4)*x*x*x/3.-((p1+p3*y)*(v5+2.*u4)
         +p2*(u2+2.*v6*y+u5*y+v3))*x*x/2.
         -(p1+p3*y)*(u2+2.*v6*y+u5*y+v3)*x);
  }

double dpdty(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-p2*(u1*y+u2*x*y+u3*y*y/2.+u4*x*x*y+u5*x*y*y/2.+u6*y*y*y/3.)
         -p3*(v1*y+v2*x*y+v3*y*y/2.+v4*x*x*y+v5*x*y*y/2.+v6*y*y*y/3.)
         -p3*(u5+2.*v6)*y*y*y/3.-((p1+p2*x)*(u5+2.*v6)
         +p3*(u2+2.*u4*x+v3+v5*x))*y*y/2.
         -(p1+p2*x)*(u2+2.*u4*x+v3+v5*x)*y);
  }

double dudt(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-(u1+u2*x+u3*y+u4*x*x+u5*x*y+u6*y*y)*(u2+2.*u4*x+u5*y)
         -(v1+v2*x+v3*y+v4*x*x+v5*x*y+v6*y*y)*(u3+u5*x+2.*u6*y)
         -p2/rho+(2.*u4+2.*u6)*mu/rho);
  }

double dudtx(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-u4*u4*x*x*x*x/2.-(u5*y+u2)*u4*x*x*x
    -(2.*(u1+u6*y*y+u3*y)*u4+pow(u5*y+u2,2.))*x*x/2.
    -(u1+u6*y*y+u3*y)*(u5*y+u2)*x-v4*u5*x*x*x*x/4.
    -((v5*y+v2)*u5+v4*(u3+2.*u6*y))*x*x*x/3.
    -((v1+v6*y*y+v3*y)*u5+(v5*y+v2)*(u3+2.*u6*y))*x*x/2.
    -(v1+v6*y*y+v3*y)*(u3+2.*u6*y)*x-p2/rho*x+(2.*u4+2.*u6)*mu/rho*x);
  }

double dudty(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-u6*u5*y*y*y*y/4.-((u3+u5*x)*u5+u6*(u2+2.*u4*x))*y*y*y/3.
    -((u1+u2*x+u4*x*x)*u5+(u3+u5*x)*(u2+2.*u4*x))*y*y/2.
    -(u1+u2*x+u4*x*x)*(u2+2.*u4*x)*y-v6*u6*y*y*y*y/2.
    -(2.*(v3+v5*x)*u6+v6*(u3+u5*x))*y*y*y/3.
    -(2.*(v1+v2*x+v4*x*x)*u6+(v3+v5*x)*(u3+u5*x))*y*y/2.
    -(v1+v2*x+v4*x*x)*(u3+u5*x)*y-p2/rho*y+(2.*u4+2.*u6)*mu/rho*y);
  }

double dvdt(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-(u1+u2*x+u3*y+u4*x*x+u5*x*y+u6*y*y)*(v2+2.*v4*x+v5*y)
         -(v1+v2*x+v3*y+v4*x*x+v5*x*y+v6*y*y)*(v3+v5*x+2.*v6*y)
         -p3/rho+(2.*v4+2.*v6)*mu/rho);
  }

double dvdtx(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-u4*v4*x*x*x*x/2.-(2.*(u5*y+u2)*v4+u4*(v5*y+v2))*x*x*x/3.
    -(2.*(u1+u6*y*y+u3*y)*v4+(u5*y+u2)*(v5*y+v2))*x*x/2.
    -(u1+u6*y*y+u3*y)*(v5*y+v2)*x-v4*v5*x*x*x*x/4.
    -((v5*y+v2)*v5+v4*(v3+2.*v6*y))*x*x*x/3.
    -((v1+v6*y*y+v3*y)*v5+(v5*y+v2)*(v3+2.*v6*y))*x*x/2.
    -(v1+v6*y*y+v3*y)*(v3+2.*v6*y)*x-p3/rho*x+(2.*v4+2.*v6)*mu/rho*x);
  }

double dvdty(double x,double y,double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return(-u6*v5*y*y*y*y/4.-((u3+u5*x)*v5+u6*(v2+2.*v4*x))*y*y*y/3.
    -((u1+u2*x+u4*x*x)*v5+(u3+u5*x)*(v2+2.*v4*x))*y*y/2.
    -(u1+u2*x+u4*x*x)*(v2+2.*v4*x)*y-v6*v6*y*y*y*y/2.
    -(v3+v5*x)*v6*y*y*y-(2.*(v1+v2*x+v4*x*x)*v6+pow(v3+v5*x,2.))*y*y/2.
    -(v1+v2*x+v4*x*x)*(v3+v5*x)*y-p3/rho*y+(2.*v4+2.*v6)*mu/rho*y);
  }

double C(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*
    (3.*p2*u4*x3*x3+p3*u5*y3*y3+p2*u6*y1*y1
    +6.*p2*u1+4.*p1*v6*y1+p3*v4*x3*x3+p3*u5*y2*y2
    +p3*v5*y1*x2+2.*p3*v5*y1*x1+2.*p3*v5*y2*x2
    +p3*v5*y2*x1+p3*v4*x2*x1+3.*p3*v6*y1*y2
    +p3*u5*y1*y2+p2*u5*y1*x2+2.*p2*u5*y1*x1
    +2.*p2*u5*y2*x2+p2*u5*y2*x1+p2*v6*y1*x2
    +2.*p2*v6*y1*x1+2.*p2*v6*y2*x2+p2*v6*y2*x1
    +p3*u4*y1*x2+2.*p3*u4*y1*x1+2.*p3*u4*y2*x2
    +p3*u4*y2*x1+3.*p2*u4*x2*x1+p2*u6*y1*y2
    +p2*v5*x2*x1+p3*v4*x2*x2+p3*v4*x1*x1
    +3.*p3*v6*y1*y1+3.*p3*v6*y2*y2+p3*u5*y1*y1
    +3.*p2*u4*x2*x2+3.*p2*u4*x1*x1+p2*u6*y2*y2
    +p2*v5*x2*x2+p2*v5*x1*x1+2.*p1*u5*y1
    +2.*p1*u5*y2+4.*p1*v6*y2+2.*p3*u2*y1
    +2.*p3*u2*y2+4.*p3*v3*y1+4.*p3*v3*y2
    +3.*p3*v6*y3*y3+p3*v5*y2*x3+2.*p3*v5*y3*x3
    +p3*v5*y3*x2+p3*v4*x3*x2+3.*p3*v6*y2*y3
    +p3*u5*y2*y3+p2*u5*y2*x3+2.*p2*u5*y3*x3
    +p2*u5*y3*x2+p2*v6*y2*x3+2.*p2*v6*y3*x3
    +p2*v6*y3*x2+p3*u4*y2*x3+2.*p3*u4*y3*x3
    +p3*u4*y3*x2+3.*p2*u4*x3*x2+p2*u6*y2*y3
    +p2*v5*x3*x2+p2*u6*y3*y3+p2*v5*x3*x3
    +2.*p1*u5*y3+4.*p1*v6*y3+2.*p3*u2*y3
    +4.*p3*v3*y3+p3*v5*y3*x1+p3*v5*y1*x3
    +p3*v4*x1*x3+3.*p3*v6*y3*y1+p3*u5*y3*y1
    +p2*u5*y3*x1+p2*u5*y1*x3+p2*v6*y3*x1
    +p2*v6*y1*x3+p3*u4*y3*x1+p3*u4*y1*x3
    +3.*p2*u4*x1*x3+p2*u6*y3*y1+p2*v5*x1*x3
    +4.*x3*p2*u2+2.*p2*v3*x2+4.*p2*u2*x2
    +2.*p1*v5*x2+4.*p1*u4*x2+2.*p2*u3*y3
    +2.*p3*v2*x2+2.*x1*p2*v3+2.*x1*p3*v2
    +4.*x1*p1*u4+2.*x1*p1*v5+4.*x1*p2*u2
    +2.*y2*p2*u3+2.*y1*p2*u3+4.*x3*p1*u4
    +2.*x3*p2*v3+2.*x3*p1*v5+2.*x3*v2*p3
    +6.*v3*p1+6.*u2*p1+6.*v1*p3)/12.);
  }

double dCdp1(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*x1*v5+2.*u5*y3
    +6.*v3+6.*u2+4.*u4*x2+2.*u5*y1+4.*x1*u4+2.*u5*y2+4.*v6*y1
    +2.*x3*v5+4.*x3*u4+2.*v5*x2+4.*v6*y3+4.*v6*y2)/12.);
  }

double dCdp2(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(6.*u1+u5*y1*x2+u5*y2*x1
    +v6*y1*x2+v6*y2*x1+u6*y1*y2+v5*x2*x1+u5*y2*x3+u5*y3*x2+v6*y2*x3+v6*y3*x2
    +u6*y2*y3+v5*x3*x2+u5*y3*x1+u5*y1*x3+v6*y3*x1+v6*y1*x3+u6*y3*y1+v5*x1*x3
    +2.*u5*y1*x1+2.*u5*y2*x2+2.*v6*y1*x1+2.*v6*y2*x2+3.*u4*x2*x1+2.*u5*y3*x3
    +2.*v6*y3*x3+3.*u4*x3*x2+3.*u4*x1*x3+4.*x1*u2+v5*x3*x3+u6*y3*y3+v5*x1*x1
    +2.*y1*u3+4.*x3*u2+u6*y1*y1+u6*y2*y2+v5*x2*x2+3.*u4*x3*x3+2.*v3*x2
    +4.*u2*x2+2.*u3*y3+2.*x1*v3+2.*y2*u3+2.*x3*v3+3.*u4*x2*x2+3.*u4*x1*x1)/12.);
  }

double dCdp3(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(v5*y2*x1+3.*v6*y3*y3
    +2.*u2*y1+4.*v3*y1+4.*v3*y2+2.*u2*y3+4.*v3*y3+2.*v2*x2+2.*x3*v2
    +3.*v6*y1*y1+2.*x1*v2+u5*y3*y3+v4*x3*x3+u5*y2*y2+v4*x2*x2+v4*x1*x1
    +u5*y1*y1+2.*u2*y2+6.*v1+v5*y2*x3+v5*y3*x2+u4*y1*x2+u4*y3*x1+u4*y1*x3
    +2.*v5*y1*x1+2.*v5*y2*x2+3.*v6*y1*y2+2.*u4*y1*x1+2.*u4*y2*x2+2.*v5*y3*x3
    +3.*v6*y2*y3+2.*u4*y3*x3+3.*v6*y3*y1+v5*y1*x3+v4*x1*x3+u5*y3*y1+v4*x3*x2
    +u5*y2*y3+u4*y2*x3+u4*y3*x2+v5*y3*x1+v5*y1*x2+v4*x2*x1+3.*v6*y2*y2
    +u4*y2*x1+u5*y1*y2)/12.);
  }

double dCdu1(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*p2/2.);
  }

double dCdu2(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*p3*y1+4.*x3*p2+6.*p1
    +4.*x1*p2+4.*p2*x2+2.*p3*y3+2.*p3*y2)/12.);
  }

double dCdu3(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*y1*p2+2.*y2*p2
    +2.*p2*y3)/12.);
  }

double dCdu4(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(4.*p1*x2+p3*y1*x2
    +p3*y2*x1+p3*y2*x3+p3*y3*x2+p3*y3*x1+p3*y1*x3+4.*x1*p1+3.*p2*x1*x1
    +3.*p2*x2*x2+4.*x3*p1+3.*p2*x3*x3+2.*p3*y1*x1+2.*p3*y2*x2+3.*p2*x2*x1
    +2.*p3*y3*x3+3.*p2*x3*x2+3.*p2*x1*x3)/12.);
  }

double dCdu5(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*p1*y3+2.*p1*y1+p3*y3*y3
    +p3*y2*y2+p3*y1*y2+p2*y1*x2+p2*y2*x1+p3*y1*y1+p3*y2*y3+p2*y2*x3
    +p2*y3*x2+p3*y3*y1+p2*y3*x1+p2*y1*x3+2.*p1*y2+2.*p2*y1*x1+2.*p2*y2*x2
    +2.*p2*y3*x3)/12.);
  }

double dCdu6(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(p2*y1*y1+p2*y1*y2
    +p2*y2*y2+p2*y2*y3+p2*y3*y3+p2*y3*y1)/12.);
  }

double dCdv1(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*p3/2.);
  }

double dCdv2(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*x1*p3+2.*x3*p3
    +2.*p3*x2)/12.);
  }

double dCdv3(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(6.*p1+2.*x3*p2+2.*x1*p2
    +2.*p2*x2+4.*p3*y3+4.*p3*y2+4.*p3*y1)/12.);
  }

double dCdv4(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(p3*x3*x3+p3*x2*x1+p3*x2*x2
    +p3*x1*x1+p3*x3*x2+p3*x1*x3)/12.);
  }

double dCdv5(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(2.*x1*p1+p3*y1*x2
    +p3*y2*x1+p2*x2*x1+p2*x2*x2+p2*x1*x1+p3*y2*x3+p3*y3*x2+p2*x3*x2
    +p2*x3*x3+p3*y3*x1+p3*y1*x3+2.*x3*p1+2.*p1*x2+p2*x1*x3+2.*p3*y1*x1
    +2.*p3*y2*x2+2.*p3*y3*x3)/12.);
  }

double dCdv6(double x1,double x2,double x3,double y1,double y2,double y3,
  double p1,double p2,double p3,
  double u1,double u2,double u3,double u4,double u5,double u6,
  double v1,double v2,double v3,double v4,double v5,double v6)
  {
  return((y3*x2-y1*x2+y1*x3+y2*x1-y2*x3-y3*x1)*(3.*p3*y3*y3+3.*p3*y2*y2
    +p2*y1*x2+p2*y2*x1+p2*y2*x3+p2*y3*x2+p2*y3*x1+p2*y1*x3+3.*p3*y1*y1
    +4.*p1*y1+4.*p1*y3+4.*p1*y2+3.*p3*y1*y2+2.*p2*y1*x1+2.*p2*y2*x2
    +3.*p3*y2*y3+2.*p2*y3*x3+3.*p3*y3*y1)/12.);
  }

int hi,le,lo,rt;
POLY poly[3],ylop[3];

double g1(double y)
  {
  double g,x;
  g=poly[rt].x;
  if((poly[0].y-y)*(poly[1].y-y)<=0.)
    {
    if(fabs(poly[0].y-poly[1].y)>DBL_EPSILON)
      x=poly[0].x+(poly[1].x-poly[0].x)*(y-poly[0].y)/(poly[1].y-poly[0].y);
    else
      x=(poly[0].x+poly[1].x)/2.;
    g=fmin(g,x);
    }
  if((poly[1].y-y)*(poly[2].y-y)<=0.)
    {
    if(fabs(poly[1].y-poly[2].y)>DBL_EPSILON)
      x=poly[1].x+(poly[2].x-poly[1].x)*(y-poly[1].y)/(poly[2].y-poly[1].y);
    else
      x=(poly[1].x+poly[2].x)/2.;
    g=fmin(g,x);
    }
  if((poly[2].y-y)*(poly[0].y-y)<=0.)
    {
    if(fabs(poly[2].y-poly[0].y)>DBL_EPSILON)
      x=poly[2].x+(poly[0].x-poly[2].x)*(y-poly[2].y)/(poly[0].y-poly[2].y);
    else
      x=(poly[2].x+poly[0].x)/2.;
    g=fmin(g,x);
    }
  return(g);
  }

double g2(double y)
  {
  double g,x;
  g=poly[le].x;
  if((poly[0].y-y)*(poly[1].y-y)<=0.)
    {
    if(fabs(poly[0].y-poly[1].y)>DBL_EPSILON)
      x=poly[0].x+(poly[1].x-poly[0].x)*(y-poly[0].y)/(poly[1].y-poly[0].y);
    else
      x=(poly[0].x+poly[1].x)/2.;
    g=fmax(g,x);
    }
  if((poly[1].y-y)*(poly[2].y-y)<=0.)
    {
    if(fabs(poly[1].y-poly[2].y)>DBL_EPSILON)
      x=poly[1].x+(poly[2].x-poly[1].x)*(y-poly[1].y)/(poly[2].y-poly[1].y);
    else
      x=(poly[1].x+poly[2].x)/2.;
    g=fmax(g,x);
    }
  if((poly[2].y-y)*(poly[0].y-y)<=0.)
    {
    if(fabs(poly[2].y-poly[0].y)>DBL_EPSILON)
      x=poly[2].x+(poly[0].x-poly[2].x)*(y-poly[2].y)/(poly[0].y-poly[2].y);
    else
      x=(poly[2].x+poly[0].x)/2.;
    g=fmax(g,x);
    }
  return(g);
  }

void PrepareQuadrature()
  {
  lo=0;
  if(poly[1].y<poly[0].y)
    {
    lo=1;
    if(poly[2].y<poly[1].y)
      lo=2;
    }
  else if(poly[2].y<poly[0].y)
    lo=2;
  hi=0;
  if(poly[1].y>poly[0].y)
    {
    hi=1;
    if(poly[2].y>poly[1].y)
      hi=2;
    }
  else if(poly[2].y>poly[0].y)
    hi=2;
  if(lo==hi)
    {
    fprintf(stderr,"y=%lG,%lG,%lG\n",poly[0].y,poly[1].y,poly[2].y);
    exit(1);
    }
  le=0;
  if(poly[1].x<poly[0].x)
    {
    le=1;
    if(poly[2].x<poly[1].x)
      le=2;
    }
  else if(poly[2].x<poly[0].x)
    le=2;
  rt=0;
  if(poly[1].x>poly[0].x)
    {
    rt=1;
    if(poly[2].x>poly[1].x)
      rt=1;
    }
  else if(poly[2].x>poly[0].x)
    rt=2;
  if(le==rt)
    {
    fprintf(stderr,"x=%lG,%lG,%lG\n",poly[0].x,poly[1].x,poly[2].x);
    exit(1);
    }
  }

int randbetween(int lo,int hi)
  {
  return(lo+rand()%(hi-lo+1));
  }

void RandomTriangle()
  {
  int i;
  poly[0].x=randbetween(-5000,5000)/1000.;
  poly[0].y=randbetween(-5000,5000)/1000.;
  poly[1].x=randbetween(-5000,5000)/1000.;
  poly[1].y=randbetween(-5000,5000)/1000.;
  poly[2].x=randbetween(-5000,5000)/1000.;
  poly[2].y=randbetween(-5000,5000)/1000.;
  if(PolygonArea(poly,3)<0.)
    {
    POLY p;
    p=poly[1];
    poly[1]=poly[2];
    poly[2]=p;
    }
  for(i=0;i<3;i++)
    {
    ylop[i].x=poly[i].y;
    ylop[i].y=poly[i].x;
    }
  }

double x1,x2,x3,y1,y2,y3,p1,p2,p3,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6;

double f(double x,double y)
  {
  return(dpdt(x,y,p1,p2,p3,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6));
  }

double fx(double x,double y)
  {
  return(dpdty(x,y,p1,p2,p3,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6));
  }

double fy(double y,double x)
  {
  return(-dpdtx(x,y,p1,p2,p3,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6));
  }

int main(int argc,char**argv,char**envp)
  {
  int i,n=10;
  double A;
  printf("testing conservation of mass\n");
  printf("  GQ2D LemmaX LemmaY Analyt\n");
  for(i=0;i<n;i++)
    {
    rho=randbetween(5000,10000)/10000.;
    mu =randbetween(5000,10000)/10000.;
    RandomTriangle();
    x1=poly[0].x;
    x2=poly[1].x;
    x3=poly[2].x;
    y1=poly[0].y;
    y2=poly[1].y;
    y3=poly[2].y;
    p1=randbetween(-5000,5000)/1000.;
    p2=randbetween(-5000,5000)/1000.;
    p3=randbetween(-5000,5000)/1000.;
    u1=randbetween(-5000,5000)/1000.;
    u2=randbetween(-5000,5000)/1000.;
    u3=randbetween(-5000,5000)/1000.;
    u4=randbetween(-5000,5000)/1000.;
    u5=randbetween(-5000,5000)/1000.;
    u6=randbetween(-5000,5000)/1000.;
    v1=randbetween(-5000,5000)/1000.;
    v2=randbetween(-5000,5000)/1000.;
    v3=randbetween(-5000,5000)/1000.;
    v4=randbetween(-5000,5000)/1000.;
    v5=randbetween(-5000,5000)/1000.;
    v6=randbetween(-5000,5000)/1000.;
    A=PolygonArea(poly,3);
    PrepareQuadrature();
    printf("%6.1lf %6.1lf %6.1lf %6.1lf\n",
      GQ2D(f,g1,g2,poly[lo].y,poly[hi].y)/A,
      GreensLemma(poly,3,fx)/A,GreensLemma(ylop,3,fy)/A,
      C(x1,x2,x3,y1,y2,y3,p1,p2,p3,u1,u2,u3,u4,u5,u6,v1,v2,v3,v4,v5,v6)/A);
    }
  return(0);
  }
