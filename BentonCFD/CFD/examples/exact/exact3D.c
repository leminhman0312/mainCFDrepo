/******************************************************************************/
/*
  Licensing:
    This code is distributed under the GNU LGPL license.
  Author:
    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.
  Reference:
    William Cody,"Rational Chebyshev approximations for the error function",Mathematics of Computation,1969,pages 631-638.
  Reference:
    Paul Bratley,Bennett Fox,Linus Schrage,A Guide to Simulation,Second Edition,Springer,1987,ISBN: 0387964673,LC: QA76.9.C65.B73.
    Bennett Fox,Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,ACM Transactions on Mathematical Software,Volume 12,Number 4,December 1986,pages 362-376.
    Pierre L'Ecuyer,Random Number Generation,in Handbook of Simulation,edited by Jerry Banks,Wiley,1998,ISBN: 0471134031,LC: T57.62.H37.
    Peter Lewis,Allen Goodman,James Miller,A Pseudo-Random Number Generator for the System/360,IBM Systems Journal,Volume 8,Number 2,1969,pages 136-143.
    C Ross Ethier,David Steinman,Exact fully 3D Navier-Stokes solutions for benchmarking,International Journal for Numerical Methods in Fluids,Volume 19,Number 5,March 1994,pages 369-375.
*/

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>


double r8_erf(double x)
  {
   double a[5]={3.16112374387056560,1.13864154151050156E+02,3.77485237685302021E+02,3.20937758913846947E+03,1.85777706184603153E-01};
   double b[4]={2.36012909523441209E+01,2.44024637934444173E+02,1.28261652607737228E+03,2.84423683343917062E+03};
   double c[9]={5.64188496988670089E-01,8.88314979438837594,6.61191906371416295E+01,2.98635138197400131E+02,8.81952221241769090E+02,1.71204761263407058E+03,2.05107837782607147E+03,1.23033935479799725E+03,2.15311535474403846E-08};
   double d[8]={1.57449261107098347E+01,1.17693950891312499E+02,5.37181101862009858E+02,1.62138957456669019E+03,3.29079923573345963E+03,4.36261909014324716E+03,3.43936767414372164E+03,1.23033935480374942E+03};
  double del;
  double erfx;
  int i;
  double p[6]={3.05326634961232344E-01,3.60344899949804439E-01,1.25781726111229246E-01,1.60837851487422766E-02,6.58749161529837803E-04,1.63153871373020978E-02};
  double q[5]={2.56852019228982242,1.87295284992346047,5.27905102951428412E-01,6.05183413124413191E-02,2.33520497626869185E-03};
  double sqrpi=0.56418958354775628695;
  double thresh=0.46875;
  double xbig=26.543;
  double xabs;
  double xden;
  double xnum;
  double xsmall=1.11E-16;
  double xsq;
  xabs=fabs(x);
/* Evaluate ERF(X) for |X|<=0.46875. */
  if(xabs<=thresh)
    {
    if(xsmall<xabs)
      xsq=xabs*xabs;
    else
      xsq=0.;
    xnum=a[4]*xsq;
    xden=xsq;
    for(i=0;i<3;i++)
      {
      xnum=(xnum+a[i])*xsq;
      xden=(xden+b[i])*xsq;
      }
    erfx=x*(xnum+a[3])/(xden+b[3]);
    }
/* Evaluate ERFC(X) for 0.46875<=|X|<=4.0. */
  else if(xabs<=4.)
    {
    xnum=c[8]*xabs;
    xden=xabs;
    for(i=0;i<7;i++)
      {
      xnum=(xnum+c[i])*xabs;
      xden=(xden+d[i])*xabs;
      }
    erfx=(xnum+c[7])/(xden+d[7]);
    xsq=(double)((int)((xabs*16.)/16.));
    del=(xabs-xsq)*(xabs+xsq);
    erfx=exp(-xsq*xsq)*exp(-del)*erfx;
    erfx=(0.5-erfx)+0.5;
    if(x<0.)
      erfx=-erfx;
    }
/* Evaluate ERFC(X) for 4.0<|X|. */
  else
    {
    if(xbig<=xabs)
      {
      if(0.0<x)
        erfx=1.;
      else
        erfx=-1.;
      }
    else
      {
      xsq=1./(xabs*xabs);
      xnum=p[5]*xsq;
      xden=xsq;
      for(i=0;i<4;i++)
        {
        xnum=(xnum+p[i])*xsq;
        xden=(xden+q[i])*xsq;
        }
      erfx=xsq*(xnum+p[4])/(xden+q[4]);
      erfx=(sqrpi-erfx)/xabs;
      xsq=(double) ((int) ((xabs*16.)/16.));
      del=(xabs-xsq)*(xabs+xsq);
      erfx=exp(-xsq*xsq)*exp(-del)*erfx;
      erfx=(0.5-erfx)+0.5;
      if(x<0.)
        erfx=-erfx;
      }
    }
  return(erfx);
  }

double r8vec_amax(int n,double a[])
  {
  double amax;
  int i;
  amax=0.;
  for(i=0;i<n;i++)
    if(amax<fabs(a[i]))
      amax=fabs(a[i]);
  return(amax);
  }

double r8vec_amin(int n,double a[])
  {
  double amin;
  int i;
  const double r8_huge=1.79769313486231571E+308;
  amin=r8_huge;
  for(i=0;i<n;i++)
    if(fabs(a[i])<amin)
      amin=fabs(a[i]);
  return(amin);
  }

double r8vec_max(int n,double r8vec[])
  {
  int i;
  double value;
  if(n<=0)
    {
    value=0.;
    return(value);
    }
  value=r8vec[0];
  for(i=1;i<n;i++)
    if(value<r8vec[i])
      value=r8vec[i];
  return(value);
  }

double r8vec_min(int n,double r8vec[])
  {
  int i;
  double value;
  value=r8vec[0];
  for(i=1;i<n;i++)
    if(r8vec[i]<value)
      value=r8vec[i];
  return(value);
  }

double *r8vec_uniform_ab_new(int n,double a,double b,int *seed)
  {
  int i;
  const int i4_huge=2147483647;
  int k;
  double *r;
  if(*seed==0)
    {
    fprintf(stderr,"\n");
    fprintf(stderr,"R8VEC_UNIFORM_AB_NEW-Fatal error!\n");
    fprintf(stderr,"  Input value of SEED=0.\n");
    exit(1);
    }
  r=(double *) malloc(n*sizeof(double));
  for(i=0;i<n;i++)
    {
    k=*seed/127773;
    *seed=16807*(*seed-k*127773)-k*2836;
    if(*seed<0)
      *seed=*seed+i4_huge;
    r[i]=a+(b-a)*(double) (*seed)*4.656612875E-10;
    }
  return(r);
  }

double *r8vec_zeros_new(int n)
  {
  int i;
  double *r;
  r=(double *) malloc(n*sizeof(double));
  for(i=0;i<n;i++)
    r[i]=0.;
  return(r);
  }

void resid_burgers(double nu,int n,double x[],double y[],double z[],double t[],double ur[],double vr[],double wr[],double pr[])
  {
  int i;
  double px;
  double py;
  double pz;
  double u;
  double ut;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double uz;
  double uzz;
  double v;
  double vt;
  double vx;
  double vxx;
  double vy;
  double vyy;
  double vz;
  double vzz;
  double w;
  double wt;
  double wx;
  double wxx;
  double wy;
  double wyy;
  double wz;
  double wzz;
  for(i=0;i<n;i++)
    {
/* Form the functions and derivatives. */
    u=2.*x[i];
    ux=2.;
    uxx=0.;
    uy=0.;
    uyy=0.;
    uz=0.;
    uzz=0.;
    ut=0.;
    v=-2.*y[i];
    vx=0.;
    vxx=0.;
    vy=-2.;
    vyy=0.;
    vz=0.;
    vzz=0.;
    vt=0.;
    w=r8_erf(y[i]/sqrt(nu));
    wx=0.;
    wxx=0.;
    wy=2.*sqrt(1./nu/M_PI)*exp(-y[i]*y[i]/nu);
    wyy=-4.*sqrt(1./nu/M_PI)*y[i]*exp(-y[i]*y[i]/nu)/nu;
    wz=0.;
    wzz=0.;
    wt=0.;
    px=-4.*x[i];
    py=-4.*y[i];
    pz=0.;
/* Evaluate the residuals. */
    ur[i]=ut+u*ux+v*uy+w*uz+px-nu*(uxx+uyy+uzz);
    vr[i]=vt+u*vx+v*vy+w*vz+py-nu*(vxx+vyy+vzz);
    wr[i]=wt+u*wx+v*wy+w*wz+pz-nu*(wxx+wyy+wzz);
    pr[i]=ux+vy+wz;
    }
  }

void resid_ethier(double a,double d,int n,double x[],double y[],double z[],double t[],double ur[],double vr[],double wr[],double pr[])
  {
  double cxy;
  double cyz;
  double czx;
  double e2t;
  double e2x;
  double e2y;
  double e2z;
  double e4t;
  double ex;
  double exy;
  double ey;
  double eyz;
  double ez;
  double ezx;
  int i;
  double px;
  double py;
  double pz;
  double sxy;
  double syz;
  double szx;
  double u;
  double ut;
  double ux;
  double uxx;
  double uy;
  double uyy;
  double uz;
  double uzz;
  double v;
  double vt;
  double vx;
  double vxx;
  double vy;
  double vyy;
  double vz;
  double vzz;
  double w;
  double wt;
  double wx;
  double wxx;
  double wy;
  double wyy;
  double wz;
  double wzz;
/* Make some temporaries. */
  for(i=0;i<n;i++)
    {
    ex=exp(a*x[i]);
    ey=exp(a*y[i]);
    ez=exp(a*z[i]);
    e2x=exp(2.*a*x[i]);
    e2y=exp(2.*a*y[i]);
    e2z=exp(2.*a*z[i]);
    e2t=exp(-d*d*t[i]);
    e4t=exp(-2.*d*d*t[i]);
    exy=exp(a*(x[i]+y[i]));
    eyz=exp(a*(y[i]+z[i]));
    ezx=exp(a*(z[i]+x[i]));
    sxy=sin(a*x[i]+d*y[i]);
    syz=sin(a*y[i]+d*z[i]);
    szx=sin(a*z[i]+d*x[i]);
    cxy=cos(a*x[i]+d*y[i]);
    cyz=cos(a*y[i]+d*z[i]);
    czx=cos(a*z[i]+d*x[i]);
/* Form the functions and derivatives. */
    u=-a*(ex*syz+ez*cxy)*e2t;
    ux=-a*(a*ex*syz-a*ez*sxy)*e2t;
    uxx=-a*(a*a*ex*syz-a*a*ez*cxy)*e2t;
    uy=-a*(a*ex*cyz-d*ez*sxy)*e2t;
    uyy=-a*(-a*a*ex*syz-d*d*ez*cxy)*e2t;
    uz=-a*(d*ex*cyz+a*ez*cxy)*e2t;
    uzz=-a*(-d*d*ex*syz+a*a*ez*cxy)*e2t;
    ut=+d*d*a*(ex*syz+ez*cxy)*e2t;
    v=-a*(ey*szx+ex*cyz)*e2t;
    vx=-a*(d*ey*czx+a*ex*cyz)*e2t;
    vxx=-a*(-d*d*ey*szx+a*a*ex*cyz)*e2t;
    vy=-a*(a*ey*szx-a*ex*syz)*e2t;
    vyy=-a*(a*a*ey*szx-a*a*ex*cyz)*e2t;
    vz=-a*(a*ey*czx-d*ex*syz)*e2t;
    vzz=-a*(-a*a*ey*szx-d*d*ex*cyz)*e2t;
    vt=+d*d*a*(ey*szx+ex*cyz)*e2t;
    w=-a*(ez*sxy+ey*czx)*e2t;
    wx=-a*(a*ez*cxy-d*ey*szx)*e2t;
    wxx=-a*(-a*a*ez*sxy-d*d*ey*czx)*e2t;
    wy=-a*(d*ez*cxy+a*ey*czx)*e2t;
    wyy=-a*(-d*d*ez*sxy+a*a*ey*czx)*e2t;
    wz=-a*(a*ez*sxy-a*ey*szx)*e2t;
    wzz=-a*(a*a*ez*sxy-a*a*ey*czx)*e2t;
    wt=+d*d*a*(ez*sxy+ey*czx)*e2t;
    px=-0.5*a*a*e4t*(+2.*a*e2x+2.*a*cxy*czx*eyz-2.*d*sxy*szx*eyz-2.*a*syz*sxy*ezx+2.*a*syz*cxy*ezx+2.*d*czx*cyz*exy+2.*a*szx*cyz*exy);
    py=-0.5*a*a*e4t*(+2.*d*cxy*czx*eyz+2.*a*sxy*czx*eyz+2.*a*e2y+2.*a*cyz*cxy*ezx-2.*d*syz*sxy*ezx-2.*a*szx*syz*exy+2.*a*szx*cyz*exy);
    pz=-0.5*a*a*e4t*(-2.*a*sxy*szx*eyz+2.*a*sxy*czx*eyz+2.*d*cyz*cxy*ezx+2.*a*syz*cxy*ezx+2.*a*e2z+2.*a*czx*cyz*exy-2.*d*szx*syz*exy);
/* Evaluate the residuals. */
    ur[i]=ut+u*ux+v*uy+w*uz+px-(uxx+uyy+uzz);
    vr[i]=vt+u*vx+v*vy+w*vz+py-(vxx+vyy+vzz);
    wr[i]=wt+u*wx+v*wy+w*wz+pz-(wxx+wyy+wzz);
    pr[i]=ux+vy+wz;
    }
  }

void uvwp_burgers(double nu,int n,double x[],double y[],double z[],double t[],double u[],double v[],double w[],double p[])
  {
  int i;
  for(i=0;i<n;i++)
    {
/* Form the functions and derivatives. */
    u[i]=2.*x[i];
    v[i]=-2.*y[i];
    w[i]=r8_erf(y[i]/sqrt(nu));
    p[i]=-2.*(x[i]*x[i]+y[i]*y[i]);
    }
  }

void uvwp_ethier(double a,double d,int n,double x[],double y[],double z[],double t[],double u[],double v[],double w[],double p[])
  {
  double cxy;
  double cyz;
  double czx;
  double e2t;
  double ex;
  double exy;
  double ey;
  double eyz;
  double ez;
  double ezx;
  int i;
  double sxy;
  double syz;
  double szx;
  for(i=0;i<n;i++)
    {
    ex=exp(a*x[i]);
    ey=exp(a*y[i]);
    ez=exp(a*z[i]);
    e2t=exp(-d*d*t[i]);
    exy=exp(a*(x[i]+y[i]));
    eyz=exp(a*(y[i]+z[i]));
    ezx=exp(a*(z[i]+x[i]));
    sxy=sin(a*x[i]+d*y[i]);
    syz=sin(a*y[i]+d*z[i]);
    szx=sin(a*z[i]+d*x[i]);
    cxy=cos(a*x[i]+d*y[i]);
    cyz=cos(a*y[i]+d*z[i]);
    czx=cos(a*z[i]+d*x[i]);
    u[i]=-a*(ex*syz+ez*cxy)*e2t;
    v[i]=-a*(ey*szx+ex*cyz)*e2t;
    w[i]=-a*(ez*sxy+ey*czx)*e2t;
    p[i]=0.5*a*a*e2t*e2t*(+ex*ex+2.*sxy*czx*eyz+ey*ey+2.*syz*cxy*ezx+ez*ez+2.*szx*cyz*exy);
    }
  }

void uvwp_burgers_test()
  {
  int n;
  double nu;
  double *p;
  int seed;
  double *t;
  double *u;
  double *v;
  double *w;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;
  nu=0.25;
  printf("\n");
  printf("UVWP_BURGERS_TEST\n");
  printf("  UVWP_BURGERS evaluates the Burgers solution.\n");
  printf("  Estimate the range of velocity and pressure\n");
  printf("  at the initial time T=0, in a region that is the\n");
  printf("  cube centered at (0,0,0) with 'radius' 1.0.\n");
  printf("  Viscosity nu=%lG\n",nu);
  n=1000;
  p=(double *) malloc(n*sizeof(double));
  u=(double *) malloc(n*sizeof(double));
  v=(double *) malloc(n*sizeof(double));
  w=(double *) malloc(n*sizeof(double));
  xyz_lo=-1.;
  xyz_hi=+1.;
  seed=123456789;
  x=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  y=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  z=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  t=r8vec_zeros_new(n);
  uvwp_burgers(nu,n,x,y,z,t,u,v,w,p);
  printf("  Minimum Maximum\n");
  printf("  U: %lG %lG\n",r8vec_min(n,u),r8vec_max(n,u));
  printf("  V: %lG %lG\n",r8vec_min(n,v),r8vec_max(n,v));
  printf("  W: %lG %lG\n",r8vec_min(n,w),r8vec_max(n,w));
  printf("  P: %lG %lG\n",r8vec_min(n,p),r8vec_max(n,p));
  free(p);
  free(t);
  free(u);
  free(v);
  free(w);
  free(x);
  free(y);
  free(z);
  }

void resid_burgers_test()
  {
  int n;
  double nu;
  double *pr;
  int seed;
  double *t;
  double *ur;
  double *vr;
  double *wr;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;
  nu=0.25;
  printf("\n");
  printf("RESID_BURGERS_TEST\n");
  printf("  RESID_BURGERS evaluates the Burgers residual.\n");
  printf("  Sample the Navier-Stokes residuals\n");
  printf("  at the initial time T=0, using a region that is\n");
  printf("  the cube centered at (0,0,0) with radius 1\n");
  printf("  Viscosity nu=%lG\n",nu);
  n=1000;
  pr=(double *) malloc(n*sizeof(double));
  ur=(double *) malloc(n*sizeof(double));
  vr=(double *) malloc(n*sizeof(double));
  wr=(double *) malloc(n*sizeof(double));
  xyz_lo=-1.;
  xyz_hi=+1.;
  seed=123456789;
  x=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  y=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  z=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  t=r8vec_zeros_new(n);
  resid_burgers(nu,n,x,y,z,t,ur,vr,wr,pr);
  printf("  Minimum Maximum\n");
  printf("  Ur: %lG %lG\n",r8vec_amin(n,ur),r8vec_amax(n,ur));
  printf("  Vr: %lG %lG\n",r8vec_amin(n,vr),r8vec_amax(n,vr));
  printf("  Wr: %lG %lG\n",r8vec_amin(n,wr),r8vec_amax(n,wr));
  printf("  Pr: %lG %lG\n",r8vec_amin(n,pr),r8vec_amax(n,pr));
  free(pr);
  free(t);
  free(ur);
  free(vr);
  free(wr);
  free(x);
  free(y);
  free(z);
  }

void uvwp_ethier_test()
  {
  double a;
  double d;
  int n;
  double *p;
  int seed;
  double *t;
  double *u;
  double *v;
  double *w;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;
  a=M_PI/4.;
  d=M_PI/2.;
  printf("\n");
  printf("UVWP_ETHIER_TEST\n");
  printf("  UVWP_ETHIER evaluates the Ethier solution.\n");
  printf("  Estimate the range of velocity and pressure\n");
  printf("  at the initial time T=0, in a region that is the\n");
  printf("  cube centered at (0,0,0) with 'radius' 1.0.\n");
  printf("  Parameter A=%lG\n",a);
  printf("  Parameter D=%lG\n",d);
  n=1000;
  p=(double *) malloc(n*sizeof(double));
  u=(double *) malloc(n*sizeof(double));
  v=(double *) malloc(n*sizeof(double));
  w=(double *) malloc(n*sizeof(double));
  xyz_lo=-1.;
  xyz_hi=+1.;
  seed=123456789;
  x=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  y=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  z=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  t=r8vec_zeros_new(n);
  uvwp_ethier(a,d,n,x,y,z,t,u,v,w,p);
  printf("  Minimum Maximum\n");
  printf("  U: %lG %lG\n",r8vec_min(n,u),r8vec_max(n,u));
  printf("  V: %lG %lG\n",r8vec_min(n,v),r8vec_max(n,v));
  printf("  W: %lG %lG\n",r8vec_min(n,w),r8vec_max(n,w));
  printf("  P: %lG %lG\n",r8vec_min(n,p),r8vec_max(n,p));
  free(p);
  free(t);
  free(u);
  free(v);
  free(w);
  free(x);
  free(y);
  free(z);
  }

void resid_ethier_test()
  {
  double a;
  double d;
  int n;
  double *pr;
  int seed;
  double *t;
  double *ur;
  double *vr;
  double *wr;
  double *x;
  double xyz_hi;
  double xyz_lo;
  double *y;
  double *z;
  a=M_PI/4.;
  d=M_PI/2.;
  printf("\n");
  printf("RESID_ETHIER_TEST\n");
  printf("  RESID_ETHIER evaluates the Ethier residual.\n");
  printf("  Sample the Navier-Stokes residuals\n");
  printf("  at the initial time T=0, using a region that is\n");
  printf("  the cube centered at (0,0,0) with radius 1\n");
  printf("  Parameter A=%lG\n",a);
  printf("  Parameter D=%lG\n",d);
  n=1000;
  pr=(double *) malloc(n*sizeof(double));
  ur=(double *) malloc(n*sizeof(double));
  vr=(double *) malloc(n*sizeof(double));
  wr=(double *) malloc(n*sizeof(double));
  xyz_lo=-1.;
  xyz_hi=+1.;
  seed=123456789;
  x=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  y=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  z=r8vec_uniform_ab_new(n,xyz_lo,xyz_hi,&seed);
  t=r8vec_zeros_new(n);
  resid_ethier(a,d,n,x,y,z,t,ur,vr,wr,pr);
  printf("  Minimum Maximum\n");
  printf("  Ur: %lG %lG\n",r8vec_amin(n,ur),r8vec_amax(n,ur));
  printf("  Vr: %lG %lG\n",r8vec_amin(n,vr),r8vec_amax(n,vr));
  printf("  Wr: %lG %lG\n",r8vec_amin(n,wr),r8vec_amax(n,wr));
  printf("  Pr: %lG %lG\n",r8vec_amin(n,pr),r8vec_amax(n,pr));
  free(pr);
  free(t);
  free(ur);
  free(vr);
  free(wr);
  free(x);
  free(y);
  free(z);
  }

int main(int argc,char**argv,char**envp)
  {
  printf("Test the NS3DE library.\n");
  uvwp_burgers_test();
  resid_burgers_test();
  uvwp_ethier_test();
  resid_ethier_test();
  printf("NS3DE_TEST\n");
  return(0);
  }
