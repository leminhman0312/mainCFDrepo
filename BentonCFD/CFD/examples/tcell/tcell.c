/* TCELL solves the steady incompressible Navier-Stokes equations
   in a 2D "T"-shaped region. The fluid flow problem is formulated
   in terms of primitive variables using finite element techniques
   with piecewise linear plus quadratic functions on triangles,
   that is, Taylor-Hood isoparametric elements.

   FORTRAN version (tcell.f) by Hyung-Chun Lee, Ajou University, Korea
   C implementation (tcell.c) by dudley.benton@gmail.com
 */
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define BC 0
#define WARP 1
double Courant=0.5;
double deltat;
double tolns=1E-6;
double velo=1.;
double visc=0.01;
double xlngth=1.;
double ylngth=1.;

int minun,nband,ncol1,nelemn,neqn1,nlband,np,nrow1,nuband,nx,ny;
int*indx,*insc,*node,*pivot;
double*a,*area,*f,*g,*h,*xc,*xm,*yc,*ym;

double sq(double x)
  {
  return(x*x);
  }

double dot(double*u,double*v,int n) /* sum(u[]*v[]) */
  {
  int i;
  double d;
  for(d=i=0;i<n;i++)
    d+=u[i]*v[i];
  return(d);
  }

int ivmax(double*v,int n) /* index of maximum */
  {
  int i,m;
  double a,x;
  m=0;
  x=fabs(v[m]);
  for(i=1;i<n;i++)
    {
    a=fabs(v[i]);
    if(a>x)
      {
      x=a;
      m=i;
      }
    }
  return(m);
  }

void scal(double*v,int n,double a) /* v[]*a */
  {
  int i;
  for(i=0;i<n;i++)
    v[i]*=a;
  }

void axpy(double*x,double*y,int n,double a) /* y[]+=a*x[] */
  {
  int i;
  for(i=0;i<n;i++)
    y[i]+=a*x[i];
  }

int Eliminate()
  {
  int i,i0,j,j0,j1,ju,jz,k,l,lm,m,mm,nm1;
  double t;
  m=nlband+nuband+1;
  j0=nuband+2;
  j1=min(neqn1,m)-1;
  if(j1>=j0)
    {
    for(jz=j0;jz<=j1;jz++)
      {
      i0=m+1-jz;
      for(i=i0;i<=nlband;i++)
        a[minun*(jz-1)+i-1]=0.;
      }
    }
  jz=j1;
  ju=0;
  nm1=neqn1-1;
  if(nm1>=1)
    {
    for(k=0;k<nm1;k++)
      {
      jz++;
      if(jz<=neqn1&&nlband>=1)
        for(i=0;i<nlband;i++)
          a[minun*(jz-1)+i]=0.;
      lm=min(nlband,neqn1-k-1);
      l=ivmax(&a[minun*k+m-1],lm+1)+m;
      pivot[k]=l+k+1-m;
      if(a[minun*k+l-1]==0.)
        return(k+1);
      if(l!=m)
        {
        t=a[minun*k+l-1];
        a[minun*k+l-1]=a[minun*k+m-1];
        a[minun*k+m-1]=t;
        }
      t=-1./a[minun*k+m-1];
      scal(&a[minun*k+m],lm,t);
      ju=min(max(ju,nuband+pivot[k]),neqn1);
      mm=m;
      if(ju<k+2)
        continue;
      for(j=k+1;j<ju;j++)
        {
        l--;
        mm--;
        t=a[minun*j+l-1];
        if(l!=mm)
          {
          a[minun*j+l-1]=a[minun*j+mm-1];
          a[minun*j+mm-1]=t;
          }
        axpy(&a[minun*k+m],&a[minun*j+mm],lm,t);
        }
      }
    }
  pivot[neqn1-1]=neqn1;
  if(a[minun*(neqn1-1)+m-1]==0.)
    return(neqn1);
  return(0);
  }

void SolveBanded()
  {
  int k,kb,l,la,lb,lm,m,nm1;
  double t;
  m=nuband+nlband+1;
  nm1=neqn1-1;
  if(nlband!=0&&nm1>=1)
    {
    for(k=0;k<nm1;k++)
      {
      lm=min(nlband,neqn1-k-1);
      l=pivot[k];
      t=f[l-1];
      if(l!=k+1)
        {
        f[l-1]=f[k];
        f[k]=t;
        }
      axpy(&a[minun*k+m],f+k+1,lm,t);
      }
    }
  for(kb=0;kb<neqn1;kb++)
    {
    k=neqn1-kb;
    f[k-1]/=a[minun*(k-1)+m-1];
    lm=min(k,m)-1;
    la=m-lm;
    lb=k-lm;
    t=-f[k-1];
    axpy(&a[minun*(k-1)+la-1],f+lb-1,lm,t);
    }
  }

void trans(int it,double xq,double yq,double*det,double*pj11,double*pj21,double*pj12,double*pj22)
  { /* transforms data between the reference and physical elements */
  int i1,i2,i3,i4,i5,i6;
  double d,f1x,f1y,f2x,f2y,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6;
  it*=6;
  i1=node[it++]-1;
  i2=node[it++]-1;
  i3=node[it++]-1;
  i4=node[it++]-1;
  i5=node[it++]-1;
  i6=node[it]-1;
  x1=xc[i1];
  y1=yc[i1];
  x2=xc[i2];
  y2=yc[i2];
  x3=xc[i3];
  y3=yc[i3];
  x4=xc[i4];
  y4=yc[i4];
  x5=xc[i5];
  y5=yc[i5];
  x6=xc[i6];
  y6=yc[i6];
  f1x=x1*(-3.+4.*xq+4.*yq)+x2*(4.*xq-1.)+x4*4.*(1.-2.*xq-yq)+x5*4.*yq-x6*4.*yq;
  f1y=x1*(-3.+4.*xq+4.*yq)+x3*(4.*yq-1.)+x6*4.*(1.-2.*yq-xq)+x5*4.*xq-x4*4.*xq;
  f2x=y1*(-3.+4.*xq+4.*yq)+y2*(4.*xq-1.)+y4*4.*(1.-2.*xq-yq)+y5*4.*yq-y6*4.*yq;
  f2y=y1*(-3.+4.*xq+4.*yq)+y3*(4.*yq-1.)+y6*4.*(1.-2.*yq-xq)+y5*4.*xq-y4*4.*xq;
  d=f1x*f2y-f1y*f2x;
  *pj11= f2y/d;
  *pj12=-f2x/d;
  *pj21=-f1y/d;
  *pj22= f1x/d;
  *det=fabs(d);
  }

double bc(int i,int j) /* boundary condition */
  {
  double x,y;
  x=xc[j];
  y=yc[j];
  if(WARP==1)
    {
    if(i==1&&y>0.6)
      return(velo);
    if(i==2&&y>0.6)
      return((x-0.5)*velo);
    return(0.);
    }
  if(WARP==2)
    {
    if(i==1&&y>0.6)
      return(velo);
    if(i==2&&y>0.6)
      return((0.5-x)*velo);
    if(i==2&&y<0.)
      return(velo);
    return(0.);
    }
  if(BC==1)
    {
    if(i==1)
      return(velo*16.*(y-0.5)*(1.-y));
    return(0.);
    }
  if(x<0.01&&i==1)
    return(velo);
  if(y<0.01&&i==2)
    return(velo);
  return(0.);
  }

double lbf(double x,double y,int iq)
  { /* linear basis function */
  if(iq==1)
    return(1.-x-y);
  if(iq==2)
    return(x);
  if(iq==3)
    return(y);
  return(0.);
  }

void qbf(double x,double y,int i,double*bb,double*bx,double*by)
  { /* quadradic basis function */
  if(i==1)
    {
    *bb=(1.-x-y)*(1.-2.*x-2.*y);
    *bx=-3.+4.*x+4.*y;
    *by=-3.+4.*x+4.*y;
    }
  else if(i==2)
    {
    *bb=x*(2.*x-1.);
    *bx=4.*x-1.;
    *by=0.;
    }
  else if(i==3)
    {
    *bb=y*(2.*y-1.);
    *bx=0.;
    *by=4.*y-1.;
    }
  else if(i==4)
    {
    *bb= 4.*x*(1.-x-y);
    *bx= 4.*(1.-2.*x-y);
    *by=-4.*x;
    }
  else if(i==5)
    {
    *bb=4.*x*y;
    *bx=4.*y;
    *by=4.*x;
    }
  else if(i==6)
    {
    *bb= 4.*y*(1.-x-y);
    *bx=-4.*y;
    *by= 4.*(1.-x-2.*y);
    }
  else
    {
    *bb=0.;
    *bx=0.;
    *by=0.;
    }
  }

void NavierStokes()
  {
  int i,ip,ipp,iq,iqq,iquad,it,iter,iuk,iukk,iun,iuse,j,kk;
  double aij,ar,arr,bb,bbb,bbbl,bbl,bbx,bby,bx,by,csim,det,diff,etax,etay,
    tbbx,tbby,tbx,tby,ubc,un[2],unx[2],uny[2],uqp,vqp,x,xix,xiy,y;
  for(iter=0;iter<10;iter++)
    {
    memcpy(g,f,neqn1*sizeof(double));
    memset(f,0,neqn1*sizeof(double));
    if(iter<3)
      csim=0.;
    else
      csim=1.;
    for(i=0;i<nrow1;i++)
      for(j=0;j<ncol1;j++)
        a[minun*j+i]=0.;
    for(it=0;it<nelemn;it++)
      {
      arr=area[it]/3.;
      for(iquad=0;iquad<3;iquad++)
        {
        y=ym[3*it+iquad];
        x=xm[3*it+iquad];
        trans(it,x,y,&det,&xix,&xiy,&etax,&etay);
        ar=arr*det;
        for(kk=0;kk<2;kk++)
          un[kk]=uny[kk]=unx[kk]=0.;
        uqp=vqp=0.;
        for(iq=0;iq<6;iq++)
          {
          qbf(x,y,iq+1,&bb,&tbx,&tby);
          bx=tbx*xix+tby*etax;
          by=tbx*xiy+tby*etay;
          ip=node[6*it+iq];
          for(iuk=0;iuk<2;iuk++)
            {
            iun=indx[2*(ip-1)+iuk];
            if(0<iun)
              {
              un[iuk]+=bb*g[iun-1];
              unx[iuk]+=bx*g[iun-1];
              uny[iuk]+=by*g[iun-1];
              if(iuk==0)
                uqp+=bb*h[iun-1];
              else
                vqp+=bb*h[iun-1];
              }
            else if(iun<0)
              {
              ubc=bc(iuk+1,ip-1);
              un[iuk]+=bb*ubc;
              unx[iuk]+=bx*ubc;
              uny[iuk]+=by*ubc;
              if(iuk==0)
                uqp+=bb*ubc;
              }
            }
          }
        for(iq=0;iq<6;iq++)
          {
          ip=node[6*it+iq];
          qbf(x,y,iq+1,&bb,&tbx,&tby);
          bx=tbx*xix+tby*etax;
          by=tbx*xiy+tby*etay;
          bbl=lbf(x,y,iq+1);
          for(iuk=0;iuk<3;iuk++)
            {
            if(iuk==2)
              i=insc[ip-1];
            else
              i=indx[2*(ip-1)+iuk];
            if(i<=0)
              continue;
            if(iuk==0)
              f[i-1]+=csim*((un[0]*unx[0]+un[1]*uny[0])*bb)*ar+uqp*bb*ar/deltat;
            else if(iuk==1)
              f[i-1]+=csim*((un[0]*unx[1]+un[1]*uny[1])*bb)*ar+vqp*bb*ar/deltat;
            for(iqq=0;iqq<6;iqq++)
              {
              ipp=node[6*it+iqq];
              qbf(x,y,iqq+1,&bbb,&tbbx,&tbby);
              bbx=tbbx*xix+tbby*etax;
              bby=tbbx*xiy+tbby*etay;
              bbbl=lbf(x,y,iqq+1);
              for(iukk=0;iukk<3;iukk++)
                {
                if(iukk<2)
                  j=indx[2*(ipp-1)+iukk];
                else
                  j=insc[ipp-1];
                if(j==0)
                  continue;
                aij=0.;
                if(i==neqn1)
                  continue;
                if(iuk==0)
                  {
                  if(iukk==0)
                    aij=visc*(by*bby+bx*bbx)+(bbb*unx[0]*bb)*csim+bb*bbx*un[0]+bb*bby*un[1]+bb*bbb/deltat;
                  else if(iukk==1)
                    aij=csim*bb*bbb*uny[0];
                  else
                    aij=-bx*bbbl;
                  }
                else if(iuk==1)
                  {
                  if(iukk==0)
                    aij=csim*bb*bbb*unx[1];
                  else if(iukk==1)
                    aij=(visc*(by*bby+bx*bbx)+(bb*bbb*uny[1])*csim+bb*bby*un[1]+bb*bbx*un[0])+bb*bbb/deltat;
                  else
                    aij=-by*bbbl;
                  }
                else
                  {
                  if(iukk==0)
                    aij=bbx*bbl;
                  else if(iukk==1)
                    aij=bby*bbl;
                  }
                if(j<0)
                  {
                  ubc=bc(iukk+1,ipp-1);
                  f[i-1]+=-ar*aij*ubc;
                  }
                else
                  {
                  iuse=i-j+nband;
                  a[minun*(j-1)+iuse-1]+=ar*aij;
                  }
                }
              }
            }
          }
        }
      }
    f[neqn1-1]=0.;
    for(j=neqn1-nlband;j<=neqn1-1;j++)
      {
      i=neqn1-j+nband;
      a[minun*(j-1)+i-1]=0.;
      }
    a[minun*(neqn1-1)+nband-1]=1.;
    if((i=Eliminate())!=0)
      {
      fprintf(stderr,"elimination error %i\n",i);
      exit(0);
      }
    SolveBanded();
    diff=0.;
    for(i=0;i<neqn1;i++)
      diff+=sq(g[i]-f[i]);
    diff=sqrt(diff);
    printf("  %lG\n",diff);
    if(diff<=tolns)
      break;
    }
  memcpy(h,f,neqn1*sizeof(double));
  }

double warpx(double x,double y)
  {
  if(WARP==1)
    return(0.1715299524+(-0.07573451603-0.00360371302*y)*y
      +(0.2260905319*y+0.7785826144-0.1789714093*x)*x);
  if(WARP==2)
    return(0.1055124062+(-0.1938031045+0.06873900085*y)*y
      +(0.3742084588*y+0.9948444536-0.2934551137*x)*x);
  return(x);
  }

double warpy(double x,double y)
  {
  if(WARP==1)
    return(0.1015474727+(1.101787627-0.2383772211*y)*y
      +(0.02883376506*y-0.5053718420+0.4726106261*x)*x);
  if(WARP==2)
    return(-0.1348679135+(1.589652436-0.4727760877*y)*y
    +(0.05600850975*y+0.502923559-0.5408293645*x)*x);
  return(y);
  }

void CreateGrid()
  {
  int i,ic,icnt,ij,ip,ip1,ip2,ipp,iq,iqq,iquater1,iquater2,it,it1,iuk,iukk,
    j,jc,jcnt,ncol,nrow,nxm1,nym1;
  double hx,hx1,hy,hy1,xx,yy;
  nym1=ny-1;
  nxm1=nx-1;
  nrow=nx+nxm1;
  ncol=ny+nym1;
  hx=xlngth/nxm1;
  hx1=hx/2.;
  hy=ylngth/nym1;
  hy1=hy/2.;
  i=0;
  ip=0;
  it1=-1;
  xx=0.-hx1;
  iquater1=(nrow-1)/4-2;
  iquater2=(3*(nrow-1))/4;
  for(ic=1;ic<=iquater1;ic++)
    {
    xx+=hx1;
    icnt=ic%2;
    for(jc=1;jc<=ny;jc++)
      {
      jcnt=jc%2;
      yy=0.5+hy1*(jc-1);
      ip++;
      xc[ip-1]=warpx(xx,yy);
      yc[ip-1]=warpy(xx,yy);
      if(icnt==1&&jcnt==1&&jc<ny)
        {
        it1+=2;
        nelemn=it1+1;
        ip1=ip+ny;
        ip2=ip+ny+ny;
        node[6*(it1-1)]=ip;
        node[6*(it1-1)+1]=ip2+2;
        node[6*(it1-1)+2]=ip2;
        node[6*(it1-1)+3]=ip1+1;
        node[6*(it1-1)+4]=ip2+1;
        node[6*(it1-1)+5]=ip1;
        node[6*(nelemn-1)]=ip;
        node[6*(nelemn-1)+1]=ip2+2;
        node[6*(nelemn-1)+2]=ip+2;
        node[6*(nelemn-1)+3]=ip1+1;
        node[6*(nelemn-1)+4]=ip1+2;
        node[6*(nelemn-1)+5]=ip+1;
        }
      if(jc!=1&&jc!=ny&&ic!=1)
        {
        i+=2;
        indx[2*(ip-1)]=i-1;
        indx[2*(ip-1)+1]=i;
        }
      else
        {
        indx[2*(ip-1)]=-1;
        indx[2*(ip-1)+1]=-1;
        }
      if(jcnt==1&&icnt==1)
        insc[ip-1]=++i;
      else
        insc[ip-1]=0;
      }
    }
  for(ic=iquater1+1;ic<=(iquater1+2);ic++)
    {
    xx+=hx1;
    icnt=ic%2;
    for(jc=1;jc<=ny;jc++)
      {
      jcnt=jc%2;
      yy=0.5+hy1*(jc-1);
      ip++;
      xc[ip-1]=warpx(xx,yy);
      yc[ip-1]=warpy(xx,yy);
      if(icnt==1&&jcnt==1)
        goto L_15;
      goto L_110;
  L_15:
      if(jc==ny)
        goto L_110;
      it1+=2;
      nelemn=it1+1;
      ip1=ip+ny;
      ip2=ip+ny+ncol;
      node[6*(it1-1)]=ip;
      node[6*(it1-1)+1]=ip2+2;
      node[6*(it1-1)+2]=ip2;
      node[6*(it1-1)+3]=ip1+1;
      node[6*(it1-1)+4]=ip2+1;
      node[6*(it1-1)+5]=ip1;
      node[6*(nelemn-1)]=ip;
      node[6*(nelemn-1)+1]=ip2+2;
      node[6*(nelemn-1)+2]=ip+2;
      node[6*(nelemn-1)+3]=ip1+1;
      node[6*(nelemn-1)+4]=ip1+2;
      node[6*(nelemn-1)+5]=ip+1;
  L_110:;
      if(jc==1||jc==ny)
        goto L_120;
      i+=2;
      indx[2*(ip-1)]=i-1;
      indx[2*(ip-1)+1]=i;
      goto L_129;
  L_120:;
      indx[2*(ip-1)]=-1;
      indx[2*(ip-1)+1]=-1;
  L_129:;
      if(jcnt==0||icnt==0)
        goto L_130;
      i++;
      insc[ip-1]=i;
      continue;
  L_130:
      insc[ip-1]=0;
      }
    }
  for(ic=iquater1+3;ic<=iquater2;ic++)
    {
    xx+=hx1;
    icnt=ic%2;
    for(jc=1;jc<=ncol;jc++)
      {
      jcnt=jc%2;
      yy=hy1*(jc-1);
      ip++;
      xc[ip-1]=warpx(xx,yy);
      yc[ip-1]=warpy(xx,yy);
      if(icnt==1&&jcnt==1)
        goto L_25;
      goto L_210;
  L_25:
      if(jc==ncol)
        goto L_210;
      it1+=2;
      nelemn=it1+1;
      ip1=ip+ncol;
      ip2=ip+ncol+ncol;
      node[6*(it1-1)]=ip;
      node[6*(it1-1)+1]=ip2+2;
      node[6*(it1-1)+2]=ip2;
      node[6*(it1-1)+3]=ip1+1;
      node[6*(it1-1)+4]=ip2+1;
      node[6*(it1-1)+5]=ip1;
      node[6*(nelemn-1)]=ip;
      node[6*(nelemn-1)+1]=ip2+2;
      node[6*(nelemn-1)+2]=ip+2;
      node[6*(nelemn-1)+3]=ip1+1;
      node[6*(nelemn-1)+4]=ip1+2;
      node[6*(nelemn-1)+5]=ip+1;
  L_210:;
      if(jc==1||jc==ncol)
        goto L_220;
      if(ic==iquater1+3&&jc<=ny)
        goto L_220;
      i+=2;
      indx[2*(ip-1)]=i-1;
      indx[2*(ip-1)+1]=i;
      goto L_229;
  L_220:;
      indx[2*(ip-1)]=-1;
      indx[2*(ip-1)+1]=-1;
  L_229:;
      if(jcnt==0||icnt==0)
        goto L_230;
      i++;
      insc[ip-1]=i;
      continue;
  L_230:
      insc[ip-1]=0;
      }
    }
  xx+=hx1;
  icnt=(iquater2+1)%2;
  for(jc=1;jc<=(ny-1);jc++)
    {
    jcnt=jc%2;
    yy=hy1*(jc-1);
    ip++;
    xc[ip-1]=warpx(xx,yy);
    yc[ip-1]=warpy(xx,yy);
    indx[2*(ip-1)]=-1;
    indx[2*(ip-1)+1]=-1;
    if(jcnt==0||icnt==0)
      goto L_930;
    i++;
    insc[ip-1]=i;
    continue;
 L_930:
    insc[ip-1]=0;
    }
  xx-=hx1;
  for(ic=iquater2+1;ic<=nrow;ic++)
    {
    xx+=hx1;
    icnt=ic%2;
    for(jc=1;jc<=ny;jc++)
      {
      jcnt=jc%2;
      yy=0.5+hy1*(jc-1);
      ip++;
      xc[ip-1]=warpx(xx,yy);
      yc[ip-1]=warpy(xx,yy);
      if(icnt==1&&jcnt==1)
        goto L_435;
      goto L_310;
  L_435:
      if(ic==nrow||jc==ny)
        goto L_310;
      it1+=2;
      nelemn=it1+1;
      ip1=ip+ny;
      ip2=ip+ny+ny;
      node[6*(it1-1)]=ip;
      node[6*(it1-1)+1]=ip2+2;
      node[6*(it1-1)+2]=ip2;
      node[6*(it1-1)+3]=ip1+1;
      node[6*(it1-1)+4]=ip2+1;
      node[6*(it1-1)+5]=ip1;
      node[6*(nelemn-1)]=ip;
      node[6*(nelemn-1)+1]=ip2+2;
      node[6*(nelemn-1)+2]=ip+2;
      node[6*(nelemn-1)+3]=ip1+1;
      node[6*(nelemn-1)+4]=ip1+2;
      node[6*(nelemn-1)+5]=ip+1;
  L_310:;
      if(jc==1||jc==ny)
        goto L_320;
      i+=2;
      indx[2*(ip-1)]=i-1;
      indx[2*(ip-1)+1]=i;
      goto L_329;
  L_320:;
      indx[2*(ip-1)]=-1;
      indx[2*(ip-1)+1]=-1;
  L_329:;
      if(jcnt==0||icnt==0)
        goto L_330;
      i++;
      insc[ip-1]=i;
      continue;
  L_330:
      insc[ip-1]=0;
      }
    }
  np=ip;
  neqn1=i;
  for(it=0;it<nelemn;it++)
    {
    xm[3*it  ]=0.5;
    xm[3*it+1]=0.5;
    xm[3*it+2]=0.0;
    ym[3*it  ]=0.0;
    ym[3*it+1]=0.5;
    ym[3*it+2]=0.5;
    area[it]=0.5;
    }
  nlband=0;
  for(it=0;it<nelemn;it++)
    {
    for(iq=0;iq<6;iq++)
      {
      ip=node[6*it+iq];
      for(iuk=0;iuk<3;iuk++)
        {
        if(iuk==2)
          i=insc[ip-1];
        else
          i=indx[2*(ip-1)+iuk];
        if(0<i)
          {
          for(iqq=0;iqq<6;iqq++)
            {
            ipp=node[6*it+iqq];
            for(iukk=0;iukk<3;iukk++)
              {
              if(iukk==2)
                j=insc[ipp-1];
              else
                j=indx[2*(ipp-1)+iukk];
              if(i<=j)
                {
                ij=j-i;
                if(ij>nlband)
                  nlband=ij;
                }
              }
            }
          }
        }
      }
    }
  nband=nlband+nlband+1;
  }

void Write2DV(char*fname)
  {
  int i;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"%i\n",np);
  for(i=0;i<np;i++)
    fprintf(fp,"%lG %lG\n",xc[i],yc[i]);
  fprintf(fp,"%i\n",nelemn);
  for(i=0;i<nelemn;i++)
    fprintf(fp,"%i %i %i\n",node[6*i],node[6*i+1],node[6*i+2]);
  fclose(fp);
  }

void WriteV2D(char*fname)
  {
  int ic,iu,iv;
  double u,v;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  for(ic=0;ic<np;ic++)
    {
    iu=indx[2*ic];
    iv=indx[2*ic+1];
    if(iu==0)
      u=0.;
    else if(iu==-1)
      u=bc(1,ic);
    else
      u=f[iu-1];
    if(iv<=0)
      v=0.;
    else
      v=f[iv-1];
    fprintf(fp,"%lG %lG %lG %lG\n",xc[ic],yc[ic],u,v);
    }
  fclose(fp);
  }

void WritePLT(char*fname)
  {
  int i,ic,ip,iu,iv;
  double p,u,v;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"TITLE=\"T-Cell\"\n");
  fprintf(fp,"VARIABLES=\"X\" \"Y\" \"P\" \"U\" \"V\"\n");
  fprintf(fp,"ZONE T=\"results\" N=%i, E=%i, F=FEPOINT, ET=TRIANGLE\n",np,nelemn);
  fprintf(fp,"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n");
  for(ic=0;ic<np;ic++)
    {
    iu=indx[2*ic];
    iv=indx[2*ic+1];
    ip=insc[ic];
    if(iu==0)
      u=0.;
    else if(iu==-1)
      u=bc(1,ic);
    else
      u=f[iu-1];
    if(iv<=0)
      v=0.;
    else
      v=f[iv-1];
    if(ip==0)
      p=0.;
    else
      p=f[ip-1];
    fprintf(fp,"%lG %lG %lG %lG %lG\n",xc[ic],yc[ic],p,u,v);
    }
  for(i=0;i<nelemn;i++)
    fprintf(fp,"%i %i %i\n",node[6*i],node[6*i+1],node[6*i+2]);
  fclose(fp);
  }

void*allocate(size_t n,size_t s)
  {
  void*ptr;
  if((ptr=calloc(n,s))==NULL)
    {
    fprintf(stderr,"can't allocate memory\n");
    exit(1);
    }
  return(ptr);
  }

int main(int argc,char**argv,char**envp)
  {
  int i,j,maxel,maxel3,maxel6,maxmin,maxnd,maxnd2,maxun,mx,my;
  double t,tend=1.;
  printf("Navier-Stokes Equations of Fluid Flow\n");
  printf("using finite elements\n");
  nx=ny=21;
  printf("grid=%ix%i\n",nx,ny);
  printf("allocating memory\n");
  mx=2*nx-1;
  my=2*ny-1;
  maxel=2*(nx-1)*(ny-1);
  maxnd=mx*my;
  maxun=2*mx*my+nx*ny;
  minun=27*ny;
  maxnd2=maxnd*2;
  maxel3=maxel*3;
  maxel6=maxel*6;
  maxmin=maxun*minun;
  indx =allocate(maxnd2,sizeof(int));
  insc =allocate(maxnd ,sizeof(int));
  node =allocate(maxel6,sizeof(int));
  pivot=allocate(maxun ,sizeof(int));
  a    =allocate(maxmin,sizeof(double));
  area =allocate(maxel ,sizeof(double));
  f    =allocate(maxun ,sizeof(double));
  g    =allocate(maxun ,sizeof(double));
  h    =allocate(maxun ,sizeof(double));
  xc   =allocate(maxnd ,sizeof(double));
  yc   =allocate(maxnd ,sizeof(double));
  xm   =allocate(maxel3,sizeof(double));
  ym   =allocate(maxel3,sizeof(double));
  CreateGrid();
  printf("number of nodes=%i\n",np);
  printf("number of elements=%i\n",nelemn);
  printf("lower bandwidth=%i\n",nlband);
  printf("total bandwidth=%i\n",nband);
  Write2DV("tcell.2dv");
  nuband=nlband;
  nrow1=nlband+nlband+nuband+1;
  ncol1=neqn1;
  for(i=0;i<neqn1;i++)
    h[i]=0.5;
  printf("Re=%lG\n",min(xlngth,ylngth)*velo/visc);
  deltat=Courant*(xlngth/nx+ylngth/ny)/velo;
  printf("deltat=%lG\n",deltat);
  t=0.;
  do{
    NavierStokes();
    t+=deltat;
    printf("t=%lG\n",t);
    if(_kbhit())
      {
      _getch();
      break;
      }
    }while(t<tend);
  WriteV2D("tcell.v2d");
  WritePLT("tcell.plt");
if(0)
  for(i=0;i<100;i++)
    {
    printf("%i %i %i",i,indx[i],insc[i]);
    for(j=0;j<6;j++)
      printf(" %i",node[6*i+j]);
    printf("\n");
    }
  return(0);
  }
