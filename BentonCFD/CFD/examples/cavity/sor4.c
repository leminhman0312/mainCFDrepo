/* Vorticity-Stream Function Method
   for the cavity flow problem using
   SOR and 4th order differences
   based on the work of
   Erturk, Corke, and Gokcol
   this implementation by
   dudley.benton@gmail.com
 */

#define _CRT_SECURE_NO_DEPRECATE
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

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

#define n 100
double a[n*n];
double b[n*n];
double c[n*n];
double d[n*n];
double e[n*n];
double f[n*n];
double So[n*n];
double Sx[n*n];
double Sy[n*n];
double S[n*n];
double U[n*n];
double V[n*n];
double Wo[n*n];
double W[n*n];
double X[n];
double Y[n];

void WriteTB2(char*fname,double*Z,double Zm,double Zx)
  {
  int i,j;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"%i\n",n);
  for(i=0;i<n;i++)
    fprintf(fp,"%lG\n",X[i]);
  fprintf(fp,"%i\n",n);
  for(j=0;j<n;j++)
    fprintf(fp,"%lG\n",Y[j]);
  fprintf(fp,"%i\n",n*n);
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      fprintf(fp,"%lG\n",fmax(Zm,fmin(Zx,Z[n*j+i])));
  fclose(fp);
  }

void WriteV2D(char*fname)
  {
  int i,j;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      fprintf(fp,"%lG %lG %lG %lG\n",X[i],Y[j],U[n*j+i],V[n*j+i]);
  fclose(fp);
  }

void WriteTP2(char*fname,double Sm,double Sx,double Re)
  {
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"TP2\n");
  fprintf(fp,"!PLOT\n");
  fprintf(fp,"TICKS\n");
  fprintf(fp,"EQUAL\n");
  fprintf(fp,"LEGENDS\n");
  fprintf(fp,"-1 5 2\n");
  fprintf(fp,"COLOR\n");
  fprintf(fp,"%lG %lG 2\n",Sm-(Sx-Sm)/1000.,Sx+(Sx-Sm)/1000.);
  fprintf(fp,"RAINBOW\n");
  fprintf(fp,"ZTABLE\n");
  fprintf(fp,"NOTITLE\n");
  fprintf(fp,"X\n");
  fprintf(fp,"%lG %lG\n",X[0],X[n-1]);
  fprintf(fp,"Y\n");
  fprintf(fp,"%lG %lG\n",Y[0],Y[n-1]);
  fprintf(fp,"PSI.TB2\n");
  fprintf(fp,"NOCOLOR\n");
  fprintf(fp,"ARROW\n");
  fprintf(fp,"0 0\n");
  fprintf(fp,"NEWFILE\n");
  fprintf(fp,"SOR4.V2D\n");
  fprintf(fp,"%lG %lG .    Re=%lG\n",X[n-1],Y[0]+(Y[n-1]-Y[0])/100.,Re);
  fprintf(fp,"!END\n");
  fclose(fp);
  }

void WritePLT(char*fname)
  {
  int i,j;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"TITLE=\"Cavity\"\n");
  fprintf(fp,"VARIABLES=\"X\" \"Y\" \"U\" \"V\" \"Psi\" \"Omega\"\n");
  fprintf(fp,"ZONE T=\"results\"  I=%i, J=%i, K=1, F=POINT\n",n,n);
  fprintf(fp,"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n");
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      fprintf(fp,"%lG %lG %lG %lG %lG %lG\n",X[i],Y[j],U[n*j+i],V[n*j+i],S[n*j+i],W[n*j+i]);
  fclose(fp);
  }

int main(int argc,char**argv,char**envp)
  {
  int done=-1,i,j;
  double beta=0.6,dh,dh2,dh3,dh4,r,Re,rr,tol=1E-6,Zm,Zx;
  double Sxx,Sxy,Syy,Sxxy,Sxyy,Sxxyy,Vxx,Vxy,Vyy,Vxxy,Vxyy,Vxxyy;
  Re=1000.;
  dh=1./(n-1);
  dh2=dh*dh;
  dh3=dh2*dh;
  dh4=dh3*dh;
  for(i=0;i<n;i++)
    X[i]=Y[i]=i*dh;
  for(i=0;i<n;i++)
    {
    for(j=0;j<n;j++)
      S[n*j+i]=W[n*j+i]=tol;
    S[n*i]=W[n*i]=0.;
    S[n*i+n-1]=W[n*i+n-1]=0.;
    S[i]=W[i]=0.;
    S[n*(n-1)+i]=W[n*(n-1)+i]=0.;
    }
  printf("solving\n");
  do{
    memcpy(So,S,sizeof(So));
    memcpy(Wo,W,sizeof(Wo));
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        Vxx=(W[n*j+i-1]-2.*W[n*j+i]+W[n*j+i+1])/dh2;
        Vyy=(W[n*(j-1)+i]-2.*W[n*j+i]+W[n*(j+1)+i])/dh2;
        Sxxyy=((S[n*(j+1)+i+1]+S[n*(j+1)+i-1]+S[n*(j-1)+i+1]
          +S[n*(j-1)+i-1])-2.*(S[n*(j+1)+i]+S[n*j+i+1]
          +S[n*(j-1)+i]+S[n*j+i-1])+4.*S[n*j+i])/dh4;
        a[n*j+i]=dh2*Vxx/12.+dh2*Vyy/12.+dh2*Sxxyy/6.;
        S[n*j+i]=beta*((S[n*j+i-1]+S[n*j+i+1]+S[n*(j-1)+i]
          +S[n*(j+1)+i]+dh2*W[n*j+i]+dh2*a[n*j+i])/4.)
          +(1.-beta)*S[n*j+i];
        }
      }
    for(i=1;i<n-1;i++)
      {
      W[i]=(-(S[n+i-1]+S[n+i]+S[n+i+1])/(3.*dh2)
        -(0.5*W[i-1]+0.5*W[i+1]+0.25*W[n+i-1]
        +W[n+i]+0.25*W[n+i+1])/9.)*9./2.;
      W[n*(n-1)+i]=(-1./dh-(S[n*(n-2)+i-1]+S[n*(n-2)+i]
        +S[n*(n-2)+i+1])/(3.*dh2)-(0.5*W[n*(n-1)+i-1]
        +0.5*W[n*(n-1)+i+1]+0.25*W[n*(n-2)+i-1]+W[n*(n-2)+i]
        +0.25*W[n*(n-2)+i+1])/9.)*9./2.;
      W[n*i]=(-(S[n*(i-1)+1]+S[n*i+1]+S[n*(i+1)+1])/(3.*dh2)
        -(0.5*W[n*(i-1)]+0.5*W[n*(i+1)]+0.25*W[n*(i-1)+1]
        +W[n*i+1]+0.25*W[n*(i+1)+1])/9.)*9./2.;
      W[n*i+n-1]=(-(S[n*(i-1)+n-2]+S[n*i+n-2]+S[n*(i+1)+n-2])/(3.*dh2)
        -(0.5*W[n*(i-1)+n-1]+0.5*W[n*(i+1)+n-1]+0.25*W[n*(i-1)+n-2]
        +W[n*i+n-2]+0.25*W[n*(i+1)+n-2])/9.)*9./2.;
      }
    W[0]=(-(S[n+1])/(3.*dh2)-(0.5*W[1]+0.5*W[n]
      +0.25*W[n+1])/9.)*9.;
    W[n-1]=(-(S[n*2-2])/(3.*dh2)-(0.5*W[n-2]+0.5*W[n*2-1]
      +0.25*W[n+n-2])/9.)*9.;
    W[n*(n-1)+n-1]=(-0.5/dh-(S[n*(n-2)+n-2])/(3.*dh2)-(0.5*W[n*(n-1)+n-2]
      +0.5*W[n*(n-2)+n-1]+0.25*W[n*(n-2)+n-2])/9.)*9.;
    W[n*(n-1)]=(-0.5/dh-(S[n*(n-2)+1])/(3.*dh2)-(0.5*W[n*(n-1)+1]
      +0.5*W[n*(n-2)]+0.25*W[n*(n-2)+1])/9.)*9.;
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        Sx[n*j+i]=Re*(S[n*j+i+1]-S[n*j+i-1])/(2.*dh);
        Sy[n*j+i]=Re*(S[n*(j+1)+i]-S[n*(j-1)+i])/(2.*dh);
        Sxx=Re*(S[n*j+i-1]-2.*S[n*j+i]+S[n*j+i+1])/dh2;
        Syy=Re*(S[n*(j-1)+i]-2.*S[n*j+i]+S[n*(j+1)+i])/dh2;
        Sxy=Re*(0.25*(S[n*(j+1)+i+1]-S[n*(j+1)+i-1]-S[n*(j-1)+i+1]+S[n*(j-1)+i-1]))/dh2;
        Sxxy=Re*(0.5*(S[n*(j+1)+i+1]+S[n*(j+1)+i-1]-S[n*(j-1)+i+1]
          -S[n*(j-1)+i-1])+S[n*(j-1)+i]-S[n*(j+1)+i])/dh3;
        Sxyy=Re*(0.5*(S[n*(j+1)+i+1]-S[n*(j+1)+i-1]+S[n*(j-1)+i+1]
          -S[n*(j-1)+i-1])+S[n*j+i-1]-S[n*j+i+1])/dh3;
        Vxy=(0.25*(W[n*(j+1)+i+1]-W[n*(j+1)+i-1]-W[n*(j-1)+i+1]+W[n*(j-1)+i-1]))/dh2;
        Vxxy=(0.5*(W[n*(j+1)+i+1]+W[n*(j+1)+i-1]-W[n*(j-1)+i+1]-W[n*(j-1)+i-1])
          +W[n*(j-1)+i]-W[n*(j+1)+i])/dh3;
        Vxyy=(0.5*(W[n*(j+1)+i+1]-W[n*(j+1)+i-1]+W[n*(j-1)+i+1]-W[n*(j-1)+i-1])
          +W[n*j+i-1]-W[n*j+i+1])/dh3;
        Vxxyy=((W[n*(j+1)+i+1]+W[n*(j+1)+i-1]+W[n*(j-1)+i+1]+W[n*(j-1)+i-1])
          -2.*(W[n*(j+1)+i]+W[n*j+i+1]+W[n*(j-1)+i]+W[n*j+i-1])+4.*W[n*j+i])/dh4;
        b[n*j+i]=dh2*Sy[n*j+i]*Sy[n*j+i]/12.-dh2*Sxy/6.;
        c[n*j+i]=dh2*Sx[n*j+i]*Sx[n*j+i]/12.+dh2*Sxy/6.;
        d[n*j+i]=dh2*Sxxy/6.-dh2*Sy[n*j+i]*Sxy/12.+dh2*Sx[n*j+i]*Syy/12.;
        e[n*j+i]=dh2*Sxyy/6.-dh2*Sy[n*j+i]*Sxx/12.+dh2*Sx[n*j+i]*Sxy/12.;
        f[n*j+i]=dh2*Vxxyy/6.+dh2*Sx[n*j+i]*Vxxy/6.-dh2*Sy[n*j+i]*Vxyy/6.
          -dh2*Sx[n*j+i]*Sy[n*j+i]*Vxy/6.+dh2*Sxx*Vxy/6.-dh2*Syy*Vxy/6.;
        W[n*j+i]=beta*(((1.+b[n*j+i])*W[n*j+i-1]+(1.+b[n*j+i])*W[n*j+i+1]
          +(1.+c[n*j+i])*W[n*(j-1)+i]+(1.+c[n*j+i])*W[n*(j+1)+i]
          -0.5*(Sy[n*j+i]+d[n*j+i])*dh*(W[n*j+i+1]-W[n*j+i-1])
          +0.5*(Sx[n*j+i]+e[n*j+i])*dh*(W[n*(j+1)+i]-W[n*(j-1)+i])
          +dh2*f[n*j+i])/(2.*(1.+b[n*j+i])+2.*(1.+c[n*j+i])))
          +(1.-beta)*W[n*j+i];
        }
      }
    for(rr=0.,i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        r=fabs((S[n*j+i-1]-2.*S[n*j+i]+S[n*j+i+1])/dh2
          +(S[n*(j-1)+i]-2.*S[n*j+i]+S[n*(j+1)+i])/dh2+W[n*j+i]+a[n*j+i]);
        if(r>rr)
          rr=r;
        r=fabs(((1.+b[n*j+i])*(W[n*j+i-1]-2.*W[n*j+i]+W[n*j+i+1])/dh2
          +(1.+c[n*j+i])*(W[n*(j-1)+i]-2.*W[n*j+i]+W[n*(j+1)+i])/dh2
          -(Sy[n*j+i]+d[n*j+i])*(W[n*j+i+1]-W[n*j+i-1])/(2.*dh)
          +(Sx[n*j+i]+e[n*j+i])*(W[n*(j+1)+i]-W[n*(j-1)+i])/(2.*dh)+f[n*j+i])/Re);
        if(r>rr)
          rr=r;
        r=fabs(S[n*j+i]-So[n*j+i]);
        if(r>rr)
          rr=r;
        r=fabs(W[n*j+i]-Wo[n*j+i]);
        if(r>rr)
          rr=r;
        r=fabs((S[n*j+i]-So[n*j+i])/So[n*j+i]);
        if(r>rr)
          rr=r;
        r=fabs((W[n*j+i]-Wo[n*j+i])/Wo[n*j+i]);
        if(r>rr)
          rr=r;
        }
      }
    i=(int)(100.*pow(tol/rr,0.2)+0.5);
    if(i>done)
      {
      done=i;
      fprintf(stderr,"\r%i%%",done);
      }
    if(_kbhit())
      {
      _getch();
      rr=0.;
      }
    }while(rr>tol);
  printf("\rdone\n");
  printf("calculating U\n");
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      if(j==0)
        U[n*j+i]=0.;
      else if(j==n-1)
        U[n*j+i]=1.;
      else
        U[n*j+i]=(S[n*(j+1)+i]-S[n*(j-1)+i])/(2.*dh);
  printf("calculating V\n");
  for(j=0;j<n;j++)
    for(i=0;i<n;i++)
      if(i==0||i==n-1)
        V[n*j+i]=0.;
      else
        V[n*j+i]=-(S[n*j+i+1]-S[n*j+i-1])/(2.*dh);
  printf("bounding psi\n");
  Zm=Zx=S[0];
  for(i=0;i<n;i++)
    {
    for(j=0;j<n;j++)
      {
      if(S[n*i+j]<Zm)
        Zm=S[n*i+j];
      if(S[n*i+j]>Zx)
        Zx=S[n*i+j];
      }
    }
  WriteTB2("psi.tb2",S,-DBL_MAX,DBL_MAX);
  WriteTB2("omega.tb2",W,-5.,5.);
  WriteV2D("sor4.v2d");
  WriteTP2("sor4.tp2",Zm,Zx,Re);
  WritePLT("sor4.plt");
  return(0);
  }
