/* Navier-Stokes by way of Vorticity-Stream Function
   for the cavity flow problem using the
   Alternating Direction Implicit (ADI) method
   4th order accuracy based on the work of Erturk
   this implementation by dudley.benton@gmail.com
 */

#define _CRT_SECURE_NO_DEPRECATE
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
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

double aa[n];
double bb[n];
double b[n*n];
double cc[n];
double c[n*n];
double dd[n];
double a[n*n];
double d[n*n];
double e[n*n];
double f[n*n];
double R1[n*n];
double R2[n*n];
double S[n*n];
double So[n*n];
double U[n*n];
double V[n*n];
double W[n*n];
double Wo[n*n];
double X[n];
double Y[n];

void WriteTB2(char*fname,double Z[n*n],double Zm,double Zx)
  {
  int i,j,k;
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
  for(k=j=0;j<n;j++)
    for(i=0;i<n;i++,k++)
      fprintf(fp,"%lG\n",fmax(Zm,fmin(Zx,Z[k])));
  fclose(fp);
  }

void WriteV2D(char*fname)
  {
  int i,j,k;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  for(k=j=0;j<n;j++)
    for(i=0;i<n;i++,k++)
      fprintf(fp,"%lG %lG %lG %lG\n",X[i],Y[j],U[k],V[k]);
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
  fprintf(fp,"ADI4.V2D\n");
  fprintf(fp,"%lG %lG .    Re=%lG\n",X[n-1],Y[0]+(Y[n-1]-Y[0])/100.,Re);
  fprintf(fp,"!END\n");
  fclose(fp);
  }

void WritePLT(char*fname)
  {
  int i,j,k;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"TITLE=\"Cavity\"\n");
  fprintf(fp,"VARIABLES=\"X\" \"Y\" \"U\" \"V\" \"Psi\" \"Omega\"\n");
  fprintf(fp,"ZONE T=\"results\"  I=%i, J=%i, K=1, F=POINT\n",n,n);
  fprintf(fp,"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n");
  for(k=j=0;j<n;j++)
    for(i=0;i<n;i++,k++)
      fprintf(fp,"%lG %lG %lG %lG %lG %lG\n",X[i],Y[j],U[k],V[k],S[k],W[k]);
  fclose(fp);
  }

int main(int argc,char**argv,char**envp)
  {
  int done=-1,i,j,k;
  double alpha=0.6,dh,dh2,dh3,dh4,dt,r,Re,rr,tol=1E-6,Zm,Zx;
  double s2x2,s2xy,s2y2,s3x2y,s3xy2,s4x2y2,sx,sy,v2x2,v2xy,v2y2,v3x2y,v3xy2,v4x2y2;
  printf("solving the cavity problem using\n"
         "the vorticity-stream function method with\n"
         "4th Order Alternating Direction Implicit Scheme\n");
  Re=1000.;
  printf("Re=%lG\n",Re);
  dh=1./(n-1);
  dh2=dh*dh;
  dh3=dh2*dh;
  dh4=dh3*dh;
  dt=alpha*dh*dh;
  printf("initializing\n");
  for(i=0;i<n;i++)
    X[i]=Y[i]=i*dh;
  for(i=0;i<n;i++)
    {
    for(j=0;j<n;j++)
      S[n*j+i]=W[n*j+i]=tol;
    S[n*i]=W[n*i]=0.;
    S[i]=W[i]=0.;
    S[n*i+n-1]=W[n*i+n-1]=0.;
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
        v2x2=(W[n*j+i-1]-2.*W[n*j+i]+W[n*j+i+1])/dh2;
        v2y2=(W[n*(j-1)+i]-2.*W[n*j+i]+W[n*(j+1)+i])/dh2;
        s4x2y2=((S[n*(j+1)+i+1]+S[n*(j+1)+i-1]+S[n*(j-1)+i+1]
          +S[n*(j-1)+i-1])-2.*(S[n*(j+1)+i]+S[n*j+i+1]
          +S[n*(j-1)+i]+S[n*j+i-1])+4.*S[n*j+i])/dh4;
        a[n*j+i]=dh2*v2x2/12.+dh2*v2y2/12.+dh2*s4x2y2/6.;
        }
      }
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        R1[n*j+i]=S[n*j+i]+0.5*dt*(S[n*(j-1)+i]
          -2.*S[n*j+i]+S[n*(j+1)+i])/dh2
          +0.5*dt*W[n*j+i]+0.5*dt*a[n*j+i];
        }
      }
    for(j=1;j<n-1;j++)
      {
      for(i=1;i<n-1;i++)
        {
        aa[i]=-0.5*dt/dh2;
        bb[i]=1.+0.5*dt*2./dh2;
        cc[i]=-0.5*dt/dh2;
        dd[i]=R1[n*j+i];
        }
      for(i=2;i<n-1;i++)
        {
        bb[i]+=-aa[i]*cc[i-1]/bb[i-1];
        dd[i]+=-aa[i]*dd[i-1]/bb[i-1];
        }
      S[n*j+n-2]=dd[n-2]/bb[n-2];
      for(i=n-3;i>=1;i--)
        S[n*j+i]=(dd[i]-S[n*j+i+1]*cc[i])/bb[i];
      }
    for(i=1;i<n-1;i++)
      for(j=1;j<n-1;j++)
        R2[n*j+i]=S[n*j+i]+0.5*dt*(S[n*j+i-1]
          -2.*S[n*j+i]+S[n*j+i+1])/dh2+0.5*dt*W[n*j+i]
          +0.5*dt*a[n*j+i];
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        aa[j]=-0.5*dt/dh2;
        bb[j]=1.+0.5*dt*2./dh2;
        cc[j]=-0.5*dt/dh2;
        dd[j]=R2[n*j+i];
        }
      for(j=2;j<n-1;j++)
        {
        bb[j]+=-aa[j]*cc[j-1]/bb[j-1];
        dd[j]+=-aa[j]*dd[j-1]/bb[j-1];
        }
      S[n*(n-2)+i]=dd[n-2]/bb[n-2];
      for(j=n-3;j>=1;j--)
        S[n*j+i]=(dd[j]-S[n*(j+1)+i]*cc[j])/bb[j];
      }
    for(i=1;i<n-1;i++)
      {
      W[i]=(-(S[n*1+i-1]+S[n*1+i]+S[n*1+i+1])/(3.*dh2)
        -(0.5*W[i-1]+0.5*W[i+1]+0.25*W[n*1+i-1]+W[n*1+i]
        +0.25*W[n*1+i+1])/9.)*9./2.;
      W[n*(n-1)+i]=(-1./dh-(S[n*(n-2)+i-1]+S[n*(n-2)+i]
        +S[n*(n-2)+i+1])/(3.*dh2)-(0.5*W[n*(n-1)+i-1]
        +0.5*W[n*(n-1)+i+1]+0.25*W[n*(n-2)+i-1]+W[n*(n-2)+i]
        +0.25*W[n*(n-2)+i+1])/9.)*9./2.;
      W[n*i]=(-(S[n*(i-1)+1]+S[n*i+1]+S[n*(i+1)+1])/(3.*dh2)
        -(0.5*W[n*(i-1)]+0.5*W[n*(i+1)]+0.25*W[n*(i-1)+1]+W[n*i+1]
        +0.25*W[n*(i+1)+1])/9.)*9./2.;
      W[n*i+n-1]=(-(S[n*(i-1)+n-2]+S[n*i+n-2]+S[n*(i+1)+n-2])/
        (3.*dh2)-(0.5*W[n*(i-1)+n-1]+0.5*W[n*(i+1)+n-1]
        +0.25*W[n*(i-1)+n-2]+W[n*i+n-2]+0.25*W[n*(i+1)+n-2])/9.)*9./2.;
      }
    W[0]=(-(S[n*1+1])/(3.*dh2)-(0.5*W[1]
      +0.5*W[n*1]+0.25*W[n*1+1])/9.)*9.;
    W[n-1]=(-(S[n*1+n-2])/(3.*dh2)-(0.5*W[n-2]
      +0.5*W[n*1+n-1]+0.25*W[n*1+n-2])/9.)*9.;
    W[n*(n-1)+n-1]=(-0.5/dh-(S[n*(n-2)+n-2])/(3.*dh2)
      -(0.5*W[n*(n-1)+n-2]+0.5*W[n*(n-2)+n-1]+0.25*W[n*(n-2)+n-2])/9.)*9.;
    W[n*(n-1)]=(-0.5/dh-(S[n*(n-2)+1])/(3.*dh2)
      -(0.5*W[n*(n-1)+1]+0.5*W[n*(n-2)]+0.25*W[n*(n-2)+1])/9.)*9.;
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        sx=Re*(S[n*j+i+1]-S[n*j+i-1])/(2.*dh);
        sy=Re*(S[n*(j+1)+i]-S[n*(j-1)+i])/(2.*dh);
        s2x2=Re*(S[n*j+i-1]-2.*S[n*j+i]+S[n*j+i+1])/dh2;
        s2y2=Re*(S[n*(j-1)+i]-2.*S[n*j+i]+S[n*(j+1)+i])/dh2;
        s2xy=Re*(0.25*(S[n*(j+1)+i+1]-S[n*(j+1)+i-1]
          -S[n*(j-1)+i+1]+S[n*(j-1)+i-1]))/dh2;
        s3x2y=Re*(0.5*(S[n*(j+1)+i+1]+S[n*(j+1)+i-1]
          -S[n*(j-1)+i+1]-S[n*(j-1)+i-1])+S[n*(j-1)+i]-S[n*(j+1)+i])/dh3;
        s3xy2=Re*(0.5*(S[n*(j+1)+i+1]-S[n*(j+1)+i-1]
          +S[n*(j-1)+i+1]-S[n*(j-1)+i-1])+S[n*j+i-1]-S[n*j+i+1])/dh3;
        b[n*j+i]=1.+dh2*sy*sy/12.-dh2*s2xy/6.;
        c[n*j+i]=1.+dh2*sx*sx/12.+dh2*s2xy/6.;
        d[n*j+i]=sy+dh2*s3x2y/6.-dh2*sy*s2xy/12.+dh2*sx*s2y2/12.;
        e[n*j+i]=sx+dh2*s3xy2/6.-dh2*sy*s2x2/12.+dh2*sx*s2xy/12.;
        v2xy=(0.25*(W[n*(j+1)+i+1]-W[n*(j+1)+i-1]-W[n*(j-1)+i+1]+W[n*(j-1)+i-1]))/dh2;
        v3x2y=(0.5*(W[n*(j+1)+i+1]+W[n*(j+1)+i-1]-W[n*(j-1)+i+1]-W[n*(j-1)+i-1])
          +W[n*(j-1)+i]-W[n*(j+1)+i])/dh3;
        v3xy2=(0.5*(W[n*(j+1)+i+1]-W[n*(j+1)+i-1]+W[n*(j-1)+i+1]-W[n*(j-1)+i-1])
          +W[n*j+i-1]-W[n*j+i+1])/dh3;
        v4x2y2=((W[n*(j+1)+i+1]+W[n*(j+1)+i-1]+W[n*(j-1)+i+1]+W[n*(j-1)+i-1])
          -2.*(W[n*(j+1)+i]+W[n*j+i+1]+W[n*(j-1)+i]+W[n*j+i-1])+4.*W[n*j+i])/dh4;
        f[n*j+i]=dh2*v4x2y2/6.+dh2*sx*v3x2y/6.-dh2*sy*v3xy2/6.
          -dh2*sx*sy*v2xy/6.+dh2*s2x2*v2xy/6.-dh2*s2y2*v2xy/6.;
        }
      }
    for(i=1;i<n-1;i++)
      for(j=1;j<n-1;j++)
        R1[n*j+i]=W[n*j+i]+0.5*dt*(c[n*j+i]*(W[n*(j-1)+i]-2.*W[n*j+i]
          +W[n*(j+1)+i])/dh2+e[n*j+i]*(W[n*(j+1)+i]-W[n*(j-1)+i])/(2.*dh)+f[n*j+i]);
    for(j=1;j<n-1;j++)
      {
      for(i=1;i<n-1;i++)
        {
        aa[i]=-0.5*dt*b[n*j+i]/dh2-0.5*dt*d[n*j+i]/(2.*dh);
        bb[i]=1.+0.5*dt*b[n*j+i]*2./dh2;
        cc[i]=-0.5*dt*b[n*j+i]/dh2+0.5*dt*d[n*j+i]/(2.*dh);
        dd[i]=R1[n*j+i];
        }
      dd[1]+=-aa[1]*W[n*j];
      dd[n-2]+=-cc[n-2]*W[n*j+n-1];
      for(i=2;i<n-1;i++)
        {
        bb[i]+=-aa[i]*cc[i-1]/bb[i-1];
        dd[i]+=-aa[i]*dd[i-1]/bb[i-1];
        }
      W[n*j+n-2]=dd[n-2]/bb[n-2];
      for(i=n-3;i>=1;i--)
        W[n*j+i]=(dd[i]-W[n*j+i+1]*cc[i])/bb[i];
      }
    for(i=1;i<n-1;i++)
      for(j=1;j<n-1;j++)
        R2[n*j+i]=W[n*j+i]+0.5*dt*(b[n*j+i]*(W[n*j+i-1]-2.*W[n*j+i]
          +W[n*j+i+1])/dh2-d[n*j+i]*(W[n*j+i+1]-W[n*j+i-1])/(2.*dh)+f[n*j+i]);
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        aa[j]=-0.5*dt*c[n*j+i]/dh2+0.5*dt*e[n*j+i]/(2.*dh);
        bb[j]=1.+0.5*dt*c[n*j+i]*2./dh2;
        cc[j]=-0.5*dt*c[n*j+i]/dh2-0.5*dt*e[n*j+i]/(2.*dh);
        dd[j]=R2[n*j+i];
        }
      dd[1]+=-aa[1]*W[i];
      dd[n-2]+=-cc[n-2]*W[n*(n-1)+i];
      for(j=2;j<n-1;j++)
        {
        bb[j]+=-aa[j]*cc[j-1]/bb[j-1];
        dd[j]+=-aa[j]*dd[j-1]/bb[j-1];
        }
      W[n*(n-2)+i]=dd[n-2]/bb[n-2];
      for(j=n-3;j>=1;j--)
        W[n*j+i]=(dd[j]-W[n*(j+1)+i]*cc[j])/bb[j];
      }
    for(rr=0.,i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        r=fabs((S[n*j+i-1]-2.*S[n*j+i]+S[n*j+i+1])/dh2+(S[n*(j-1)+i]
          -2.*S[n*j+i]+S[n*(j+1)+i])/dh2+W[n*j+i]+a[n*j+i]);
        if(r>rr)
          rr=r;
        r=fabs((b[n*j+i]*(W[n*j+i-1]-2.*W[n*j+i]+W[n*j+i+1])/dh2
          +c[n*j+i]*(W[n*(j-1)+i]-2.*W[n*j+i]+W[n*(j+1)+i])/dh2
          -d[n*j+i]*(W[n*j+i+1]-W[n*j+i-1])/(2.*dh)
          +e[n*j+i]*(W[n*(j+1)+i]-W[n*(j-1)+i])/(2.*dh)+f[n*j+i])/Re);
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
    {
    for(i=0;i<n;i++)
      {
      k=n*i+j;
      if(j==0)
        U[k]=0.;
      else if(j==n-1)
        U[k]=1.;
      else
        U[k]=(S[k+n]-S[k-n])/(2.*dh);
      }
    }
  printf("calculating V\n");
  for(j=0;j<n;j++)
    {
    for(i=0;i<n;i++)
      {
      k=n*i+j;
      if(i==0||i==n-1)
        V[k]=0.;
      else
        V[k]=-(S[k+1]-S[k-1])/(2.*dh);
      }
    }
  printf("bounding psi\n");
  Zm=Zx=S[0];
  for(k=i=0;i<n;i++)
    {
    for(j=0;j<n;j++,k++)
      {
      if(S[k]<Zm)
        Zm=S[k];
      if(S[k]>Zx)
        Zx=S[k];
      }
    }
  WriteTB2("psi.tb2",S,-DBL_MAX,DBL_MAX);
  WriteTB2("omega.tb2",W,-5.,5.);
  WriteV2D("adi4.v2d");
  WriteTP2("adi4.tp2",Zm,Zx,Re);
  WritePLT("adi4.plt");
  return(0);
  }
