/* Vorticity-Stream Function Method
   for the cavity flow problem using
   SOR and 2nd order differences
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
double s[n*n];
double S[n*n];
double U[n*n];
double V[n*n];
double w[n*n];
double W[n*n];
double X[n];
double Y[n];

void WriteTB2(char*fname,double*Z,double Zm,double Zx)
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
  for(j=0;j<n;j++)
    {
    for(i=0;i<n;i++)
      {
      k=n*i+j;
      fprintf(fp,"%lG\n",fmax(Zm,fmin(Zx,Z[k])));
      }
    }
  fclose(fp);
  }

void WriteV2D(char*fname)
  {
  int i,j,k;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  for(j=0;j<n;j++)
    {
    for(i=0;i<n;i++)
      {
      k=n*i+j;
      fprintf(fp,"%lG %lG %lG %lG\n",X[i],Y[j],U[k],V[k]);
      }
    }
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
  fprintf(fp,"SOR2.V2D\n");
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
  for(j=0;j<n;j++)
    {
    for(i=0;i<n;i++)
      {
      k=n*i+j;
      fprintf(fp,"%lG %lG %lG %lG %lG %lG\n",X[i],Y[j],U[k],V[k],S[k],W[k]);
      }
    }
  fclose(fp);
  }

int main(int argc,char**argv,char**envp)
  {
  int done=-1,i,j,k;
  double beta=0.6,dh,dhdh,r,rr,Re,Sm,Sn,So,Sx,Sy,tol=1E-6,Wn,Wo,Wx,Wy;
  printf("solving the cavity problem using\n"
         "the vorticity-stream function method\n");
  Re=100.;
  dh=1./(n-1);
  dhdh=dh*dh;
  for(k=0;k<n;k++)
    X[k]=Y[k]=k/(n-1.);
  printf("initializing\n");
  for(k=i=0;i<n;i++)
    {
    for(j=0;j<n;j++,k++)
      S[k]=W[k]=tol;
    S[i]=W[i]=0.;
    S[n*i]=W[n*i]=0.;
    S[n*i+n-1]=W[n*i+n-1]=0.;
    S[n*(n-1)+i]=W[n*(n-1)+i]=0.;
    }
  printf("solving\n");
  do{
    memcpy(s,S,sizeof(s));
    memcpy(w,W,sizeof(w));
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        So=S[n*i+j];
        Sn=S[n*(i-1)+j]+S[n*(i+1)+j]+S[n*i+j-1]+S[n*i+j+1];
        S[n*i+j]=beta*(Sn+dhdh*W[n*i+j])/4.+(1.-beta)*So;
        }
      }
    for(k=1;k<n-1;k++)
      {
      Sn=S[n*(k-1)+1]+S[n*k+1]+S[n*(k+1)+1];
      Wn=2.*W[n*(k-1)]+2.*W[n*(k+1)]+W[n*(k-1)+1]+4.*W[n*k+1]+W[n*(k+1)+1];
      W[n*k]=-3.*Sn/(2.*dhdh)-Wn/8.;
      Sn=S[n*(k-1)+n-2]+S[n*k+n-2]+S[n*(k+1)+n-2];
      Wn=2.*W[n*(k-1)+n-1]+2.*W[n*(k+1)+n-1]+W[n*(k-1)+n-2]+4.*W[n*k+n-2]+W[n*(k+1)+n-2];
      W[n*k+n-1]=-4.5/dh-1.5*Sn/dhdh-Wn/8.;
      Sn=S[n*1+k-1]+S[n*1+k]+S[n*1+k+1];
      Wn=2.*W[k-1]+2.*W[k+1]+W[n*1+k-1]+4.*W[n*1+k]+W[n*1+k+1];
      W[k]=-1.5*Sn/dhdh-Wn/8.;
      Sn=S[n*(n-2)+k-1]+S[n*(n-2)+k]+S[n*(n-2)+k+1];
      Wn=2.*W[n*(n-1)+k-1]+2.*W[n*(n-1)+k+1]+W[n*(n-2)+k-1]+4.*W[n*(n-2)+k]+W[n*(n-2)+k+1];
      W[n*(n-1)+k]=-1.5*Sn/dhdh-Wn/8.;
      }
    Sn=S[n*1+1];
    Wn=2.*W[n*1]+2.*W[1]+W[n*1+1];
    W[0]=-3.*Sn/dhdh-Wn/4.;
    Sn=S[n*(n-2)+1];
    Wn=2.*W[n*(n-2)]+2.*W[n*(n-1)+1]+W[n*(n-2)+1];
    W[n*(n-1)]=-3.*Sn/dhdh-Wn/4.;
    Sn=S[n*(n-2)+n-2];
    Wn=2.*W[n*(n-2)+n-1]+2.*W[n*(n-1)+n-2]+W[n*(n-2)+n-2];
    W[n*(n-1)+n-1]=-4.5/dh-3.*Sn/dhdh-Wn/4.;
    Sn=S[n*1+n-2];
    Wn=2.*W[n*1+n-1]+2.*W[n-2]+W[n*1+n-2];
    W[n-1]=-4.5/dh-3.*Sn/dhdh-Wn/4.;
    for(i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        Wo=W[n*i+j];
        Wn=W[n*(i-1)+j]+W[n*(i+1)+j]+W[n*i+j-1]+W[n*i+j+1];
        Sx=S[n*i+j+1]-S[n*i+j-1];
        Wy=W[n*(i+1)+j]-W[n*(i-1)+j];
        Sy=S[n*(i+1)+j]-S[n*(i-1)+j];
        Wx=W[n*i+j+1]-W[n*i+j-1];
        W[n*i+j]=(4.*Wn-Re*(Sx*Wy-Sy*Wx))*beta/16.+(1.-beta)*Wo;
        }
      }
    for(rr=0.,i=1;i<n-1;i++)
      {
      for(j=1;j<n-1;j++)
        {
        Sn=S[n*(i-1)+j]+S[n*(i+1)+j]+S[n*i+j-1]+S[n*i+j+1];
        So=S[n*i+j];
        r=fabs(W[n*i+j]+(Sn-4.*So)/dhdh);
        if(r>rr)
          rr=r;
        Wn=W[n*(i-1)+j]+W[n*(i+1)+j]+W[n*i+j-1]+W[n*i+j+1];
        Wo=W[n*i+j];
        Sy=S[n*(i+1)+j]-S[n*(i-1)+j];
        Wx=W[n*i+j+1]-W[n*i+j-1];
        Sx=S[n*i+j+1]-S[n*i+j-1];
        Wy=W[n*(i+1)+j]-W[n*(i-1)+j];
        r=fabs(((Wn-4.*Wo)/Re+(Sy*Wx-Sx*Wy)/4.)/dhdh);
        if(r>rr)
          rr=r;
        r=fabs(S[n*i+j]-s[n*i+j]);
        if(r>rr)
          rr=r;
        r=fabs(W[n*i+j]-w[n*i+j]);
        if(r>rr)
          rr=r;
        r=fabs((S[n*i+j]-s[n*i+j])/s[n*i+j]);
        if(r>rr)
          rr=r;
        r=fabs((W[n*i+j]-w[n*i+j])/w[n*i+j]);
        if(r>rr)
          rr=r;
        }
      }
    i=(int)(100.*pow(tol/rr,0.1)+0.5);
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
        U[k]=(S[k+1]-S[k-1])/(2.*dh);
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
        V[k]=-(S[k+n]-S[k-n])/(2.*dh);
      }
    }
  printf("bounding psi\n");
  Sm=Sx=S[0];
  for(k=i=0;i<n;i++)
    {
    for(j=0;j<n;j++,k++)
      {
      if(S[k]<Sm)
        Sm=S[k];
      if(S[k]>Sx)
        Sx=S[k];
      }
    }
  WriteTB2("psi.tb2",S,-DBL_MAX,DBL_MAX);
  WriteTB2("omega.tb2",W,-5.,5.);
  WriteV2D("sor2.v2d");
  WriteTP2("sor2.tp2",Sm,Sx,Re);
  WritePLT("sor2.plt");
  return(0);
  }
