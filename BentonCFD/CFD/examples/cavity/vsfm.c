/* Vorticity-Stream Function Method
   for the cavity flow problem
   minimalist low order accuracy
   dudley.benton@gmail.com
 */

#define _CRT_SECURE_NO_DEPRECATE
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

#define Nx 60
#define Ny 60

int b[Nx*Ny];
double S[Nx*Ny];
double U[Nx*Ny];
double W[Nx*Ny];
double V[Nx*Ny];
double X[Nx];
double Y[Ny];

void WriteTB2(char*fname,double*Z,double Zm,double Zx)
  {
  int i,j,k;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"%i\n",Nx);
  for(j=0;j<Nx;j++)
    fprintf(fp,"%lG\n",X[j]);
  fprintf(fp,"%i\n",Ny);
  for(i=0;i<Ny;i++)
    fprintf(fp,"%lG\n",Y[i]);
  fprintf(fp,"%i\n",Nx*Ny);
  for(k=i=0;i<Ny;i++)
    for(j=0;j<Nx;j++,k++)
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
  for(k=i=0;i<Ny;i++)
    for(j=0;j<Nx;j++,k++)
      fprintf(fp,"%lG %lG %lG %lG\n",X[j],Y[i],U[k],V[k]);
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
  fprintf(fp,"%lG %lG\n",X[0],X[Nx-1]);
  fprintf(fp,"Y\n");
  fprintf(fp,"%lG %lG\n",Y[0],Y[Ny-1]);
  fprintf(fp,"PSI.TB2\n");
  fprintf(fp,"NOCOLOR\n");
  fprintf(fp,"ARROW\n");
  fprintf(fp,"0 0\n");
  fprintf(fp,"NEWFILE\n");
  fprintf(fp,"VSFM.V2D\n");
  fprintf(fp,"%lG %lG .    Re=%lG\n",X[Nx-1],Y[0]+(Y[Ny-1]-Y[0])/100.,Re);
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
  fprintf(fp,"ZONE T=\"results\"  I=%i, J=%i, K=1, F=POINT\n",Nx,Ny);
  fprintf(fp,"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n");
  for(k=i=0;i<Ny;i++)
    for(j=0;j<Nx;j++,k++)
      fprintf(fp,"%lG %lG %lG %lG %lG %lG\n",X[j],Y[i],U[k],V[k],S[k],W[k]);
  fclose(fp);
  }

int main(int argc,char**argv,char**envp)
  {
  int i,j,k,l,m;
  double dt,dX,dY,H,L,nu,Re,Sb[3],Sm,Sx,Sy,t,T,Uo,Wx,Wxx,Wy,Wyy;
  printf("solving the cavity problem using\n"
         "the vorticity-stream function method\n");
  H=L=1.;
  Uo=1.;
  nu=0.1;
  Re=Uo*H/nu;
  printf("Re=%lG\n",Re);
  if(Re>75.)
    dt=0.0001;
  else
    dt=0.001;
  T=5.;
  printf("initializing\n");
  dX=L/(Nx-1);
  dY=H/(Ny-1);
  for(j=0;j<Nx;j++)
    X[j]=j*dX;
  for(i=0;i<Ny;i++)
    Y[i]=i*dY;
  for(k=i=0;i<Ny;i++)
    {
    for(j=0;j<Nx;j++,k++)
      {
      W[k]=0.1;
#if(defined(CASE12)|defined(CASE13))
      if(hypot(j-Nx/4.,i-Ny/4.)<=7.)
        b[k]=1;
      else if(hypot(j-3.*Nx/4.,i-Ny/4.)<=7.)
        b[k]=2;
      else if(hypot(j-Nx/4.,i-3.*Ny/4.)<=7.)
        b[k]=3;
      else if(hypot(j-3.*Nx/4.,i-3.*Ny/4.)<=7.)
        b[k]=4;
#else
      if(hypot(j-Nx/4.,i-Ny/4.)<=7.)
        b[k]=1;
      else if(hypot(j-3.*Nx/4.,i-Ny/4.)<=7.)
        b[k]=2;
      else if(hypot(j-Nx/2.,i-3.*Ny/4.)<=7.)
        b[k]=3;
      }
#endif
    }
  printf("time-stepping\n");
  for(t=0.;t<T;t+=dt)
    {
    for(k=i=0;i<Ny;i++)
      {
      for(j=0;j<Nx;j++,k++)
        {
        if(i!=0&&i!=Ny-1&&j!=0&&j!=Nx-1)
          {
          Sx=(S[k+1]-S[k-1])/(2.*dX);
          Wx=(W[k+1]-W[k-1])/(2.*dX);
          Sy=(S[k+Nx]-S[k-Nx])/(2.*dY);
          Wy=(W[k+Nx]-W[k-Nx])/(2.*dY);
          Wxx=(W[k+1]-2.*W[k]+W[k-1])/(dX*dX);
          Wyy=(W[k+Nx]-2.*W[k]+W[k-Nx])/(dY*dY);
          W[k]+=dt*(Sx*Wy-Sy*Wx+nu*(Wxx+Wyy));
          }
        else if(i==0&&j!=0&&j!=Nx-1)            /* bottom */
#if(defined(CASE12))
          W[k]=2.*(S[k]-S[k+Nx])/(dY*dY)+2.*Uo/dY;
#elif(defined(CASE13))
          W[k]=2.*(S[k]-S[k+Nx])/(dY*dY)-2.*Uo/dY;
#else
          W[k]=2.*(S[k]-S[k+Nx])/(dY*dY);
#endif
        else if(i==Ny-1&&j!=0&&j!=Nx-1)         /* top */
          W[k]=2.*(S[k]-S[k-Nx])/(dY*dY)-2.*Uo/dY;
        else if(j==0&&i!=0&&i!=Ny-1)            /* left */
          W[k]=2.*(S[k]-S[k+1])/(dX*dX);
        else if(j==Nx-1&&i!=0&&i!=Ny-1)         /* right */
          W[k]=2.*(S[k]-S[k-1])/(dX*dX);
        else if(j==0&&i==Ny-1)                  /* top left */
          W[k]=W[k+1];
        else if(j==Nx-1&&i==Ny-1)               /* top right */
          W[k]=W[k-1];
        }
      }
    for(k=i=0;i<Ny;i++)
      for(j=0;j<Nx;j++,k++)
        if(i!=0&&i!=Ny-1&&j!=0&&j!=Nx-1)
          if(b[k]==0)
            S[k]=(2.*dX*dX*dY*dY*W[k]+dY*dY*(S[k+1]+S[k-1])
              +dX*dX*(S[k+Nx]+S[k-Nx]))/2./(dX*dX+dY*dY);
    for(l=0;l<sizeof(Sb)/sizeof(Sb[0]);l++)
      {
      for(Sb[l]=m=0,i=1;i<Ny-1;i++)
        {
        for(j=1;j<Nx-1;j++)
          {
          k=Nx*i+j;
          if(b[k]!=l+1)
            continue;
          if(b[k-1]==0)
            {
            Sb[l]+=S[k-1];
            m++;
            }
          if(b[k+1]==0)
            {
            Sb[l]+=S[k+1];
            m++;
            }
          if(b[k-Nx]==0)
            {
            Sb[l]+=S[k-Nx];
            m++;
            }
          if(b[k+Nx]==0)
            {
            Sb[l]+=S[k+Nx];
            m++;
            }
          }
        }
      if(m==0)
        {
        printf("boundary %i is isolated\n",l+1);
        exit(1);
        }
      Sb[l]/=m;
      }
    for(k=i=0;i<Ny;i++)
      for(j=0;j<Nx;j++)
        if(b[k]==l+1)
          S[k]=Sb[l];
    fprintf(stderr,"\r%.0lf%%",t*100./T);
    }
  printf("\rdone\n");
  printf("Sb=%lG",Sb[0]);
  for(l=1;l<sizeof(Sb)/sizeof(Sb[0]);l++)
    printf(",%lG",Sb[l]);
  printf("\n");
  printf("calculating U\n");
  for(k=i=0;i<Ny;i++)
    for(j=0;j<Nx;j++,k++)
      if(i==0)
        U[k]=0.;
      else if(i==Ny-1)
        U[k]=Uo;
      else
        U[k]=(S[k+Nx]-S[k-Nx])/(2.*dY);
  printf("calculating V\n");
  for(k=i=0;i<Ny;i++)
    for(j=0;j<Nx;j++,k++)
      if(j==0||j==Nx-1||i==Ny-1)
        V[k]=0.;
      else
        V[k]=-(S[k+1]-S[k-1])/(2.*dX);
  printf("bounding psi\n");
  Sm=Sx=S[0];
  for(k=i=0;i<Ny;i++)
    {
    for(j=0;j<Nx;j++,k++)
      {
      if(S[k]<Sm)
        Sm=S[k];
      if(S[k]>Sx)
        Sx=S[k];
      }
    }
  WriteTB2("psi.tb2",S,-DBL_MAX,DBL_MAX);
  WriteTB2("omega.tb2",W,-5.,5.);
  WriteV2D("vsfm.v2d");
  WriteTP2("vsfm.tp2",Sm,Sx,Re);
  WritePLT("vsfm.plt");
  return(0);
  }
