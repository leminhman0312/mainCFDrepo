/* Navier-Stokes Equations in 2D
   Pressure from Poisson with SOR
   dudley.benton@gmail.com
 */

#define _CRT_SECURE_NO_DEPRECATE
#include <conio.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#define CASE 7
#undef DUCT

#ifdef SCROLL
  #define Nx 80
  #define Ny 80
  double lX=6.;
  double lY=6.;
#elif(defined(DUCT))
  #define Nx 200
  #define Ny  80
  double lX=10.;
  double lY= 4.;
#else
  #define Nx  80
  #define Ny 120
  double lX=4.;    /* X length */
  double lY=6.;    /* Y length */
#endif
double dT=0.003; /* time step */
double Re=100.;  /* Reynolds number */
double Ub=1.;    /* X velocity along boundary */
double Vb=0.;    /* Y velocity along boundary */

typedef unsigned char byte;

byte   B[Nx+2][Ny+2];
double F[Nx+2][Ny+2];
double G[Nx+2][Ny+2];
double P[Nx+2][Ny+2];
double R[Nx+2][Ny+2];
double S[Nx+2][Ny+2];
double U[Nx+2][Ny+2];
double V[Nx+2][Ny+2];
double W[Nx+2][Ny+2];
double dX,dXdX,dY,dYdY;

byte B_O=0x00; /* obstacle/boundary */
byte B_N=0x01; /* active north */
byte B_S=0x02; /* active south */
byte B_W=0x04; /* active west */
byte B_E=0x08; /* active east */
byte B_F=0x10; /* fluid cell */
byte B_NW;     /* north or west */
byte B_SW;     /* south or west */
byte B_NE;     /* north or east */
byte B_SE;     /* south or east */
byte B_NSEW;   /* north, south, east, or west */

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

double sq(double x)
  {
  return(x*x);
  }

void BoundaryConditions()
  {
  int i,j;
  for(j=0;j<=Ny+1;j++)
    {
    U[0][j]=U[1][j];
    V[0][j]=V[1][j];
    U[Nx][j]=U[Nx-1][j];
    V[Nx+1][j]=V[Nx][j];
    }
  for(i=0;i<=Nx+1;i++)
    {
    V[i][0]=V[i][Ny]=0.;
    U[i][Ny+1]=U[i][Ny];
    U[i][0]=U[i][1];
    }
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(B[i][j]&B_NSEW)
        {
        if(B[i][j]==B_N)
          U[i][j]=-U[i][j+1];
        else if(B[i][j]==B_E)
          U[i][j]=0.;
        else if(B[i][j]==B_NE)
          U[i][j]=0.;
        else if(B[i][j]==B_SE)
          U[i][j]=0.;
        else if(B[i][j]==B_NW)
          U[i][j]=-U[i][j+1];
        else if(B[i][j]==B_S)
          U[i][j]=-U[i][j-1];
        else if(B[i][j]==B_SW)
          U[i][j]=-U[i][j-1];
        }
      }
    }
  for(i=0;i<=(Nx-1);i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(B[i+1][j]&B_NSEW)
        {
        if(B[i+1][j]==B_N)
          U[i][j]=-U[i][j+1];
        else if(B[i+1][j]==B_W)
          U[i][j]=0.;
        else if(B[i+1][j]==B_NE)
          U[i][j]=-U[i][j+1];
        else if(B[i+1][j]==B_SW)
          U[i][j]=0.;
        else if(B[i+1][j]==B_NW)
          U[i][j]=0.;
        else if(B[i+1][j]==B_S)
          U[i][j]=-U[i][j-1];
        else if(B[i+1][j]==B_SE)
          U[i][j]=-U[i][j-1];
        }
      }
    }
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(B[i][j]&B_NSEW)
        {
        if(B[i][j]==B_N)
          V[i][j]=0.;
        else if(B[i][j]==B_E)
          V[i][j]=-V[i+1][j];
        else if(B[i][j]==B_NE)
          V[i][j]=0.;
        else if(B[i][j]==B_SE)
          V[i][j]=-V[i+1][j];
        else if(B[i][j]==B_NW)
          V[i][j]=0.;
        else if(B[i][j]==B_W)
          V[i][j]=-V[i-1][j];
        else if(B[i][j]==B_SW)
          V[i][j]=-V[i-1][j];
        }
      }
    }
  for(i=1;i<=Nx;i++)
    {
    for(j=0;j<=(Ny-1);j++)
      {
      if(B[i][j+1]&B_NSEW)
        {
        if(B[i][j+1]==B_E)
          V[i][j]=-V[i+1][j];
        else if(B[i][j+1]==B_S)
          V[i][j]=0.;
        else if(B[i][j+1]==B_NE)
          V[i][j]=-V[i+1][j];
        else if(B[i][j+1]==B_SE)
          V[i][j]=0.;
        else if(B[i][j+1]==B_SW)
          V[i][j]=0.;
        else if(B[i][j+1]==B_W)
          V[i][j]=-V[i-1][j];
        else if(B[i][j+1]==B_NW)
          V[i][j]=-V[i-1][j];
        }
      }
    }
#if(defined(DUCT))
  for(j=0;j<=Ny;j++)
    U[0][j]=U[1][j]=U[Nx-1][j]=U[Nx][j]=Ub;
#else
  V[0][0]=2.*Vb-V[1][0];
  for(j=1;j<=Ny;j++)
    {
    U[0][j]=Ub;
    V[0][j]=2.*Vb-V[1][j];
    }
#endif
  }

void Vorticity()
  {
  int i,j;
  printf("computing vorticity\n");
  for(j=1;j<=Ny-1;j++)
    for(i=1;i<=Nx-1;i++)
      if((B[i][j]&B_F)&&(B[i+1][j]&B_F)&&(B[i][j+1]&B_F)&&(B[i+1][j+1]&B_F))
        W[i][j]=(U[i][j+1]-U[i][j])/dY-(V[i+1][j]-V[i][j])/dX;
  }

void StreamFunction()
  {
  int i,j;
  printf("computing stream funciton\n");
  for(j=1;j<=Ny;j++)
    {
    S[1][j]=S[1][j-1]+U[1][j]*dY;
    for(i=2;i<=Nx;i++)
      S[i][j]=S[i-1][j]-V[i][j]*dX;
    }
  }

void ProvisionalVelocity()
  {
  int i,j;
  double gamma=0.9,UUx,UVx,UVy,Uxx,VVy,Vyy;
  for(i=1;i<=Nx-1;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if((B[i][j]&B_F)&&(B[i+1][j]&B_F))
        {
        UUx=((U[i][j]+U[i+1][j])*(U[i][j]+U[i+1][j])
          +gamma*fabs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])
          -(U[i-1][j]+U[i][j])*(U[i-1][j]+U[i][j])
          -gamma*fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))/(4.*dX);
        UVy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])
          +gamma*fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])
          -(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])
          -gamma*fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))/(4.*dY);
        Uxx=(U[i+1][j]-2.*U[i][j]+U[i-1][j])/dX/dX
          +(U[i][j+1]-2.*U[i][j]+U[i][j-1])/dY/dY;
        F[i][j]=U[i][j]+dT*(Uxx/Re-UUx-UVy);
        }
      else
        F[i][j]=U[i][j];
      }
    }
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny-1;j++)
      {
      if((B[i][j]&B_F)&&(B[i][j+1]&B_F))
        {
        UVx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])
          +gamma*fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])
          -(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])
          -gamma*fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))/(4.*dX);
        VVy=((V[i][j]+V[i][j+1])*(V[i][j]+V[i][j+1])
          +gamma*fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])
          -(V[i][j-1]+V[i][j])*(V[i][j-1]+V[i][j])
          -gamma*fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))/(4.*dY);
        Vyy=(V[i+1][j]-2.*V[i][j]+V[i-1][j])/dX/dX
          +(V[i][j+1]-2.*V[i][j]+V[i][j-1])/dY/dY;
        G[i][j]=V[i][j]+dT*(Vyy/Re-UVx-VVy);
        }
      else
        G[i][j]=V[i][j];
      }
    }
  for(j=1;j<=Ny;j++)
    {
    F[0][j]=U[0][j];
    F[Nx][j]=U[Nx][j];
    }
  for(i=1;i<=Nx;i++)
    {
    G[i][0]=V[i][0];
    G[i][Ny]=V[i][Ny];
    }
  }

void PressureEquation()
  {
  int i,j;
  for(i=1;i<=Nx;i++)
    for(j=1;j<=Ny;j++)
      if(B[i][j]&B_F)
        R[i][j]=((F[i][j]-F[i-1][j])/dX+(G[i][j]-G[i][j-1])/dY)/dT;
  }

void PoissonSolver()
  {
  int i,iter,j,k;
  double a,b,eE,eN,eS,eW,res,Rij,sor=5./3.,tol=0.001;
  a=-sor*(dXdX*dYdY)/(dXdX+dYdY)/2.;
  for(iter=0;iter<100;iter++)
    {
    for(k=0;k<2;k++)
      {
      for(i=1;i<=Nx;i++)
        {
        for(j=1;j<=Ny;j++)
          {
          if((i+j)%2!=k)
            continue;
          if(B[i][j]==(B_F|B_NSEW))
            P[i][j]=(1.-sor)*P[i][j]-a*((P[i+1][j]+P[i-1][j])/dXdX
              +(P[i][j+1]+P[i][j-1])/dYdY-R[i][j]);
          else if(B[i][j]&B_F)
            {
            eE=(B[i+1][j]&B_F)?1.:0.;
            eW=(B[i-1][j]&B_F)?1.:0.;
            eN=(B[i][j+1]&B_F)?1.:0.;
            eS=(B[i][j-1]&B_F)?1.:0.;
            b=sor/((eE+eW)/dXdX+(eN+eS)/dYdY);
            P[i][j]=b*((eE*P[i+1][j]+eW*P[i-1][j])/dXdX
                      +(eN*P[i][j+1]+eS*P[i][j-1])/dYdY
                      -R[i][j])+(1.-sor)*P[i][j];
            }
          }
        }
      }
    for(res=0.,i=1;i<=Nx;i++)
      {
      for(j=1;j<=Ny;j++)
        {
        if(B[i][j]&B_F)
          {
          eE=(B[i+1][j]&B_F)?1.:0.;
          eW=(B[i-1][j]&B_F)?1.:0.;
          eN=(B[i][j+1]&B_F)?1.:0.;
          eS=(B[i][j-1]&B_F)?1.:0.;
          Rij=(eE*(P[i+1][j]-P[i][j])-eW*(P[i][j]-P[i-1][j]))/dXdX
             +(eN*(P[i][j+1]-P[i][j])-eS*(P[i][j]-P[i][j-1]))/dYdY;
          res+=sq(Rij-R[i][j]);
          }
        }
      }
    res=sqrt(res/(Nx*Ny));
    if(res<tol)
      break;
    }
  }

void UpdateVelocity()
  {
  int i,j;
  for(i=1;i<=Nx-1;i++)
    for(j=1;j<=Ny;j++)
      if((B[i][j]&B_F)&&(B[i+1][j]&B_F))
        U[i][j]=F[i][j]-(P[i+1][j]-P[i][j])*dT/dX;
  for(i=1;i<=Nx;i++)
    for(j=1;j<=Ny-1;j++)
      if((B[i][j]&B_F)&&(B[i][j+1]&B_F))
        V[i][j]=G[i][j]-(P[i][j+1]-P[i][j])*dT/dY;
  }

void Courant()
  {
  int i,j;
  double dU,dV,dW,Ux,Vx;
  Ux=Vx=FLT_EPSILON;
  for(i=0;i<=Nx+1;i++)
    for(j=1;j<=Ny+1;j++)
      Ux=fmax(fabs(U[i][j]),Ux);
  for(i=1;i<=Nx+1;i++)
    for(j=0;j<=Ny+1;j++)
      Vx=fmax(fabs(V[i][j]),Vx);
  dU=dX/Ux;
  dV=dY/Vx;
  dW=1./(1./(dX*dX)+1./(dY*dY))*Re/2.;
  if(dU<dV)
    dT=fmin(dU,dW);
  else
    dT=fmin(dV,dW);
  dT/=2.;
  }

void BoundaryFlags()
  {
  int i,j;
  B_NW=B_N|B_W;
  B_SW=B_S|B_W;
  B_NE=B_N|B_E;
  B_SE=B_S|B_E;
  B_NSEW=B_N|B_S|B_E|B_W;
  for(i=0;i<=Nx+1;i++)
    {
    B[i][0]=B_O;
    B[i][Ny+1]=B_O;
    }
  for(j=1;j<=Ny;j++)
    {
    B[0][j]=B_O;
    B[Nx+1][j]=B_O;
    }
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(!(B[i][j]&B_F))
        {
        if(B[i-1][j]&B_F)
          B[i][j]|=B_W;
        if(B[i+1][j]&B_F)
          B[i][j]|=B_E;
        if(B[i][j-1]&B_F)
          B[i][j]|=B_S;
        if(B[i][j+1]&B_F)
          B[i][j]|=B_N;
        }
      }
    }
  }

#if(defined(SCROLL)||defined(DUCT))
typedef struct {double x,y;}POLY;
int InsidePolygon(POLY*P,int n,double X,double Y)
  {
  int above1,above2,i,right=0;
  double X1,X2,Y1,Y2;
  X2=P[n-1].x;
  Y2=P[n-1].y;
  above2=Y2>Y?1:0;
  for(i=0;i<n;i++)
    {
    X1=X2;
    Y1=Y2;
    X2=P[i].x;
    Y2=P[i].y;
    above1=above2;
    above2=Y2>Y?1:0;
    if(above1==above2)
      continue;
    if(X1>X&&X2>X)
      right++;
    else if(Y1<Y2)
      {
      if((X-X1)*(Y2-Y1)<(X2-X1)*(Y-Y1))
        right++;
      }
    else if(Y1>Y2)
      {
      if((X-X1)*(Y2-Y1)>(X2-X1)*(Y-Y1))
        right++;
      }
    }
  return(right&1);
  }
#endif

#ifdef SCROLL
POLY bndy[]={ {0.00,0.00},{0.00,1.52},{1.88,1.52},{1.57,1.66},{1.30,1.85},
  {1.06,2.09},{0.87,2.37},{0.72,2.67},{0.64,3.00},{0.61,3.33},{0.64,3.67},
  {0.72,4.00},{0.87,4.30},{1.06,4.58},{1.30,4.81},{1.57,5.01},{1.88,5.15},
  {2.20,5.24},{2.54,5.27},{2.88,5.24},{3.20,5.15},{3.51,5.01},{3.78,4.81},
  {4.02,4.58},{4.21,4.30},{6.02,4.30},{6.02,3.26},{4.31,3.26},{4.27,2.86},
  {4.21,2.46},{4.10,2.07},{3.94,1.69},{3.72,1.35},{3.45,1.05},{3.14,0.80},
  {2.80,0.57},{2.45,0.36},{2.09,0.18},{1.72,0.00}};
#elif(defined(DUCT))
POLY bndy[]={{0,2},{2,2},{2,4},{8,4},{8,2},{10,2},{10,0},{6,0},{6,2},{4,2},{4,0},{0,0}};
#endif

void Obstacles()
  {
  int i,j;
#if(defined(SCROLL)||defined(DUCT))
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(InsidePolygon(bndy,sizeof(bndy)/sizeof(bndy[0]),(i-1)*dX,(j-1)*dY))
        B[i][j]=B_F;
      else
        B[i][j]=B_O;
      }
    }
#elif CASE==1
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(hypot(i-Nx/4.,j-Ny/2.)<Ny/8.)
        B[i][j]=B_O;
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==2
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(hypot(i-Nx/4.,j-Ny/2.)<Ny/8.)
        B[i][j]=B_O;
      else if(hypot(i-3.*Nx/4.,j-Ny/2.)<Ny/8.)
        B[i][j]=B_O;
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==3
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(hypot((i-Nx/2.)/3.,(j-Ny/2.)*i*4./Nx)<Ny/8.)
        B[i][j]=B_O;
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==4
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(abs(i-Nx/3)<2&&abs(j-Ny/2)<Ny/6)
        B[i][j]=B_O;
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==5
  for(i=1;i<=Nx;i++)
    {
    for(j=1;j<=Ny;j++)
      {
      if(i>=Nx/4&&i<=Nx/2&&abs(j-Ny/2)<=abs(i-Nx/4))
        B[i][j]=B_O;
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==6
  for(j=1;j<=Ny;j++)
    {
    for(i=1;i<=Nx;i++)
      {
      if(i>Nx/2)
        B[i][j]=B[Nx-i][j];
      else if(i>=Nx/4&&i<=Nx/2&&abs(j-Ny/4)<=abs(i-Nx/4))
        B[i][j]=B_O;
      else if(j>Ny/2)
        B[i][j]=B[i][Ny-j];
      else
        B[i][j]=B_F;
      }
    }
#elif CASE==7
  for(j=1;j<=Ny;j++)
    {
    for(i=1;i<=Nx;i++)
      {
      if(hypot(i-Nx/2.,j-Ny/4.)<Nx/4.)
        B[i][j]=B_O;
      else if(j>Ny/2)
        B[i][j]=B[i][Ny-j];
      else
        B[i][j]=B_F;
      }
    }
#else
#error unexpected case CASE
#endif
  }

void Initialize()
  {
  int i,j;
  printf("initializing\n");
  dX=lX/(Nx-1);
  dXdX=dX*dX;
  dY=lY/(Ny-1);
  dYdY=dY*dY;
  for(i=0;i<=Nx+1;i++)
    {
    for(j=0;j<=Ny+1;j++)
      {
      U[i][j]=Ub;
      V[i][j]=Vb;
      P[i][j]=0.;
      }
    }
  printf("Re=%lG\n",Re);
  }

void WriteTB2(char*fname,double Z[Nx+2][Ny+2],double Zo)
  {
  int i,j;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"%i\n",Nx);
  for(i=1;i<Nx+1;i++)
    fprintf(fp,"%lG\n",(i-1)*dX);
  fprintf(fp,"%i\n",Ny);
  for(j=1;j<Ny+1;j++)
    fprintf(fp,"%lG\n",(j-1)*dY);
  fprintf(fp,"%i\n",Nx*Ny);
  for(j=1;j<Ny+1;j++)
    for(i=1;i<Nx+1;i++)
      fprintf(fp,"%lG\n",Z[i][j]-Zo);
  fclose(fp);
  }

void WriteV2D(char*fname)
  {
  int i,j;
  FILE*fp;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  for(j=1;j<Ny+1;j++)
    for(i=1;i<Nx+1;i++)
      if(B[i][j]&B_F)
        fprintf(fp,"%lG %lG %lG %lG\n",i*dX,j*dY,U[i][j],V[i][j]);
  fclose(fp);
  }

void WriteTP2(char*fname,char*dname,double Z[Nx+2][Ny+2])
  {
  int i,j;
  double Zm,Zo,Zx;
  FILE*fp;
  Zm=Zx=Z[1][1];
  for(j=1;j<Ny+1;j++)
    {
    for(i=1;i<Nx+1;i++)
      {
      if(B[i][j]&B_F)
        {
        if(Z[i][j]<Zm)
          Zm=Z[i][j];
        if(Z[i][j]>Zx)
          Zx=Z[i][j];
        }
      }
    }
  Zo=Zm;
  WriteTB2(dname,Z,Zo);
  Zx-=Zo;
  Zm-=Zo;
  printf("writing: %s\n",fname);
  if((fp=fopen(fname,"wt"))==NULL)
    return;
  fprintf(fp,"TP2\n");
  fprintf(fp,"!PLOT\n");
  fprintf(fp,"TICKS\n");
  fprintf(fp,"EQUAL\n");
  fprintf(fp,"LEGENDS\n");
  fprintf(fp,"-1 3 1\n");
  fprintf(fp,"COLOR\n");
  fprintf(fp,"%lG %lG 2\n",Zm,Zx);
  fprintf(fp,"RAINBOW\n");
  fprintf(fp,"ZTABLE\n");
  fprintf(fp,"NOTITLE\n");
  fprintf(fp,"X\n");
  fprintf(fp,"0 1 %lG 1 0\n",lX);
  fprintf(fp,"Y\n");
  fprintf(fp,"0 1 %lG 1 0\n",lY);
  fprintf(fp,"%s\n",dname);
  fprintf(fp,"NOCOLOR\n");
  fprintf(fp,"ARROW\n");
  fprintf(fp,"0 0\n");
  fprintf(fp,"NEWFILE\n");
  fprintf(fp,"CART2D.V2D\n");
  fprintf(fp,"%lG %lG .Re=%lG\n",lX,lY/100.,Re);
  fprintf(fp,"!END\n");
  fclose(fp);
  }

void WritePLT(char*fname)
  {
  int i,j;
  FILE*fp;
  printf("output file: %s\n",fname);
  if((fp=fopen(fname,"wb"))==NULL)
    return;
  fprintf(fp,"TITLE=\"2D Fluid Flow\"\n");
  fprintf(fp,"VARIABLES=\"X\" \"Y\" \"P\" \"U\" \"V\" \"psi\" \"omega\"\n");
  fprintf(fp,"ZONE T=\"boundary\"  I=%i, J=%i, K=1, F=POINT\n",Ny,Nx);
  fprintf(fp,"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )\n");
  for(j=1;j<Ny+1;j++)
    for(i=1;i<Nx+1;i++)
      fprintf(fp,"%lG %lG %lG %lG %lG %lG %lG\n",(i-1)*dX,(j-1)*dY,P[i][j],U[i][j],V[i][j],S[i][j],W[i][j]);
  fclose(fp);
  }

int main(int argc,char**argv,char**envp)
  {
  double T,Tend=40.;
  printf("2D Navier-Stokes Equations\n");
  printf("structured grid finite difference method\n");
  Initialize();
  Obstacles();
  BoundaryFlags();
  BoundaryConditions();
  printf("solving\n");
  for(T=0.;T<Tend;T+=dT)
    {
    Courant();
    ProvisionalVelocity();
    PressureEquation();
    PoissonSolver();
    UpdateVelocity();
    BoundaryConditions();
    fprintf(stderr,"\r%.0lf%%",T*100./Tend);
    if(_kbhit())
      {
      _getch();
      break;
      }
    }
  fprintf(stderr,"\rdone\n");
  Vorticity();
  StreamFunction();
  WriteTB2("psi.tb2",S,0.);
  WriteTB2("omega.tb2",W,0.);
  WriteV2D("cart2d.v2d");
  WriteTP2("cart2d.tp2","p.tb2",P);
  WritePLT("cart2d.plt");
  return(0);
  }
