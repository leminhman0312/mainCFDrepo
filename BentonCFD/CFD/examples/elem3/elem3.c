/*******************************************************************************

                   EEEEEEE LLLL    EEEEEEE MM   MM  33333
                    EE   E  LL      EE   E MMM MMM 33   33
                    EE      LL      EE     MMMMMMM      33
                    EEEE    LL      EEEE   MM M MM    333
                    EE      LL      EE     MM   MM      33
                    EE   E  LL  LL  EE   E MM   MM 33   33
                   EEEEEEE LLLLLLL EEEEEEE MM   MM  33333
                   ---------------------------------------
                        Triangular Element Generator
                            by Dudley J. Benton

  This program reads a data file containing a sequence of nodes which
  defing a polygon. This polygon is subdivided into triangular regions
  such that the longest side is no longer than a specified length. The
  subdivision is such that the all of the triangles are as near to
  equilateral as possible. An output file is created which contains the
  elements. The optional commands are as follows:

             ELEM3 input_file [element_size]

  V1.00 (10/07/1993): original
  V1.10 (01/19/1995): create output files *.NOD and *.ELM
  V1.20 (06/04/1996): clean up and optimize some loops
  V1.21 (04/18/1997): modify .PLT file for TP2 and change .NOD->.NDE
  V1.22 (10/22/1999): modify font output format
  V2.00 (01/23/2001): eliminate graphics and font file and convert to 32-bit

*******************************************************************************/

char*Version="2.00";

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <malloc.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

void Abort(int code,char*format,...) /* fatal error exit */
  {
  va_list arg_marker;
  fprintf(stderr,"\nFatal Error Code: %i\n",code);
  va_start(arg_marker,format);
  vfprintf(stderr,format,arg_marker);
  fprintf(stderr,"\n");
  _fcloseall();
  exit(1);
  }

void*allocate(int line,size_t count,size_t size)
  {
  void*ptr;
  if((ptr=calloc(count,size))==NULL)
    Abort(line,"can't allocate memory [%u*%u=%u]",count,size,count*size);
  return(ptr);
  }

double fmax(double x,double y)
  {
  if(x>y)
    return(x);
  return(y);
  }

double fmin(double x,double y)
  {
  if(x<y)
    return(x);
  return(y);
  }

char cbuf[256];
FILE*Finp;
FILE*Fout;
int*Ie;
int*Je;
int*Ke;
int Kp;
int*Lp;
int Mp;
int Nb;
int Ne;
int Nn;
int Np;
double*Xb;
double Xx;
double Xm;
double*Xn;
double*Xp;
double*Yb;
double Yx;
double Ym;
double*Yn;
double*Yp;

double PolygonArea(double*Xp,double*Yp,int Np)
  {
  int i,j;
  double A;
  for(i=Np-1,A=j=0;j<Np;i=j++)
    A+=(Xp[j]-Xp[i])*(Yp[i]+Yp[j]);
  return(fabs(A)/2.);
  }

double DistanceToLine(double X1,double Y1,double X2,double Y2,double X3,double Y3)
  {
  double B,R1,R2;
  B=hypot(X2-X1,Y2-Y1);
  R1=hypot(X3-X1,Y3-Y1);
  R2=hypot(X3-X2,Y3-Y2);
  if(((X3-X1)*(X3-X2)<=0.||(Y3-Y1)*(Y3-Y2)<=0.)&&B>DBL_EPSILON)
    return(fabs((X3-X1)*(Y2-Y1)-(X2-X1)*(Y3-Y1))/B);
  return(fmin(R1,R2));
  }

int AddNode(double X,double Y)
  {
  int i;
  double*xn,*yn;
  for(i=0;i<Nn;i++)
    if(fabs(X-Xn[i])<(Xx-Xm)/1E3&&fabs(Y-Yn[i])<(Yx-Ym)/1E3)
      return(i);
  xn=Xn;
  yn=Yn;
  Xn=allocate(__LINE__,Nn+1,sizeof(double));
  Yn=allocate(__LINE__,Nn+1,sizeof(double));
  if(Nn)
    {
    memcpy(Xn,xn,Nn*sizeof(double));
    memcpy(Yn,yn,Nn*sizeof(double));
    free(xn);
    free(yn);
    }
  Xn[Nn]=X;
  Yn[Nn]=Y;
  Nn++;
  fprintf(stderr,"\r%i nodes, %i elements",Nn,Ne);
  return(Nn-1);
  }

void AddElement(double X1,double Y1,double X2,double Y2,double X3,double Y3)
  {
  double X,Xe[3],Y,Ye[3];
  int i,*ie,j,*je,k,*ke;
  if((X3-X1)*(Y2-Y1)-(X2-X1)*(Y3-Y1)<0)
    {
    X=X1;
    X1=X2;
    X2=X;
    Y=Y1;
    Y1=Y2;
    Y2=Y;
    }
  Xe[0]=X1;
  Ye[0]=Y1;
  Xe[1]=X2;
  Ye[1]=Y2;
  Xe[2]=X3;
  Ye[2]=Y3;
  i=AddNode(X1,Y1);
  j=AddNode(X2,Y2);
  k=AddNode(X3,Y3);
  ie=Ie;
  je=Je;
  ke=Ke;
  Ie=allocate(__LINE__,Ne+1,sizeof(int));
  Je=allocate(__LINE__,Ne+1,sizeof(int));
  Ke=allocate(__LINE__,Ne+1,sizeof(int));
  if(Ne)
    {
    memcpy(Ie,ie,Ne*sizeof(int));
    memcpy(Je,je,Ne*sizeof(int));
    memcpy(Ke,ke,Ne*sizeof(int));
    free(ie);
    free(je);
    free(ke);
    }
  Ie[Ne]=i;
  Je[Ne]=j;
  Ke[Ne]=k;
  Ne++;
  fprintf(stderr,"\r%i nodes, %i elements",Nn,Ne);
  }

double InsidePolygon(double*Xb,double*Yb,int Nb,double Xp,double Yp)
  {
  int above1,above2,i,right=0;
  double X1,X2,Y1,Y2;
  X2=Xb[Nb-1];
  Y2=Yb[Nb-1];
  above2=Y2>Yp?1:0;
  for(i=0;i<Nb;i++)
    {
    X1=X2;
    Y1=Y2;
    X2=Xb[i];
    Y2=Yb[i];
    above1=above2;
    above2=Y2>Yp?1:0;
    if(above1==above2)
      continue;
    if(X1>Xp&&X2>Xp)
      right++;
    else if(Y1<Y2)
      {
      if((Xp-X1)*(Y2-Y1)<(X2-X1)*(Yp-Y1))
        right++;
      }
    else if(Y1>Y2)
      {
      if((Xp-X1)*(Y2-Y1)>(X2-X1)*(Yp-Y1))
        right++;
      }
    }
  return(right&1);
  }

int LinesIntersect(double X1,double Y1,double X2,double Y2,double X3,double Y3,double X4,double Y4)
  {
  double A,B,C;
  A=(Y2-Y1)*(X4-X3)-(X2-X1)*(Y4-Y3);
  B=(X1-X3)*(Y4-Y3)-(Y1-Y3)*(X4-X3);
  C=(X1-X3)*(Y2-Y1)-(Y1-Y3)*(X2-X1);
  if(A==0)
    {
    if(B==0)
      return(1);
    return(0);
    }
  else
    {
    if(A>0)
      {
      if(B>=0&&B<=A&&C>=0&&C<=A)
        return(1);
      return(0);
      }
    else
      {
      if(B<=0&&B>=A&&C<=0&&C>=A)
        return(1);
      return(0);
      }
    }
  }

int PolygonIntersections(double*Xb,double*Yb,int Nb,double X1,double Y1,double X2,double Y2)
  {
  int i,j,n;
  for(n=i=0,j=Nb-1;i<Nb;j=i++)
    n+=LinesIntersect(X1,Y1,X2,Y2,Xb[i],Yb[i],Xb[j],Yb[j]);
  return(n);
  }

double LawOfCosines(double X1,double Y1,double X2,double Y2,double X3,double Y3)
  {
  double S1,S2,S3;
  S1=hypot(X2-X1,Y2-Y1);
  S2=hypot(X3-X2,Y3-Y2);
  S3=hypot(X1-X3,Y1-Y3);
  if(S1<DBL_EPSILON)
    return(-4);
  if(S2<DBL_EPSILON)
    return(-4);
  if(S3<DBL_EPSILON)
    return(-4);
  return((S1*S1+S2*S2-S3*S3)/2/S1/S2);
  }

void FindElements(double S) /* subdivide polygon into triangles */
  {
  int i,i1,i2,i3,j,j1,j2,j3,k,l,*lp,mp,n,np;
  double A,a,dX,dY,s,X1,X2,X3,X4,X5,X6,*xp,Y1,Y2,Y3,Y4,Y5,Y6,*yp;
  while(Np>3)
    {
    A=-2;
    for(i2=0;i2<Np;i2++)
      {
      i1=(Np+i2-1)%Np;
      i3=(i2+1)%Np;
      X1=Xp[i1];
      Y1=Yp[i1];
      X2=Xp[i2];
      Y2=Yp[i2];
      X3=Xp[i3];
      Y3=Yp[i3];
      if(Np>4)
        {
        if(PolygonIntersections(Xb,Yb,Nb,X1,Y1,X3,Y3)>4)
          continue;
        if(!InsidePolygon(Xb,Yb,Nb,(X1+X3)/2,(Y1+Y3)/2))
          continue;
        }
      a=LawOfCosines(X1,Y1,X2,Y2,X3,Y3);
      if(a>A)
        {
        A=a;
        k=i2;
        }
      }
    if(A<0) /* if no acute angles, bisect region */
      goto bisect;
    i=(Np+k-1)%Np;
    X1=Xp[i];
    Y1=Yp[i];
    i=(Np+k)%Np;
    X2=Xp[i];
    Y2=Yp[i];
    i=(Np+k+1)%Np;
    X3=Xp[i];
    Y3=Yp[i];
    s=hypot(X3-X1,Y3-Y1);
    n=(int)(s/S/1.5);
    if(n) /* add 2 elements and change 1 boundary point */
      {
      X4=(X1+X3)/2;
      Y4=(Y1+Y3)/2;
      AddElement(X1,Y1,X2,Y2,X4,Y4);
      AddElement(X2,Y2,X3,Y3,X4,Y4);
      i=(Np+k)%Np;
      Xp[i]=X4;
      Yp[i]=Y4;
      }
    else /* add 1 element and delete 1 boundary point */
      {
      AddElement(X1,Y1,X2,Y2,X3,Y3);
      Np--;
      Mp--;
      if(k<Mp)
        {
        memcpy(Xp+k,Xp+k+1,(Mp-k)*sizeof(double));
        memcpy(Yp+k,Yp+k+1,(Mp-k)*sizeof(double));
        }
      }
    continue;
    bisect: /* if no acute angles, bisect region */
    A=6;
    i=-1;
    j=-1;
    i1=Np-1;
    for(i2=0;i2<Np-1;i2++)
      {
      i3=i2+1;
      X1=Xp[i1];
      Y1=Yp[i1];
      X2=Xp[i2];
      Y2=Yp[i2];
      X3=Xp[i3];
      Y3=Yp[i3];
      for(j2=i2+2;j2<Np-1;j2++)
        {
        j1=j2-1;
        j3=j2+1;
        X4=Xp[j1];
        Y4=Yp[j1];
        X5=Xp[j2];
        Y5=Yp[j2];
        X6=Xp[j3];
        Y6=Yp[j3];
        if(Np>4)
          {
          if(PolygonIntersections(Xb,Yb,Nb,X3,Y3,X5,Y5)>4)
            continue;
          if(!InsidePolygon(Xb,Yb,Nb,(X2+X5)/2,(Y2+Y5)/2))
            continue;
          dX=hypot(X2-X1,X2-Y1);
          dX=fmin(dX,hypot(X3-X2,Y3-Y2));
          dX=fmin(dX,hypot(X5-X4,Y5-Y4));
          dX=fmin(dX,hypot(X6-X5,Y6-Y5));
          dY=dX;
          for(k=0;k<Np;k++)
            {
            if(k==i1)
              continue;
            if(k==i2)
              continue;
            if(k==i3)
              continue;
            if(k==j1)
              continue;
            if(k==j2)
              continue;
            if(k==j3)
              continue;
            dY=fmin(dY,DistanceToLine(X2,Y2,X5,Y5,Xp[k],Yp[k]));
            }
          if(dY<dX)
            continue;
          }
        a=0;
        a=fmax(a,LawOfCosines(X1,Y1,X2,Y2,X5,Y5));
        a=fmax(a,LawOfCosines(X3,Y3,X2,Y2,X5,Y5));
        a=fmax(a,LawOfCosines(X4,Y4,X5,Y5,X2,Y2));
        a=fmax(a,LawOfCosines(X6,Y6,X5,Y5,X2,Y2));
        if(a>=A)
          continue;
        A=a;
        i=i2;
        j=j2;
        }
      i1=i2;
      }
    if(i<0||j<0)
      Abort(__LINE__,"polygon bisection failed");
    X1=Xp[i]; /* create new polygon */
    Y1=Yp[i];
    X2=Xp[j];
    Y2=Yp[j];
    dX=X2-X1;
    dY=Y2-Y1;
    s=sqrt(dX*dX+dY*dY);
    n=(int)(s/S/1.5);
    l=j-i+1;
    np=l+n;
    mp=Mp;
    Mp+=np;
    lp=allocate(__LINE__,Kp+1,sizeof(int));
    if(Kp)
      {
      memcpy(lp+1,Lp,Kp*sizeof(int));
      free(Lp);
      }
    Lp=lp;
    Lp[0]=np;
    Kp++;
    xp=allocate(__LINE__,Mp,sizeof(double));
    memcpy(xp+np,Xp,mp*sizeof(double));
    memcpy(xp,Xp,Np*sizeof(double));
    memcpy(xp+Np,Xp+i,l*sizeof(double));
    free(Xp);
    Xp=xp;
    yp=allocate(__LINE__,Mp,sizeof(double));
    memcpy(yp+np,Yp,mp*sizeof(double));
    memcpy(yp,Yp,Np*sizeof(double));
    memcpy(yp+Np,Yp+i,l*sizeof(double));
    free(Yp);
    Yp=yp;
    dX/=n+1;
    dY/=n+1;
    for(k=0;k<n;k++)
      {
      X2-=dX;
      Y2-=dY;
      Xp[Np+l]=X2;
      Yp[Np+l]=Y2;
      l++;
      }
    l=i+1; /* bisect old polygon */
    for(k=0;k<n;k++)
      {
      X1+=dX;
      Y1+=dY;
      Xp[l]=X1;
      Yp[l]=Y1;
      l++;
      }
    if(l<j)
      {
      memcpy(Xp+l,Xp+j,(Mp-j)*sizeof(double));
      memcpy(Yp+l,Yp+j,(Mp-j)*sizeof(double));
      Np-=j-l;
      Mp-=j-l;
      }
    }
  if(Np<3)
    return;
  AddElement(Xp[0],Yp[0],Xp[1],Yp[1],Xp[2],Yp[2]);
  Np--;
  Mp--;
  memcpy(Xp,Xp+1,Mp*sizeof(double));
  memcpy(Yp,Yp+1,Mp*sizeof(double));
  }

int ShortenSides(double S) /* no side inter than S */
  {
  int i,n,Ns;
  double dX,dY,s,X1,X2,*xp,Y1,Y2,*yp;
  printf("  shortening sides\n");
  Ns=0;
  while(1)
    {
    i=Np-1;
    X2=Xp[i];
    Y2=Yp[i];
    for(i=0;i<Np;i++)
      {
      X1=X2;
      Y1=Y2;
      X2=Xp[i];
      Y2=Yp[i];
      dX=X2-X1;
      dY=Y2-Y1;
      s=sqrt(dX*dX+dY*dY);
      n=(int)(s/S/1.5);
      if(n)
        break;
      }
    if(!n) /* no, they are all shorter */
      break;
    xp=allocate(__LINE__,Np+n,sizeof(double));
    if(i)
      memcpy(xp,Xp,i*sizeof(double));
    memcpy(xp+i+n,Xp+i,(Np-i)*sizeof(double));
    free(Xp);
    Xp=xp;
    yp=allocate(__LINE__,Np+n,sizeof(double));
    if(i)
      memcpy(yp,Yp,i*sizeof(double));
    memcpy(yp+i+n,Yp+i,(Np-i)*sizeof(double));
    free(Yp);
    Yp=yp;
    dX/=n+1; /* insert closer points */
    dY/=n+1;
    while(n)
      {
      X1+=dX;
      Y1+=dY;
      Xp[i]=X1;
      Yp[i]=Y1;
      i++;
      n--;
      Np++;
      Ns++;
      }
    }
  return(Ns);
  }

int main(int argc,char**argv,char**envp)
  {
  int i,line;
  double R,S=0,X,X1,X2,Y,Y1,Y2;
  printf("ELEM3/V%s: Triangularization of a Closed Polygon\n",Version);
  if(argc<2)
    Abort(__LINE__,"no file specified");
  if(argc>3)
    Abort(__LINE__,"wrong number of parameters");
  printf("input file: %s\n",argv[1]);
  if((Finp=fopen(argv[1],"rt"))==NULL)
    Abort(__LINE__,"can't open input file");
  printf("  scanning input file\n");
  line=Nb=0;
  while(fgets(cbuf,sizeof(cbuf),Finp))
    {
    line++;
    if(sscanf(cbuf,"%lf%*[ ,]%lf%",&X,&Y)==2)
      Nb++;
    }
  printf("  %li lines read\n",line);
  printf("  %li boundary points found\n",Nb);
  if(Nb<3)
    Abort(__LINE__,"<3 boundary points");
  printf("  allocating memory\n");
  Xb=allocate(__LINE__,Nb,sizeof(double));
  Yb=allocate(__LINE__,Nb,sizeof(double));
  printf("  reading input file\n");
  rewind(Finp);
  Nb=0;
  while(fgets(cbuf,sizeof(cbuf),Finp))
    if(sscanf(cbuf,"%lf%*[ ,]%lf%",Xb+Nb,Yb+Nb)==2)
      Nb++;
  fclose(Finp);
  printf("analyzing boundary polygon\n");
  Xm=Ym=DBL_MAX;
  Xx=Yx=-DBL_MAX;
  X2=Xb[Nb-1];
  Y2=Yb[Nb-1];
  for(i=0;i<Nb;i++)
    {
    X1=X2;
    Y1=Y2;
    X2=Xb[i];
    Y2=Yb[i];
    Xm=fmin(Xm,X2);
    Xx=fmax(Xx,X2);
    Ym=fmin(Ym,Y2);
    Yx=fmax(Yx,Y2);
    R=hypot(X2-X1,Y2-Y1);
    if(R<DBL_EPSILON)
      Abort(__LINE__,"polygon has zero-length side");
    S=fmax(S,R);
    }
  if(Xm>=Xx||Ym>=Yx)
    Abort(__LINE__,"polygon is degenerate");
  if(argc>2)
    if(sscanf(argv[2],"%lf",&R)!=1)
      Abort(__LINE__,"scan error on parameter 2 (element size)");
    else if(R>DBL_EPSILON)
      S=R;
  if(S<fmin((Xx-Xm)/1E3,(Yx-Ym)/1E3))
    Abort(__LINE__,"element size too small");
  if(PolygonArea(Xb,Yb,Nb)<(Xx-Xm)*(Yx-Ym)/1E2)
    Abort(__LINE__,"polygon is too acute");
  printf("  initializing points to boundary\n");
  Mp=Np=Nb;
  Xp=allocate(__LINE__,Mp,sizeof(double));
  Yp=allocate(__LINE__,Mp,sizeof(double));
  memcpy(Xp,Xb,Np*sizeof(double));
  memcpy(Yp,Yb,Np*sizeof(double));
  ShortenSides(S);
  Mp=Np;
  Nn=0;
  for(i=0;i<Np;i++)
    AddNode(Xp[i],Yp[i]);
  Ne=0;
  while(1)
    {
    FindElements(S);
    if(Np==Mp)
      break;
    Mp-=2;
    if(Kp==0)
      break;
    Kp--;
    Np=Lp[0];
    if(Kp)
      memcpy(Lp,Lp+1,Kp*sizeof(int));
    memcpy(Xp,Xp+2,Mp*sizeof(double));
    memcpy(Yp,Yp+2,Mp*sizeof(double));
    }
  fprintf(stderr,"\r%li nodes, %li elements\n",Nn,Ne);
  strcpy(cbuf,argv[1]);
  if(strchr(cbuf,'.'))
    strcpy(strchr(cbuf,'.'),".2dv");
  else
    strcat(cbuf,".2dv");
  fprintf(stderr,"output file: %s\n",cbuf);
  if((Fout=fopen(cbuf,"wt"))==NULL)
    Abort(__LINE__,"can't create output file");
  fprintf(Fout,"%li nodes (below X,Y)\n",Nn);
  for(i=0;i<Nn;i++)
    fprintf(Fout,"%lG %lG\n",Xn[i],Yn[i]);
  fprintf(Fout,"%li elements (below n1,n2,n3)\n",Ne);
  for(i=0;i<Ne;i++)
    fprintf(Fout,"%li %li %li\n",Ie[i]+1,Je[i]+1,Ke[i]+1);
  fclose(Fout);
  exit(0);
  }
