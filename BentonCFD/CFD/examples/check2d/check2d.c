/***********************************************************************
          CCCC  HH  HH  EEEEEEE   CCCC  KKK  KK  2222   DDDDD
         CC  CC HH  HH   EE   E  CC  CC  KK  KK 22  22   DD DD
        CC      HH  HH   EE E   CC       KK KK      22   DD  DD
        CC      HHHHHH   EEEE   CC       KKKK     222    DD  DD
        CC      HH  HH   EE E   CC       KK KK   22      DD  DD
         CC  CC HH  HH   EE   E  CC  CC  KK  KK 22  22   DD DD
          CCCC  HH  HH  EEEEEEE   CCCC  KKK  KK 222222  DDDDD
        -------------------------------------------------------
    The purpose of this program is to check a file containing nodes
    and triangular or quadlateral elements for inconsistencies.
       developed by Dudley J. Benton, Knoxville, Tennessee
***********************************************************************/

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

double PolygonArea(double*Xp,double*Yp,int Np)
  {
  int i,j;
  double A;
  for(i=Np-1,A=j=0;j<Np;i=j++)
    A+=(Xp[j]-Xp[i])*(Yp[i]+Yp[j]);
  return(A/2.);
  }

int InsidePolygon(double*Xb,double*Yb,int Nb,double Xp,double Yp)
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

char bufr[256];
FILE*file;

double Xmin=DBL_MAX,Xmax=-DBL_MAX;
double Ymin=DBL_MAX,Ymax=-DBL_MAX;

int*Ie,*In,Nn,Ne;
double Scale,Xoff,*Xn,Yoff,*Yn;

int main(int argc,char**argv,char**envp)
  {
  int clockwise,e,i,j,k,l,lines,m,n,n1,n2,n3,n4,node[8];
  double A,Asmall,X,Xi[4],Xp[4],Xsmall,Xx,Y,Yi[4],Yp[4],Ysmall,Yy;

  printf("CHECK2D: Check 2D Elements\n");

  if(argc<2)
    {
    printf("no file specified... try something like\n"
           "CHECK2D star.2dv\n");
    return(__LINE__);
    }
  if(argc>2)
    {
    printf("too many parameters... try something like\n"
           "CHECK2D tree.2dv\n");
    return(__LINE__);
    }

  printf("input file=%s\n",argv[1]);
  if((file=fopen(argv[1],"rt"))==NULL)
    {
    printf("can't open file\n");
    return(__LINE__);
    }

  printf("reading nodes\n");
  lines=0;
  while(1)
    {
    if(!fgets(bufr,sizeof(bufr),file))
      {
      printf("can't read file\n");
      return(__LINE__);
      }
    lines++;
    if(sscanf(bufr,"%i",&Nn))
      break;
    }
  printf("  expecting %i nodes\n",Nn);
  if(Nn<3)
    {
    printf("less than 3 nodes\n");
    return(__LINE__);
    }

  printf("  allocating memory for nodes\n");

  if((Xn=calloc(Nn,sizeof(double)))==NULL)
    {
    printf("can't allocate memory\n");
    return(__LINE__);
    }
  if((Yn=calloc(Nn,sizeof(double)))==NULL)
    {
    printf("can't allocate memory\n");
    return(__LINE__);
    }
  if((In=calloc(Nn,sizeof(int)))==NULL)
    {
    printf("can't allocate memory\n");
    return(__LINE__);
    }

  n=0;
  while(fgets(bufr,sizeof(bufr),file))
    {
    lines++;
    i=sscanf(bufr,"%lf%*[ ,]%lf",&X,&Y);
    if(i==0)
      continue;
    if(i<2)
      {
      printf("\n");
      printf(bufr);
      printf("\nX but not Y found on line %i\n",lines);
      return(__LINE__);
      }
    Xmin=min(Xmin,X);
    Xmax=max(Xmax,X);
    Ymin=min(Ymin,Y);
    Ymax=max(Ymax,Y);
    Xn[n]=X;
    Yn[n]=Y;
    if(++n==Nn)
      break;
    }
  printf("  %i nodes found\n",n);
  if(Nn!=n)
    {
    printf("nodes found does not match nodes expected\n");
    return(__LINE__);
    }

  printf("  %lGóXó%lG\n",Xmin,Xmax);
  if(Xmin>=Xmax)
    {
    printf("data is degenerate in X\n");
    return(__LINE__);
    }
  printf("  %lGóYó%lG\n",Ymin,Ymax);
  if(Ymin>=Ymax)
    {
    printf("data is degenerate in Y\n");
    return(__LINE__);
    }

  printf("reading elements\n");
  while(1)
    {
    if(!fgets(bufr,sizeof(bufr),file))
      {
      printf("can't read file\n");
      return(__LINE__);
      }
    lines++;
    if(sscanf(bufr,"%i",&Ne))
      break;
    }
  printf("  expecting %i elements\n",Ne);
  if(Ne<1)
    {
    printf("no elements\n");
    return(__LINE__);
    }

  printf("  allocating memory for elements\n");

  if((Ie=calloc(4*Ne,sizeof(int)))==NULL)
    {
    printf("can't allocate memory\n");
    return(__LINE__);
    }

  e=0;
  while(fgets(bufr,sizeof(bufr),file))
    {
    lines++;
    node[3]=-1;
    i=sscanf(bufr,"%i%*[ ,]%i%*[ ,]%i%*[ ,]%i",node+0,node+1,node+2,node+3);
    if(i==0)
      continue;
    if(i<3)
      {
      printf("\n%s",bufr);
      printf("at least 1, but less than 3 nodes found on line %i\n",lines);
      return(__LINE__);
      }
    for(j=0;j<i;j++)
      {
      if(node[j]<1||node[j]>Nn)
        {
        printf("\n%s",bufr);
        printf("no such node %i found on line %i\n",node[j],lines);
        return(__LINE__);
        }
      }
    for(j=0;j<4;j++)
      Ie[4*e+j]=node[j]-1;
    if(++e==Ne)
      break;
    }
  printf("  %i elements found\n",e);
  if(Ne!=e)
    {
    printf("elements found does not match elements expected\n");
    return(__LINE__);
    }

  fclose(file);

  printf("checking for errors\n");
  printf("  unused nodes...");
  for(n=0;n<Nn;n++)
    In[n]=0;
  for(e=0;e<Ne;e++)
    {
    for(i=0;i<4;i++)
      {
      n=Ie[4*e+i];
      if(n>=0)
        In[n]++;
      }
    }
  for(n=0;n<Nn;n++)
    {
    if(In>0)
      continue;
    printf("\nnode %i does not appear in any element\n",n+1);
    return(__LINE__);
    }
  printf(" none\n");

  printf("  coincident boundaries...");
  for(e=0;e<Ne-1;e++)
    {
    for(j=0;j<4;j++)
      node[j]=Ie[4*e+j];
    if(node[3]<0)
      n2=node[2];
    else
      n2=node[3];
    for(i=0;i<4;i++)
      {
      m=0;
      n1=n2;
      n2=node[i];
      if(n2<0)
        break;
      for(l=e+1;l<Ne;l++)
        {
        for(j=0;j<4;j++)
          node[4+j]=Ie[4*l+j];
        if(node[7]<0)
          n4=node[6];
        else
          n4=node[7];
        for(k=0;k<4;k++)
          {
          n3=n4;
          n4=node[k+4];
          if(n4<0)
            break;
          if(!((n1==n3&&n2==n4)||(n1==n4&&n2==n3)))
            continue;
          if(m==0)
            {
            m=l+1;
            continue;
            }
          printf("\nnodes %i and %i appear together in at least 3 elements %i %i %i\n",n1,n2,e+1,m,l+1);
          return(__LINE__);
          }
        }
      }
    }
  printf(" none\n");

  printf("  overlapping boundaries...");
  for(e=0;e<Ne-1;e++)
    {
    for(j=0;j<4;j++)
      {
      n=Ie[4*e+j];
      node[j]=n;
      Xi[j]=(Xn[n]-Xmin)/(Xmax-Xmin);
      Yi[j]=(Yn[n]-Ymin)/(Ymax-Ymin);
      }
    if(n<0)
      m=3;
    else
      m=4;
    for(n=0;n<Nn;n++)
      {
      if(n==node[0])
        continue;
      if(n==node[1])
        continue;
      if(n==node[2])
        continue;
      if(n==node[3])
        continue;
      Xx=(Xn[n]-Xmin)/(Xmax-Xmin);
      Yy=(Yn[n]-Ymin)/(Ymax-Ymin);
      if(!InsidePolygon(Xi,Yi,m,Xx,Yy))
        continue;
      printf("\nnode %i lies inside element %i\n",n+1,e+1);
      return(__LINE__);
      }
    }
  printf(" none\n");

  printf("  coincident nodes...");
  Xsmall=(Xmax-Xmin)*sqrt(DBL_EPSILON);
  Ysmall=(Ymax-Ymin)*sqrt(DBL_EPSILON);
  for(n=0;n<Nn-1;n++)
    {
    for(i=n+1;i<Nn;i++)
      {
      if(fabs(Xn[n]-Xn[i])>Xsmall)
        continue;
      if(fabs(Yn[n]-Yn[i])>Ysmall)
        continue;
      printf("\nnodes %i and %i are coincident\n",n+1,i+1);
      return(__LINE__);
      }
    }
  printf(" none\n");

  printf("  elements having duplicate nodes...");
  for(e=0;e<Ne;e++)
    {
    for(i=0;i<4;i++)
      node[i]=Ie[4*e+i];
    for(i=0;i<3;i++)
      {
      if(node[i]<1)
        continue;
      for(j=i+1;j<4;j++)
        {
        if(node[j]<1)
          continue;
        if(node[i]!=node[j])
          continue;
        printf("\nelement %i contains duplicate node %i\n",e+1,node[i]+1);
        return(__LINE__);
        }
      }
    }
  printf(" none\n");

  printf("  degenerate elements...");
  Asmall=(Xmax-Xmin)*(Ymax-Ymin)*sqrt(DBL_EPSILON);
  clockwise=0;
  for(e=0;e<Ne;e++)
    {
    for(i=0;i<4;i++)
      {
      n=Ie[4*e+i];
      if(n<0)
        break;
      Xp[i]=Xn[n];
      Yp[i]=Yn[n];
      }
    if(n<0)
      A=PolygonArea(Xp,Yp,3);
    else
      A=PolygonArea(Xp,Yp,4);
    if(fabs(A)<Asmall)
      {
      printf("\nelement %i is degenerate\n",e+1);
      return(__LINE__);
      }
    if(clockwise==0)
      {
      if(A>0.)
        clockwise=1;
      else
        clockwise=-1;
      }
    else if(clockwise>0)
      {
      if(A<0.)
        {
        printf("\nelement %i is not clockwise\n",e+1);
        return(__LINE__);
        }
      }
    else
      {
      if(A>0.)
        {
        printf("\nelement %i is not counter-clockwise\n",e+1);
        return(__LINE__);
        }
      }
    }
  printf(" none\n");

  if(clockwise>0)
    printf("elements are clockwise\n");
  else
    printf("elements are counter-clockwise\n");

  return(0);
  }
