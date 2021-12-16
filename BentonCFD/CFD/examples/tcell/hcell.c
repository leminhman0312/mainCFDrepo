#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <malloc.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define ADR(t,x) ( t=(x), &t )

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

int idamax(int*n,double dx[],int*incx)
  {
  int i,idamax_v,ix;
  double dmax;
  idamax_v=0;
  if(*n>=1)
    {
    idamax_v=1;
    if(*n!=1)
      {
      if(*incx==1)
        {
        dmax=fabs(dx[0]);
        for(i=2;i<=*n;i++)
          {
          if(dmax<fabs(dx[i-1]))
            {
            idamax_v=i;
            dmax=fabs(dx[i-1]);
            }
          }
        }
      else
        {
        ix=1;
        dmax=fabs(dx[0]);
        ix+=*incx;
        for(i=2;i<=*n;i++)
          {
          if(dmax<fabs(dx[ix-1]))
            {
            idamax_v=i;
            dmax=fabs(dx[ix-1]);
            ix+=*incx;
            }
          }
        }
      }
    }
  return(idamax_v);
  }

void dscal(int*n,double*da,double dx[],int*incx)
  {
  int i,m,nincx;
  if(*n>0)
    {
    if(*incx==1)
      {
      m=*n%5;
      for(i=1;i<=m;i++)
        {
        dx[i-1]*=*da;
        }
      for(i=m+1;i<=*n;i+=5)
        {
        dx[i-1]*=*da;
        dx[i]*=*da;
        dx[i+1]*=*da;
        dx[i+2]*=*da;
        dx[i+3]*=*da;
        }
      }
    else
      {
      nincx=*n ** incx;
      for(i=1;i<=nincx;i+=*incx)
        dx[i-1]*=*da;
      }
    }
  }

void daxpy(int*n,double*da,double dx[],int*incx,double dy[],int*incy)
  {
  int i,ix,iy,m;
  if(*n>0)
    {
    if(*da!=0.)
      {
      if(*incx==1&&*incy==1)
        {
        m=*n%4;
        for(i=1;i<=m;i++)
          dy[i-1]+=*da*dx[i-1];
        for(i=m+1;i<=*n;i+=4)
          {
          dy[i-1]+=*da*dx[i-1];
          dy[i]+=*da*dx[i];
          dy[i+1]+=*da*dx[i+1];
          dy[i+2]+=*da*dx[i+2];
          }
        }
      else
        {
        ix=1;
        iy=1;
        if(*incx<0)
          ix=(-*n+1) ** incx+1;
        if(*incy<0)
          iy=(-*n+1) ** incy+1;
        for(i=1;i<=*n;i++)
          {
          dy[iy-1]+=*da*dx[ix-1];
          ix+=*incx;
          iy+=*incy;
          }
        }
      }
    }
  }

#define MAXEL (2*(NX-1)*(NY-1))
#define MAXND (MX*MY)
#define MAXUN (2*MX*MY+NX*NY)
#define MINUN (27*NY)
#define MX (2*NX-1)
#define MY (2*NY-1)
#define N_TIME 10
#define NNODES 6
#define NQUAD 3
#define NUK 3
#define NX1 10
#define NX2 10
#define NX3 10
#define NX (NX1+NX2+NX3+1)
#define NY1 10
#define NY2 10
#define NY3 10
#define NY (NY1+NY2+NY3+1)

int indx[NUK][MAXND];
int node[NNODES][MAXEL];
double a[MAXUN][MINUN];
double xm[NQUAD][MAXEL];
double ym[NQUAD][MAXEL];

#undef MAXEL
#undef MAXND
#undef MAXUN
#undef MINUN
#undef MX
#undef MY
#undef N_TIME
#undef NNODES
#undef NQUAD
#undef NUK
#undef NX
#undef NX1
#undef NX2
#undef NX3
#undef NY
#undef NY1
#undef NY2
#undef NY3

void hcell_element_count(int ireg_density_x[],int ireg_density_y[],int*element_num)
  {
  *element_num=2*ireg_density_x[0]*ireg_density_y[0]+2*ireg_density_x[0] *
  ireg_density_y[2]+2*ireg_density_x[1]*ireg_density_y[0]+2 *
  ireg_density_x[1]*ireg_density_y[1]+2*ireg_density_x[1]*ireg_density_y[2] +
  2*ireg_density_x[2]*ireg_density_y[0]+2*ireg_density_x[2]*ireg_density_y[2];
  printf("number of elements=%i\n",*element_num);
  }

void hcell_element_node(int ireg_density_x[],int ireg_density_y[],int*maxel,int*element_num,int*nnodes)
  {
  int col,col2,element,inc1,inc2,node_sw,row,row2;
  element=0;
  for(col=1;col<=3;col++)
    {
    for(col2=1;col2<=ireg_density_x[col-1];col2++)
      {
      for(row=1;row<=3;row++)
        {
        if(row!=2||col==2)
          {
          if(col==1)
            {
            if(col2<ireg_density_x[0])
              {
              if(row!=1)
                node_sw+=1;
              else if(col2==1)
                node_sw=1;
              else
                node_sw+=inc1+1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              }
            else if(row==1)
              {
              node_sw+=inc1+1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              }
            else if(row==3)
              {
              node_sw+=1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[1]-1)+(2*ireg_density_y[2] +
               1);
              }
            }
          else if(col==2)
            {
            if(row==1)
              node_sw+=inc1+1;
            inc1=(2*ireg_density_y[0]+1)+(2*ireg_density_y[1] -
             1)+(2*ireg_density_y[2]+1);
            inc2=(2*ireg_density_y[0]+1)+(2*ireg_density_y[1] -
             1)+(2*ireg_density_y[2]+1);
            }
          else if(col==3)
            {
            if(1<col2)
              {
              if(row==1)
                node_sw+=inc1+1;
              else
                node_sw+=1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              }
            else if(row==1)
              {
              node_sw+=inc1+1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[1]-1)+(2*ireg_density_y[2] +
               1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              }
            else if(row==3)
              {
              node_sw+=(2*ireg_density_y[1]-1) +
               1;
              inc1=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              inc2=(2*ireg_density_y[0]+1)+(2 *
               ireg_density_y[2]+1);
              }
            }
          for(row2=1;row2<=ireg_density_y[row-1];row2++)
            {
            element+=1;
            node[0][element-1]=node_sw;
            node[1][element-1]=node_sw+inc1+inc2;
            node[2][element-1]=node_sw+2;
            node[3][element-1]=node_sw+inc1;
            node[4][element-1]=node_sw+inc1+1;
            node[5][element-1]=node_sw+1;
            element+=1;
            node[0][element-1]=node_sw+inc1+inc2+2;
            node[1][element-1]=node_sw+2;
            node[2][element-1]=node_sw+inc1+inc2;
            node[3][element-1]=node_sw+inc1+2;
            node[4][element-1]=node_sw+inc1+1;
            node[5][element-1]=node_sw+inc1+inc2+1;
            node_sw+=2;
            }
          }
        }
      }
    }
  }

void hcell_node_count(int ireg_density_x[],int ireg_density_y[],int*node_num)
  {
  *node_num=(2*ireg_density_x[0]+1)*(2*ireg_density_y[0]+1) +
  (2*ireg_density_x[0]+1)*(2*ireg_density_y[2]+1)+(2*ireg_density_x[1] -
   1)*(2*ireg_density_y[0]+1)+(2*ireg_density_x[1]+1)*(2*ireg_density_y[1] -
   1)+(2*ireg_density_x[1]-1)*(2*ireg_density_y[2]+1)+(2 *
   ireg_density_x[2]+1)*(2*ireg_density_y[0]+1)+(2*ireg_density_x[2] +
   1)*(2*ireg_density_y[2]+1);
  printf("number of nodes=%i\n",*node_num);
  }

void hcell_node_xy(int ireg_density_x[],int ireg_density_y[],int*node_num,double region_x[],double region_y[],double x_node[],double y_node[])
  {
  int col,i,j,node,row;
  node=0;
  j=0;
  for(col=1;col<=(2*ireg_density_x[0]);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[0]+1 -
        col)*region_x[0]+(double) (-1+col)*region_x[1]) /
       (double) (2*ireg_density_x[0]);
      y_node[node-1]=((double) (2*ireg_density_y[0]+1 -
        row)*region_y[0]+(double) (-1+row)*region_y[1]) /
       (double) (2*ireg_density_y[0]);
      }
    i+=2*ireg_density_y[1]-1;
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[0]+1 -
        col)*region_x[0]+(double) (-1+col)*region_x[1]) /
       (double) (2*ireg_density_x[0]);
      y_node[node-1]=((double) (2*ireg_density_y[2]+1 -
        row)*region_y[2]+(double) (-1+row)*region_y[3]) /
       (double) (2*ireg_density_y[2]);
      }
    }
  for(col=1;col<=(2*ireg_density_x[1]+1);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[1]+1 -
        col)*region_x[1]+(double) (-1+col)*region_x[2]) /
       (double) (2*ireg_density_x[1]);
      y_node[node-1]=((double) (2*ireg_density_y[0]+1 -
        row)*region_y[0]+(double) (-1+row)*region_y[1]) /
       (double) (2*ireg_density_y[0]);
      }
    for(row=2;row<=(2*ireg_density_y[1]);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[1]+1 -
        col)*region_x[1]+(double) (-1+col)*region_x[2]) /
       (double) (2*ireg_density_x[1]);
      y_node[node-1]=((double) (2*ireg_density_y[1]+1 -
        row)*region_y[1]+(double) (-1+row)*region_y[2]) /
       (double) (2*ireg_density_y[1]);
      }
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[1]+1 -
        col)*region_x[1]+(double) (-1+col)*region_x[2]) /
       (double) (2*ireg_density_x[1]);
      y_node[node-1]=((double) (2*ireg_density_y[2]+1 -
        row)*region_y[2]+(double) (-1+row)*region_y[3]) /
       (double) (2*ireg_density_y[2]);
      }
    }
  for(col=2;col<=(2*ireg_density_x[2]+1);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[2]+1 -
        col)*region_x[2]+(double) (-1+col)*region_x[3]) /
       (double) (2*ireg_density_x[2]);
      y_node[node-1]=((double) (2*ireg_density_y[0]+1 -
        row)*region_y[0]+(double) (-1+row)*region_y[1]) /
       (double) (2*ireg_density_y[0]);
      }
    i+=2*ireg_density_y[1]-1;
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      x_node[node-1]=((double) (2*ireg_density_x[2]+1 -
        col)*region_x[2]+(double) (-1+col)*region_x[3]) /
       (double) (2*ireg_density_x[2]);
      y_node[node-1]=((double) (2*ireg_density_y[2]+1 -
        row)*region_y[2]+(double) (-1+row)*region_y[3]) /
       (double) (2*ireg_density_y[2]);
      }
    }
  }

void hcell_dof_count(int ireg_density_x[],int ireg_density_y[],int*dof_num)
  {
  int node_num_p,node_num_u;
  node_num_u=(2*ireg_density_x[0]+1)*(2*ireg_density_y[0] +
   1)+(2*ireg_density_x[0]+1)*(2*ireg_density_y[2]+1)+(2 *
   ireg_density_x[1]-1)*(2*ireg_density_y[0]+1)+(2*ireg_density_x[1] +
   1)*(2*ireg_density_y[1]-1)+(2*ireg_density_x[1]-1)*(2*ireg_density_y[2] +
   1)+(2*ireg_density_x[2]+1)*(2*ireg_density_y[0]+1)+(2 *
   ireg_density_x[2]+1)*(2*ireg_density_y[2]+1);
  node_num_p=(ireg_density_x[0]+1)*(ireg_density_y[0]+1) +
   (ireg_density_x[0]+1)*(ireg_density_y[2]+1)+(ireg_density_x[1] -
   1)*(ireg_density_y[0]+1)+(ireg_density_x[1]+1)*(ireg_density_y[1] -
   1)+(ireg_density_x[1]-1)*(ireg_density_y[2]+1)+(ireg_density_x[2] +
   1)*(ireg_density_y[0]+1)+(ireg_density_x[2]+1)*(ireg_density_y[2] +
   1);
  *dof_num=2*node_num_u+node_num_p;
  printf("number of degrees of freedom=%i\n",*dof_num);
  }

void hcell_dof_set(int ireg_density_x[],int ireg_density_y[],int*maxnd,int*node_num,int*dof)
  {
  int col,i,j,node,row;
  node=0;
  *dof=0;
  j=0;
  for(col=1;col<=(2*ireg_density_x[0]);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      if(col==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(row>=(2*2*ireg_density_y[0])/5+1&&row <=
         (3*2*ireg_density_y[0])/5+1)
          indx[0][node-1]=-3;
        }
      else if(row==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==2*ireg_density_y[0]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[0])/45+1&&col <=
         (11*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-3;
        if(col>=(22*2*ireg_density_x[0])/45+1&&col <=
         (23*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-3;
        if(col>=(34*2*ireg_density_x[0])/45+1&&col <=
         (35*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-3;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    i+=2*ireg_density_y[1]-1;
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      if(col==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(row>=(2*2*ireg_density_y[2])/5+1&&row <=
         (3*2*ireg_density_y[2])/5+1)
          indx[0][node-1]=-1;
        }
      else if(row==2*ireg_density_y[2]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[0])/45+1&&col <=
         (11*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-5;
        if(col>=(22*2*ireg_density_x[0])/45+1&&col <=
         (23*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-5;
        if(col>=(34*2*ireg_density_x[0])/45+1&&col <=
         (35*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-5;
        }
      else if(row==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[0])/45+1&&col <=
         (11*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-1;
        if(col>=(22*2*ireg_density_x[0])/45+1&&col <=
         (23*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-1;
        if(col>=(34*2*ireg_density_x[0])/45+1&&col <=
         (35*2*ireg_density_x[0])/45+1)
          indx[1][node-1]=-1;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    }
  for(col=1;col<=(2*ireg_density_x[1]+1);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      if(row==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==2*ireg_density_y[0]+1&&col==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==2*ireg_density_y[0]+1&&col==2*ireg_density_x[1] +
       1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    for(row=2;row<=(2*ireg_density_y[1]);row++)
      {
      i+=1;
      node+=1;
      if(col==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(col==2*ireg_density_x[1]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      if(row==1&&col==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==1&&col==2*ireg_density_x[1]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==2*ireg_density_y[2]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(5*2*ireg_density_x[1])/15&&col<=(10 *
          2*ireg_density_x[1])/15)
          {
          *dof+=1;
          indx[1][node-1]=*dof;
          }
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    }
  for(col=2;col<=(2*ireg_density_x[2]+1);col++)
    {
    j+=1;
    i=0;
    for(row=1;row<=(2*ireg_density_y[0]+1);row++)
      {
      i+=1;
      node+=1;
      if(col==2*ireg_density_x[2]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(row>=(2*2*ireg_density_y[0])/5+1&&row <=
         (3*2*ireg_density_y[0])/5+1)
          indx[0][node-1]=-4;
        }
      else if(row==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        }
      else if(row==2*ireg_density_y[0]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[2])/45+1&&col <=
         (11*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-4;
        if(col>=(22*2*ireg_density_x[2])/45+1&&col <=
         (23*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-4;
        if(col>=(34*2*ireg_density_x[2])/45+1&&col <=
         (35*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-4;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    i+=2*ireg_density_y[1]-1;
    for(row=1;row<=(2*ireg_density_y[2]+1);row++)
      {
      i+=1;
      node+=1;
      if(col==2*ireg_density_x[2]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(row>=(2*2*ireg_density_y[2])/5+1&&row <=
         (3*2*ireg_density_y[2])/5+1)
          indx[0][node-1]=-2;
        }
      else if(row==2*ireg_density_y[2]+1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[2])/45+1&&col <=
         (11*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-6;
        if(col>=(22*2*ireg_density_x[2])/45+1&&col <=
         (23*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-6;
        if(col>=(34*2*ireg_density_x[2])/45+1&&col <=
         (35*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-6;
        }
      else if(row==1)
        {
        indx[0][node-1]=0;
        indx[1][node-1]=0;
        if(col>=(10*2*ireg_density_x[2])/45+1&&col <=
         (11*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-2;
        if(col>=(22*2*ireg_density_x[2])/45+1&&col <=
         (23*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-2;
        if(col>=(34*2*ireg_density_x[2])/45+1&&col <=
         (35*2*ireg_density_x[2])/45+1)
          indx[1][node-1]=-2;
        }
      else
        {
        *dof+=1;
        indx[0][node-1]=*dof;
        *dof+=1;
        indx[1][node-1]=*dof;
        }
      if((j%2)==1&&(i%2)==1)
        {
        *dof+=1;
        indx[2][node-1]=*dof;
        }
      else
        indx[2][node-1]=0;
      }
    }
  }

void quad_a_set(int*maxel,int*nelemn,int*nquad,double area[])
  {
  int element;
  for(element=1;element<=*nelemn;element++)
    {
    xm[0][element-1]=0.5;
    xm[1][element-1]=0.5;
    xm[2][element-1]=0.0;
    }
  for(element=1;element<=*nelemn;element++)
    {
    ym[0][element-1]=0.0;
    ym[1][element-1]=0.5;
    ym[2][element-1]=0.5;
    }
  for(element=1;element<=*nelemn;element++)
    area[element-1]=0.5;
  }

void element_node_bandwidth(int*maxnd,int*maxel,int*nnodes,int*nelemn,int*neqn,int*np,int*nlband,int*nuband,int*nband)
  {
  int element,ieqn,l1,l2,n1,n2,p1,p2,u1,u2,v1,v2;
  int*dof_max,*dof_min;
  dof_max=calloc(*neqn,sizeof(int));
  dof_min=calloc(*neqn,sizeof(int));
  for(ieqn=1;ieqn<=*neqn;ieqn++)
    {
    dof_min[ieqn-1]=ieqn;
    dof_max[ieqn-1]=ieqn;
    }
  for(element=1;element<=*nelemn;element++)
    {
    for(l1=1;l1<=*nnodes;l1++)
      {
      n1=node[l1-1][element-1];
      u1=indx[0][n1-1];
      v1=indx[1][n1-1];
      p1=indx[2][n1-1];
      for(l2=1;l2<=*nnodes;l2++)
        {
        n2=node[l2-1][element-1];
        u2=indx[0][n2-1];
        v2=indx[1][n2-1];
        p2=indx[2][n2-1];
        if(1<=u1&&1<=u2)
          {
          dof_min[u1-1]=min(dof_min[u1-1],u2);
          dof_max[u1-1]=max(dof_max[u1-1],u2);
          }
        if(1<=u1&&1<=v2)
          {
          dof_min[u1-1]=min(dof_min[u1-1],v2);
          dof_max[u1-1]=max(dof_max[u1-1],v2);
          }
        if(1<=u1&&1<=p2)
          {
          dof_min[u1-1]=min(dof_min[u1-1],p2);
          dof_max[u1-1]=max(dof_max[u1-1],p2);
          }
        if(1<=v1&&1<=u2)
          {
          dof_min[v1-1]=min(dof_min[v1-1],u2);
          dof_max[v1-1]=max(dof_max[v1-1],u2);
          }
        if(1<=v1&&1<=v2)
          {
          dof_min[v1-1]=min(dof_min[v1-1],v2);
          dof_max[v1-1]=max(dof_max[v1-1],v2);
          }
        if(1<=v1&&1<=p2)
          {
          dof_min[v1-1]=min(dof_min[v1-1],p2);
          dof_max[v1-1]=max(dof_max[v1-1],p2);
          }
        if(1<=p1&&1<=u2)
          {
          dof_min[p1-1]=min(dof_min[p1-1],u2);
          dof_max[p1-1]=max(dof_max[p1-1],u2);
          }
        if(1<=p1&&1<=v2)
          {
          dof_min[p1-1]=min(dof_min[p1-1],v2);
          dof_max[p1-1]=max(dof_max[p1-1],v2);
          }
        }
      }
    }
  *nlband=0;
  *nuband=0;
  for(ieqn=1;ieqn<=*neqn;ieqn++)
    {
    *nlband=max(*nlband,ieqn-dof_min[ieqn-1]);
    *nuband=max(*nuband,dof_max[ieqn-1]-ieqn);
    }
  *nband=*nlband+*nuband+1;
  printf("lower half bandwidth=%i\n",*nlband);
  printf("upper half bandwidth=%i\n",*nuband);
  printf("total bandwidth=%i\n",*nband);
  free(dof_max);
  free(dof_min);
  }

double ubdry(int*iuk,int*ip,int*iun,double xc[],double yc[])
  {
  double ubdry_v,x,y;
  ubdry_v=0e0;
  x=xc[*ip-1];
  y=yc[*ip-1];
  if(*iuk==1)
    {
    if(y<=3e0)
      {
      if(x<5e0)
        ubdry_v=1e0;
      else
        ubdry_v=-1e0;
      }
    else if(y>=6e0)
      {
      if(x<5e0)
        ubdry_v=-1e0;
      else
        ubdry_v=1e0;
      }
    }
  return(ubdry_v);
  }

double refbsp(double*x,double*y,int*iq)
  {
  double refbsp_v;
  if(*iq==1)
    refbsp_v=1.-*x-*y;
  else if(*iq==2)
    refbsp_v=*x;
  else if(*iq==3)
    refbsp_v=*y;
  else
    {
    puts("REFBSP-FATAL ERROR");
    exit(0);
    }
  return(refbsp_v);
  }

void refqbf(double*x,double*y,int*in,double*bb,double*bx,double*by)
  {
  if(*in==1)
    {
    *bb=(1.-*x-*y)*(1.-2. ** x-2. ** y);
    *bx=-3.+4. ** x+4. ** y;
    *by=-3.+4. ** x+4. ** y;
    }
  else if(*in==2)
    {
    *bb=*x*(2. ** x-1.);
    *bx=4. ** x-1.;
    *by=0.;
    }
  else if(*in==3)
    {
    *bb=*y*(2. ** y-1.);
    *bx=0.;
    *by=4. ** y-1.;
    }
  else if(*in==4)
    {
    *bb=4. ** x*(1.-*x-*y);
    *bx=4.*(1.-2. ** x-*y);
    *by=-4. ** x;
    }
  else if(*in==5)
    {
    *bb=4. ** x ** y;
    *bx=4. ** y;
    *by=4. ** x;
    }
  else if(*in==6)
    {
    *bb=4. ** y*(1.-*x-*y);
    *bx=-4. ** y;
    *by=4.*(1.-*x-2. ** y);
    }
  else
    {
    *bb=0.;
    *bx=0.;
    *by=0.;
    }
  }

void trans(int*it,double*xq,double*yq,double*det,double*pj11,double*pj21,double*pj12,double*pj22,double xc[],double yc[],int*maxel)
  {
  int i1,i2,i3,i4,i5,i6;
  double f1x,f1y,f2x,f2y,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6;
  i1=node[0][*it-1];
  i2=node[1][*it-1];
  i3=node[2][*it-1];
  i4=node[3][*it-1];
  i5=node[4][*it-1];
  i6=node[5][*it-1];
  x1=xc[i1-1];
  y1=yc[i1-1];
  x2=xc[i2-1];
  y2=yc[i2-1];
  x3=xc[i3-1];
  y3=yc[i3-1];
  x4=xc[i4-1];
  y4=yc[i4-1];
  x5=xc[i5-1];
  y5=yc[i5-1];
  x6=xc[i6-1];
  y6=yc[i6-1];
  f1x=x1*(-3.+4. ** xq+4. ** yq)+x2*(4. ** xq-1.) +
   x4*4.*(1.-2. ** xq-*yq)+x5*4. ** yq+x6*4.*(-*yq);
  f1y=x1*(-3.+4. ** xq+4. ** yq)+x3*(4. ** yq-1.) +
   x4*4.*(-*xq)+x5*4. ** xq+x6*4.*(1.-*xq-2. ** yq);
  f2x=y1*(-3.+4. ** xq+4. ** yq)+y2*(4. ** xq-1.) +
   y4*4.*(1.-2. ** xq-*yq)+y5*4. ** yq+y6*4.*(-*yq);
  f2y=y1*(-3.+4. ** xq+4. ** yq)+y3*(4. ** yq-1.) +
   y4*4.*(-*xq)+y5*4. ** xq+y6*4.*(1.-*xq-2. ** yq);
  *det=f1x*f2y-f1y*f2x;
  *pj11=f2y/ *det;
  *pj12=-f2x/ *det;
  *pj21=-f1y/ *det;
  *pj22=f1x/ *det;
  *det=fabs(*det);
  }

void dgbfa(int*lda,int*n,int*ml,int*mu,int ipvt[],int*info)
  {
  int _i0,_i1,i,i0,j,j0,j1,ju,jz,k,kp1,l,lm,m,mm,nm1;
  double t;
  m=*ml+*mu+1;
  *info=0;
  j0=*mu+2;
  j1=min(*n,m)-1;
  for(jz=j0;jz<=j1;jz++)
    {
    i0=m+1-jz;
    for(i=i0;i<=*ml;i++)
      a[jz-1][i-1]=0.;
    }
  jz=j1;
  ju=0;
  nm1=*n-1;
  for(k=1;k<=nm1;k++)
    {
    kp1=k+1;
    jz+=1;
    if(jz<=*n)
      for(i=1;i<=*ml;i++)
        a[jz-1][i-1]=0.;
    lm=min(*ml,*n-k);
    l=idamax(ADR(_i0,lm+1),&a[k-1][m-1],ADR(_i1,1)) +
     m-1;
    ipvt[k-1]=l+k-m;
    if(a[k-1][l-1]==0.)
      *info=k;
    else
      {
      if(l!=m)
        {
        t=a[k-1][l-1];
        a[k-1][l-1]=a[k-1][m-1];
        a[k-1][m-1]=t;
        }
      t=-1./a[k-1][m-1];
      dscal(&lm,&t,&a[k-1][m],ADR(_i0,1));
      ju=min(max(ju,*mu+ipvt[k-1]),*n);
      mm=m;
      for(j=kp1;j<=ju;j++)
        {
        l-=1;
        mm-=1;
        t=a[j-1][l-1];
        if(l!=mm)
          {
          a[j-1][l-1]=a[j-1][mm-1];
          a[j-1][mm-1]=t;
          }
        daxpy(&lm,&t,&a[k-1][m],ADR(_i0,1),&a[j-1][mm],ADR(_i1,1));
        }
      }
    }
  ipvt[*n-1]=*n;
  if(a[*n-1][m-1]==0.)
    *info=*n;
  }

void dgbsl(int*lda,int*n,int*ml,int*mu,int ipvt[],double b[],int*job)
  {
  int _i0,_i1,k,kb,l,la,lb,lm,m,nm1;
  double t;
  m=*mu+*ml+1;
  nm1=*n-1;
  if(*job==0)
    {
    if(*ml!=0)
      {
      for(k=1;k<=nm1;k++)
        {
        lm=min(*ml,*n-k);
        l=ipvt[k-1];
        t=b[l-1];
        if(l!=k)
          {
          b[l-1]=b[k-1];
          b[k-1]=t;
          }
        daxpy(&lm,&t,&a[k-1][m],ADR(_i0,1),&b[k],ADR(_i1,1));
        }
      }
    for(kb=1;kb<=*n;kb++)
      {
      k=*n+1-kb;
      b[k-1]/=a[k-1][m-1];
      lm=min(k,m)-1;
      la=m-lm;
      lb=k-lm;
      t=-b[k-1];
      daxpy(&lm,&t,&a[k-1][la-1],ADR(_i0,1),&b[lb-1],ADR(_i1,1));
      }
    }
  else
    {
    for(k=1;k<=*n;k++)
      {
      lm=min(k,m)-1;
      la=m-lm;
      lb=k-lm;
      t=dot(&a[k-1][la-1],&b[lb-1],lm);
      b[k-1]=(b[k-1]-t)/a[k-1][m-1];
      }
    if(*ml!=0)
      {
      for(kb=1;kb<=nm1;kb++)
        {
        k=*n-kb;
        lm=min(*ml,*n-k);
        b[k-1]+=dot(&a[k-1][m],&b[k],lm);
        l=ipvt[k-1];
        if(l!=k)
          {
          t=b[l-1];
          b[l-1]=b[k-1];
          b[k-1]=t;
          }
        }
      }
    }
  }

void nstoke(double xc[],double yc[],double area[],double f[],double g[],double uold[],double*reynld,double*tolns,double*xlngth,double*ylngth,int ipivot[],int*mrow1,int*nlband,int*nuband,int*nband,int*nrow1,int*ncol1,int*nelemn,int*np,int*nnodes,int*nuk,int*nquad,int*neqn,int*nsteps,int*nsim,int*maxnd,int*maxel,double*rdel,double*alpha)
  {
  int _i0,i,info,ip,ipp,iq,iqq,iquad,it,iter,iuk,iukk,iun,iuse,j,job,kk,niter,unk_u,unk_v;
  double aij,ar,arr,bb,bbb,bbbl,bbl,bbx,bby,bx,by,csim,det,diff,etax,etay,tbbx,tbby,tbx,tby,ubc,un[2],unx[2],uny[2],uold_qp,visc,vold_qp,x,xix,xiy,y;
  visc=1./ *reynld;
  for(iter=1;iter<=*nsteps;iter++)
    {
    niter=iter;
    if(iter<*nsim)
      csim=0.;
    else
      csim=1.;
    for(i=1;i<=*nrow1;i++)
      for(j=1;j<=*ncol1;j++)
        a[j-1][i-1]=0.;
    for(it=1;it<=*nelemn;it++)
      {
      arr=area[it-1]/3.;
      for(iquad=1;iquad<=*nquad;iquad++)
        {
        y=ym[iquad-1][it-1];
        x=xm[iquad-1][it-1];
        trans(&it,&x,&y,&det,&xix,&xiy,&etax,&etay,xc,yc,maxel);
        ar=arr*det;
        for(kk=1;kk<=2;kk++)
          {
          un[kk-1]=0.;
          uny[kk-1]=0.;
          unx[kk-1]=0.;
          }
        uold_qp=0.;
        vold_qp=0.;
        for(iq=1;iq<=*nnodes;iq++)
          {
          refqbf(&x,&y,&iq,&bb,&tbx,&tby);
          bx=tbx*xix+tby*etax;
          by=tbx*xiy+tby*etay;
          ip=node[iq-1][it-1];
          for(iuk=1;iuk<=2;iuk++)
            {
            iun=indx[iuk-1][ip-1];
            if(0<iun)
              {
              un[iuk-1]+=bb*g[iun-1];
              unx[iuk-1]+=bx*g[iun-1];
              uny[iuk-1]+=by*g[iun-1];
              if(iuk==1)
                uold_qp+=bb*uold[iun-1];
              if(iuk==2)
                vold_qp+=bb*uold[iun-1];
              }
            else if(iun<0)
              {
              ubc=*alpha*ubdry(ADR(_i0,1),&iuk,&ip,xc,yc);
              un[iuk-1]+=bb*ubc;
              unx[iuk-1]+=bx*ubc;
              uny[iuk-1]+=by*ubc;
              if(iuk==1)
                uold_qp+=bb*ubc;
              if(iuk==2)
                vold_qp+=bb*ubc;
              }
            }
          }
        for(iq=1;iq<=*nnodes;iq++)
          {
          ip=node[iq-1][it-1];
          refqbf(&x,&y,&iq,&bb,&tbx,&tby);
          bx=tbx*xix+tby*etax;
          by=tbx*xiy+tby*etay;
          if(iq<=3)
            bbl=refbsp(&x,&y,&iq);
          for(iuk=1;iuk<=*nuk;iuk++)
            {
            i=indx[iuk-1][ip-1];
            if(i>0)
              {
              if(iuk==1)
                {
                f[i-1]+=csim*((un[0]*unx[0] +
                  un[1]*uny[0])*bb)*ar+*rdel*uold_qp *
                 bb*ar;
                }
              else if(iuk==2)
                {
                f[i-1]+=csim*((un[0]*unx[1] +
                  un[1]*uny[1])*bb)*ar+*rdel*vold_qp *
                 bb*ar;
                }
              for(iqq=1;iqq<=*nnodes;iqq++)
                {
                ipp=node[iqq-1][it-1];
                refqbf(&x,&y,&iqq,&bbb,&tbbx,&tbby);
                bbx=tbbx*xix+tbby*etax;
                bby=tbbx*xiy+tbby*etay;
                if(iqq<=3)
                  bbbl=refbsp(&x,&y,&iqq);
                for(iukk=1;iukk<=*nuk;iukk++)
                  {
                  j=indx[iukk-1][ipp-1];
                  if(j!=0)
                    {
                    aij=0.;
                    if(i!=*neqn)
                      {
                      if(iuk==1)
                        {
                        if(iukk==1)
                          aij=visc*(by*bby+bx*bbx)+(bbb*unx[0]*bb)*csim+bb*bbx*un[0]+bb*bby*un[1]+*rdel*(bb*bbb);
                        else if(iukk==2)
                          aij=csim*(bb*bbb*uny[0]);
                        else if(iukk==3)
                          aij=-bx*bbbl;
                        }
                      else if(iuk==2)
                        {
                        if(iukk==1)
                          aij=csim*(bb*bbb*unx[1]);
                        else if(iukk==2)
                          aij=(visc*(by*bby+bx*bbx)+(bb*bbb*uny[1])*csim+bb*bby*un[1]+bb*bbx*un[0])+*rdel*(bb*bbb);
                        else if(iukk==3)
                          aij=-by*bbbl;
                        }
                      else if(iuk==3)
                        {
                        if(iukk==1)
                          aij=bbx*bbl;
                        else if(iukk==2)
                          aij=bby*bbl;
                        }
                      if(0<j)
                        {
                        iuse=i-j+*nband;
                        a[j-1][iuse-1]+=aij*ar;
                        }
                      else
                        f[i-1]+=-ar ** alpha*ubdry(&iukk,&ipp,&j,xc,yc)*aij;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    f[*neqn-1]=0.;
    for(j=*neqn-*nlband;j<=(*neqn-1);j++)
      {
      i=*neqn-j+*nband;
      a[j-1][i-1]=0.;
      }
    a[*neqn-1][*nband-1]=1.;
    job=0;
    dgbfa(mrow1,neqn,nlband,nuband,ipivot,&info);
    dgbsl(mrow1,neqn,nlband,nuband,ipivot,f,&job);
    diff=0.;
    for(i=1;i<=*np;i++)
      {
      unk_u=indx[0][i-1];
      if(0<unk_u)
        diff+=sq(g[unk_u-1]-f[unk_u-1]);
      unk_v=indx[1][i-1];
      if(0<unk_v)
        diff+=sq(g[unk_v-1]-f[unk_v-1]);
      }
    diff=sqrt(diff);
    printf("iteration %i difference is %lG\n",iter,diff);
    if(diff<=*tolns)
      break;
    for(i=1;i<=*neqn;i++)
      {
      g[i-1]=f[i-1];
      f[i-1]=0.;
      }
    }
  }

void setgrd(double xc[],double yc[],double area[],int ireg_density_x[],int ireg_density_y[],double region_x[],double region_y[],int*nlband,int*nuband,int*nband,int*nelemn,int*np,int*nnodes,int*nuk,int*nquad,int*neqn,int*maxnd,int*maxel)
  {
  hcell_node_count(ireg_density_x,ireg_density_y,np);
  hcell_node_xy(ireg_density_x,ireg_density_y,np,region_x,region_y,xc,yc);
  hcell_element_count(ireg_density_x,ireg_density_y,nelemn);
  hcell_element_node(ireg_density_x,ireg_density_y,maxel,nelemn,nnodes);
  hcell_dof_set(ireg_density_x,ireg_density_y,maxnd,np,neqn);
  quad_a_set(maxel,nelemn,nquad,area);
  element_node_bandwidth(maxnd,maxel,nnodes,nelemn,neqn,np,nlband,nuband,nband);
  }

#define MAXEL (2*(NX-1)*(NY-1))
#define MAXND (MX*MY)
#define MAXUN (2*MX*MY+NX*NY)
#define MINUN (27*NY)
#define MX (2*NX-1)
#define MY (2*NY-1)
#define N_TIME 10
#define NNODES 6
#define NQUAD 3
#define NUK 3
#define NX1 10
#define NX2 10
#define NX3 10
#define NX (NX1+NX2+NX3+1)
#define NY1 10
#define NY2 10
#define NY3 10
#define NY (NY1+NY2+NY3+1)

int ipivot[MAXUN];
double area[MAXEL];
double f[MAXUN];
double g[MAXUN];
double uold[MAXUN];
double xc[MAXND];
double yc[MAXND];

int main(int argc,char**argv,char**envp)
  {
  int ireg_density_x[3],ireg_density_y[3];
  int _i0,_i1,_i2,_i3,_i4,i,ibc_type,iter,mrow1,n,nband,ncol1,nelemn,
    neqn,nlband,np,nrow1,nsim,nsteps,nuband;
  double region_x[4],region_y[4];
  double alpha,deltat,p,rdel,reynld,tolns,tolopt,u,v,xlngth,ylngth;
  FILE*fp,*fq;
  printf("solve the navier stokes fluid flow\n");
  printf("equations in an h-shaped region,\n");
  printf("using finite elements.\n");
  ibc_type=1;
  ireg_density_x[0]=NX1;
  ireg_density_x[1]=NX2;
  ireg_density_x[2]=NX3;
  ireg_density_y[0]=NY1;
  ireg_density_y[1]=NY2;
  ireg_density_y[2]=NY3;
  region_x[0]=0.;
  region_x[1]=3.;
  region_x[2]=6.;
  region_x[3]=9.;
  region_y[0]=0.;
  region_y[1]=3.;
  region_y[2]=6.;
  region_y[3]=9.;
  reynld=1.;
  mrow1=MINUN;
  nsim=3;
  nsteps=15;
  tolns=1E-6;
  tolopt=1E-6;
  printf("maximum number of nodes=%i\n",MAXND);
  printf("maximum number of elements=%i\n",MAXEL);
  printf("maximum number of unknowns=%i\n",MAXUN);
  printf("maximum matrix dimension=%i\n",MINUN);
  setgrd(xc,yc,area,ireg_density_x,ireg_density_y,region_x,region_y,&nlband,&nuband,&nband,&nelemn,&np,ADR(_i0,NNODES),ADR(_i1,NUK),ADR(_i2,NQUAD),&neqn,ADR(_i3,MAXND),ADR(_i4,MAXEL));
  printf("number of nodes=%i\n",np);
  printf("number of elements=%i\n",nelemn);
  printf("number of unknowns=%i\n",neqn);
  printf("writing: hcell.2dv\n");
  fp=fopen("hcell.2dv","wt");
  fprintf(fp,"%i\n",np);
  for(i=1;i<=np;i++)
    fprintf(fp,"%lG %lG\n",xc[i-1],yc[i-1]);
  fprintf(fp,"%i\n",nelemn);
  for(i=1;i<=nelemn;i++)
    fprintf(fp,"%i %i %i\n",node[0][i-1],node[1][i-1],node[2][i-1]);
  fclose(fp);
  nrow1=nlband+nlband+nuband+1;
  ncol1=neqn;
  for(i=1;i<=neqn;i++)
    f[i-1]=0.;
  for(i=1;i<=neqn;i++)
    uold[i-1]=0.5;
  deltat=0.0002;
  rdel=1./deltat;
  printf("delta t=%lG\n",deltat);
  for(iter=1;iter<=N_TIME;iter++)
    {
    printf("time step %i\n",iter);
    if(ibc_type==1)
      {
      if(iter<=250)
        alpha=5.;
      else
        alpha=1.;
      }
    else if(ibc_type==2)
      {
      if(iter<=250)
        alpha=80.*iter*deltat+1.;
      else
        alpha=-80.*iter*deltat+9.;
      }
    else if(ibc_type==3)
      alpha=2.*sin(0.01*iter*M_PI);
    for(i=1;i<=neqn;i++)
      g[i-1]=f[i-1];
    for(i=1;i<=neqn;i++)
      f[i-1]=0.;
    nstoke(xc,yc,area,f,g,uold,&reynld,&tolns,&xlngth,&ylngth,ipivot,&mrow1,&nlband,&nuband,&nband,&nrow1,&ncol1,&nelemn,&np,ADR(_i0,NNODES),ADR(_i1,NUK),ADR(_i2,NQUAD),&neqn,&nsteps,&nsim,ADR(_i3,MAXND),ADR(_i4,MAXEL),&rdel,&alpha);
    for(i=1;i<=neqn;i++)
      uold[i-1]=f[i-1];
    }
  printf("writing: hcell.v2d\n");
  fp=fopen("hcell.v2d","wt");
  printf("writing: hcell.plt\n");
  fq=fopen("hcell.plt","wt");
  fprintf(fq,"TITLE=\"H-CELL\"\n");
  fprintf(fq,"VARIABLES=\"X\" \"Y\" \"P\" \"U\" \"V\"\n");
  fprintf(fq,"ZONE T=\"RESULTS\" N=%i, E=%i, F=FEPOINT, ET=TRIANGLE\n",np,nelemn);
  fprintf(fq,"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n");
  for(n=1;n<=np;n++)
    {
    i=indx[0][n-1];
    if(0<i)
      u=f[i-1];
    else if(i<0)
      u=ubdry(ADR(_i0,1),&n,&i,xc,yc);
    else
      u=0.;
    i=indx[1][n-1];
    if(0<i)
      v=f[i-1];
    else if(i<0)
      v=ubdry(ADR(_i0,2),&n,&i,xc,yc);
    else
      v=0.;
    i=indx[2][n-1];
    if(0<i)
      p=f[i-1];
    else
      p=0.;
    fprintf(fp,"%lG %lG %lG %lG\n",xc[n-1],yc[n-1],u,v);
    fprintf(fq,"%lG %lG %lG %lG %lG\n",xc[n-1],yc[n-1],p,u,v);
    }
  fclose(fp);
  for(i=1;i<=nelemn;i++)
    fprintf(fq,"%i %i %i\n",node[0][i-1],node[1][i-1],node[2][i-1]);
  fclose(fq);
  return(0);
  }
