/* Modified slightly by D. Orchard (2010) from the classic code from:
   Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer, Numerical
   Simulation in Fluid Dynamics, SIAM, 1998.
   http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

   V1.00 (12/05/2019)
   This code has also been streamlined by dudley.benton@gmail.com
   and modified to compile with Visual Studio and run on Windows.
   The .PPM file output has been converted to .BMP for compatibility.

   V1.10 (03/01/2021)
   reduce 24-bit to 8-bit
   add rainbow option
   add split 3 views
   add embedded text

   V1.20 (03/29/2021)
   eliminate 2D matrices
   elevate arrays and parameters to global level
   reduce argument passing
   streamline calculations
   enable full optimizations
*/
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include "12x08.h"
#define Font C12x08

#define rainbow
#undef  CONTRACTION
#undef  MAKE_TWO
#undef  MAKE_THREE
#undef  MAKE_FOUR
#undef  MAKE_OVAL
#undef  MINIVAN
#undef  AIRFOIL
#define SUPERGIRL

typedef struct{double x,y;}XY;

XY minivan[]={
  {-0.607,-0.382},{-0.691,-0.350},{-0.754,-0.285},{-0.830,-0.252},
  {-0.919,-0.233},{-1.000,-0.198},{-0.997,-0.108},{-0.979,-0.020},
  {-0.918, 0.045},{-0.836, 0.083},{-0.751, 0.114},{-0.661, 0.129},
  {-0.572, 0.144},{-0.483, 0.160},{-0.396, 0.185},{-0.318, 0.230},
  {-0.239, 0.274},{-0.156, 0.311},{-0.071, 0.341},{ 0.017, 0.363},
  { 0.107, 0.377},{ 0.197, 0.382},{ 0.288, 0.382},{ 0.378, 0.376},
  { 0.468, 0.368},{ 0.558, 0.357},{ 0.648, 0.346},{ 0.738, 0.335},
  { 0.811, 0.306},{ 0.862, 0.231},{ 0.929, 0.172},{ 0.972, 0.098},
  { 0.984, 0.008},{ 0.995,-0.081},{ 1.000,-0.172},{ 0.923,-0.196},
  { 0.841,-0.233},{ 0.800,-0.313},{ 0.728,-0.366},{ 0.683,-0.382},
  { 0.638,-0.370},{ 0.559,-0.329},{ 0.503,-0.258},{ 0.428,-0.235},
  { 0.337,-0.238},{ 0.247,-0.240},{ 0.156,-0.242},{ 0.066,-0.244},
  {-0.025,-0.247},{-0.115,-0.249},{-0.206,-0.252},{-0.296,-0.254},
  {-0.387,-0.256},{-0.455,-0.297},{-0.521,-0.357}};

XY airfoil[]={
  { 2.049, 0.192},{ 1.955, 0.225},{ 1.858, 0.252},{ 1.761, 0.277},
  { 1.665, 0.304},{ 1.567, 0.327},{ 1.469, 0.349},{ 1.371, 0.371},
  { 1.273, 0.392},{ 1.175, 0.413},{ 1.076, 0.431},{ 0.978, 0.447},
  { 0.879, 0.464},{ 0.779, 0.478},{ 0.680, 0.490},{ 0.580, 0.502},
  { 0.481, 0.514},{ 0.381, 0.521},{ 0.281, 0.528},{ 0.180, 0.530},
  { 0.080, 0.530},{-0.020, 0.530},{-0.120, 0.530},{-0.221, 0.526},
  {-0.320, 0.515},{-0.420, 0.501},{-0.518, 0.484},{-0.617, 0.463},
  {-0.713, 0.437},{-0.807, 0.401},{-0.897, 0.358},{-0.982, 0.304},
  {-1.043, 0.226},{-1.037, 0.134},{-0.949, 0.088},{-0.852, 0.064},
  {-0.752, 0.052},{-0.652, 0.047},{-0.552, 0.041},{-0.452, 0.042},
  {-0.351, 0.048},{-0.251, 0.052},{-0.151, 0.055},{-0.051, 0.061},
  { 0.049, 0.068},{ 0.149, 0.069},{ 0.250, 0.075},{ 0.350, 0.082},
  { 0.450, 0.087},{ 0.550, 0.096},{ 0.650, 0.098},{ 0.750, 0.103},
  { 0.850, 0.110},{ 0.950, 0.117},{ 1.050, 0.124},{ 1.150, 0.127},
  { 1.250, 0.133},{ 1.350, 0.140},{ 1.451, 0.143},{ 1.551, 0.149},
  { 1.651, 0.153},{ 1.751, 0.155},{ 1.851, 0.159},{ 1.952, 0.162},
  { 2.052, 0.162}};

XY supergirl[]={
  { 0.034, 0.408},{ 0.079, 0.408},{ 0.123, 0.407},{ 0.167, 0.401},
  { 0.210, 0.390},{ 0.251, 0.373},{ 0.291, 0.354},{ 0.332, 0.336},
  { 0.373, 0.318},{ 0.414, 0.301},{ 0.455, 0.284},{ 0.503, 0.293},
  { 0.551, 0.305},{ 0.600, 0.318},{ 0.648, 0.330},{ 0.697, 0.339},
  { 0.747, 0.345},{ 0.796, 0.346},{ 0.846, 0.340},{ 0.894, 0.317},
  { 0.929, 0.279},{ 0.974, 0.245},{ 1.031, 0.234},{ 1.083, 0.230},
  { 1.135, 0.226},{ 1.187, 0.223},{ 1.239, 0.219},{ 1.292, 0.215},
  { 1.344, 0.211},{ 1.396, 0.207},{ 1.448, 0.202},{ 1.500, 0.196},
  { 1.552, 0.189},{ 1.603, 0.181},{ 1.655, 0.173},{ 1.706, 0.163},
  { 1.758, 0.154},{ 1.809, 0.144},{ 1.830, 0.183},{ 1.877, 0.188},
  { 1.918, 0.161},{ 1.989, 0.114},{ 2.036, 0.077},{ 2.040, 0.050},
  { 1.983, 0.058},{ 1.936, 0.071},{ 1.884, 0.080},{ 1.837, 0.080},
  { 1.788, 0.080},{ 1.738, 0.080},{ 1.689, 0.081},{ 1.639, 0.083},
  { 1.590, 0.086},{ 1.540, 0.088},{ 1.491, 0.089},{ 1.441, 0.089},
  { 1.392, 0.088},{ 1.342, 0.085},{ 1.293, 0.082},{ 1.244, 0.078},
  { 1.194, 0.073},{ 1.145, 0.069},{ 1.096, 0.065},{ 1.046, 0.062},
  { 0.997, 0.059},{ 0.947, 0.056},{ 0.898, 0.053},{ 0.848, 0.051},
  { 0.799, 0.049},{ 0.749, 0.049},{ 0.700, 0.050},{ 0.650, 0.053},
  { 0.601, 0.059},{ 0.552, 0.066},{ 0.503, 0.075},{ 0.455, 0.085},
  { 0.407, 0.097},{ 0.359, 0.110},{ 0.311, 0.123},{ 0.272, 0.125},
  { 0.214, 0.125},{ 0.165, 0.118},{ 0.107, 0.102},{ 0.050, 0.097},
  { 0.011, 0.123},{-0.017, 0.164},{-0.055, 0.185},{-0.104, 0.176},
  {-0.152, 0.163},{-0.199, 0.149},{-0.246, 0.132},{-0.293, 0.114},
  {-0.338, 0.094},{-0.383, 0.072},{-0.426, 0.047},{-0.462, 0.033},
  {-0.512, 0.028},{-0.561, 0.037},{-0.610, 0.044},{-0.660, 0.051},
  {-0.709, 0.058},{-0.759, 0.064},{-0.808, 0.070},{-0.857, 0.077},
  {-0.907, 0.083},{-0.956, 0.089},{-1.006, 0.095},{-1.020, 0.140},
  {-0.970, 0.146},{-0.921, 0.149},{-0.871, 0.154},{-0.821, 0.151},
  {-0.772, 0.147},{-0.724, 0.147},{-0.685, 0.143},{-0.638, 0.142},
  {-0.590, 0.136},{-0.531, 0.140},{-0.477, 0.157},{-0.431, 0.177},
  {-0.386, 0.197},{-0.340, 0.217},{-0.295, 0.237},{-0.249, 0.257},
  {-0.224, 0.278},{-0.265, 0.306},{-0.307, 0.334},{-0.347, 0.363},
  {-0.386, 0.395},{-0.419, 0.431},{-0.440, 0.476},{-0.434, 0.525},
  {-0.407, 0.567},{-0.372, 0.602},{-0.332, 0.632},{-0.287, 0.646},
  {-0.238, 0.644},{-0.188, 0.638},{-0.141, 0.624},{-0.084, 0.589},
  {-0.047, 0.528},{-0.028, 0.453},{-0.010, 0.408}};

int InsidePolygon(XY*p,int n,double x,double y)
  {
  int above1,above2,i,right=0;
  double x1,x2,y1,y2;
  x2=p[n-1].x;
  y2=p[n-1].y;
  above2=y2>y?1:0;
  for(i=0;i<n;i++)
    {
    x1=x2;
    y1=y2;
    x2=p[i].x;
    y2=p[i].y;
    above1=above2;
    above2=y2>y?1:0;
    if(above1==above2)
      continue;
    if(x1>x&&x2>x)
      right++;
    else if(y1<y2)
      {
      if((x-x1)*(y2-y1)<(x2-x1)*(y-y1))
        right++;
      }
    else if(y1>y2)
      {
      if((x-x1)*(y2-y1)>(x2-x1)*(y-y1))
        right++;
      }
    }
  return(right&1);
  }

#define C_B      0x0000 /* This cell is an obstacle/boundary cell */
#define B_N      0x0001 /* This obstacle cell has a fluid cell to the north */
#define B_S      0x0002 /* This obstacle cell has a fluid cell to the south */
#define B_W      0x0004 /* This obstacle cell has a fluid cell to the west */
#define B_E      0x0008 /* This obstacle cell has a fluid cell to the east */
#define B_NW     (B_N|B_W)
#define B_SW     (B_S|B_W)
#define B_NE     (B_N|B_E)
#define B_SE     (B_S|B_E)
#define B_NSEW   (B_N|B_S|B_E|B_W)
#define C_F      0x0010 /* This cell is a fluid cell */

/* Macros for poisson(),denoting whether there is an obstacle cell
   adjacent to some direction */

#define eps_E ((flag[i+1][j]&C_F)?1:0)
#define eps_W ((flag[i-1][j]&C_F)?1:0)
#define eps_N ((flag[i][j+1]&C_F)?1:0)
#define eps_S ((flag[i][j-1]&C_F)?1:0)

DWORD palette[]={
#ifdef rainbow
  0x000000,0xFFFFFF,0x0000FC,0x0008FC,0x0010FC,0x0018FC,0x0020FC,0x0028FC,
  0x0030FC,0x0038FC,0x0040FC,0x0048FC,0x0050FC,0x0058FC,0x0060FC,0x0068FC,
  0x0070FC,0x0078FC,0x0080FC,0x0088FC,0x0090FC,0x0098FC,0x00A0FC,0x00A8FC,
  0x00B0FC,0x00B8FC,0x00C0FC,0x00C8FC,0x00D0FC,0x00D8FC,0x00E0FC,0x00E8FC,
  0x00F0FC,0x00FCFC,0x00FCF0,0x00FCE8,0x00FCE0,0x00FCD8,0x00FCD0,0x00FCC8,
  0x00FCC0,0x00FCB8,0x00FCB0,0x00FCA8,0x00FCA0,0x00FC98,0x00FC90,0x00FC88,
  0x00FC80,0x00FC78,0x00FC70,0x00FC68,0x00FC60,0x00FC58,0x00FC50,0x00FC48,
  0x00FC40,0x00FC38,0x00FC30,0x00FC28,0x00FC20,0x00FC18,0x00FC10,0x00FC08,
  0x00FC00,0x08FC00,0x10FC00,0x18FC00,0x20FC00,0x28FC00,0x30FC00,0x38FC00,
  0x40FC00,0x48FC00,0x50FC00,0x58FC00,0x60FC00,0x68FC00,0x70FC00,0x78FC00,
  0x80FC00,0x88FC00,0x90FC00,0x98FC00,0xA0FC00,0xA8FC00,0xB0FC00,0xB8FC00,
  0xC0FC00,0xC8FC00,0xD0FC00,0xD8FC00,0xE0FC00,0xE8FC00,0xF0FC00,0xFCFC00,
  0xFCF000,0xFCE800,0xFCE000,0xFCD800,0xFCD000,0xFCC800,0xFCC000,0xFCB800,
  0xFCB000,0xFCA800,0xFCA000,0xFC9800,0xFC9000,0xFC8800,0xFC8000,0xFC7800,
  0xFC7000,0xFC6800,0xFC6000,0xFC5800,0xFC5000,0xFC4800,0xFC4000,0xFC3800,
  0xFC3000,0xFC2800,0xFC2000,0xFC1800,0xFC1000,0xFC0800,0xFC0000};
#else
  0xFF0000,0x00FF00,0x000000,0x020202,0x040404,0x060606,0x080808,0x0A0A0A,
  0x0C0C0C,0x0E0E0E,0x101010,0x121212,0x141414,0x161616,0x181818,0x1A1A1A,
  0x1C1C1C,0x1E1E1E,0x202020,0x222222,0x252525,0x272727,0x292929,0x2B2B2B,
  0x2D2D2D,0x2F2F2F,0x313131,0x333333,0x353535,0x373737,0x393939,0x3B3B3B,
  0x3D3D3D,0x3F3F3F,0x414141,0x434343,0x454545,0x474747,0x4A4A4A,0x4C4C4C,
  0x4E4E4E,0x505050,0x525252,0x545454,0x565656,0x585858,0x5A5A5A,0x5C5C5C,
  0x5E5E5E,0x606060,0x626262,0x646464,0x666666,0x686868,0x6A6A6A,0x6C6C6C,
  0x6F6F6F,0x717171,0x737373,0x757575,0x777777,0x797979,0x7B7B7B,0x7D7D7D,
  0x7F7F7F,0x818181,0x838383,0x858585,0x878787,0x898989,0x8B8B8B,0x8D8D8D,
  0x8F8F8F,0x929292,0x949494,0x969696,0x989898,0x9A9A9A,0x9C9C9C,0x9E9E9E,
  0xA0A0A0,0xA2A2A2,0xA4A4A4,0xA6A6A6,0xA8A8A8,0xAAAAAA,0xACACAC,0xAEAEAE,
  0xB0B0B0,0xB2B2B2,0xB4B4B4,0xB7B7B7,0xB9B9B9,0xBBBBBB,0xBDBDBD,0xBFBFBF,
  0xC1C1C1,0xC3C3C3,0xC5C5C5,0xC7C7C7,0xC9C9C9,0xCBCBCB,0xCDCDCD,0xCFCFCF,
  0xD1D1D1,0xD3D3D3,0xD5D5D5,0xD7D7D7,0xD9D9D9,0xDCDCDC,0xDEDEDE,0xE0E0E0,
  0xE2E2E2,0xE4E4E4,0xE6E6E6,0xE8E8E8,0xEAEAEA,0xECECEC,0xEEEEEE,0xF0F0F0,
  0xF2F2F2,0xF4F4F4,0xF6F6F6,0xF8F8F8,0xFAFAFA,0xFCFCFC,0xFFFFFF};
#endif

#define imax 660 /* Number of cells horizontally */
#define jmax 120 /* Number of cells vertically */

double u[imax+2][jmax+2];
double v[imax+2][jmax+2];
double f[imax+2][jmax+2];
double g[imax+2][jmax+2];
double p[imax+2][jmax+2];
double rhs[imax+2][jmax+2];
double psi[imax+2][jmax+2];
double zeta[imax+2][jmax+2];
char flag[imax+2][jmax+2];
double delx,dely;

double delt=0.003; /* Duration of each timestep */
double eps=0.001;  /* Stopping error threshold for SOR */
double gamma=0.9;  /* Upwind differencing factor in PDE discretisation */
double omega=1.7;  /* Relaxation parameter for SOR */
double Re=150.;    /* Reynolds number */
double tau=0.5;    /* Safety factor for timestep control */
double t_end=40.;  /* Simulation runtime */
double ui=1.;      /* Initial X velocity */
double vi=0.;      /* Initial Y velocity */
double xlength=22.;/* Width of simulated domain */
double ylength=4.1;/* Height of simulated domain */

/* Given the boundary conditions defined by the flag matrix,update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix. */
void ApplyBoundaryConditions(double ui,double vi)
  {
  int i,j;
  for(j=0;j<=jmax+1;j++)
    {
    /* Fluid freely flows in from the west */
    u[0][j]=u[1][j];
    v[0][j]=v[1][j];
    /* Fluid freely flows out to the east */
    u[imax][j]=u[imax-1][j];
    v[imax+1][j]=v[imax][j];
    }
  for(i=0;i<=imax+1;i++)
    {
    /* The vertical velocity approaches 0 at the north and south boundaries,but fluid flows freely in the horizontal direction */
    v[i][jmax]=0.;
    u[i][jmax+1]=u[i][jmax];
    v[i][0]=0.;
    u[i][0]=u[i][1];
    }
  /* Apply no-slip boundary conditions to cells that are adjacent to internal obstacle cells. This forces the u and v velocity
  *to tend towards zero in these cells. */
  for(i=1;i<=imax;i++)
    {
    for(j=1;j<=jmax;j++)
      {
      if(flag[i][j]&B_NSEW)
        {
        switch (flag[i][j])
          {
         case B_N:
          u[i][j]=-u[i][j+1];
          break;
         case B_E:
          u[i][j]=0.;
          break;
         case B_NE:
          u[i][j]=0.;
          break;
         case B_SE:
          u[i][j]=0.;
          break;
         case B_NW:
          u[i][j]=-u[i][j+1];
          break;
         case B_S:
          u[i][j]=-u[i][j-1];
          break;
         case B_SW:
          u[i][j]=-u[i][j-1];
          break;
          }
        }
      }
    }
  for(i=0;i<=(imax-1);i++)
    {
    for(j=1;j<=jmax;j++)
      {
      if(flag[i+1][j]&B_NSEW)
        {
        switch (flag[i+1][j])
          {
         case B_N:
          u[i][j]=-u[i][j+1];
          break;
         case B_W:
          u[i][j]=0.;
          break;
         case B_NE:
          u[i][j]=-u[i][j+1];
          break;
         case B_SW:
          u[i][j]=0.;
          break;
         case B_NW:
          u[i][j]=0.;
          break;
         case B_S:
          u[i][j]=-u[i][j-1];
          break;
         case B_SE:
          u[i][j]=-u[i][j-1];
          break;
          }
        }
      }
    }
  for(i=1;i<=imax;i++)
    {
    for(j=1;j<=jmax;j++)
      {
      if(flag[i][j]&B_NSEW)
        {
        switch (flag[i][j])
          {
         case B_N:
          v[i][j]=0.;
          break;
         case B_E:
          v[i][j]=-v[i+1][j];
          break;
         case B_NE:
          v[i][j]=0.;
          break;
         case B_SE:
          v[i][j]=-v[i+1][j];
          break;
         case B_NW:
          v[i][j]=0.;
          break;
         case B_W:
          v[i][j]=-v[i-1][j];
          break;
         case B_SW:
          v[i][j]=-v[i-1][j];
          break;
          }
        }
      }
    }
  for(i=1;i<=imax;i++)
    {
    for(j=0;j<=(jmax-1);j++)
      {
      if(flag[i][j+1]&B_NSEW)
        {
        switch (flag[i][j+1])
          {
         case B_E:
          v[i][j]=-v[i+1][j];
          break;
         case B_S:
          v[i][j]=0.;
          break;
         case B_NE:
          v[i][j]=-v[i+1][j];
          break;
         case B_SE:
          v[i][j]=0.;
          break;
         case B_SW:
          v[i][j]=0.;
          break;
         case B_W:
          v[i][j]=-v[i-1][j];
          break;
         case B_NW:
          v[i][j]=-v[i-1][j];
          break;
          }
        }
      }
    }
  /* Finally,fix the horizontal velocity at the  western edge to have a continual flow of fluid into the simulation. */
  v[0][0]=2*vi-v[1][0];
  for(j=1;j<=jmax;j++)
    {
    u[0][j]=ui;
    v[0][j]=2*vi-v[1][j];
    }
  }

/* Initialize the flag array,marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too. */
void InitFlag(int*ibound)
  {
  int i,j;
  double mx,my,x,y,rad1;
  /* Mask of a circular obstacle */
  mx=20.0/41.0*jmax*dely;
  my=mx;
  rad1=5.0/41.0*jmax*dely;
#if(defined(CONTRACTION))
  for(i=1;i<=imax;i++)
    for(j=1;j<=jmax;j++)
      flag[i][j]=C_F;
  for(i=imax/4,jj=1;jj<jmax/4;i++,jj++)
    for(j=0;j<jj;j++)
      flag[i][j+1]=flag[i][jmax-j]=C_B;
  for(;jj>0;i++,jj--)
    for(j=0;j<jj;j++)
      flag[i][j+1]=flag[i][jmax-j]=C_B;
#else
  for(i=1;i<=imax;i++)
    {
    for(j=1;j<=jmax;j++)
      {
#if(defined(MINIVAN))
      x=(i-imax/5.0)*delx/3.;
      y=(j-jmax/3.5)*delx/3.;
      if(InsidePolygon(minivan,sizeof(minivan)/sizeof(minivan[0]),x,y))
        flag[i][j]=C_B;
      else
        flag[i][j]=C_F;
#elif(defined(AIRFOIL))
      x=(i-imax/5.0)*delx/3.;
      y=(j-jmax/3.5)*delx/3.;
      if(InsidePolygon(airfoil,sizeof(airfoil)/sizeof(airfoil[0]),x,y))
        flag[i][j]=C_B;
      else
        flag[i][j]=C_F;
#elif(defined(SUPERGIRL))
      x=(i-imax/5.0)*delx/3.;
      y=(j-jmax/3.5)*delx/3.;
      if(InsidePolygon(supergirl,sizeof(supergirl)/sizeof(supergirl[0]),x,y))
        flag[i][j]=C_B;
      else
        flag[i][j]=C_F;
#elif(defined(MAKE_OVAL))
      x=((i-0.5)*delx-2*mx)/4.;
      y=(j-0.5)*dely-my;
      r=hypot(x,y);
      flag[i][j]=(r<=rad1)?C_B:C_F;
#elif(defined(MAKE_TWO))
      x=(i-0.5)*delx-mx;
      y=(j-0.5)*dely-my;
      r=hypot(x,y);
      x=(i-0.5)*delx-(mx+9.*rad1);
      r=min(r,hypot(x,y));
      flag[i][j]=(r<=rad1)?C_B:C_F;
#elif(defined(MAKE_THREE))
      x=(i-0.5)*delx-mx;
      y=(j-0.5)*dely-my;
      r=hypot(x,y);
      x=(i-0.5)*delx-(mx+9.*rad1);
      r=min(r,hypot(x,y));
      x=(i-0.5)*delx-(mx+18.*rad1);
      r=min(r,hypot(x,y));
      flag[i][j]=(r<=rad1)?C_B:C_F;
#elif(defined(MAKE_FOUR))
      x=(i-0.5)*delx-mx;
      y=(j-0.5)*dely-my;
      r=hypot(x,y);
      x=(i-0.5)*delx-(mx+9.*rad1);
      r=min(r,hypot(x,y));
      x=(i-0.5)*delx-(mx+18.*rad1);
      r=min(r,hypot(x,y));
      x=(i-0.5)*delx-(mx+27.*rad1);
      r=min(r,hypot(x,y));
      flag[i][j]=(r<=rad1)?C_B:C_F;
#else
      x=(i-0.5)*delx-mx;
      y=(j-0.5)*dely-my;
      r=hypot(x,y);
      flag[i][j]=(r<=rad1)?C_B:C_F;
#endif
      }
    }
#endif
  /* Mark the north&south boundary cells */
  for(i=0;i<=imax+1;i++)
    {
    flag[i][0]=C_B;
    flag[i][jmax+1]=C_B;
    }
  /* Mark the east and west boundary cells */
  for(j=1;j<=jmax;j++)
    {
    flag[0][j]=C_B;
    flag[imax+1][j]=C_B;
    }
  /* flags for boundary cells */
  *ibound=0;
  for(i=1;i<=imax;i++)
    {
    for(j=1;j<=jmax;j++)
      {
      if(!(flag[i][j]&C_F))
        {
        (*ibound)++;
        if(flag[i-1][j]&C_F)
          flag[i][j] |= B_W;
        if(flag[i+1][j]&C_F)
          flag[i][j] |= B_E;
        if(flag[i][j-1]&C_F)
          flag[i][j] |= B_S;
        if(flag[i][j+1]&C_F)
          flag[i][j] |= B_N;
        }
      }
    }
  }

/* Computation of stream function and vorticity */
/* Computation of the stream function at the upper right corner */
/* of cell (i,j) (only if bother lower cells are fluid cells) */
void CalcZeta()
  {
  int i,j;
  /* Computation of the vorticity zeta at the upper right corner */
  /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
  for(i=1;i<=imax-1;i++)
    {
    for(j=1;j<=jmax-1;j++)
      {
      if((flag[i][j]&C_F)&&(flag[i+1][j]&C_F)&&
       (flag[i][j+1]&C_F)&&(flag[i+1][j+1]&C_F))
        zeta[i][j]=(u[i][j+1]-u[i][j])/dely-(v[i+1][j]-v[i][j])/delx;
      else
        zeta[i][j]=0.;
      }
    }
  }

void CalcPsi()
  {
  int i,j;
/* Computation of the stream function at the upper right corner */
/* of cell (i,j) (only if bother lower cells are fluid cells) */
  for(i=0;i<=imax;i++)
    {
    psi[i][0]=0.0;
    for(j=1;j<=jmax;j++)
      {
      psi[i][j]=psi[i][j-1];
      if((flag[i][j]&C_F)||(flag[i+1][j]&C_F))
        psi[i][j]+=u[i][j]*dely;
      }
    }
  }

void EmbedText(int col,int row,char*string,BYTE fore,BYTE back,int transparent,BYTE*bits,int wide)
  {
  BYTE*font,mask;
  int h,i,j,w;
  while(*string)
    {
    i=Font[0]*(DWORD)Font[1]*(DWORD)(*string++);
    font=Font+2+i/8;
    mask=0x80>>(i%8);
    for(h=0;h<Font[0];h++)
      {
      j=wide*(row-h)+col;
      for(w=0;w<Font[1];w++,j++)
        {
        if(mask&*font)
          bits[j]=fore;
        else if(!transparent)
          bits[j]=back;
        mask>>=1;
        if(mask==0)
          {
          mask=0x80;
          font++;
          }
        }
      }
    col+=Font[1];
    }
  }

void FillBitmap(double xlength,double ylength,BYTE*bits,int wide,int outmode)
  {
  int i,j,k;
  double pmax=-1E10,pmin=1E10;
  if(outmode==0)
    CalcZeta(imax,jmax,xlength/imax,ylength/jmax);
  else if(outmode==1)
    CalcPsi(imax,jmax,xlength/imax,ylength/jmax);
  else if(outmode==2)
    {
    for(j=1;j<jmax+1;j++)
      {
      for(i=1;i<imax+1;i++)
        {
        if(flag[i][j]&C_F)
          {
          pmax=max(pmax,p[i][j]);
          pmin=min(pmin,p[i][j]);
          }
        }
      }
    }
  else
    {
    printf("outmode=%i?\n",outmode);
    exit(1);
    }
  for(j=1;j<jmax+1;j++)
    {
    for(i=1;i<imax+1;i++)
      {
      if(!(flag[i][j]&C_F))
        k=0;
      else
        {
        if(outmode==0)
          {
          double z=(i<imax&&j<jmax)?zeta[i][j]:0.;
          k=2+(int)(pow(fabs(z/12.6),0.4)*124.);
          }
        else if(outmode==1)
          {
          double p=(i<imax&&j<jmax)?psi[i][j]:0.;
          k=2+(int)((p+3.0)/7.5*125.);
          }
        else if(outmode==2)
          k=2+(int)((p[i][j]-pmin)/(pmax-pmin)*125.);
        if(k<2)
          k=2;
        else if(k>126)
          k=126;
        }
      bits[wide*(j-1)+i-1]=k;
      }
    }
#ifdef MINIVAN
  i=124;
  j=48;
#elif defined(SUPERGIRL)
  i=116;
  j=66;
#else
  i=104;
  j=64;
#endif
  if(outmode==0)
    EmbedText(i,j,"vort",1,0,TRUE,bits,wide);
  else if(outmode==1)
    EmbedText(i,j,"strm",1,0,TRUE,bits,wide);
  else if(outmode==2)
    EmbedText(i,j,"pres",1,0,TRUE,bits,wide);
  }

/* Computation of tentative velocity field (f,g) */
void ComputeTentativeVelocity()
  {
  int i,j;
  double du2dx,duvdy,duvdx,dv2dy,laplu,laplv;
  for(i=1;i<=imax-1;i++)
    {
    for(j=1;j<=jmax;j++)
      {
      /* only if both adjacent cells are fluid cells */
      if((flag[i][j]&C_F)&&(flag[i+1][j]&C_F))
        {
        du2dx=((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])
          +gamma*fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])
          -(u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j])
          -gamma*fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j]))
          /(4.0*delx);
        duvdy=((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])
          +gamma*fabs(v[i][j]+v[i+1][j])*(u[i][j]-u[i][j+1])
          -(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j])
          -gamma*fabs(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]-u[i][j]))
          /(4.0*dely);
        laplu=(u[i+1][j]-2.0*u[i][j]+u[i-1][j])/delx/delx
          +(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dely/dely;
        f[i][j]=u[i][j]+delt*(laplu/Re-du2dx-duvdy);
        }
      else
        f[i][j]=u[i][j];
      }
    }
  for(i=1;i<=imax;i++)
    {
    for(j=1;j<=jmax-1;j++)
      {
      /* only if both adjacent cells are fluid cells */
      if((flag[i][j]&C_F)&&(flag[i][j+1]&C_F))
        {
        duvdx=((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])
          +gamma*fabs(u[i][j]+u[i][j+1])*(v[i][j]-v[i+1][j])
          -(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j])
          -gamma*fabs(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]-v[i][j]))
          /(4.0*delx);
        dv2dy=((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])
          +gamma*fabs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])
          -(v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j])
          -gamma*fabs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j]))
          /(4.0*dely);
        laplv=(v[i+1][j]-2.0*v[i][j]+v[i-1][j])/delx/delx
          +(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dely/dely;
        g[i][j]=v[i][j]+delt*(laplv/Re-duvdx-dv2dy);
        }
      else
        g[i][j]=v[i][j];
      }
    }
  /* f&g at external boundaries */
  for(j=1;j<=jmax;j++)
    {
    f[0][j]=u[0][j];
    f[imax][j]=u[imax][j];
    }
  for(i=1;i<=imax;i++)
    {
    g[i][0]=v[i][0];
    g[i][jmax]=v[i][jmax];
    }
  }

/* Calculate the right hand side of the pressure equation */
void ComputeRHS()
  {
  int i,j;
  for(i=1;i<=imax;i++)
    for(j=1;j<=jmax;j++)
      if(flag[i][j]&C_F)
        rhs[i][j]=((f[i][j]-f[i-1][j])/delx+(g[i][j]-g[i][j-1])/dely)/delt;
  }

/* Red/Black SOR to solve the poisson equation */
int Poisson(double eps,int itermax,double omega,double*res,int ifull)
  {
  int i,j,iter;
  double add,beta_2,beta_mod;
  double p0=0.;
  int rb;/* Red-black value. */
  double rdx2=1.0/(delx*delx);
  double rdy2=1.0/(dely*dely);
  beta_2=-omega/(2.0*(rdx2+rdy2));
  /* Calculate sum of squares */
  for(i=1;i<=imax;i++)
    for(j=1;j<=jmax;j++)
      if(flag[i][j]&C_F)
        p0+=p[i][j]*p[i][j];
  p0=sqrt(p0/ifull);
  if(p0<0.0001)
    p0=1.;
  /* Red/Black SOR-iteration */
  for(iter=0;iter<itermax;iter++)
    {
    for(rb=0;rb<=1;rb++)
      {
      for(i=1;i<=imax;i++)
        {
        for(j=1;j<=jmax;j++)
          {
          if((i+j)%2!=rb)
            continue;
          if(flag[i][j]==(C_F|B_NSEW))
            {
            /* five point star for interior fluid cells */
            p[i][j]=(1.-omega)*p[i][j]-beta_2*((p[i+1][j]+p[i-1][j])*rdx2
              +(p[i][j+1]+p[i][j-1])*rdy2-rhs[i][j]);
            }
          else if(flag[i][j]&C_F)
            {
            /* modified star near boundary */
            beta_mod=-omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
            p[i][j]=(1.-omega)*p[i][j]-beta_mod*
              ((eps_E*p[i+1][j]+eps_W*p[i-1][j])*rdx2
              +(eps_N*p[i][j+1]+eps_S*p[i][j-1])*rdy2-rhs[i][j]);
            }
          }
        }
      }
    /* Partial computation of residual */
    *res=0.;
    for(i=1;i<=imax;i++)
      {
      for(j=1;j<=jmax;j++)
        {
        if(flag[i][j]&C_F)
          {
          /* only fluid cells */
          add=(eps_E*(p[i+1][j]-p[i][j])-eps_W*(p[i][j]-p[i-1][j]))*rdx2
            +(eps_N*(p[i][j+1]-p[i][j])-eps_S*(p[i][j]-p[i][j-1]))*rdy2
            -rhs[i][j];
          *res+=add*add;
          }
        }
      }
    *res=sqrt((*res)/ifull)/p0;
    /* convergence? */
    if(*res<eps)
      break;
    }
  return(iter);
  }

/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix */
void UpdateVelocity()
  {
  int i,j;
  for(i=1;i<=imax-1;i++)
    for(j=1;j<=jmax;j++)
      /* only if both adjacent cells are fluid cells */
      if((flag[i][j]&C_F)&&(flag[i+1][j]&C_F))
        u[i][j]=f[i][j]-(p[i+1][j]-p[i][j])*delt/delx;
  for(i=1;i<=imax;i++)
    for(j=1;j<=jmax-1;j++)
      /* only if both adjacent cells are fluid cells */
      if((flag[i][j]&C_F)&&(flag[i][j+1]&C_F))
        v[i][j]=g[i][j]-(p[i][j+1]-p[i][j])*delt/dely;
  }

/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable. */
void SetTimeStepInterval()
  {
  int i,j;
  double umax,vmax,deltu,deltv,deltRe;
  /* delt satisfying CFL conditions */
  if(tau>=1E-10)
    {  /* else no time stepsize control */
    umax=1E-10;
    vmax=1E-10;
    for(i=0;i<=imax+1;i++)
      for(j=1;j<=jmax+1;j++)
        umax=max(fabs(u[i][j]),umax);
    for(i=1;i<=imax+1;i++)
      for(j=0;j<=jmax+1;j++)
        vmax=max(fabs(v[i][j]),vmax);
    deltu=delx/umax;
    deltv=dely/vmax;
    deltRe=1/(1/(delx*delx)+1/(dely*dely))*Re/2.;
    if(deltu<deltv)
      delt=min(deltu,deltRe);
    else
      delt=min(deltv,deltRe);
    delt=tau*delt; /* multiply by safety factor */
    }
  }

int main(int argc,char**argv,char**envp)
  {
  int outmode;
  int output_frequency=40;
  int itermax=100;   /* Maximum number of iterations in SOR */
  int verbose=1;     /* Verbosity level */
  unsigned checker=0;
  double checker1=0.;
  int i,j,itersor=0,ifluid=0,ibound=0;
  double res,t;
  int iters=0,wide;
  int show_help=0,show_usage=0,show_version=0;
  BITMAPFILEHEADER bf;
  BITMAPINFOHEADER bm;
  BYTE*bits;
  memset(&bm,0,sizeof(bm));
  memset(&bf,0,sizeof(bf));
  bm.biSize=sizeof(BITMAPINFOHEADER);
  bm.biPlanes=1;
  bm.biBitCount=8;
  bm.biWidth=imax;
  bm.biHeight=jmax*3+2;
  bm.biClrUsed=sizeof(palette)/sizeof(DWORD);
  wide=4*((bm.biWidth*bm.biBitCount+31)/32);
  bm.biSizeImage=wide*bm.biHeight;
  bits=calloc(bm.biSizeImage,sizeof(BYTE));
  bf.bfType=0x4D42;
  bf.bfOffBits=sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+sizeof(palette);
  bf.bfSize=bf.bfOffBits+bm.biSizeImage;
  printf("NAST2D by Griebel, Dornseifer, and Neunhoeffer\n");
  printf("modified by Orchard and also Benton\n");
  if(argc>1)
    output_frequency=atoi(argv[1]);
  output_frequency=max(1,min(1000,output_frequency));
  delx=xlength/imax;
  dely=ylength/jmax;
  for(i=0;i<=imax+1;i++)
    {
    for(j=0;j<=jmax+1;j++)
      {
      checker+=(i*jmax)+j+1;
      checker1+=(i*jmax)+j+1.;
      u[i][j]=ui;
      v[i][j]=vi;
      p[i][j]=0.;
      }
    }
  InitFlag(&ibound);
  ApplyBoundaryConditions(ui,vi);
  for(t=0.;t<t_end;t+=delt,iters++)
    {
    SetTimeStepInterval();
    ifluid=(imax*jmax)-ibound;
    ComputeTentativeVelocity();
    ComputeRHS(delt);
    if(ifluid>0)
      itersor=Poisson(eps,itermax,omega,&res,ifluid);
    else
      itersor=0;
    printf("%d t:%lG,delt:%lG,SOR iters:%3d,res:%lG,bcells:%d\n",iters,t+delt,delt,itersor,res,ibound);
    UpdateVelocity();
    ApplyBoundaryConditions(ui,vi);
    if(iters%output_frequency==0)
      {
      char outfile[13];
      FILE*fp;
      for(outmode=0;outmode<3;outmode++)
        FillBitmap(xlength,ylength,bits+((jmax+1)*wide)*outmode,wide,outmode);
      sprintf(outfile,"%06d.BMP",iters/output_frequency);
      printf("output file: %s\n",outfile);
      if((fp=fopen(outfile,"wb"))==NULL)
        {
        printf("can't create file\n");
        exit(1);
        }
      fwrite(&bf,sizeof(bf),1,fp);
      fwrite(&bm,sizeof(bm),1,fp);
      fwrite(palette,sizeof(palette)/sizeof(DWORD),sizeof(DWORD),fp);
      fwrite(bits,1,bm.biSizeImage,fp);
      fclose(fp);
      }
    if(_kbhit())
      break;
    }
  free(bits);
  return(0);
  }
