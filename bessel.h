#include <cmath>

double bespr(int n, double x)
{
  double bes(int n, double x);
  int n1;
  double b,b1,bpr;
  n1 = n+1;
  b1 = bes(n1,x);
  b  = bes(n,x);
  bpr=-b1+n/x*b;
  return bpr;
}

///////////////////////////////////////////////////////////////////////////////

double bes(int n, double x)
{
  double b=0.0;
  double bessj0(double x);
  double bessj1(double x);
  double bessj(int n, double x);
 
  if (n == 0) b=bessj0(x);
  if (n == 1) b=bessj1(x);
  if (n >= 2) b=bessj(n,x);
  return b; 
}
///////////////////////////////////////////////////////////////////////////////

double bessj0(double x)
{
  double y,bessj0,z,xx,ax;

  double p[5] = {1.e0,-0.1098628627e-2,0.2734510407e-4,
      -0.2073370639e-5,0.2093887211e-6};
  double q[5] = {-0.1562499995e-1,
      0.1430488765e-3,-0.6911147651e-5,0.7621095161e-6,-0.934945152e-7};
  double r[6] = {57568490574.0e0,-13362590354.e0,651619640.7e0,
      -11214424.18e0,77392.33017e0,-184.9052456e0};
  double s[6] = {57568490411.e0,1029532985.0e0,
      9494680.718e0,59272.64853e0,267.8532712e0,1.0e0};

  if(fabs(x) < 8.0)
    {
      y=x*x;
      bessj0=(r[0]+y*(r[1]+y*(r[2]+y*(r[3]+y*(r[4]+y*r[5])))))
	/(s[0]+y*(s[1]+y*(s[2]+y*(s[3]+y*(s[4]+y*s[5])))));
    }
  else
    {
      ax=fabs(x);
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      bessj0=sqrt(0.636619772/ax)*(cos(xx)*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*p[4]))))
				  -z*sin(xx)*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*q[4])))));
    }
  return bessj0;
}

///////////////////////////////////////////////////////////////////////////////

double bessj1(double x)
{

  double y,z,ax,xx,bessj1;

  double r[6] = {72362614232.0e0,-7895059235.0e0,242396853.1e0,
       -2972611.439e0,15704.48260e0,-30.16036606e0};
  
  double s[6] = {144725228442.0e0,2300535178.0e0,
       18583304.74e0,99447.43394e0,376.9991397e0,1.0e0};
  
  double p[5] = {1.e0,0.183105e-2,-0.3516396496e-4,0.2457520174e-5,
      -0.240337019e-6};

  double q[5] = {.04687499995e0,-0.2002690873e-3,
      0.8449199096e-5,-0.88228987e-6,0.105787412e-6};

    if (fabs(x) < 8.0)
    {
      y=x*x;
      bessj1=x*(r[0]+y*(r[1]+y*(r[2]+y*(r[3]+y*(r[4]+y*r[5])))))
	/(s[0]+y*(s[1]+y*(s[2]+y*(s[3]+y*(s[4]+y*s[5])))));
    }
    else
      {
	ax=fabs(x);
	z=8.0/ax;
	y=z*z;
	xx=ax-2.356194491;
	bessj1=sqrt(0.636619772/ax)*(cos(xx)*(p[0]+y*(p[1]+y*(p[2]+y*(p[3]+y*p[4]))))
				     -z*sin(xx)*(q[0]+y*(q[1]+y*(q[2]+y*(q[3]+y*q[4])))))*(x/ax);//*sign(1.,x);
	  }
  return bessj1;
}

///////////////////////////////////////////////////////////////////////////////

double bessj(int n, double x)
{

  const int iacc=40;
  int m;
  const double bigno=1.0e10,bigni=1.0e-10;
  double tox,bjm,bj,bjp,bessj,jsum,sum;
    
    //if(n < 2) return // pause 'bad argument n in bessj'
  tox=2.0/x;
  if (x > float(n))
    {
      bjm=bessj0(x);
      bj=bessj1(x);
      for (int j=1;j < n; j++)
	{
	  bjp=j*tox*bj-bjm;
	  bjm=bj;
	  bj=bjp;
	}
      bessj=bj;
    } 
  else
    {
      m=2*((n+int(sqrt(float(iacc*n))))/2);
      bessj=0.;
      jsum=0;
      sum=0.;
      bjp=0.;
      bj=1.;
      for(int j=m;j >= 1;j--)
	{
	  bjm=j*tox*bj-bjp;
	  bjp=bj;
	  bj=bjm;
	  if(fabs(bj) > bigno)
	    {
	      bj=bj*bigni;
	      bjp=bjp*bigni;
	      bessj=bessj*bigni;
	      sum=sum*bigni;
	    }
	  if(jsum != 0)
	    sum=sum+bj;
	  jsum=1-jsum;
	  if(j == n)
	    bessj=bjp;
	}
      sum=2.0*sum-bj;
      bessj=bessj/sum;
    }
  return bessj;
}


void beselj_wh(double x, int s, double Jn[2])
{
  /*
    b = beselj_wh(x,n,db)
    Bessel function approximations by Wild and Hill, 1971.
    Use this function for gyrosynchrotron integration
    x = xs
    n = s
    note that the input x is actually s*x, while in the approximation
    x is x/s.
  */
  
#define A1 0.503297
#define B1 1.193000
#define PI 3.1415926535897931160
  
  double a,b,z;

  x = x/s;

  a = pow(pow(1.0-x*x,1.5) + A1/s,(1.0/6.0));
  b = pow(pow(1.0-x*x,1.5) + B1/s,(1.0/6.0)) * (1.0-1.0/(5.0*pow(s,2.0/3.0)));

  z = (x * exp(sqrt(1.0-x*x))) / (1.0+sqrt(1.0-x*x));

  Jn[0] = 1.0/sqrt(2.0*PI*s)*pow(z,s)/a;
  Jn[1] = a*b*Jn[0]/x;
 
}
