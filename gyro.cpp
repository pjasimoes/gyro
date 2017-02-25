#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "bessel.h"
#include "gauleg.h"
#include "gyro.h"

using namespace std;

/*
  program: 
  gyro

  purpose:
  calculation of gyrosynchrotron emission and self-absorption coefficients,
  based on the full equations of Ramaty, ApJ, 1969.
  See gyro.pro (IDL) for input/output information.

  author: 
  Paulo Jose de Aguiar Simoes (PJAS)

  history:
  Ramaty FORTRAN, 199? (isotropic)
  TY     IDL, 199? (isotropic)
  PJAS   IDL, 2003-2009 (anisotropy terms, Gauss-Legendre integration, vectorization)
         several inputs from JER Costa and CG Gimenez de Castro
  PJAS   C++, 2011 (Wild-Hill Bessel approx., OpenMP)
*/

/* global variables */
int m;
int ned, nfreq, npd;
double cs,ss,ffp;
double betaMax;
double bnor=0.0,p0,p1;

/* prototypes */
/******************************************************************************/
void refr(double &an1, double &an2, double &ath1, double &ath2, double &ffb);

void getAnor(double delta[], double anor[], double energy[], int ned);

void gsycore(double g[2], double an, double ath, double &ffb, const double delta[]
	     ,const double anor[],const double energy[], const int ned, double phi[]
	     ,double gphi[], double dgdphi[], int npd);

double linterpol(double x, double xm[], double ym[], int n);

void deriv(double x[], double y[], int n, double d[]);






/* main program */
/******************************************************************************/
int main(int argc, char *argv[]) 
{
  bool output = 1;
  bool info   = 0;
 
  string outputFile="gyro.out";
  string inputFile ="gyro.in";

  int optind=1;
  /* decode arguments */
  while ((optind < argc) && (argv[optind][0]=='-')) 
    {
      string sw = argv[optind];
     
      if (sw=="-o") { /* use switch '-o filename' to name the output file */
	optind++;
	outputFile = argv[optind];
      }
      if (sw=="-i") { /* use switch '-i filename' to name the input file */
	optind++;
	inputFile = argv[optind];
      }
    }

  /* handling input file: */
  ifstream fin;
  fin.open(inputFile.c_str());
  if (!fin.is_open())
    {
      cerr << "ERROR: Could not open " << inputFile << " for reading." << endl;
      fin.clear();
      fin.close();
      exit(EXIT_FAILURE);
    }

  /* handling output file: */
  ofstream fout;
  fout.open(outputFile.c_str());
   if (!fout.is_open())
    {
      cerr << "ERROR: Could not open " << inputFile << " for writing." <<endl;
      fout.clear();
      fout.close();
      exit(EXIT_FAILURE);
    }


   /* scalar variables */
  double vp,vb,ffb;               // plasma freq., gyro-freq., gyro-harmonic
  double bmag,np,nel,cte_j,cte_k; // magnetic field, thermal and non-thermal density, ctes.
  double angle;//  view-angle
  double an1,an2,ath1,ath2;       // refraction and polarization index

  /* reading input */
  fin >> bmag;  // magnetic field [Gauss]
  fin >> angle; // viewing angle [degree]
  fin >> np;    // plasma density [cm-3]
  fin >> nel;   // non-thermal electron density [cm-3]

  fin >> ned;   // energy grid size
  fin >> m;     // type of anisotropy
  fin >> npd;   // grid size for pitch-angle  (m=2 or else) 
  fin >> nfreq; // frequency grid size
  /*
    m defines the type of anisotropy:
    m = 0 -> isotropic
    m = 1 -> anisotropic, gaussian distribution
    m = 2 -> input array distribution
    else  -> f(E,u)
  */

  /* prepare to read or create freq array */
  int createFreq=0;

  if (nfreq == 0) 
    {
      createFreq = 1;
      nfreq = 111;
    }
  
  double * freq = new double [nfreq];

  if (createFreq == 1) 
    {   
      /* setting mw frequency array (if not passed as an input) */
      float fi;
      for (int k=0; k < nfreq; k++)
	{
	  if (k <= 100) fi=pow(10.0,(0.01*(k-01)));
	  if (k >  100) fi=pow(10.0,(0.10*(k-91)));
	  freq[k]=0.5*E/(PI*C*M0)*bmag*fi;
	}
    }
  else 
    {
      for (int i=0; i < nfreq; i++)
	fin >> freq[i];
    }

  /* allocate arrays for energy distribution */
  double * energy = new double [ned];
  double * delta  = new double [ned-1];
  double * anor   = new double [ned-1];

  /* setting energy array (input) [keV] */
  for (int i=0; i < ned; i++)
      fin >> energy[i];
  
  /* setting delta array (input) */
  for (int i=0; i < ned-1; i++)
    fin >> delta[i];
  

  int parray;
  if (m < 2) 
    parray=npd;
  else
    parray=ned*npd;
  
  /* allocate arrays for pitchangle distribution */
  double * phi     = new double [npd];
  double * gphi    = new double [parray];
  double * dgdphi  = new double [parray];

  switch (m)
    {
    case 0:  /* isotropic case */
      bnor = 1.0/4.0/PI;
      break;
    case 1:
      {  /* anisotropic: gaussian distribution */
	fin >> p0;    // center (m=1)
	fin >> p1;    // width  (m=1)
	int nd = ND;
	double wb[nd];
	double t_phi[nd];
	double t_gphi[nd],pp0,sum=0;
	gauleg(0.0,PI,t_phi,wb,nd);
	for (int i=0; i < nd; i++)
	  {
	    pp0=(t_phi[i]-p0);
	    t_gphi[i] = exp(-0.5*(pp0*pp0)/(p1*p1));
	    sum += t_gphi[i]*wb[i];
	  }
	bnor = 1.0/(sum*4.0);
	break;
      }
    case 2:
      { 
	/* anisotropic case: input distribution array.
	   input gphi is interpolated with gauleg to find bnor,
	   derivative is calculated */
	// read phi array
	for (int i=0;i<npd;i++)
	  fin >> phi[i];
	// read gphi array
	for (int i=0;i<npd;i++)
	  fin >> gphi[i];

	int nd = ND; 
	double wb[nd];
	double t_phi[nd];
	double t_gphi[nd],sum=0.0;
	gauleg(0.0,PI,t_phi,wb,nd);
	for (int i=0;i<nd;i++)
	  {
	    t_gphi[i] = linterpol(t_phi[i],phi,gphi,nd);
	    sum+=t_gphi[i]*wb[i];
	  }
	bnor = 1.0/(sum*4.0);
	deriv(phi,gphi,npd,dgdphi);
	break;
      }
    default:
      {
	/* here we will handle a f(E,u) distribution with 
	   arbitrary grid size   */
	int nd = ND, index; 
	double wb[nd];
	double t_phi[nd];
	double t_gphi[nd],sum=0.0;
	//	double * tempgphi = new double [npd];
	//double * tempdgdphi = new double [npd];
	double tempgphi[200];
	double tempdgdphi[200];
	
	gauleg(0.0,PI,t_phi,wb,nd);
	
	/* read phi array */
	for (int i=0;i<npd;i++)
	  fin >> phi[i];
	  	
	for (int j=0;j<ned;j++)
	  {
	    /* read gphi array for each e */
	    for (int i=0;i<npd;i++)
	      fin >> tempgphi[i];
	
	    /* then normalize it */
	    for (int i=0;i<nd;i++)
	      {
		t_gphi[i] = linterpol(t_phi[i],phi,tempgphi,nd);
		sum+=t_gphi[i]*wb[i];
	      }
	    bnor = 1.0/(sum*4.0);
	    /* finally, do the derivative */
	    deriv(phi,tempgphi,npd,tempdgdphi);
	    /* and put into larger array */
	    for (int k=0;k<npd;k++)
	      {
		index=j*npd+k;
		gphi[index]=tempgphi[k]*bnor;
		dgdphi[index]=tempdgdphi[k]*bnor;
	      }
	    sum=0.0; /* clear sum */
	  }
	bnor=1.0;
	//	delete [] tempgphi;
	//delete [] tempdgdphi;
	break;
      }
    }
 
  // closing input file
  fin.clear();
  fin.close();

  // allocate arrays for gyro-coefficients (output)
  double * jo   = new double [nfreq];
  double * jx   = new double [nfreq];
  double * ko   = new double [nfreq];
  double * kx   = new double [nfreq];
  
  /* initialize some constants */
  cs = abs(cos(angle*DTOR));            // cosine of viewing angle
  //cs = (cos(angle*DTOR));            // cosine of viewing angle
  ss = sqrt(1.0-cs*cs);                 // sine
  // signal_cs = cos(angle*DTOR)/cs;       // sign of cos(angle)
  vb = 0.5*E/(PI*C*M0)*bmag;            // gyrofrequency
  vp = E*sqrt(np/PI/M0);                // plasma frequency
  ffp = vp / vb;                        // ratio plasma/gyro frequency
  // alpha = 1.5 * vb / vp;                // Razin parameter
  //vr = 2.0/3.0*vp*vp/vb/ss;             // Razin effect cutoff frequency

  /* normalizaton of energy distribution */
  getAnor(delta, anor, energy, ned);

  betaMax = sqrt(energy[ned-1]*energy[ned-1]-1.0)/energy[ned-1];

  /* cte = e^3d / m0 / c^2d * BMag * Nel */
  cte_j = pow(E,3.0) / M0 / (C*C) * bmag * nel;
  cte_k = 4.0*PI*PI * E / bmag * nel;
  
  double g1[2];
  double g2[2];

  /* openMP thread number */
  #ifdef _OPENMP
    int tid[nfreq];
  #endif

  /* OpenMP directives */
#pragma omp parallel shared(jo,ko,jx,kx,tid) private(g1,g2,an1,an2,ath1,ath2,ffb)
#pragma omp for schedule(dynamic)
  for (int i=0; i < nfreq; i++)
    {

      /* Obtain thread number */
      #ifdef _OPENMP
            tid[i] = omp_get_thread_num();
      #endif

      /* gyro-harmonic */
      ffb = freq[i] / vb;

      /* refraction index */
      refr(an1,an2,ath1,ath2,ffb);

      /* ordinary mode */
      if ((ffb > ffp) && (an1 > 0.0))
	{
	  gsycore(g1,an1,ath1,ffb,delta,anor,energy,ned,phi,gphi,dgdphi,npd);
	  jo[i] = g1[0] * cte_j;
	  ko[i] = g1[1] * cte_k;
	} 
      else 
	{
	  jo[i] = 0.0;
	  ko[i] = 0.0;
	}
     
      /* extraordinary mode */
      if ((ffb > (sqrt(ffp*ffp + 0.25) + 0.5)) && (an2 > 0.0))
	{
	  gsycore(g2,an2,ath2,ffb,delta,anor,energy,ned,phi,gphi,dgdphi,npd);
	  jx[i] = g2[0] * cte_j;
	  kx[i] = g2[1] * cte_k;
	}
      else
	{
	  jx[i] = 0.0;
	  kx[i] = 0.0;
	}
    }

  /* output to file (ascii) */
  if (output) 
    {
      fout << nfreq << endl;
      for (int i=0; i < nfreq; i++)
	{

	  fout << freq[i] << "\t"
	       << jo[i]   << "\t" 
	       << jx[i]   << "\t" 
	       << ko[i]   << "\t" 
	       << kx[i]   << endl; 
	}
    }
  else 
    if (info) 
      {
	{
	  cout << "--------------------------------" << endl;
	  cout << "bmag=" << bmag << endl;
	  cout << "angle=" << angle << endl;
	  cout << "np="   << np << endl;
	  cout << "nel=" << nel << endl;
	  cout << "ned=" << ned << endl;
	  cout << "m=" << m << endl;
	  cout << "p0=" << p0 << endl;
	  cout << "p1=" << p1 << endl;
	  cout << "npd=" << npd << endl;
	  cout << "cs=" << cs << endl;
	  cout << "ss=" << ss << endl;
	  cout << "vb=" << vb << endl;  
	  cout << "vp=" << vp << endl;
	  cout << "ffp=" << ffp << endl; 

	  cout << "--------------------------------" << endl;
	  for (int i=0; i<ned; i++) 
	    cout << "energy=" << (energy[i]-1.0)*511. << endl;
	  cout << "--------------------------------" << endl;
	  for (int i=0; i<ned-1; i++) 
	    cout << "delta=" << delta[i] << endl;
	  cout << "--------------------------------" << endl;
	  cout << "nfreq=" << nfreq << endl;
	  //  for (int i=0; i < nfreq-1; i++)
	  //  cout << "freq=" << freq[i] << endl;	  
	  cout << "--------------------------------" << endl;
	  /*	  cout << scientific;
	  for (int i=0; i < ned-1; i++)
	    cout << "anor=" << anor[i] << endl;
	  cout << "--------------------------------" << endl;
	  */
	  /* testing OpenMP */
	  /*#ifdef _OPENMP
	  for (int i=0;i<nfreq;i++)
	    cout << "i="<<i << "  tid[i]="<< tid[i]<< " ffb[i]="<< ar_ffb[i]*vb/1e9 
		 << "  an1[i]  " << ar_an1[i] 
		 << "  an2[i]  " << ar_an2[i] 
		 << "  ath1[i] " << ar_ath1[i] 		 
		 << "  ath2[i] " << ar_ath2[i] 
		 << endl;
#endif 
	  */
	}
      }

  /* closing output file */
  fout.clear();
  fout.close();

  /* free allocated memory */
  delete [] phi;
  delete [] gphi;
  delete [] dgdphi;
  delete [] energy;
  delete [] delta;
  delete [] anor;
  delete [] jo;
  delete [] ko;
  delete [] jx;
  delete [] kx;
  delete [] freq;

  /*   
#ifdef _OPENMP
  for (int i=0;i<nfreq;i++)
    cout << i << "\t" << tid[i]<<endl;
#endif
  */

return 0;

}

/*
*******************************************************************************
*/

void gsycore(double g[2], double an, double ath, double &ffb, const double delta[], 
	     const double anor[], const double energy[], const int ned, double phi[], 
	     double gphi[], double dgdphi[], int npd)
{
  /* core function to calculate the coefficients */ 

  int j; 
  int s,s1,s2;
  double tmp1,tmp2,tmp3,b1;
  double ginf,gsup,gammai;
  double gamm[NPTS]; 
  double weight[NPTS];
  double gx[NPTS];
  double hx[NPTS];
  double Jb[2];
  double beta,cphis,sphis,tphis,phis;
  double fe,zs,xs,b,bpr,h3a,h3b;
  double gphis=1.0;
  double dgdphis=0.0;
  double g12 = 0.0, h12 = 0.0;
  double gErr, hErr, old_g12, old_h12, pp0;
  bool gFlag = 0, hFlag = 0;  

  // double * tempgphi = new double [npd];  
  //double * tempdgdphi = new double [npd];
  double tempgphi[200];
  double tempdgdphi[200];
		  
  /* s harmonic limits */
  s2 = floor(ffb * energy[ned-1] * (1.0 +  an * betaMax * cs));
  s1 = floor(ffb * sqrt(1.0-an*an*cs*cs))+1;
  
// cout << "s2=" << s2 << endl;

  tmp3 = (1.0-an*an*cs*cs);

  for (s = s1; s < s2; s++)
    {
      
      /* Energy limits to integrate in each s */
      tmp1 = s/ffb;
      tmp2 = an*abs(cs)*sqrt(tmp1*tmp1+an*an*cs*cs-1.0);
      ginf = (tmp1-tmp2)/tmp3;
      if (ginf < energy[0]) ginf = energy[0];
      gsup = (tmp1+tmp2)/tmp3; 
      if (gsup > energy[ned-1]) gsup = energy[ned-1];

      if (ginf >= gsup) continue;

      gauleg(ginf, gsup, gamm, weight, NPTS);
     
      /* loop thru gamma points */
      for (int i=0; i < NPTS; i++) 
	{
	  gammai = gamm[i];

	  beta = sqrt(gammai*gammai-1.0)/gammai;
	  cphis = (1.0-s/ffb/gammai)/beta/cs/an;
	  sphis = sqrt(1.0-cphis*cphis);
	  phis = acos(cphis);  /* phis: pitchangle (!pi-0) [rad] */
	  xs = s*an*beta*ss*sphis/(1.0-an*beta*cs*cphis);

	  /* where function (IDL): j is the index of energy array for each loop */
	  for (j=0;j<ned-2;j++)    
	    if ((gammai > energy[j]) && (gammai <= energy[j+1])) 
	      break;
	  
	  switch (m) /* anisotropy type */
	    {

	    case 0: /* isotropic */
	      gphis = 1.0;
	      break;

	    case 1: /* anisotropic gaussian */
	      pp0=phis-p0;
	      gphis = exp(-0.5*(pp0*pp0) / (p1*p1));
	      dgdphis = (p0-phis)/(p1*p1)*exp(-0.5*((p0-phis)*(p0-phis))/(p1*p1));
	      break;

	    case 2: /* array distribution: y = linterpol(x,x_arr[],y_arr[],n); */
	      gphis   = linterpol(phis,phi,gphi,npd);
	      dgdphis = linterpol(phis,phi,dgdphi,npd);
	      break;

	    default: /* f(E,u) (gphi and dgdphi are interpolated in phis) y = linterpol(x,xm,ym,n); */
	      /* get gphi and dgdphi from larger array for this j energy*/
	      for (int k=0;k<npd;k++)
		{
		  tempgphi[k]=gphi[j*npd+k];
		  tempdgdphi[k]=dgdphi[j*npd+k];
		}
	      gphis   = linterpol(phis,phi,tempgphi,npd);
	      dgdphis = linterpol(phis,phi,tempdgdphi,npd);
	      break;
	    }
	  
	  /* energy distribution function: power-law */
	  fe = anor[j] * pow(gammai-1.0,-delta[j]);
	  
	  if (ffb < BESSEL_LIMIT)
	    {
	      /* full Bessel function and derivative */
	      b  = bes(s,xs);
	      b1 = bes(s+1,xs);
	      bpr = -b1 + s / xs * b;      
	    } 
	  else 
	    {
	      /* Wild & Hill (1971) approximation. */
	      beselj_wh(xs,s,Jb);
	      b = Jb[0];
	      bpr = Jb[1];
	    }
	  
	  /* Function Zs (Klein,1984) 2o.term */
	  zs=-beta*sphis*bpr+ath*(cs/ss/an-beta*cphis/ss)*b;
	  gx[i] = bnor*zs*zs*gphis*(2.0*PI/cs/(1.0+ath*ath))*ffb/beta*fe;
	  hx[i] = 0.0;

	  if (hFlag == 0)
	    {
	      h3a = (delta[j]*gammai*(gammai+1.0)+2.0*gammai*gammai-1.0)
		/gammai/(gammai*gammai-1.0);
	      tphis=(an*beta*cs-cphis)/(gammai*beta*beta*sphis);
	      h3b = tphis / gphis * dgdphis;
	      h3b = (m != 0) ? h3b : 0.0;
	      hx[i] = (h3a+h3b)*gx[i]/ffb/ffb/an/an; /* Function H */
	    }
	} /* endfor (i) */

      old_g12 = g12;
      for (int i=0;i<NPTS;i++)
	g12+=gx[i]*weight[i];
      gErr = fabs((g12-old_g12)/old_g12);
      if (gErr < RERR) gFlag=1;
      if (hFlag == 0) 
	{
	  old_h12 = h12;
	  for (int i=0;i<NPTS;i++)
	    h12+=hx[i]*weight[i];
	  hErr = fabs((h12-old_h12)/old_h12);
	  if (hErr < RERR) hFlag=1;
	}

      if ((gFlag == 1) && (hFlag == 1)) break; 
    }
  
  g[0]=g12;
  g[1]=h12;

  // delete [] tempgphi;
  //delete [] tempdgdphi;
}

/*
*******************************************************************************
*/

void refr(double &an1, double &an2, double &ath1, double &ath2, double &ffb)

/*   calculation of indices of refraction (an) and polarization 
     coefficient (ath) of magnetoionic modes ordinary (1) and 
     extraordinary (2) (Appleton-Hartree equation) */
{
  
  double anum, dnum1, dnum2, aknum, dknum1, dknum2;
  double ar,br,cr;

  br = (ffp*ffp-ffb*ffb);  
  ar = sqrt(pow(ffb,4.0)*pow(ss,4.0)+4.0*ffb*ffb*pow(br,2.0)*cs*cs);
  cr = ffb*ffb*ss*ss;

  anum = 2.0*ffp*ffp*br;
  dnum1 = +ar - 2.0*ffb*ffb*br - cr;
  dnum2 = -ar - 2.0*ffb*ffb*br - cr;

  an1 = sqrt(1.0 + anum / dnum1);
  an2 = sqrt(1.0 + anum / dnum2);

  aknum = 2.0*ffb*br*cs;
  dknum1 = +ar-cr;
  dknum2 = -ar-cr;
  ath1 = -aknum / dknum1;
  ath2 = -aknum / dknum2;
}

/* 
*******************************************************************************
*/

void getAnor(double delta[], double anor[], double energy[], int ned)
{
  /* normalization for single or multi power-law
     original IDL routine by Guigue */

  double * ce = new double [ned-1];
  double * f  = new double [ned-1];
  double ce1,ce2,v;
  
  /* converting energy [keV] to Lorentz factor energy */
  for (int i=0; i < ned; i++) /* loop n_elements(energy) */
    energy[i] = energy[i]/E0+1.0;

  for (int i=0; i < (ned-1); i++) /* loop n_elements(delta) */
    {
      ce1 = pow(energy[i]  -1.0, 1.0-delta[i]);
      ce2 = pow(energy[i+1]-1.0, 1.0-delta[i]);
      ce[i]  = (ce1-ce2) / (delta[i]-1.0);
    }

  f[0] = 1.0;
  for (int i=0; i < (ned-2); i++) /* loop n_elements(delta)-1 */
      f[i+1] = pow(energy[i+1]-1.0, delta[i+1]-delta[i]);

  v = 0.0;
  for (int i=0; i < (ned-1); i++)
    {
      f[i] *= ce[i];
      v += f[i];
    }
  v = 1.0 / v;

  anor[0] = v;
  for (int i=0; i < (ned-2); i++)
    anor[i+1] = anor[i]*pow(energy[i+1]-1.0,delta[i+1]-delta[i]);

  delete [] f;
  delete [] ce;

}

/* 
*******************************************************************************
*/

void deriv(double x[], double y[], int n, double d[])
/* deriv: 3-point, Lagrangian interpolation. */
{

  double x12[n];
  double x01[n];
  double x02[n];

  for (int i=1; i<n-1; i++)
    {
      x12[i] = x[i] - x[i+1];   /* x1 - x2 */
      x01[i] = x[i-1] - x[i];   /* x0 - x1 */
      x02[i] = x[i-1] - x[i+1]; /* x0 - x2 */
      
      d[i] = y[i-1] *(x12[i] / (x01[i]*x02[i])) +
	y[i] * (1./x12[i] - 1./x01[i]) - 
	y[i+1]* (x01[i] / (x02[i] * x12[i]));
    }
 
  /* Formulae for the first and last points */
 
  d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - /* first point */
    y[1] * x02[1]/(x01[1]*x12[1]) + 
    y[2] * x01[1]/(x02[1]*x12[1]);

  int n2 = n-2;
  d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) + /* last point */
    y[n-2] * x02[n2]/(x01[n2]*x12[n2]) - 
    y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2]);
  
}

/*
*******************************************************************************
*/

double linterpol(double x, double xm[], double ym[], int n)
{
  /* linear interpolation */
  double r;

  for (int i=0; i < n; i++)
    {
      if ((x >= xm[i]) && (x <= xm[i+1]))
	{
	  r = ym[i]+(x-xm[i])*(ym[i+1]-ym[i])/(xm[i+1]-xm[i]); 
	  return r;
	}
    }
 
  r = ym[0];
  if (x >= xm[n-1])
    r = ym[n-1];

  return r;
}  
               
/* end of file */
/******************************************************************************/

/* NOTES FOR LATER

 
   FIND powerlaw in arbitrary f(E) distr.
   log_e=alog10(e)
   log_fe=alog10(fe)
   a=fltarr(n_elements(e)-1)
   b=fltarr(n_elements(e)-1)
   i=indgen(n_elements(e)-1)
   ;; a=delta, b=k f(E)=10.^k * E^(-delta)
   a[i] = (log_fe[1+i]-log_fe[i])/(log_e[1+i]-log_e[i])
   b[i] = log_fe[i]-a[i]*log_e[i]

*/
