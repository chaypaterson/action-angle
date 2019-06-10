#include <stdio.h>
#include <math.h>
#include <complex.h>

/*		Action Angle Calculator
 *   This program calculates periods of oscillation from a Hamiltonian
 *   and two formulas:
 *   	T = dJ/dE
 *   	J = i*dp*dx (prev. used Pick's theorem but it was inappropriate)
 *   	    here, i is the number of squares inside an implicit surface
 * */

double T(double p) 
{	// kinetic energy
	double T = p*p/2;	// unit mass m=1
	return T;
}

double V(double x)
{	// this is the potential
	double V;
	V = x*x/2;	// harmonic oscillator, k=1
//	V = x*x*x*x/4 - x*x/2; // double dip well
//	V = 1-1/(1+x*x/2);	// lorentzian
	return V;
}

double H(double x, double p)
{
	double E = T(p)+V(x);
	return E;
}

double Pick(int i, int b)
{
    // Area by Pick's theorem
	double Area;
	Area += (double)i;
	Area += 0.5*(double)b;
	Area = Area - 1;

	return Area;
}

int DumbSrc(double x, double p, double deltax, double deltap, double E)
{
    // find edges of surface defined by H(x,p) <= E
    // old linear search
    int m = 0;
	while (H(x,p)<E)
	{
		x += deltax;
		p += deltap;
		m++;
	}

    return m;
}

int BinSrc(double x, double p, double deltax, double deltap, double E)
{
    // find edges of surface defined by H(x,p) <= E
	// NEW: binary search
	int m = 0;
	int n = 0;

    // Binary search:
	while(H(x+m*deltax,p+m*deltap)<E)
	{
		n = 0;
		while(H(x+m*deltax+deltax*(2<<(n+1)),p+m*deltap+deltap*(2<<(n+1)))<E)
		{
			n++;
		}
		m += 2<<n;
	}

	return m;
}

double J(double E, double dx, double dp)
{	// Calculate J for a given E.
	// TODO This algorithm is dumb! 
	// Find a more efficient way to calculate areas.
	double dJ = dx * dp; // area unit

	int i =0;
	// We will start at the origin:
	// Find the edges on the x-axis:
	// NEW: binary search
	double Apos = dx*(BinSrc(0,0,dx,0,E));
	double Aneg = -dx*BinSrc(0,0,-dx,0,E);

	// Now we scan off the X axis to find other edges:
	double x;	
	for (x=Aneg; x<Apos; x+=dx)
	{	// scanning between Aneg<x<Apos:
		i += BinSrc(x,0,0,dp,E); // up from x axis
		i += BinSrc(x,-dp,0,-dp,E); // down from below x
	}

	return (double)i*dJ; // cruder area
}

int main()
{	// Constant increments:
	double dx = 0.0001;
	double dp = 0.0001;
    double errT = 1E-6; // absolute error in period
    // accuracy in errT = dx*dp/dE, so set dE to:
	double dE = dp*dx/errT;

	double Emin = 0+dE;
	double Emax = 1;
	double E;

	for (E = Emin; E<Emax; E+=dE)
	{
		double dJ = J(E+dE,dx,dp)-J(E,dx,dp);
		double T = dJ/dE; // derivative dJ/dE
        // return energy level and period:
		printf("%lf, %lf\n",E,T);
	}

	return 0;
}
