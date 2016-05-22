#include <stdio.h>
#include <math.h>

/*		Action Angle Calculator
 *   This program calculates periods of oscillation from a Hamiltonian
 *   and two formulas:
 *   	T = dJ/dE
 *   	J = i*dp*dx (prev. used Pick's theorem but it was inappropriate)
 * */

double T(double p) 
{	// kinetic energy
	double T = p*p/2;	// unit mass m=1
	return T;
}

double V(double x)
{	// this is the potential
	double V;
//	V = x*x/2;	// harmonic oscillator, k=1
//	V = x*x*x*x/4 - x*x/2; // double dip well
	V = 1-1/(1+x*x/2);	// lorentzian
	return V;
}

double H(double x, double p)
{
	double E = T(p)+V(x);
	return E;
}

double Pick(int i, int b)
{
	double Area;
	Area += (double)i;
	Area += 0.5*(double)b;
	Area = Area - 1;

	return Area;
}

double J(double E, double dx, double dp)
{	// Calculate J for a given E.
	// TODO This algorithm is dumb! 
	// Find a more efficient way to calculate areas.
	double dJ = dx * dp; // area unit

	// initialise counters:
	int i=0; // the origin is counted in loop

	// We will start at the origin:
	double x = 0;
	double p = 0;

	// Find the edges on the x-axis:
	while (H(x,p)<E)
	{	// first scan across x:
		i++;	// if the loop has not exited then we are inside
		x += dx;
	}

	double Apos = x; // note the +ve amplitude

	x = 0; // back to the origin

	while (H(x,p)<E)
	{	// now we sweep left
		i++;
		x -= dx;
	}

	double Aneg = x;

	// Now we scan off the X axis to find other edges:
	
	for (x=Aneg; x<Apos; x+=dx)
	{	// scanning between Aneg<x<Apos:
		p=dp; // start above the x axis.
		while (H(x,p)<E)
		{	// sweep UP p...
			i++;
			p += dp;
		}

		p=-dp; // start below the x axis.
		while (H(x,p)<E)
		{	// sweep DOWN p..
			i++;
			p -= dp;
		}

	}

	return (double)i*dJ; // cruder area
}

int main()
{	// Constant increments:
	double dx = 0.001;
	double dp = 0.001;
	double dE = 0.001;

	double Emin = 0+dE;
	double Emax = 1;
	double E;

	for (E = Emin; E<Emax; E+=dE)
	{
		double T = J(E+dE,dx,dp)-J(E,dx,dp);
		T = T/dE; // derivative dJ/dE
		printf("%lf, %lf\n",E,T);
	}

	return 0;
}
