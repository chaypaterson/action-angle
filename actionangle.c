#include <stdio.h>
#include <math.h>
#include <complex.h>

/*		Action Angle Calculator
 *   This program calculates periods of oscillation from a Hamiltonian
 *   and two formulas:
 *   	T = dJ/dE
 *   	J = i*dp*dx (prev. used Pick's theorem but it was inappropriate)
 *   		here, i is the number of squares inside an implicit surface
 * */

typedef double (*real_func)(double);

typedef struct {
    double p;
    double x;
} point;

typedef double (*hamiltonian)(point);

double T(double p) {
    // kinetic energy
	double T = p*p/2;	// unit mass m=1
	return T;
}

/* Different functions for harmonic, quartic, and lorentzian well:
 */

double harmonic_V(double x) {
    // harmonic oscillator, k=1
    return x*x/2;
}

double quartic(double x) {
    // quartic well with nontrivial dip
    return x*x*x*x/4 - x*x/2;
}

double lorentzian(double x) {
    // well with finite escape energy
    return 1-1/(1+x*x/2);
}

/* double H(double x, double p) {
    // H should be an argument of type hamiltonian
	double E = T(p)+V(x);
	return E;
}
*/

double H_harmonic(point q) {
    return T(q.p) + harmonic_V(q.x);
}

double Pick(int i, int b) {
	// Area by Pick's theorem
	double Area;
	Area += (double)i;
	Area += 0.5*(double)b;
	Area = Area - 1;

	return Area;
}

int BinSrc(hamiltonian H, point q, point deltaq, double E) {
	// find edges of surface defined by H(x,p) <= E
	// NEW: binary search
	int m = 0;

    double p = q.p;
    double x = q.x;
    double deltap = deltaq.p;
    double deltax = deltaq.x;
    point edge;
    edge.p = q.p, edge.x = q.x;

	// Binary search: move out in a straight line in direction deltaq
	while(H(edge) < E) {
        point beyond = edge;
        int n;
        for (n = 0; H(beyond) < E; n++) {
            beyond.p += deltap * (2<<(n+1));
            beyond.x += deltax * (2<<(n+1));
        }
		m += 2<<n;
        edge.p += m * deltap;
        edge.x += m * deltax;
	}

	return m;
}

double J(hamiltonian H, double E, point deltaq) {
    // TODO the "walker" algorithm from my other project might be more efficient
    // here. combined with the shoelace equation/Green's theorem, I should be
    // able to calculate J in O(sqrt(J)) time. This would also avoid assuming
    // that the area to be computed is convex.
    // Calculate J for a given E.
    double dp = deltaq.p;
    double dx = deltaq.x;
	double dJ = dx * dp; // area unit

	int i = 0;
	// We will start at the origin:
	// Find the edges on the x-axis:
	// NEW: binary search
    point origin;
    origin.p = 0, origin.x = 0;
    point delta;
    delta.p = 0, delta.x = dx;
	double Apos = +dx * BinSrc(H, origin, delta, E);
    delta.x = -dx;
	double Aneg = -dx * BinSrc(H, origin, delta, E);

	// Now we scan off the X axis to find other edges:
	double x;	
	for (x = Aneg; x < Apos; x += dx) {	// scanning between Aneg<x<Apos:
        origin.x = x;
        delta.p = dp, delta.x = 0;
		i += BinSrc(H, origin, delta, E); // up from x axis
        delta.p = -dp, delta.x = 0;
		i += BinSrc(H, origin, delta, E); // down from below x
	}

	return dJ * i;
}

int main() {
    // Constant increments:
    point dq;
    dq.x = 0.01;
    dq.p = 0.01;

	double errT = 1E-3; // absolute error in period
	// accuracy in errT = dx*dp/dE, so set dE to:
	double dE = dq.p * dq.x / errT;

	double Emin = 0+dE;
	double Emax = 1;
	double E;

    hamiltonian H = H_harmonic;

	for (E = Emin; E<Emax; E+=dE) {
		double dJ = J(H, E + dE / 2,dq)-J(H, E - dE / 2,dq);
		double T = dJ / dE;
		// return energy level and period:
		printf("%lf, %lf\n",E,T);
		// TODO return amplitude and period
	}

	return 0;
}
