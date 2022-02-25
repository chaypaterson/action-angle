#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "walker.h"

/*        Action Angle Calculator
 *   This program calculates periods of oscillation from a Hamiltonian
 *   and two formulas:
 *       T = dJ/dE
 *       J = i*dp*dx (prev. used Pick's theorem but it was inappropriate)
 *           here, i is the number of squares inside an implicit surface
 * */

typedef double (*real_func)(double);
typedef landscape hamiltonian;

real T(real p) {
    // kinetic energy
    real T = p*p/2;    // unit mass m=1
    return T;
}

/* Different functions for harmonic, quartic, and lorentzian well:
 */

real harmonic_V(real x) {
    // harmonic oscillator, k=1
    return x*x/2;
}

real quartic(real x) {
    // quartic well with nontrivial dip
    return x*x*x*x/4 - x*x/2;
}

real lorentzian(real x) {
    // well with finite escape energy
    return 1-1/(1+x*x/2);
}

real H_harmonic(point q) {
    return T(q.y) + harmonic_V(q.x);
}

real J(hamiltonian H, point q, real deltax) {
    return area_inside(H, q, deltax, 0.01);
}

real diffx(hamiltonian H, point* q, real dx) {
    point qp = {
        .x = q->x + 0.5 * dx,
        .y = q->y
    };
    point qm = {
        .x = q->x - 0.5 * dx,
        .y = q->y
    };
    real dH = H(qp) - H(qm);
    return dH / dx;
}

int main() {
    // Constant increments:
    point dq;
    dq.x = 0.01;
    dq.y = 0.01;

    real dA = dq.x;

    // demonstrating that we can choose the hamiltonian programmatically at
    // runtime:
    hamiltonian H = H_harmonic;

    for (real Amp = 0; Amp < 4; Amp += dA) {
        point q = {
            .x = Amp,
            .y = 0
        };
        real E = H(q);
        real dx = 0.01 * dA;
        real dE = diffx(H, &q, dx);
        printf("%f \n", dE);
        // compute dJ:
        q.x += 0.5 * dx;
        real dJ = J(H, q, dA);
        printf("%f \n", dJ);
        q.x -= 1.0 * dx;
        dJ -= J(H, q, dA);
        printf("%f \n", dJ);
        dJ /= dx;
        printf("%f \n", dJ);

        // calculate T = dJ / dE
        real T = dJ / dE;
        // return energy level and period:
        printf("%lf, %lf\n",E,T);
        // TODO return amplitude and period
    }

    return 0;
}
