#include <math.h>

typedef double angle;
typedef double real;

typedef struct {
    real x;
    real y;
} point;

typedef real (*landscape)(point);

// the walker struct may as well be private (in this file only, invisible
// outside)
typedef struct {
    const point start;
    const real height; // walkers try to keep a constant height
    const real stride;
    const angle delta;
    point posn;
    angle dirn;
    real tot_area;
} walker;

walker walker_construct(landscape H, point start, real stride, real delta) {
    walker iter = {
        .start = start,
        .height = H(start),
        .stride = stride,
        .delta = delta,
        .posn = start,
        .dirn = 0,
        .tot_area = 0
    };

    return iter;
}

void revolve(point *edge, point *centre, real step, angle theta) {
    // rotate the point "edge" around "centre"
    edge->x = centre->x + step * cos(theta);
    edge->y = centre->y + step * sin(theta);
}

void walker_step(landscape H, walker* iter) {
    // make a note of initial position and angle for the area measurements
    point last = iter->posn;
    angle theta0 = iter->dirn;
    // first seek out a new position to step to
    point next = iter->posn;
    next.x += iter->stride * cos(iter->dirn);
    next.y += iter->stride * sin(iter->dirn);

    angle dq = iter->delta;

    for(int i=0; i<1; i++) {
        while(H(next) > iter->height) {
            iter->dirn += dq;
            revolve(&next, &last, iter->stride, iter->dirn);
        }
        dq *= 0.5;

        while (H(next) < iter->height) {
            iter->dirn -= dq;
            revolve(&next, &last, iter->stride, iter->dirn);
        }
        dq *= 0.5;
    }

    // now update the area

    iter->tot_area += 0.5 * (
                     +(last.x - iter->start.x) * (next.y - last.y)
                     -(last.y - iter->start.y) * (next.x - last.x));

    // add a small correction for error in height?
    iter->tot_area -= 0.5 * iter->delta * iter->stride * iter->stride;

    // and a very small correction for path curvature
    angle dtheta = iter->dirn - theta0; // curvature = dtheta / ds
    iter->tot_area += dtheta * (iter->stride) * (iter->stride) / 6 ;

    // finally, move to the next position

    iter->posn = next;
}

real distance(point p, point q) {
    return sqrt( (p.x - q.x) * (p.x - q.x)
                +(p.y - q.y) * (p.y - q.y));
}

real area_inside(landscape H, point start, real ds, real dq) {
    // a mountain goat tries to keep its feet level as it navigates a basin:
    walker goat = walker_construct(H, start, ds, dq);

    // take a couple of steps to get away from the start
    while (distance(goat.posn, goat.start) < ds) {
        walker_step(H, &goat);
    }
    // keep going until we loop back to the start
    while (distance(goat.posn, goat.start) > ds) {
        walker_step(H, &goat);
    }
    // we don't need to add the last little bit of area on because it is
    // guaranteed to be zero.

    return goat.tot_area;
}
