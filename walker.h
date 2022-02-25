#ifndef WALKER
#define WALKER

typedef double angle;
typedef double real;

typedef struct {
    real x;
    real y;
} point;

typedef real (*landscape)(point);

real distance(point, point);
real area_inside(landscape, point, real, real);

#endif
