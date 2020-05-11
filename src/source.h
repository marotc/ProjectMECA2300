#ifndef SOURCE_H
#define SOURCE_H

#include "BOV.h"
#include "./neighborhood_search.h"
#include "./point.h"
#include "kernel.h"
#include <time.h>
#include <math.h>

extern int NPTS;
extern double MASSE;
extern double SOURCE_TEMP;
extern double INIT_TEMP;
extern double dt;
extern double vmax;

// see stringification process
#define xstr(s) str(s)
#define str(s) #s
#define M_PI 3.14159265358979323846

// Structure to represent a point

void updateDensity(point* points, double kh);

void setDensity(point* points, double kh);

void updatePressure(point* points, double kh);

void updateVelocity(point* points, double kh, float vx[400], float vy[400]);

void updatePosition(point* points);

#endif