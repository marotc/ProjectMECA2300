#include "point.h"
#include "neighborhood_search.h"
#include "kernel.h"
#include <math.h>

void fillPointsRand(point* points) {
	for (int i = 0; i < NPTS; i++) {
		points[i].x = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		points[i].y = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		points[i].vx = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		points[i].vy = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		points[i].type = 1;
	}
}

void fillPointsGrid(point* points) {
	for (int i = 0; i < NPTS; i++) {
		double npt = (double)NPTS;
		points[i].x = (i % (int)sqrt(npt))*200/ NPTS -100;
		points[i].y = (i / (int)sqrt(npt))* 200 / NPTS -100;
		points[i].vx = 0 ;
		points[i].vy = 0 ;
		points[i].type = 1;
	}
}

void updateData(point* points, GLfloat(*data)[8]) {
	float rmax = 141.42;
	for (int i = 0; i < NPTS; i++) {
		data[i][0] = points[i].x;
		data[i][1] = points[i].y;
		data[i][2] = points[i].vx;
		data[i][3] = points[i].vy;
		float r = pow(data[i][0] * data[i][0] + data[i][1] * data[i][1], 0.5);
		colormap(r/rmax, &data[i][4]);
		data[i][7] = 0.8;
	}
}

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3])
{
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

void updateDensity(point* points, neighborhood* nh, double kh) {
	for (int i = 0; i < NPTS; i++) {
		points[i].density = nh[i].nNeighbours * MASSE / (kh * kh * M_PI / 4);
	}
}