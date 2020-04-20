#include "point.h"
#include "neighborhood_search.h"
#include "kernel.h"

void fillPointsRand(point* points) {
	for (int i = 0; i < NPTS; i++) {
		points[i].x = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		points[i].y = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		points[i].vx = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		points[i].vy = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		if (points[i].x > -10 && points[i].x<10 && points[i].y>-10 && points[i].y < 10) {
			points[i].val = SOURCE_TEMP;
		}
		else {
			points[i].val = INIT_TEMP;
		}
	}
}

void fillPointsGrid(point* points) {
	for (int i = 0; i < NPTS; i++) {
		double npt = (double)NPTS;
		points[i].x = (i % (int)sqrt(npt))*200/ NPTS -100;
		points[i].y = (i % (int)sqrt(npt))* 200 / NPTS -100;
		points[i].vx = 0 ;
		points[i].vy = 0 ;
		if (points[i].x > -10 && points[i].x<10 && points[i].y>-10 && points[i].y < 10) {
			points[i].val = SOURCE_TEMP;
		}
		else {
			points[i].val = INIT_TEMP;
		}
	}
}

// funtion to compute the radius kh of the circle of influence of a particle
// nPoints : number of particles in the simulation
// RA : int used as a boolean to choose over the algorithm of the radius choice; 
//      0 means the dummy algorithm, 1 means the more sophisticated algorithm expalined at the seminar ( (intersection between areaGrid and 1/4 areaCircle)/areaGrid = 21/nPoints)
double compute_kh(int RA) {
	double target = 21.0 / NPTS;
	if (NPTS < 21)
		return sqrt(2.0);
	else if (!RA)
		return sqrt(2.0) * target;
	else if (target <= M_PI / 4)
		return sqrt(4 * target / M_PI);
	double tolerance = 0.0000000001;
	double kh_min = 1.0;
	double kh_max = sqrt(2.0);
	while (fabs(kh_max - kh_min) >= tolerance) {
		kh_max = kh_min;
		kh_min = sqrt((target - sin(acos(1.0 / kh_min)) * kh_min) / ((M_PI / 4 - acos(1.0 / kh_min))));
	}
	return kh_min;
}

void updateData(point* points, GLfloat(*data)[8]) {
	for (int i = 0; i < NPTS; i++) {
		data[i][0] = points[i].x;
		data[i][1] = points[i].y;
		data[i][2] = points[i].vx;
		data[i][3] = points[i].vy;
		colormap((points[i].val-INIT_TEMP)/(SOURCE_TEMP-INIT_TEMP), &data[i][4]);
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