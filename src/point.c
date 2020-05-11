#include "point.h"
#include "neighborhood_search.h"
#include "kernel.h"
#include <math.h>

void fillPointsRand(point* points) {
	for (int i = 0; i < NPTS; i++)
	{
		if (i < NPTS - 4 * NPTS_boundaries) // normal particles
		{
			//data[i][0] = (i / (int)sqrt(NPTS)) * 200 / sqrt(NPTS) - 100;
			//data[i][1] = (i % (int)sqrt(NPTS)) * 200 / sqrt(NPTS) - 100;
			//points[i].x = (i / (int)sqrt(NPTS - 4 * NPTS_boundaries)) * 100 / sqrt(NPTS - 4 * NPTS_boundaries);
			//points[i].y = (i % (int)sqrt(NPTS - 4 * NPTS_boundaries)) * 100 / sqrt(NPTS - 4 * NPTS_boundaries);
			points[i].x = rand() * 100 / RAND_MAX; // x (rand between -100 and 100)
			points[i].y = rand() * 100 / RAND_MAX; // y (rand between -100 and 100)
			points[i].vx = 0;
			points[i].vy = 0;
			if (points[i].x < 3)
			{
				//points[i].vy = 40;
			}
		}
		// Boundary 
		else if (i < NPTS - 3 * NPTS_boundaries) // Left wall
		{
			int j = i - (NPTS - 4 * NPTS_boundaries);
			int row = (int)floor(j/26);
			points[i].x = -5*(row+1); // x
			points[i].y =(j%26)*5-25;// y
			points[i].vx = 0;
			points[i].vy = 0;
		}
		else if (i < NPTS - 2 * NPTS_boundaries) // right wall
		{
			int j = i - (NPTS - 3 * NPTS_boundaries);
			int row = (int)floor(j / 26);
			points[i].x = 100+5 * (row + 1); // x
			points[i].y = (j % 26) * 5;// y
			points[i].vx = 0;
			points[i].vy = 0;
		}
		else if (i < NPTS - NPTS_boundaries) // upper wall
		{
			int j = i - (NPTS - 2 * NPTS_boundaries);
			int row = (int)floor(j / 26);
			points[i].x = (j % 26) * 5-25;// y 
			points[i].y = 100 + 5 * (row + 1); // x
			points[i].vx = 0;
			points[i].vy = 0;
		}
		else
		{
			int j = i - (NPTS - NPTS_boundaries);
			int row = (int)floor(j / 26);
			points[i].x = (j % 26) * 5;// y 
			points[i].y =  -5 * (row + 1); // x
			points[i].vx = 0;
			points[i].vy = 0;

		}
		
	}
}


void updateData(point* points, GLfloat(*data)[8]) {
	for (int i = 0; i < NPTS; i++) 
	{
		data[i][0] = points[i].x;
		data[i][1] = points[i].y;
		data[i][2] = points[i].vx;
		data[i][3] = points[i].vy;
		if (i < NPTS - 4 * NPTS_boundaries)
		{
			float v = pow(data[i][2] * data[i][2] + data[i][3] * data[i][3], 0.5);
			colormap(v/500, &data[i][4]); // fill color v/4 pas mal
			data[i][7] = 0.8;
		}
		else
		{
			float r = pow(data[i][0] * data[i][0] + data[i][1] * data[i][1], 0.5);
			colormap(10, &data[i][4]); // fill color
			data[i][7] = 1;
		}
	}
}

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
static void colormap(float v, float color[3])
{
	if (v == 9999)
	{
		color[0] = 0;
		color[1] = 0;
		color[2] = 0;
	}
	else
	{
		float v1 = 3.5 * (v - 0.7);
		float v2 = 1.25 * v;
		float v3 = fminf(0.5, v) * 2.0;

		//color[0] = -v1 * v1 + 1.0f;
		//color[1] = 6.0f * v2 * v2 * (1.0f - v2);
		//color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

		// alternative: classical jet colormap
		color[0] = 1.5 - 4.0 * fabs(v - 0.75);
		color[1] = 1.5 - 4.0 * fabs(v - 0.5);
		color[2] = 1.5 - 4.0 * fabs(v - 0.25);
	}
}
