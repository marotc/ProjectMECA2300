#include "neighborhood_search.h"
#include "point.h"
#include "kernel.h"

int NPTS = 100;
double SOURCE_TEMP = 1000;
double INIT_TEMP = 273;
int DOM = 200; //Utile ?
int NUMBER_ITERATIONS = 10;

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]

/*void fillData(GLfloat(* data)[8])
{
	float rmax = 100.0 * sqrtf(2.0f);
	for (int i = 0; i < NPTS; i++) {
		data[i][0] = rand() * 200.0 / RAND_MAX - 100.0; // x (rand between -100 and 100)
		data[i][1] = rand() * 200.0 / RAND_MAX - 100.0; // y (rand between -100 and 100)
		double r = sqrt(data[i][0] * data[i][0] + data[i][1] * data[i][1]);
		data[i][2] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		data[i][3] = rand() * 2.0 / RAND_MAX - 1.0; //Random starting speed
		colormap(r / rmax, &data[i][4]); // fill color
		data[i][7] = 0.8f; // transparency
	}
}*/

int main()
{
	point* points = malloc(sizeof(point) * NPTS);
	fillPointsGrid(points);
	//kernel(points,neighbors,kh);
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);
	updateData(points, data);
	//fillData(data);


	double timestep = 0.5;
	double maxspeed = 1;
	neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	neighborhood* nh = options->nh;
	
	for (int iterations = 0; iterations < NUMBER_ITERATIONS;iterations++) {
		if(iterations)
			bouncyrandomupdate(points, timestep, options->half_length, maxspeed); 
		neighborhood_update(options, nh, points, iterations);
		//kernel(data, nh, kh);
	}
	printNeighborhood(nh, points);

	neighborhood_options_delete(options,nh);

	free(data);
	free(points);
	return EXIT_SUCCESS;
}
