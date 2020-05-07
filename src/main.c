#include "point.h"
#include "neighborhood_search.h"
#include "kernel.h"
#include <math.h>

int NPTS = 2500;
double MASSE = 1;
int DOM = 200; //Utile ?
int NUMBER_ITERATIONS = 10;

int main(){
	printf("pouit\n");
	point* points = malloc(sizeof(point) * NPTS);
	fillPointsRand(points);
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);
	updateData(points, data);


	double timestep = 0.5;
	double maxspeed = 1;
	//neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	//neighborhood* nh = options->nh;
	double kh = compute_KH()*200;
	neighborhood* nh = createNeighborhood(kh, points);
	updateDensity(points, nh, kh);
	
	for (int i = 0; i < NPTS; i++) {
		printf("%d: %d\n", nh[i].index, nh[i].nNeighbours);
	}
	deleteNeighborhood(nh);
	for (int iterations = 0; iterations < NUMBER_ITERATIONS; iterations++) {
		if (iterations) {
			//bouncyrandomupdate(points, timestep, options->half_length, maxspeed);
			//neighborhood_update(options, nh, points, iterations);
		}
	}
	//printNeighborhood(nh, points);

	//neighborhood_options_delete(options,nh);
	
	free(data);
	free(points);
	printf("pouet\n");
	return EXIT_SUCCESS;
}
