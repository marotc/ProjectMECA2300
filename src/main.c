#include "point.h"
#include "neighborhood_search.h"
#include "kernel.h"
#include "source.h"
#include <math.h>

int NPTS = 400+25*5*4+20;
int NPTS_boundaries = 5*25+5;// 24;
double vmax = 1;
double MASSE = 2.5;
int NUMBER_ITERATIONS = 10;
double dt = 0.001;

int main(){
	bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
	bov_window_set_color(window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 0.0f });
	point* points = malloc(sizeof(point) * NPTS);
	fillPointsRand(points);
	GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
	CHECK_MALLOC(data);

	updateData(points, data);
	/* send data to GPU, and receive reference to those data in a points object */
	bov_points_t* particles = bov_particles_new(data, NPTS, GL_STATIC_DRAW);

	/* setting particles appearance */
	bov_points_set_width(particles, 0.02);
	bov_points_set_outline_width(particles, 0.0025);

	/* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
	bov_points_scale(particles, (GLfloat[2]) { 0.008, 0.008 });
	bov_points_set_pos(particles, (GLfloat[2]) { 0.0, -0.1 });
	bov_text_t* msg = bov_text_new(
		(GLubyte[]) {
		"Rendering " xstr(NPTS) " particles"
	},
		GL_STATIC_DRAW);
	bov_text_set_pos(msg, (GLfloat[2]) { -0.95, 0.82 });
	bov_text_set_fontsize(msg, 0.1);
	// Seed the random
	time_t seed = time(NULL);
	//printf(" %u \n", seed);
	srand(seed);
	//neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
	//neighborhood* nh = options->nh;
	double kh = compute_KH();
	printf("%f", kh);
	float vx [400];
	float vy [400];
	float New_rho[400];
	//neighborhood* nh = createNeighborhood(kh, points);
	while (!bov_window_should_close(window)) 
	{
		if (GetAsyncKeyState(VK_RETURN) & 0x8000) {
			fillPointsRand(points);
			updateData(points, data);
		}
		setDensity(points, kh, New_rho);
		//updateDensity(points, kh);
		updatePressure(points, kh);
		updateVelocity(points, kh, vx, vy);
		updatePosition(points);
		updateData(points, data);
	
		bov_text_draw(window, msg);
		bov_particles_draw(window, particles, 0, BOV_TILL_END);
		bov_particles_update(particles, data, NPTS);
		char s[60];
		sprintf(s, " ", time);
		bov_text_update(msg, s);
		bov_window_update(window); //don't wait for events => bov_window_update(window)
	}
	//deleteNeighborhood(nh);
	// Creation des points
	free(data);
	free(points);
	free(New_rho);
	free(vx);
	free(vy);
	return EXIT_SUCCESS;
}
