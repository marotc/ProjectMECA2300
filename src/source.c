#include "./source.h"
#include <math.h>
void setDensity(point* points, double kh)
{
	for (int i = 0; i < NPTS; i++)
	{
		double density = 0;
		// double density = MASS * w_lucy(0, kh);
		//neighbours* point = nh[i].list;
		//for (int j = 0; j < (int)nh[i].nNeighbours; j++)
		for (int j = 0;j<NPTS;j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x)/100.0;
			double dy = ((double)points[i].y - (double)points[j].y)/100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh)
			{
				density += MASSE * w_lucy(distance, kh);
			}
			//point = point->next;
		}
		points[i].old_density = density;
	}
	for (int i = 0; i < NPTS; i++)
	{
		double num = 1;
		double denom = 1 / points[i].old_density;
		//neighbours* point = nh[i].list;
		for (int j = 0; j < NPTS; j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
			double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh)
			{
				num += w_lucy(distance, kh);
				denom += w_lucy(distance, kh) / points[j].old_density;
			}

		}
		points[i].density = num / denom;
		if (points[i].density < 1000)
		{
			points[i].density = 1000;
		}
	}
}
// DENSITY SET OK
void updateDensity(point* points, double kh)
{
	for (int i = 0; i < NPTS - 4 * NPTS_boundaries; i++)
	{
		double der_density = 0;
		double vxi = (double)points[i].vx / 100.0;
		double vyi = (double)points[i].vy / 100.0;
		double speed_i = sqrt(vxi * vxi + vyi * vyi);
		for (int j = 0; j < NPTS-4*NPTS_boundaries; j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
			double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh)
			{
				double vxj = (double)points[j].vx / 100.0;
				double vyj = (double)points[j].vy / 100.0;
				double speed_j = sqrt(vxj * vxj + vyj * vyj);
				der_density += MASSE * (speed_i - speed_j) * grad_w_lucy(distance, kh);
			}
		}
		points[i].density = points[i].density + dt * der_density;
		//Euler pour le moment
	}
}

void updatePressure(point* points, double kh)
{
	double rho_0= 1000; // kg/m^3 => minimum value with 400 particles
	//rho_0 = 100;
	double B = 10130;//10130; // Pa
	double gamma = 7;
	for (int i = 0; i < NPTS; i++)
	{
		points[i].pressure = B * (pow(points[i].density / rho_0, gamma)-1);
	} 
	for (int i = 0; i < NPTS-4*NPTS_boundaries; i++)
	{
		points[i].old_dpdx = 0;
		points[i].old_dpdy = 0; // Need to initialise the variable first
		for (int j = 0; j < NPTS ; j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
			double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh && i != j)
			{
				points[i].old_dpdx -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(dx, distance, kh);
				points[i].old_dpdy -= MASSE * (points[i].pressure / (points[i].density * points[i].density) + points[j].pressure / (points[j].density * points[j].density)) * grad_w_lucy1D(dy, distance, kh);
			}
		}
	}
	for (int i = 0; i < NPTS-4*NPTS_boundaries; i++)
	{
		double numx = points[i].old_dpdx;
		double numy = points[i].old_dpdy;
		double denomx = 0;
		double denomy = 0;
		for (int j = 0; j < NPTS - 4 * NPTS_boundaries; j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
			double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh && i!=j )
			{
				denomx -= grad_w_lucy1D(dx , distance, kh) * dx * MASSE / points[j].density;
				denomy -= grad_w_lucy1D(dy, distance, kh) * dy * MASSE / points[j].density;
			}
		}
		points[i].dpdx = numx / denomx;
		points[i].dpdy = numy / denomy;
	}
}

void updateVelocity(point* points, double kh, float vx[400], float vy[400])
{
	double kin_v = 0.000001;// kinematic viscosity m^2/s
	double g = -9.81;// 0.005; // gravity m/s^2
	for (int i = 0; i < NPTS-4*NPTS_boundaries; i++)
	{
		double cut_off = kh;//0.05;kh =0.25
		float dvx = points[i].dpdx;
		float dvy = points[i].dpdy+g;///g
		for (int j = 0; j < NPTS-4*NPTS_boundaries; j++)
		{
			double dx = ((double)points[i].x - (double)points[j].x) / 100.0;
			double dy = ((double)points[i].y - (double)points[j].y) / 100.0;
			double distance = sqrt(dx * dx + dy * dy);
			if (distance < kh && i!=j)
			{
				// PRESSURE PART
				//dvx -= 1 / points[i].density * MASSE * ((double)points[i].pressure / ((double)points[i].density * (double)points[i].density) + (double)points[j].pressure / ((double)points[j].density * (double)points[j].density)) * grad_w_lucy(distance, kh) * dx / distance;
				//dvy -= 1 / points[i].density * MASSE * ((double)points[i].pressure / ((double)points[i].density * (double)points[i].density) + (double)points[j].pressure / ((double)points[j].density * (double)points[j].density)) * grad_w_lucy(distance, kh) * dy / distance;
				// VELOCITY PART
				dvx += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vx - (double)points[j].vx) * grad_w_lucy1D(dy,distance, kh) /(dy+0.0001);// *dx / (dx * dx);
				dvy += 2 * kin_v * MASSE / (points[j].density) * ((double)points[i].vy - (double)points[j].vy) * grad_w_lucy1D(dx,distance, kh) /(dx+0.0001); //*dy / (dy * dy);
				// type Lennard-Jones
				if (distance < cut_off)
				{
					//double damper = 0.00000005;
					//dvx += damper*0.01 * dx / distance * (pow(cut_off / distance, 12) - pow(cut_off / distance, 4)) / MASSE;
					//dvy += damper*0.01 * dy / distance * (pow(cut_off / distance, 12) - pow(cut_off / distance, 4)) / MASSE;
				}
				// type coulomb
				//double k = 100;
				//dvx += k / points[i].density * (1/(double)points[i].density + 1/(double)points[i].density) / distance * dx / distance;
				//dvy += k / points[i].density * (1/(double)points[i].pressure + 1/(double)points[i].pressure) / distance * dy / distance;
				// loi parabolique
				//dvx += k*points[j].density*1/distance*dx / distance;
				//dvy +=k* points[j].density * 1 / distance * dy / distance;
				// d rho => f =>
				//double k = 0.000000001;
				//dvx += k * points[j].density * dx / (distance+0.001) * ((double)points[i].density - (double)points[j].density);
				//dvy += k * points[j].density * dy / (distance+0.001) * ((double)points[i].density - (double)points[j].density);
			}
		}
		vx[i] = points[i].vx + dvx * dt*100; // euler method
		vy[i] = points[i].vy + dvy * dt*100;
	}
	for (int i = 0; i < NPTS-4*NPTS_boundaries; i++)
	{
		
		points[i].vx = vx[i];
		points[i].vy = vy[i];
	}
}

void updatePosition(point* points)
{
	for (int i = 0; i < NPTS-4*NPTS_boundaries; i++)
	{
		//points[i].x += points[i].vx * dt; //euler
		//points[i].y += points[i].vy * dt;
		double cut_off = 5;
		if (points[i].x > 100)
		{
			points[i].vy = 0;
		}
		else if (points[i].x < 0)
		{
			points[i].vy = 0;
		}
		if (points[i].y > 100)
		{
			points[i].vx = 0;
		}
		else if (points[i].y<0)
		{
			points[i].vx = 500;
		}
		if (points[i].x > 100)
		{
			points[i].vx -= (pow(cut_off / (cut_off - ((double)points[i].x - 100)), 12) - pow(cut_off / (cut_off - ((double)points[i].x - 100)), 4)) / MASSE;
			/*points[i].x = 100 - (points[i].x - 100);
			points[i].vx = - 0.8*points[i].vx;
			//points[i].vy = - 1;
			if (points[i].x < 0)
			{
				points[i].x = i/2;
				points[i].vx = 0;
				points[i].vy = 0;
			}*/
		}
		else if (points[i].x < 0)
		{
			points[i].vx += (pow(cut_off / (cut_off+points[i].x), 12) - pow(cut_off / (cut_off + points[i].x), 4)) / MASSE;
			/*points[i].x = -5-(points[i].x+5);
			points[i].vx = -0.8*points[i].vx;
			points[i].vy = 1;
			if (points[i].x > 100)
			{
				points[i].x = i/2;
				points[i].vx = 0;
				points[i].vy = 0;
			}*/
		}
		if (points[i].y > 100)
		{
			points[i].vy -= (pow(cut_off / (cut_off - ((double)points[i].y - 100)), 12) - pow(cut_off / (cut_off - ((double)points[i].y - 100)), 4)) / MASSE;
			/*	
			points[i].y = 100 - (points[i].y - 100);
			points[i].vy = -0.8*points[i].vy;
			points[i].vx =  1;
			if (points[i].y < -0)
			{
				points[i].y = i/2;
				points[i].vx = 0;
				points[i].vy = 0;
			}*/

		}
		else if (points[i].y < 0)
		{
			points[i].vy += (pow(cut_off / (cut_off + points[i].y), 12) - pow(cut_off / (cut_off + points[i].y), 4)) / MASSE;
			/*
			points[i].y = -5 - (points[i].y + 5);
			points[i].vy = -0.8*points[i].vy;
			points[i].vx = - 1;
			if (points[i].y > 100)
			{
				points[i].y = i/2;
				points[i].vx = 0;
				points[i].vy = 0;
			}*/
		}
		points[i].x += points[i].vx * dt; //euler
		points[i].y += points[i].vy * dt;
	}
}
