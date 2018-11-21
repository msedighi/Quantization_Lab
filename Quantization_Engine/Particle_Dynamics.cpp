
#include "stdafx.h"
#include <iostream>
#include "Particle_Dynamics.h"
#include "Data_Structures.h"
// Particle Dynamics functions 

double* Force12(double* r1, double* r2, double(*Force_func)(double), int dim)
{
	double* force = new double[dim];

	double distance = Euclidean_Distance(r1, r2, dim);
	for (int i = 0; i < dim; i++)
		force[i] = (Force_func(distance) / distance) * (r2[i] - r1[i]);
	return force;
}

double Energy12(double* r1, double* r2, double(*Energy_func)(double), int dim)
{
	double distance = Euclidean_Distance(r1, r2, dim);
	return Energy_func(distance);
}

double*** Force_Operator(double(*Force_func)(double), double** positions, int num_points, int dim, double** force_vector)
{
	double*** force_array = new double**[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		force_array[index_r] = new double*[num_points];
		force_vector[index_r] = new double[dim];
		for (int i = 0; i < dim; i++) { force_vector[index_r][i] = 0; }

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			force_array[index_r][index_c] = new double[dim];

			if (index_r == index_c)
			{
				for (int i = 0; i < dim; i++) { force_array[index_r][index_c][i] = 0; }
			}
			else
			{
				force_array[index_r][index_c] = Force12(positions[index_c], positions[index_r], Force_func, dim);
				for (int i = 0; i < dim; i++) { force_vector[index_r][i] += force_array[index_r][index_c][i]; }
			}

		}
	}
	return force_array;
}

double* Momentum(double** velocities, double* masses, int num_points, int dim)
{
	double* momentum = new double[dim];
	for (int j = 0; j < dim; j++)
	{
		momentum[j] = 0;
		for (int i = 0; i < num_points; i++)
		{
			momentum[j] += masses[i] * velocities[i][j];
		}
	}
	return momentum;
}

double Kinetic_Energy(double** velocities, double* masses, int num_points, int dim)
{
	double energy = 0;
	for (int i = 0; i < num_points; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			energy += masses[i] * velocities[i][j] * velocities[i][j] / 2;
		}
	}
	return energy;
}

double Potential_Energy(double(*energy_func)(double), double** positions, int num_points, int dim)
{
	double energy = 0;
	for (int i = 0; i < num_points; i++)
	{
		for (int j = 0; j < i; j++)
		{
			energy += Energy12(positions[i], positions[j], energy_func, dim);
		}
	}
	return energy;
}


// RadialPower Force Class
double RadialPower_Force::Coefficient = 100.0;
template<int p>
double RadialPower_Force::Force(double r)
{
	//double p = (double)n / (double)d;
	return (-Coefficient * pow(r, p));
}
template<int p>
double RadialPower_Force::Energy(double r)
{
	//double p = (double)n / (double)d;
	return (Coefficient * pow(r, p+1) / (p+1));
}

// Gravitation Class
double Gravitation::Coefficient = 100.0;
double Gravitation::Force(double r)
{
	return (-Coefficient / sqrt(r));
}
double Gravitation::Energy(double r)
{
	return (-Coefficient / r);
}

// Spring Class
double Spring::Coefficient = 1.0;
double Spring::Force(double r)
	{
		return (-Coefficient * r);
	}
double Spring::Energy(double r)
	{
		return (Coefficient * r * r / 2.0);
	}
//

// Particle Dynamics Algorithms

void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps)
{
	double dt = t / num_steps;
	double** force = new double*[num_points];
	double** a = new double*[num_points];
	double** a_plus = new double*[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		a[index_r] = new double[dim];
		a_plus[index_r] = new double[dim];
	}

	Force_Operator(Force_func, x, num_points, dim, force);
	for (long i_step = 0; i_step < num_steps; i_step++)
	{
		for (int index_r = 0; index_r < num_points; index_r++)
		{
			for (int index_c = 0; index_c < dim; index_c++)
			{
				if (i_step == 0)
					a[index_r][index_c] = force[index_r][index_c] / m[index_r];
				x[index_r][index_c] += dt * v[index_r][index_c] + dt * dt * a[index_r][index_c] / 2.0;
			}
		}

		Force_Operator(Force_func, x, num_points, dim, force);
		for (int index_r = 0; index_r < num_points; index_r++)
		{
			for (int index_c = 0; index_c < dim; index_c++)
			{
				a_plus[index_r][index_c] = force[index_r][index_c] / m[index_r];
				v[index_r][index_c] += dt * (a_plus[index_r][index_c] + a[index_r][index_c]) / 2.0;
				a[index_r][index_c] = a_plus[index_r][index_c];
			}
		}
	}
}

void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double dt, int num_points, int dim)
{
	double** force = new double*[num_points];
	double** a = new double*[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		a[index_r] = new double[dim];
	}

	Force_Operator(Force_func, x, num_points, dim, force);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		for (int index_c = 0; index_c < dim; index_c++)
		{
			a[index_r][index_c] = force[index_r][index_c] / m[index_r];
			x[index_r][index_c] += dt * v[index_r][index_c] + dt * dt * a[index_r][index_c] / 2.0;
		}
	}

	Force_Operator(Force_func, x, num_points, dim, force);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		for (int index_c = 0; index_c < dim; index_c++)
		{
			v[index_r][index_c] += dt * (force[index_r][index_c] / m[index_r] + a[index_r][index_c]) / 2.0;
		}
	}
}

void Beeman(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps)
{
	double dt = t / num_steps;
	double** force = new double*[num_points];
	double*** force_operator, **a, **a_plus, **a_minus;
	// Step zero: It is a Verlet algorithm. Since Beeman is not self-starting	
	a_minus = new double*[num_points];
	a = new double*[num_points];
	a_plus = new double*[num_points];
	force_operator = Force_Operator(Force_func, x, num_points, dim, force);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		a_minus[index_r] = new double[dim];
		a[index_r] = new double[dim];
		a_plus[index_r] = new double[dim];
		for (int index_c = 0; index_c < dim; index_c++)
		{
			a_minus[index_r][index_c] = force[index_r][index_c] / m[index_r];
		}
	}
	
	Verlet(Force_func, x, v, m, dt, num_points, dim, 1);
	//	
	for (long i_step = 1; i_step < num_steps; i_step++)
	{
		for (int index_r = 0; index_r < num_points; index_r++)
		{
			for (int index_c = 0; index_c < dim; index_c++)
			{
				if (i_step == 1)
					a[index_r][index_c] = force[index_r][index_c] / m[index_r];
				x[index_r][index_c] += dt * v[index_r][index_c] + dt * dt * (4.0 * a[index_r][index_c] - a_minus[index_r][index_c]) / 6.0;
			}
		}

		force_operator = Force_Operator(Force_func, x, num_points, dim, force);
		for (int index_r = 0; index_r < num_points; index_r++)
		{
			for (int index_c = 0; index_c < dim; index_c++)
			{
				a_plus[index_r][index_c] = force[index_r][index_c] / m[index_r];
				v[index_r][index_c] += dt * (5.0 * a_plus[index_r][index_c] + 8.0 * a[index_r][index_c] - a_minus[index_r][index_c]) / 12.0;
				a_minus[index_r][index_c] = a[index_r][index_c];
				a[index_r][index_c] = a_plus[index_r][index_c];
			}
		}
	}
}

void Runge_Kutta_5th(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps, double** e_x, double** e_v)
{
	int num_RK_steps = 6;
	double** c = new double*[num_RK_steps];
	double* b = new double[num_RK_steps + 1];
	double* _b = new double[num_RK_steps + 1];
	for (int RK_step = 0; RK_step < num_RK_steps; RK_step++)
		c[RK_step] = new double[RK_step + 1];

	c[0][0] = 1.0 / 5.0;
	c[1][0] = 3.0 / 40.0; c[1][1] = 9.0 / 40.0;
	c[2][0] = 44.0 / 45.0; c[2][1] = -56.0 / 15.0; c[2][2] = 32.0 / 9.0;
	c[3][0] = 19372.0 / 6561.0; c[3][1] = -25360.0 / 2187.0; c[3][2] = 64448.0 / 6561.0; c[3][3] = -212.0 / 729.0;
	c[4][0] = 9017.0 / 3168.0; c[4][1] = -355.0 / 33.0; c[4][2] = 46732.0 / 5247.0; c[4][3] = 49.0 / 176.0; c[4][4] = -5103.0 / 18656.0;
	c[5][0] = 35.0 / 384.0; c[5][1] = 0.0; c[5][2] = 500.0 / 1113.0; c[5][3] = 125.0 / 192.0; c[5][4] = -2187.0 / 6784.0; c[5][5] = 11.0 / 84.0;
	b[0] = c[5][0]; b[1] = c[5][1]; b[2] = c[5][2]; b[3] = c[5][3]; b[4] = c[5][4]; b[5] = c[5][5]; b[6] = 0.0;
	_b[0] = 5179.0 / 57600.0; _b[1] = 0.0; _b[2] = 7571.0 / 16695.0; _b[3] = 393.0 / 640.0; _b[4] = -92097.0 / 339200.0; _b[5] = 187.0 / 2100.0; _b[6] = 1.0 / 40.0;

	double dt = t / num_steps;
	double ***force_operator, **force, **a_temp, **v_temp, **x_temp, ***kv, ***kx;
	force = new double*[num_points];
	a_temp = new double*[num_points];
	v_temp = new double*[num_points];
	x_temp = new double*[num_points];
	kv = new double**[num_RK_steps + 1];
	kx = new double**[num_RK_steps + 1];

	for (long i_step = 0; i_step < num_steps; i_step++)
	{
		for (int RK_step = 0; RK_step < (num_RK_steps + 1); RK_step++)
		{
			kv[RK_step] = new double*[num_points];
			kx[RK_step] = new double*[num_points];

			for (int index_r = 0; index_r < num_points; index_r++)
			{
				a_temp[index_r] = new double[num_points];
				v_temp[index_r] = new double[num_points];
				x_temp[index_r] = new double[num_points];
				kv[RK_step][index_r] = new double[num_points];
				kx[RK_step][index_r] = new double[num_points];
				for (int index_c = 0; index_c < dim; index_c++)
				{
					x_temp[index_r][index_c] = x[index_r][index_c];
					v_temp[index_r][index_c] = v[index_r][index_c];
					for (int i_RK_step = 0; i_RK_step < RK_step; i_RK_step++)
					{
						x_temp[index_r][index_c] += c[RK_step - 1][i_RK_step] * kx[i_RK_step][index_r][index_c];
						v_temp[index_r][index_c] += c[RK_step - 1][i_RK_step] * kv[i_RK_step][index_r][index_c];
					}
				}
			}

			force_operator = Force_Operator(Force_func, x_temp, num_points, dim, force);
			for (int index_r = 0; index_r < num_points; index_r++)
			{
				for (int index_c = 0; index_c < dim; index_c++)
				{
					a_temp[index_r][index_c] = force[index_r][index_c] / m[index_r];
					kv[RK_step][index_r][index_c] = dt * a_temp[index_r][index_c];
					kx[RK_step][index_r][index_c] = dt * v_temp[index_r][index_c];
				}
			}
		}

		for (int index_r = 0; index_r < num_points; index_r++)
		{
			for (int index_c = 0; index_c < dim; index_c++)
			{
				for (int RK_step = 0; RK_step < (num_RK_steps + 1); RK_step++)
				{
					v[index_r][index_c] += dt * b[RK_step] * kv[RK_step][index_r][index_c];
					x[index_r][index_c] += dt * b[RK_step] * kx[RK_step][index_r][index_c];
					e_v[index_r][index_c] += dt * (b[RK_step] - _b[RK_step]) * kv[RK_step][index_r][index_c];
					e_x[index_r][index_c] += dt * (b[RK_step] - _b[RK_step]) * kx[RK_step][index_r][index_c];
				}
			}
		}
	}
}

void Simulation_Test_Spring1(double** x0, double** v0, double** x, double** v, double* m, double t, int dim, double** e_x, double** e_v)
{
	int num_points = 2;
	double w = sqrt((m[0] + m[1]) / (m[0] * m[1]));
	double** xt = new double*[num_points];
	double** vt = new double*[num_points];
	xt[0] = new double[dim];
	xt[1] = new double[dim];
	vt[0] = new double[dim];
	vt[1] = new double[dim];
	for (int j = 0; j < dim; j++)
	{
		double p0 = m[0] * v0[0][j] + m[1] * v0[1][j];
		vt[0][j] = (p0 + m[1] * (-w * (x0[0][j] - x0[1][j]) * sin(w*t) + (v0[0][j] - v0[1][j]) * cos(w*t))) / (m[0] + m[1]);
		vt[1][j] = (p0 - m[0] * (-w * (x0[0][j] - x0[1][j]) * sin(w*t) + (v0[0][j] - v0[1][j]) * cos(w*t))) / (m[0] + m[1]);
		e_v[0][j] = abs(vt[0][j] - v[0][j]);
		e_v[1][j] = abs(vt[1][j] - v[1][j]);
		double p0t = m[0] * (x0[0][j] + t * v0[0][j]) + m[1] * (x0[1][j] + t * v0[1][j]);
		xt[0][j] = (p0t + m[1] * ((x0[0][j] - x0[1][j]) * cos(w*t) + (v0[0][j] - v0[1][j]) * sin(w*t) / w)) / (m[0] + m[1]);
		xt[1][j] = (p0t - m[0] * ((x0[0][j] - x0[1][j]) * cos(w*t) + (v0[0][j] - v0[1][j]) * sin(w*t) / w)) / (m[0] + m[1]);
		e_x[0][j] = abs(xt[0][j] - x[0][j]);
		e_x[1][j] = abs(xt[1][j] - x[1][j]);
	}
}