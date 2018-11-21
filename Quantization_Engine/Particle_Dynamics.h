#pragma once
#ifndef PARTICLE_DYNAMICS_H
#define PARTICLE_DYNAMICS_H
#endif // !1

#include "iostream"

// Particle Dynamics functions 

enum Simulation_Method { V, B, RK45, V_B, B_RK, V_RK, All };

double* Force12(double* r1, double* r2, double(*Force_func)(double), int dim);

double Energy12(double* r1, double* r2, double(*Energy_func)(double), int dim);

double*** Force_Operator(double(*Force_func)(double), double** positions, int num_points, int dim, double** force_vector);

double* Momentum(double** velocities, double* masses, int num_points, int dim);

double Kinetic_Energy(double** velocities, double* masses, int num_points, int dim);

double Potential_Energy(double(*energy_func)(double), double** positions, int num_points, int dim);


class RadialPower_Force
{
public:
	static double Coefficient;
	template<int p>	static double Force(double r);
	template<int p> static double Energy(double r);
};

class Gravitation
{
public:
	static double Coefficient;
	static double Force(double r);
	static double Energy(double r);
};

class Spring
{
public:
	static double Coefficient;
	static double Force(double r);
	static double Energy(double r);
};
//

// Particle Dynamics Algorithms

void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps);

void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double dt, int num_points, int dim);

void Beeman(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps);

void Runge_Kutta_5th(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps, double** e_x, double** e_v);

void Simulation_Test_Spring1(double** x0, double** v0, double** x, double** v, double* m, double t, int dim, double** e_x, double** e_v);
