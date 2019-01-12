#pragma once
#ifndef PARTICLE_DYNAMICS_H
#define PARTICLE_DYNAMICS_H
#endif // !1

#include "iostream"
#include "Eigen/Dense"

// Particle Dynamics functions 

enum Simulation_Method { V, B, RK45, V_B, B_RK, V_RK, All };

struct Interaction
{
public:
	//static double Coefficient;
	virtual double Force(double r) = 0;
	virtual double Energy(double r) = 0;
};

struct Global_Variables
{
public:
	Eigen::VectorXd Global_KineticEnergy;
	Eigen::VectorXd Global_EnergyExchange;
	Eigen::MatrixXd Global_Positions;
	Eigen::MatrixXd Global_Velocities;
	Eigen::MatrixXd Global_Momentum;

	Global_Variables(int num_points, int dimension)
	{
		Global_KineticEnergy = Eigen::VectorXd::Constant(num_points, 0);
		Global_EnergyExchange = Eigen::VectorXd::Constant(num_points, 0);
		Global_Positions = Eigen::MatrixXd::Constant(num_points, dimension, 0);
		Global_Velocities = Eigen::MatrixXd::Constant(num_points, dimension, 0);
		Global_Momentum = Eigen::MatrixXd::Constant(num_points, dimension, 0);
	}
};

double* Force12(double* r1, double* r2, Interaction* interaction, int dim);

double* Force12(double* r1, double* r2, double(*Force_func)(double), int dim);

double Energy12(double* r1, double* r2, Interaction* interaction, int dim);

double Energy12(double* r1, double* r2, double(*Energy_func)(double), int dim);

double*** Force_Operator(Interaction* interaction, double** positions, int num_points, int dim, double** force_vector);

double*** Force_Operator(double(*Force_func)(double), double** positions, int num_points, int dim, double** force_vector);

double* Momentum(double** velocities, double* masses, int num_points, int dim);

double* Energy_Exchange(Interaction* interaction, double** positions, double** velocities, int num_points, int dim);

Eigen::VectorXd CollectiveEnergy_Exchange(Interaction* interaction, double** positions, double** velocities, int num_points, int dim, Eigen::MatrixXd orthogonal_transformation);

Global_Variables* Collective_Variables(Interaction* interaction, double** positions, double** velocities, double* masses, int num_points, int dim, Eigen::MatrixXd orthogonal_transformation);

double* Kinetic_Energy(double** velocities, double* masses, int num_points, int dim, double total_energy);

double* Kinetic_Energy(double** velocities, double* masses, int num_points, int dim);

double Potential_Energy(Interaction* interaction, double** positions, int num_points, int dim);

double Potential_Energy(double(*energy_func)(double), double** positions, int num_points, int dim);


struct RadialPower_Force
{
public:
	static double Coefficient;
	template<int p>	static double Force(double r);
	template<int p> static double Energy(double r);
};

struct Gravitation : public Interaction
{
public:
	static double Coefficient;
	double Force(double r);
	double Energy(double r);
};

struct Logarithmic : public Interaction
{
public:
	static double Coefficient;
	double Force(double r);
	double Energy(double r);
};

struct Spring : public Interaction
{
public:
	static double Coefficient;
	double Force(double r);
	double Energy(double r);
};

struct InverseRoot : public Interaction
{
public:
	static double Coefficient;
	double Force(double r);
	double Energy(double r);
};

struct Lennard_Jones : public Interaction
{
public:
	static double Coefficient;
	static double MinPotential_Radius;
	double Force(double r);
	double Energy(double r);
};

//

// Particle Dynamics Algorithms
void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps);

void Verlet(Interaction* interaction, double** x, double** v, double* m, double dt, int num_points, int dim);

void Verlet(double(*Force_func)(double), double** x, double** v, double* m, double dt, int num_points, int dim);

void Beeman(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps);

void Runge_Kutta_5th(double(*Force_func)(double), double** x, double** v, double* m, double t, int num_points, int dim, long num_steps, double** e_x, double** e_v);

void Simulation_Test_Spring1(double** x0, double** v0, double** x, double** v, double* m, double t, int dim, double** e_x, double** e_v);

