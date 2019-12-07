
#include "stdafx.h"
#include <iostream>
#include "Particle_Dynamics.h"
#include "Data_Structures.h"
#include <cmath>

using namespace Eigen;
// Particle Dynamics functions 

double* Force12(double* r1, double* r2, Interaction* interaction, int dim)
{
	double* force = new double[dim];

	double distance = Euclidean_Distance(r1, r2, dim);
	for (int i = 0; i < dim; i++)
		force[i] = (interaction->Force(distance) / distance) * (r2[i] - r1[i]);
	return force;
}

double* Force12(double* r1, double* r2, double(*Force_func)(double), int dim)
{
	double* force = new double[dim];

	double distance = Euclidean_Distance(r1, r2, dim);
	for (int i = 0; i < dim; i++)
		force[i] = (Force_func(distance) / distance) * (r2[i] - r1[i]);
	return force;
}

double Energy12(double* r1, double* r2, Interaction* interaction, int dim)
{
	double distance = Euclidean_Distance(r1, r2, dim);
	return interaction->Energy(distance);
}

double Energy12(double* r1, double* r2, double(*Energy_func)(double), int dim)
{
	double distance = Euclidean_Distance(r1, r2, dim);
	return Energy_func(distance);
}

double*** Force_Operator(Interaction* interaction, double** positions, int num_points, int dim, double** force_vector)
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
				force_array[index_r][index_c] = Force12(positions[index_c], positions[index_r], interaction, dim);
				for (int i = 0; i < dim; i++) { force_vector[index_r][i] += force_array[index_r][index_c][i]; }
			}

		}
	}
	return force_array;
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

double*** VelocityDiff_Operator(double ** velocities, int num_points, int dim)
{
	double*** velocity_array = new double**[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		velocity_array[index_r] = new double*[num_points];

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			velocity_array[index_r][index_c] = new double[dim];

			if (index_r == index_c)
			{
				for (int i = 0; i < dim; i++) 
					velocity_array[index_r][index_c][i] = 0;
			}
			else
			{
				for (int i = 0; i < dim; i++) 
					velocity_array[index_r][index_c][i] = velocities[index_c][i] - velocities[index_r][i];
			}

		}
	}
	return velocity_array;
}

double*** Vertex2Edge_VectorField_Operator(double** vectorfield, int num_points, int dim)
{
	double*** vectorfield_operator = new double**[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		vectorfield_operator[index_r] = new double*[num_points];

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			vectorfield_operator[index_r][index_c] = new double[dim];

			if (index_r == index_c)
			{
				for (int i = 0; i < dim; i++)
					vectorfield_operator[index_r][index_c][i] = 0;
			}
			else
			{
				for (int i = 0; i < dim; i++)
					vectorfield_operator[index_r][index_c][i] = vectorfield[index_c][i] - vectorfield[index_r][i];
			}

		}
	}
	return vectorfield_operator;
}

double *** MomentumDiff_Operator(double ** velocities, double * masses, int num_points, int dim)
{
	double** massdiff_operator = MassDiff_Operator(masses, num_points);
	double*** velocitydiff_operator = VelocityDiff_Operator(velocities, num_points, dim);

	return MomentumDiff_Operator(velocitydiff_operator, massdiff_operator, num_points, dim);
}

double *** MomentumDiff_Operator(double *** velocitydiff_operator, double ** massdiff_operator, int num_points, int dim)
{
	double*** momentum_array = new double**[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		momentum_array[index_r] = new double*[num_points];

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			momentum_array[index_r][index_c] = new double[dim];

			for (int i = 0; i < dim; i++)
				momentum_array[index_r][index_c][i] = massdiff_operator[index_r][index_c] * velocitydiff_operator[index_r][index_c][i];

		}
	}
	return momentum_array;
}

double ** KEnergyDiff_Operator(double *** velocitydiff_operator, double ** massdiff_operator, int num_points, int dim)
{
	double** kenergydiff_array = new double*[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		kenergydiff_array[index_r] = new double[num_points];

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			kenergydiff_array[index_r][index_c] = 0;

			if (index_r != index_c)
			{
				for (int i = 0; i < dim; i++)
					kenergydiff_array[index_r][index_c] += velocitydiff_operator[index_r][index_c][i] * velocitydiff_operator[index_r][index_c][i];

				kenergydiff_array[index_r][index_c] = massdiff_operator[index_r][index_c] * kenergydiff_array[index_r][index_c] / 2; // Kinetic Energy between Particles
			}
		}
	}
	return kenergydiff_array;
}

Eigen::MatrixXd KEnergyDiff_Operator(double ** velocities, double * masses, int num_points, int dim)
{
	double** massdiff_operator = MassDiff_Operator(masses, num_points);
	double*** velocitydiff_operator = VelocityDiff_Operator(velocities, num_points, dim);

	Eigen::MatrixXd kenergydiff_matrix = MatrixXd::Constant(num_points, num_points, 0);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		for (int index_c = 0; index_c < num_points; index_c++)
		{
			if (index_r != index_c)
			{
				for (int i = 0; i < dim; i++)
					kenergydiff_matrix(index_r,index_c) += velocitydiff_operator[index_r][index_c][i] * velocitydiff_operator[index_r][index_c][i];

				kenergydiff_matrix(index_r, index_c) = massdiff_operator[index_r][index_c] * kenergydiff_matrix(index_r, index_c) / 2.0; // Kinetic Energy between Particles
			}
		}
	}
	return kenergydiff_matrix;
}

Eigen::MatrixXd Compute_VectorField_Commutator(double** first_vectorfield, double** second_vectorfield, int num_points, int dim)
{
	double*** First_VectorField_operator = Vertex2Edge_VectorField_Operator(first_vectorfield, num_points, dim);
	double*** Second_VectorField_operator = Vertex2Edge_VectorField_Operator(second_vectorfield, num_points, dim);

	Eigen::MatrixXd commutator = MatrixXd::Constant(num_points, num_points, 0);
	for (int index_i = 0; index_i < num_points; index_i++)
	{
		for (int index_j = 0; index_j < num_points; index_j++)
		{
			commutator(index_i, index_j) = 0;
			for (int index_k = 0; index_k < num_points; index_k++)
			{
				for (int i = 0; i < dim; i++)
					commutator(index_i, index_j) += First_VectorField_operator[index_i][index_k][i] * Second_VectorField_operator[index_k][index_j][i];
			}

			for (int index_k = 0; index_k < num_points; index_k++)
			{
				for (int i = 0; i < dim; i++)
					commutator(index_i, index_j) -= Second_VectorField_operator[index_i][index_k][i] * First_VectorField_operator[index_k][index_j][i];
			}
		}
	}
	return commutator;
}

double ** MassDiff_Operator(double * masses, int num_points)
{
	double total_mass = 0;
	for (int index = 0; index < num_points; index++)
	{
		total_mass += masses[index];
	}
	
	double** mass_array = new double*[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		mass_array[index_r] = new double[num_points];

		for (int index_c = 0; index_c < num_points; index_c++)
		{
			if (index_r == index_c)
			{
				mass_array[index_r][index_c] = 0;
			}
			else
			{
				mass_array[index_r][index_c] = masses[index_c] * masses[index_r] / total_mass;
			}

		}
	}
	return mass_array;
}

double* Energy_Exchange(Interaction* interaction, double** positions, double** velocities, int num_points, int dim)
{
	double* energy_exchange = new double[num_points];
	double** local_force = new double*[num_points];
	Force_Operator(interaction, positions, num_points, dim, local_force);
	for (int i = 0; i < num_points; i++)
	{
		energy_exchange[i] = 0;
		for (int j = 0; j < dim; j++)
		{
			energy_exchange[i] += local_force[i][j] * velocities[i][j];
		}
	}
	return energy_exchange;
}

VectorXd CollectiveEnergy_Exchange(Interaction* interaction, double** positions, double** velocities, int num_points, int dim, MatrixXd orthogonal_transformation)
{
	MatrixXd velocities_matrix = MatrixXd::Constant(num_points, dim, 0);
	MatrixXd local_force_matrix = MatrixXd::Constant(num_points, dim, 0);

	double** local_force = new double*[num_points];
	Force_Operator(interaction, positions, num_points, dim, local_force);
	for (int i = 0; i < num_points; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			velocities_matrix(i,j) = velocities[i][j];
			local_force_matrix(i, j) = local_force[i][j];
		}
	}

	MatrixXd global_velocities = orthogonal_transformation.transpose() * velocities_matrix;
	MatrixXd global_forces = orthogonal_transformation.transpose() * local_force_matrix;

	VectorXd collectiveenergy_exchange = VectorXd::Constant(num_points, 0);
	for (int i = 0; i < num_points; i++)
	{
		collectiveenergy_exchange[i] = global_forces.row(i).dot(global_velocities.row(i));
	}

	if (false)
	{
		std::cout << "Hierarchical Transformation: " << std::endl;
		std::cout << orthogonal_transformation << std::endl;
		std::cout << std::endl;
		std::cout << "Local Forces: " << std::endl;
		std::cout << local_force_matrix << std::endl;
		std::cout << std::endl;
		std::cout << "Global Forces: " << std::endl;
		std::cout << global_forces << std::endl;
		std::cout << std::endl;
	}

	return collectiveenergy_exchange;
}

Global_Variables* Collective_Variables(Interaction* interaction, double** positions, double** velocities, double* masses, int num_points, int dim, Eigen::MatrixXd orthogonal_transformation)
{
	Global_Variables* collective_variables = new Global_Variables(num_points, dim);

	MatrixXd positions_matrix = MatrixXd::Constant(num_points, dim, 0);
	MatrixXd velocities_matrix = MatrixXd::Constant(num_points, dim, 0);
	MatrixXd momentum_matrix = MatrixXd::Constant(num_points, dim, 0);
	VectorXd mass_vector = VectorXd::Constant(num_points, 0);

	MatrixXd local_force_matrix = MatrixXd::Constant(num_points, dim, 0);

	double** local_force = new double*[num_points];
	Force_Operator(interaction, positions, num_points, dim, local_force);
	for (int i = 0; i < num_points; i++)
	{
		mass_vector[i] = masses[i];
		for (int j = 0; j < dim; j++)
		{
			positions_matrix(i, j) = positions[i][j];
			velocities_matrix(i, j) = velocities[i][j];
			momentum_matrix(i, j) = masses[i] * velocities[i][j];

			local_force_matrix(i, j) = local_force[i][j];
		}
	}

	collective_variables->Global_Positions = orthogonal_transformation.transpose() * positions_matrix;
	collective_variables->Global_Velocities = orthogonal_transformation.transpose() * velocities_matrix;
	collective_variables->Global_Momentum = orthogonal_transformation.transpose() * momentum_matrix;

	VectorXd global_mass_vector = orthogonal_transformation.transpose() * mass_vector;
	MatrixXd global_forces = orthogonal_transformation.transpose() * local_force_matrix;

	for (int i = 0; i < num_points; i++)
	{
		collective_variables->Global_EnergyExchange[i] = global_forces.row(i).dot((collective_variables->Global_Velocities).row(i));
		collective_variables->Global_KineticEnergy[i] = (collective_variables->Global_Velocities).row(i).dot((collective_variables->Global_Momentum).row(i)) / 2;
	}

	if (false)
	{
		std::cout << "Hierarchical Transformation: " << std::endl;
		std::cout << orthogonal_transformation << std::endl;
		std::cout << std::endl;
		std::cout << "Global Kinetic Energy: " << std::endl;
		std::cout << collective_variables->Global_KineticEnergy << std::endl;
		std::cout << std::endl;
	}

	return collective_variables;
}

double* Kinetic_Energy(double** velocities, double* masses, int num_points, int dim, double total_energy)
{
	double* energy_vector = new double[num_points];
	total_energy = 0;
	for (int i = 0; i < num_points; i++)
	{
		energy_vector[i] = 0;
		for (int j = 0; j < dim; j++)
		{
			energy_vector[i] += masses[i] * velocities[i][j] * velocities[i][j];
		}
		total_energy += energy_vector[i] / 2.0;
	}

	return energy_vector;
}

double* Kinetic_Energy(double** velocities, double* masses, int num_points, int dim)
{
	double total_energy = 0;
	double* energy_vector = Kinetic_Energy(velocities, masses, num_points, dim, total_energy);

	return energy_vector;
}

double Potential_Energy(Interaction* interaction, double** positions, int num_points, int dim)
{
	double energy = 0;
	for (int i = 0; i < num_points; i++)
	{
		for (int j = 0; j < i; j++)
		{
			energy += Energy12(positions[i], positions[j], interaction, dim);
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
double Gravitation::Coefficient = 10.0;
double Gravitation::Force(double r)
{
	return (-Coefficient / (r*r));
}
double Gravitation::Energy(double r)
{
	return (-Coefficient / r);
}
// Logarithmic Class
double Logarithmic::Coefficient = 100.0;
double Logarithmic::Force(double r)
{
	return (-Coefficient / r);
}
double Logarithmic::Energy(double r)
{
	return (Coefficient * log(r));
}
// Inverse Root
double InverseRoot::Coefficient = 100.0;
double InverseRoot::Force(double r)
{
	return (-Coefficient / sqrt(r));
}
double InverseRoot::Energy(double r)
{
	return (2 * Coefficient * sqrt(r));
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
// Lennard Jones Class
double Lennard_Jones::Coefficient = 100.0;
double Lennard_Jones::MinPotential_Radius = 10;
double Lennard_Jones::Force(double r)
{
	return (- 12 * (Coefficient / MinPotential_Radius) * (pow(MinPotential_Radius / r, 13) - pow(MinPotential_Radius / r, 7)));
}
double Lennard_Jones::Energy(double r)
{
	return (Coefficient * (pow(MinPotential_Radius / r, 12) - 2 * pow(MinPotential_Radius / r, 6)));
}
// Constant Class
double Constant::Coefficient = 10.0;
double Constant::Force(double r)
{
	return (-Coefficient);
}
double Constant::Energy(double r)
{
	return (Coefficient * abs(r) / 2.0);
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

void Verlet(Interaction* interaction, double** x, double** v, double* m, double dt, int num_points, int dim)
{
	double** force = new double*[num_points];
	double** a = new double*[num_points];
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		a[index_r] = new double[dim];
	}

	Force_Operator(interaction, x, num_points, dim, force);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		for (int index_c = 0; index_c < dim; index_c++)
		{
			a[index_r][index_c] = force[index_r][index_c] / m[index_r];
			x[index_r][index_c] += dt * v[index_r][index_c] + dt * dt * a[index_r][index_c] / 2.0;
		}
	}

	Force_Operator(interaction, x, num_points, dim, force);
	for (int index_r = 0; index_r < num_points; index_r++)
	{
		for (int index_c = 0; index_c < dim; index_c++)
		{
			v[index_r][index_c] += dt * (force[index_r][index_c] / m[index_r] + a[index_r][index_c]) / 2.0;
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