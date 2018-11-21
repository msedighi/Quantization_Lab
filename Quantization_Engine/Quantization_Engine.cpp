// Quantization_Engine.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>

#include "Quantization_Engine.h"

using namespace std;
using namespace Eigen;
using namespace chrono;

void Compute::Perturb(Eigen::MatrixXd dH, int num_points, long num_scale_bins, Eigen::VectorXd energy_vector, Eigen::MatrixXd orthonormal_transformation)
{
	MatrixXd dH_prime0, dH_prime;
	MatrixXd dU_prime0, dU_prime;
	VectorXd dE;

	dU_prime = MatrixXd::Constant(num_points, num_points, 0);

	dH_prime0 = orthonormal_transformation.transpose() * dH * orthonormal_transformation;
	dH_prime = dH_prime0;

	for (int i_o = 0; i_o < Peturb_Order; i_o++)
	{
		dE = dH_prime.diagonal();
		for (int i = 0; i < num_points; i++)
		{
			for (int j = 0; j < num_points; j++)
			{
				if (i == j)
				{
					if (i_o > 0)
					{
						dU_prime(i, i) = -0.5*(dU_prime0 * dU_prime0.transpose())(i, i);
					}
				}
				else
				{
					dU_prime(i, j) = dH_prime(i, j) / (energy_vector[j] - energy_vector[i]);
				}
			}
		}

		dU_prime0 = dU_prime;
		dH_prime = dH_prime0 + dH_prime0 * dU_prime - dU_prime * dE.asDiagonal();

	}

	energy_vector = energy_vector + dE;
	orthonormal_transformation = orthonormal_transformation + orthonormal_transformation.transpose() * dU_prime;
}

void Compute::Perturb(double** dH, int num_points, long num_scale_bins, double* energy_vector, double** orthonormal_transformation)
{
	double** dH_prime0 = new double*[num_points];
	double** dH_prime = new double*[num_points];
	double** dU_prime0 = new double*[num_points];
	double** dU_prime = new double*[num_points];;
	double* dE = new double[num_points];
	// This is to apply the same algorithm but without the Eigen library
}

void Compute::Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, bool eigenvectors_flag, bool perturb_flag)
{
	long num_scale_bins = num_points * num_points;
	Compute::Run(positions, velocities, masses, num_points, dimension, dt, num_scale_bins, eigenvectors_flag, perturb_flag);
}

void Compute::Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag)
{	
	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, positions, num_points, dimension);
	Min_Distance = Distances.Vector[0].Distance;
	Max_Distance = Distances.Vector[Number_Pairs - 1].Distance;

	// Hierarchical Clustering
	Hierarchical_Clusters = new Clusters(Distances);
	//

	double Scale_Distance;
	MatrixXcd Commutator = MatrixXcd::Constant(num_points, num_points, 0);
	for (long Scale_Counter = 0; Scale_Counter < (num_scale_bins + 2); Scale_Counter++)
	{
		Scale_Distance = Scale_Counter * Max_Distance / num_scale_bins; 
		//Scale_Distance = (Hierarchical_Clusters.Max_Vacuum_Scales + Hierarchical_Clusters.Min_Saturation_Scales) / 2;

		Correlation_Operator = (Distances.Operator.array() < Scale_Distance).cast<double>();
		//Correlation_Operator[Scale_Counter] = (Distances.Operator.array() / Scale_Distance).cast<double>();

		Mass_Vector[Scale_Counter] = Correlation_Operator * Vac;
		Mass_Operator = (Mass_Vector[Scale_Counter]).asDiagonal();

		MatrixXd Laplacian_New = Mass_Operator - Correlation_Operator;
		complex<double> I(0, 1);
		Commutator = I * (Mass_Operator * Laplacian_New - Laplacian_New * Mass_Operator);
		//std::cout << "Commutator : " << endl;
		//std::cout << Commutator << endl;
		//std::cout << endl;
		//std::cout << "Commutator (Imag): " << endl;
		//std::cout << Commutator.imag() << endl;
		//std::cout << endl;

		if (perturb_flag)
			Perturb(Laplacian_New - Laplacian[Scale_Counter], num_points, num_scale_bins, Energy_Vector[Scale_Counter], Orthonormal_Transformation[Scale_Counter]);
		else
		{
			if (eigenvectors_flag)
			{
				Laplacian_Eigenstructure.compute(Laplacian_New);
				Orthonormal_Transformation[Scale_Counter] = Laplacian_Eigenstructure.eigenvectors();
			}
			else
			{
				Laplacian_Eigenstructure.compute(Laplacian_New, EigenvaluesOnly);
			}

			MutualInteraction_Eigenstructure.compute(Correlation_Operator, EigenvaluesOnly);
			Commutator_Eigenstructure.compute(Commutator, EigenvaluesOnly);

			Laplacian_Energy[Scale_Counter] = Laplacian_Eigenstructure.eigenvalues().cwiseAbs();
			Commutator_Energy[Scale_Counter] = Commutator_Eigenstructure.eigenvalues().real();
			Energy_Vector[Scale_Counter] = MutualInteraction_Eigenstructure.eigenvalues();
		}

		Laplacian[Scale_Counter] = Laplacian_New;
	}

	//double xx = RadialPower_Force::Force<-1>(10.0);

	// Particle Dynamics :
	// Verlet Method	
	Verlet(Gravitation::Force, positions, velocities, masses, dt, num_points, dimension);
	//
}

Compute::Compute(int num_points, long num_scale_bins, int perturb_order)
{
	Peturb_Order = perturb_order;
	Number_Pairs = num_points * (num_points - 1) / 2;

	Laplacian_Energy = new VectorXd[num_scale_bins + 2];
	Commutator_Energy = new VectorXd[num_scale_bins + 2];
	Mass_Vector = new VectorXd[num_scale_bins + 2];
	Energy_Vector = new VectorXd[num_scale_bins + 2];
	Orthonormal_Transformation = new MatrixXd[num_scale_bins + 2];

	Vac = VectorXd::Constant(num_points, 1);
	Laplacian = new MatrixXd[num_scale_bins + 2];
}

Compute::Compute(int num_points, long num_scale_bins):Compute::Compute(num_points, num_scale_bins, 2) {}

Compute::~Compute()
{
	delete[] Laplacian;
	delete[] Laplacian_Energy;
	delete[] Commutator_Energy;
	delete[] Mass_Vector;
	delete[] Energy_Vector;
	delete[] Orthonormal_Transformation;
	delete Hierarchical_Clusters;
}

int main()
{
	srand(time(NULL));

	bool error_calculation = false;

	const long Number_Points = 7;
	const long Dimension = 2;

	const long Number_Scale_Bins = Number_Points * Number_Points;
	const long Total_Time = 1;
	double dt = 0.1;

	double* Masses = new double[Number_Points];
	for (int i = 0; i < Number_Points; i++)
	{
		Masses[i] = 1.0;
	}
	double** initial_Positions = new double*[Number_Points];
	double** initial_Velocities = new double*[Number_Points];

	double** Positions = new double*[Number_Points];
	double** Positions_Error = new double*[Number_Points];
	double** Positions_Error_Actual = new double*[Number_Points];
	for (int i = 0; i < Number_Points; i++)
	{
		initial_Positions[i] = new double[Dimension];

		Positions[i] = new double[Dimension];
		if (error_calculation)
		{
			Positions_Error[i] = new double[Dimension];
			Positions_Error_Actual[i] = new double[Dimension];
		}
		for (int j = 0; j < Dimension; j++)
		{
			initial_Positions[i][j] = 10.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 5.0;
			Positions[i][j] = initial_Positions[i][j];
		}
	}

	double** Velocities = new double*[Number_Points];
	double** Velocity_Error = new double*[Number_Points];
	double** Velocity_Error_Actual = new double*[Number_Points];
	for (int i = 0; i < Number_Points; i++)
	{
		initial_Velocities[i] = new double[Dimension];

		Velocities[i] = new double[Dimension];
		if (error_calculation)
		{
			Velocity_Error[i] = new double[Dimension];
			Velocity_Error_Actual[i] = new double[Dimension];
		}
		for (int j = 0; j < Dimension; j++)
		{
			initial_Velocities[i][j] = 2.0 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX) - 1.0;
			Velocities[i][j] = initial_Velocities[i][j];
		}
	}

	double *Total_Momentum;
	double Total_Energy;

	// Calculating Forces & Power
	double** Force_Vector = new double*[Number_Points];
	double*** Gravitational_Force_Operator = Force_Operator(Spring::Force, initial_Positions, Number_Points, Dimension, Force_Vector);

	double* Power_Vector = Vector_Product(initial_Velocities, Force_Vector, Number_Points, Dimension);
	// Calculating Particle Dynamics


	// Total Energy & Momentum Calculation
	double* initial_Total_Momentum = Momentum(initial_Velocities, Masses, Number_Points, Dimension);

	double initial_Total_Energy = Potential_Energy(Spring::Energy, initial_Positions, Number_Points, Dimension)
		+ Kinetic_Energy(initial_Velocities, Masses, Number_Points, Dimension);
	//


	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, Positions, Number_Points, Dimension);
	// Hierarchical Clustering
	Clusters Hierarchical_Clusters = Clusters(Distances);

	MatrixXd Scale_Operator = MatrixXd::Constant(Number_Points, Number_Points, 0);
	MatrixXd Scale_Operator_Dual = MatrixXd::Constant(Number_Points, Number_Points, 0);
	//VectorXd Vac = VectorXd::Constant(Number_Points, 1);
	for (int i_p = 0; i_p < (Number_Points - 1); i_p++)
	{
		Scale_Operator(i_p, i_p) = Hierarchical_Clusters.Dendogram->Scale[i_p + 1];
		Scale_Operator_Dual(i_p, i_p) = Hierarchical_Clusters.Dendogram_Dual->Scale[i_p + 1];
	}
	MatrixXd ScaleInv_Distance_Matrix = -Hierarchical_Clusters.Orthogonal_Transformation_Dual.transpose() * Scale_Operator_Dual * Hierarchical_Clusters.Orthogonal_Transformation_Dual;
	for (int i_p = 0; i_p < Number_Points; i_p++)
	{
		ScaleInv_Distance_Matrix(i_p, i_p) = 0;
	}
	//
	std::cout << Scale_Operator << endl;
	std::cout << endl;
	std::cout << Hierarchical_Clusters.Orthogonal_Transformation << endl;
	std::cout << endl;
	std::cout << Scale_Operator_Dual << endl;
	std::cout << endl;
	std::cout << Hierarchical_Clusters.Orthogonal_Transformation_Dual << endl;
	std::cout << endl;
	std::cout << "ScaleInv Distance Matrix : " << endl;
	std::cout << ScaleInv_Distance_Matrix << endl;
	std::cout << endl;
	for (int i_p = 0; i_p < Number_Points; i_p++)
	{
		for (int j_p = 0; j_p < Number_Points; j_p++)
		{
			std::cout << Hierarchical_Clusters.Dendogram->Structure[i_p][j_p] << " ";
		}
		std::cout << ", " << Hierarchical_Clusters.Dendogram->Scale[i_p] << endl;
	}
	std::cout << endl;
	for (int i_p = 0; i_p < Number_Points; i_p++)
	{
		for (int j_p = 0; j_p < Number_Points; j_p++)
		{
			std::cout << Hierarchical_Clusters.Dendogram_Dual->Structure[i_p][j_p] << " ";
		}
		std::cout << ", " << Hierarchical_Clusters.Dendogram_Dual->Scale[i_p] << endl;
	}
	std::cout << endl;

	double** x = new double*[Number_Points];
	int Calculated_Dimension = Embedding_Dimension(ScaleInv_Distance_Matrix, Number_Points, x);

	//
	std::cout << "X : " << std::endl;
	for (int t1 = 0; t1 < Number_Points; t1++)
	{
		for (int t2 = 0; t2 < (Number_Points - 1); t2++)
		{
			std::cout << x[t1][t2] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	//


	Compute Q_Compute = Compute(Number_Points, Number_Scale_Bins, 10);
	bool eigenvectors_flag = false;

	long long computation_time1 = 0;
	auto start_time1 = high_resolution_clock::now();
	
	Q_Compute.Run(Positions, Velocities, Masses, Number_Points, Dimension, dt, Number_Scale_Bins);
	Verlet(Spring::Force, Positions, Velocities, Masses, dt, Number_Points, Dimension);

	auto elapsed_time1 = high_resolution_clock::now() - start_time1;
	computation_time1 = duration_cast<microseconds>(elapsed_time1).count();

	//int median_scale_index = (int)((Q_Compute.Hierarchical_Clusters->Min_Saturation_Scale + Q_Compute.Hierarchical_Clusters->Max_Vacuum_Scale) *  Number_Scale_Bins / Q_Compute.Max_Distance / 2);
	int median_scale_index = Number_Scale_Bins;

	if (eigenvectors_flag)
	{
		std::cout << "Eigenvectors of median scale at t=0 : " << endl;
		for (int i_points = 0; i_points < Number_Points; i_points++)
		{
			std::cout << Q_Compute.Orthonormal_Transformation[median_scale_index].col(i_points).transpose() << " ,  " << Q_Compute.Laplacian_Energy[median_scale_index](i_points) << endl;
		}
	}
	else
	{
		std::cout << "Laplacian Energy of median scale at t=0 : " << endl;
		std::cout << Q_Compute.Laplacian_Energy[median_scale_index] << endl;
	}
	std::cout << endl;
	std::cout << endl;
	std::cout << "Mass Vector of median scale at t=0 : " << endl;
	std::cout << Q_Compute.Mass_Vector[median_scale_index] << endl;
	std::cout << endl;
	std::cout << "Energy Vector of median scale at t=0 : " << endl;
	std::cout << Q_Compute.Energy_Vector[median_scale_index] << endl;
	std::cout << endl;
	std::cout << "Commutator Energy of median scale at t=0 : " << endl;
	std::cout << Q_Compute.Commutator_Energy[median_scale_index] << endl;
	std::cout << endl;

	if (false)
	{
		long long computation_time2 = 0;
		auto start_time2 = high_resolution_clock::now();

		Q_Compute.Run(Positions, Velocities, Masses, Number_Points, Dimension, dt, Number_Scale_Bins, true);
		Verlet(Spring::Force, Positions, Velocities, Masses, dt, Number_Points, Dimension);

		auto elapsed_time2 = high_resolution_clock::now() - start_time2;
		computation_time2 = duration_cast<microseconds>(elapsed_time2).count();

		std::cout << "Eigenvectors of median scale at t=1 : " << endl;
		for (int i_points = 0; i_points < Number_Points; i_points++)
		{
			std::cout << Q_Compute.Orthonormal_Transformation[median_scale_index].col(i_points).transpose() << " ,  " << Q_Compute.Laplacian_Energy[median_scale_index](i_points) << endl;
		}
		std::cout << endl;

		//for (int i_perturb = 0; i_perturb < Q_Compute.Peturb_Order; i_perturb++)
		//{
		//	std::cout << "Eigenvectors of median scale at t=1 (perturbed: " << i_perturb << " ): " << endl;
		//	for (int i_points = 0; i_points < Number_Points; i_points++)
		//	{
		//		double mm = min((Q_Compute.Orthonormal_Transformation_P[i_perturb][median_scale_index].col(i_points).transpose() - Q_Compute.Orthonormal_Transformation[median_scale_index].col(i_points).transpose()).norm(), (Q_Compute.Orthonormal_Transformation_P[i_perturb][median_scale_index].col(i_points).transpose() + Q_Compute.Orthonormal_Transformation[median_scale_index].col(i_points).transpose()).norm());
		//		std::cout << mm << " ,  " << abs(Q_Compute.Energy_Vector_P[i_perturb][median_scale_index](i_points)- Q_Compute.Energy_Vector[median_scale_index](i_points)) / Q_Compute.Energy_Vector[median_scale_index](i_points) << endl;
		//	}
		//	std::cout << endl;
		//}

	}

	// Hexagonal Dome Wave States 
	double Scale_Distance = 1.5;
	Wave_State Hexagon_Waves = Hexagon5(Scale_Distance);
	//

	std::cout << "Min Scale : " << endl;
	for (int t = 0; t < Total_Time; t++)
	{
		std::cout << Q_Compute.Min_Distance << "  ";
	}
	std::cout << endl;
	
	std::cout << "Min NonVacuum Scale : " << endl;
	for (int t = 0; t < Total_Time; t++)
	{
		std::cout << Q_Compute.Hierarchical_Clusters->Max_Vacuum_Scale << "  ";
	}
	std::cout << endl;

	std::cout << "Max NonVacuum Scale : " << endl;
	for (int t = 0; t < Total_Time; t++)
	{
		std::cout << Q_Compute.Hierarchical_Clusters->Min_Saturation_Scale << "  ";
	}
	std::cout << endl;

	std::cout << "Max Scale : " << endl;
	for (int t = 0; t < Total_Time; t++)
	{
		std::cout << Q_Compute.Max_Distance << "  ";
	}
	std::cout << endl;

	std::cout << endl;
	std::cout << "Calculated Dimension : " << Calculated_Dimension << endl;
	std::cout << endl;

	std::cout << endl;
	std::cout << endl;
	std::cout << "Computation Time 1 : " << computation_time1/1e6 << endl;
	std::cout << endl;
	//std::cout << "Computation Time 2 : " << computation_time2 / 1e6 << endl;
	//std::cout << endl;
	//std::cout << "Energies : " << Energy_Vector[Number_Scale_Bins / 2] << endl;

	return 0;
}

