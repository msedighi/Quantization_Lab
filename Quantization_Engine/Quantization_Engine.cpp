// Quantization_Engine.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>

#include "Quantization_Engine.h"
#include "Particle_to_Wave.h"

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

void Compute::Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag, bool smooth_flag)
{	
	//InverseRoot* interaction = new InverseRoot();
	//Spring* interaction = new Spring();
	//Gravitation* interaction = new Gravitation();
	//Lennard_Jones* interaction = new Lennard_Jones();
	Logarithmic* interaction = new Logarithmic();
	Compute::Run(interaction, positions, velocities, masses, num_points, dimension, dt, num_scale_bins, eigenvectors_flag, perturb_flag, smooth_flag);
}

void Compute::Run(Interaction* interaction, double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag, bool smooth_flag)
{	
	Distance_Struct Distances = Distance_Matrix(Euclidean_Distance, positions, num_points, dimension);
	Min_Distance = Distances.Vector[0].Distance;
	Max_Distance = Distances.Vector[Number_Pairs - 1].Distance;

	// Hierarchical Clustering
	Hierarchical_Clusters = new Clusters(Distances);
	//

	ClassicalKineticEnergy = 0;
	KineticEnergy_Vector = Kinetic_Energy(velocities, masses, num_points, dimension);
	ClassicalEnergy_Exchange = Energy_Exchange(interaction, positions, velocities, num_points, dimension);
	for (int i_p = 0; i_p < num_points; i_p++)
	{
		ClassicalKineticEnergy += KineticEnergy_Vector[i_p] / 2.0;
		ClassicalEnergy_Hamiltonian(i_p, i_p) = KineticEnergy_Vector[i_p];
		// Off-diagonal
		for (int j_p = 0; j_p < num_points; j_p++)
		{
			if (i_p != j_p)
			{
				ClassicalEnergy_Hamiltonian(i_p, j_p) = interaction->Energy(Distances.Operator(i_p, j_p));
				PotentialEnergy_Operator(i_p, j_p) = ClassicalEnergy_Hamiltonian(i_p, j_p);
			}
		}
	}
	PotentialEnergy_Vector = PotentialEnergy_Operator * Vac;
	PotentialEnergy_Laplacian = (MatrixXd)PotentialEnergy_Vector.asDiagonal() - PotentialEnergy_Operator;

	//
	ClassicalEnergy = ClassicalEnergy_Hamiltonian.sum() / 2.0;
	ClassicalPotentialEnergy = PotentialEnergy_Operator.sum() / 2.0;
	// Computing Classical Eigen-Structure
	ClassicalLaplacian_Eigenstructure.compute(PotentialEnergy_Laplacian);
	ClassicalLaplacian_Energy = ClassicalLaplacian_Eigenstructure.eigenvalues();
	ClassicalLaplacian_EigenStates = ClassicalLaplacian_Eigenstructure.eigenvectors();

	ClassicalHamiltonian_Eigenstructure.compute(ClassicalEnergy_Hamiltonian);
	ClassicalHamiltonian_Energy = ClassicalHamiltonian_Eigenstructure.eigenvalues();
	ClassicalHamiltonian_EigenStates = ClassicalHamiltonian_Eigenstructure.eigenvectors();
	//

	// Global Energy Exchange 
	GlobalEnergy_Exchange = CollectiveEnergy_Exchange(interaction, positions, velocities, num_points, dimension, Hierarchical_Clusters->Orthogonal_Transformation);

	Tile_Waves = Tile_Interface->Run(GlobalEnergy_Exchange);

	// NEED TO INCLUDE TRIANGLE TILE IN OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

	//

	if (debug_flag)
	{
		std::cout << "Hamiltonian: " << endl;
		std::cout << ClassicalEnergy_Hamiltonian << endl;
		std::cout << endl;
		std::cout << "ClASSICAL ENERGY : " << endl;
		std::cout << ClassicalEnergy << endl;
		std::cout << endl;
		std::cout << "ClASSICAL ENERGY Hamiltonian: " << endl;
		std::cout << ClassicalHamiltonian_Energy << endl;
		std::cout << endl;
		std::cout << "ClASSICAL ENERGY STATES: " << endl;
		std::cout << ClassicalHamiltonian_EigenStates << endl;
		std::cout << endl;
	}

	//

	//#pragma loop(hint_parallel(6))
	for (long Scale_Counter = 0; Scale_Counter < num_scale_bins; Scale_Counter++)
	{
		//Scale_Distance = Scale_Counter * Max_Distance / num_scale_bins; 
		//Scale_Distance = (Hierarchical_Clusters.Max_Vacuum_Scales + Hierarchical_Clusters.Min_Saturation_Scales) / 2;

		for (int i_p=0; i_p < num_points; i_p++)
			for (int j_p = 0; j_p < num_points; j_p++)
			{
				Correlation_Operator(i_p, j_p) = StepFunc_3((double)Scale_Counter / (double)num_scale_bins, Distances.Operator(i_p, j_p) / (Max_Distance + Min_Distance));
			}
		//Correlation_Operator = (Distances.Operator.array() <= Scale_Distance).cast<double>();
		//Correlation_Operator[Scale_Counter] = (Distances.Operator.array() / Scale_Distance).cast<double>();

		Mass_Vector[Scale_Counter] = Correlation_Operator * Vac;
		Mass_Operator = (Mass_Vector[Scale_Counter]).asDiagonal();

		Laplacian_New = Mass_Operator - Correlation_Operator;

		complex<double> I(0, 1);
		Commutator = I * (Mass_Operator * Laplacian_New - Laplacian_New * Mass_Operator);

		// Debug!
		if (debug_flag)
		{
			//std::cout << "Laplacian Matrix : " << endl;
			//std::cout << Laplacian_New << endl;
			//std::cout << endl;
		}

		if (perturb_flag)
			Perturb(Laplacian_New - Laplacian[Scale_Counter], num_points, num_scale_bins, Energy_Vector[Scale_Counter], Laplacian_Orthonormal_Transformation[Scale_Counter]);
		else
		{
			if (eigenvectors_flag)
			{
				Laplacian_Eigenstructure.compute(Laplacian_New);
				MutualInteraction_Eigenstructure.compute(Correlation_Operator);
				Commutator_Eigenstructure.compute(Commutator);

				Laplacian_Energy[Scale_Counter] = Laplacian_Eigenstructure.eigenvalues().cwiseAbs();
				Commutator_Energy[Scale_Counter] = Commutator_Eigenstructure.eigenvalues().real();
				Energy_Vector[Scale_Counter] = MutualInteraction_Eigenstructure.eigenvalues();

				Laplacian_Multiplicity[Scale_Counter] = Compute_Multiplicity(Laplacian_Energy[Scale_Counter], num_points, 0.04);
				// Debug!
				if (debug_flag)
				{
					if (Laplacian_Multiplicity[Scale_Counter].Num_Degenerate_States > 0)
					{
						std::cout << "Laplacian Multiplicity @ " << Scale_Counter << " : " << endl;
						std::cout << Laplacian_Energy[Scale_Counter] << endl;
						std::cout << endl;
						std::cout << "Num Degenerate States @ " << Scale_Counter << " : " << endl;
						std::cout << Laplacian_Multiplicity[Scale_Counter].Num_Degenerate_States << endl;
						std::cout << endl;
						for (int i_m = 0; i_m < Laplacian_Multiplicity[Scale_Counter].Num_Degenerate_States; i_m++)
						{
							std::cout << "index : " << endl;
							std::cout << Laplacian_Multiplicity[Scale_Counter].initial_index[i_m] << " to " << Laplacian_Multiplicity[Scale_Counter].final_index[i_m] << " : " << Laplacian_Multiplicity[Scale_Counter].degeneracy[i_m] << endl;
							std::cout << endl;
						}
					}
				}
				//

				Laplacian_Orthonormal_Transformation[Scale_Counter] = Laplacian_Eigenstructure.eigenvectors();
				Energy_Orthonormal_Transformation[Scale_Counter] = MutualInteraction_Eigenstructure.eigenvectors();
				Commutator_Orthonormal_Transformation[Scale_Counter] = Commutator_Eigenstructure.eigenvectors();

				Commutator_Orthonormal_Transformation_Real[Scale_Counter] = Commutator_Orthonormal_Transformation[Scale_Counter].real();
				Commutator_Orthonormal_Transformation_Imag[Scale_Counter] = Commutator_Orthonormal_Transformation[Scale_Counter].imag();
			}
			else
			{
				Laplacian_Eigenstructure.compute(Laplacian_New, EigenvaluesOnly);
				MutualInteraction_Eigenstructure.compute(Correlation_Operator, EigenvaluesOnly);
				Commutator_Eigenstructure.compute(Commutator, EigenvaluesOnly);

				Laplacian_Energy[Scale_Counter] = Laplacian_Eigenstructure.eigenvalues().cwiseAbs();
				Commutator_Energy[Scale_Counter] = Commutator_Eigenstructure.eigenvalues().real();
				Energy_Vector[Scale_Counter] = MutualInteraction_Eigenstructure.eigenvalues();

				Laplacian_Multiplicity[Scale_Counter] = Compute_Multiplicity(Laplacian_Energy[Scale_Counter], num_points, 0.04);
			}
		}


 		Laplacian[Scale_Counter] = Laplacian_New;
	}

	//double xx = RadialPower_Force::Force<-1>(10.0);

	// Smoothing
	if (smooth_flag)
	{
		Compute::Smooth(Laplacian_Energy, Laplacian_Orthonormal_Transformation, Laplacian_Multiplicity, num_points, num_scale_bins);
		//Compute::Smooth(Energy_Vector, Energy_Orthonormal_Transformation, num_points, num_scale_bins);
		//Compute::Smooth(Commutator_Energy, Commutator_Orthonormal_Transformation, num_points, num_scale_bins);

		// Derivative Computation
		for (long Scale_Counter = 1; Scale_Counter < num_scale_bins; Scale_Counter++)
		{
			Laplacian_Energy_Derivative[Scale_Counter - 1] = Laplacian_Energy[Scale_Counter] - Laplacian_Energy[Scale_Counter - 1];
			Laplacian_Energy_Derivative_smoothed[Scale_Counter - 1] = (Laplacian_Orthonormal_Transformation[Scale_Counter - 1].transpose().conjugate() * (Laplacian[Scale_Counter] - Laplacian[Scale_Counter - 1]) * Laplacian_Orthonormal_Transformation[Scale_Counter - 1]).diagonal();
			// Debug!
			//std::cout << "Laplacian Energy Derivative : " << endl;
			//std::cout << Laplacian_Energy_Derivative[Scale_Counter - 1] << endl;
			//std::cout << endl;
			//std::cout << "Laplacian Energy Derivative Error : " << endl;
			//std::cout << (Laplacian_Energy_Derivative_smoothed[Scale_Counter - 1] - Laplacian_Energy_Derivative[Scale_Counter - 1]).cwiseAbs() << endl;
			//std::cout << endl;
		}

	}


	// Particle Dynamics :
	// Verlet Method	
	Verlet(interaction, positions, velocities, masses, dt, num_points, dimension);
	//
}

void Compute::Smooth(Eigen::VectorXd* values, Eigen::MatrixXd* vectors, Multiplicity* values_multiplicity, int num_points, int num_scale_bins)
{
	int start_index = num_scale_bins / 2;
	int off_set = 0;
	while (values_multiplicity[start_index + off_set].Num_Degenerate_States > 0)
	{
		off_set++;
		if (start_index >= off_set)
		{
			if (values_multiplicity[start_index - off_set].Num_Degenerate_States == 0)
			{
				start_index = start_index - off_set;
				break;
			}
		}

		if (start_index + off_set >= num_scale_bins)
			break;
	}
	if (start_index == num_scale_bins / 2)
		start_index = start_index + off_set;

	bool ColSumTest_flag = false;
	bool RowSumTest_flag = false;

	MatrixXd Identity_Close_Prod = MatrixXd::Identity(num_points, num_points);
	for (int i_s = start_index; i_s < num_scale_bins; i_s++)
	{
		MatrixXd Identity_Close = MatrixXd::Constant(num_points, num_points, 0);

		for (int i_d = 0; i_d < values_multiplicity[i_s].Num_Degenerate_States; i_d++)
		{
			for (int i = values_multiplicity[i_s].initial_index[i_d]; i <= values_multiplicity[i_s].final_index[i_d]; i++)
			{
				vectors[i_s].col(i) = vectors[i_s - 1].col(i);
			}
		}

		if (i_s > start_index)
		{
			MatrixXd Id = (vectors[i_s - 1].transpose() * vectors[i_s]);
			MatrixXd Id2 = Id.cwiseProduct(Id);

			for (int i = 0; i < num_points; i++)
				for (int j = 0; j < num_points; j++)
					if (Id2(i, j) == Id2.row(i).maxCoeff())
					{
						Identity_Close(i, j) = Id(i, j) / abs(Id(i, j));
					}
				
			//Debug!
			std::cout << " INEEEEEEEEEEEEEE @ Scale_Counter = " << i_s - 1 << endl;
			std::cout << (10 * Id).array().round() / 10 << endl;
			std::cout << endl;
			std::cout << (10 * Id2).array().round() / 10 << endl;
			std::cout << endl;
			std::cout << Identity_Close << endl;
			std::cout << endl;
			std::cout << Identity_Close_Prod << endl;
			std::cout << endl;
			std::cout << endl;
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s - 2].col(i_points).transpose() << " ,  " << values[i_s - 2](i_points) << endl;
			}
			std::cout << endl;
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s - 1].col(i_points).transpose() << " ,  " << values[i_s - 1](i_points) << endl;
			}
			std::cout << endl;
			//


			ArrayXd col_test = Identity_Close.cwiseAbs().colwise().sum().array();
			ArrayXd row_test = Identity_Close.cwiseAbs().rowwise().sum().array();

			ColSumTest_flag = (col_test.minCoeff() == 1) && (col_test.maxCoeff() == 1);
			RowSumTest_flag = (row_test.minCoeff() == 1) && (row_test.maxCoeff() == 1);

			if (!ColSumTest_flag || !RowSumTest_flag)
			{
				std::cout << " Error @ Scale_Counter = " << i_s << endl;
				for (int i_points = 0; i_points < num_points; i_points++)
				{
					std::cout << vectors[i_s].col(i_points).transpose() << " ,  " << values[i_s](i_points) << endl;
				}
				std::cout << endl;
			}

			vectors[i_s - 1] = (vectors[i_s - 1] * Identity_Close_Prod.transpose()).eval();
			values[i_s - 1] = (Identity_Close_Prod.cwiseAbs() * values[i_s - 1]).eval();

			Identity_Close_Prod = (Identity_Close_Prod * Identity_Close).eval();

			//
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s - 1].col(i_points).transpose() << " ,  " << values[i_s - 1](i_points) << endl;
			}
			std::cout << endl;
			//

			// TEST!!
			MatrixXd Id_Test = (vectors[i_s - 2].transpose() * vectors[i_s - 1]);
			MatrixXd Id2_Test = Id_Test.cwiseProduct(Id_Test);
			MatrixXd Identity_Close_Test = MatrixXd::Constant(num_points, num_points, 0);
			for (int i = 0; i < num_points; i++)
				for (int j = 0; j < num_points; j++)
					if (Id2_Test(i, j) == Id2_Test.row(i).maxCoeff())
					{
						Identity_Close_Test(i, j) = Id_Test(i, j) / abs(Id_Test(i, j));
					}

			if (!Identity_Close_Test.isIdentity() && (i_s > (start_index + 1)))
			{
				std::cout << Identity_Close_Test << endl;
				std::cout << endl;
			}
			//

		}
		
	}

	Identity_Close_Prod = MatrixXd::Identity(num_points, num_points);
	for (int i_s = start_index; i_s >= 0; i_s--)
	{
		MatrixXd Identity_Close = MatrixXd::Constant(num_points, num_points, 0);

		for (int i_d = 0; i_d < values_multiplicity[i_s].Num_Degenerate_States; i_d++)
		{
			for (int i = values_multiplicity[i_s].initial_index[i_d]; i <= values_multiplicity[i_s].final_index[i_d]; i++)
			{
				vectors[i_s].col(i) = vectors[i_s + 1].col(i);
			}
		}

		if (i_s < start_index)
		{
			MatrixXd Id = (vectors[i_s + 1].transpose() * vectors[i_s]);
			MatrixXd Id2 = Id.cwiseProduct(Id);

			for (int i = 0; i < num_points; i++)
				for (int j = 0; j < num_points; j++)
					if (Id2(i, j) == Id2.row(i).maxCoeff())
					{
						Identity_Close(i, j) = Id(i, j) / abs(Id(i, j));
						// take a note of non deg states on the other side i.e. collect j's here!!
					}


			//Debug!
			std::cout << " INEEEEEEEEEEEEEE @Scale_Counter = " << i_s << endl;
			std::cout << (10 * Id).array().round() / 10 << endl;
			std::cout << endl;
			std::cout << (10 * Id2).array().round() / 10 << endl;
			std::cout << endl;
			std::cout << Identity_Close << endl;
			std::cout << endl;
			std::cout << Identity_Close_Prod << endl;
			std::cout << endl;
			std::cout << endl;
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s + 2].col(i_points).transpose() << " ,  " << values[i_s + 2](i_points) << endl;
			}
			std::cout << endl;
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s + 1].col(i_points).transpose() << " ,  " << values[i_s + 1](i_points) << endl;
			}
			std::cout << endl;
			//


			ArrayXd col_test = Identity_Close.cwiseAbs().colwise().sum().array();
			ArrayXd row_test = Identity_Close.cwiseAbs().rowwise().sum().array();

			ColSumTest_flag = (col_test.minCoeff() == 1) && (col_test.maxCoeff() == 1);
			RowSumTest_flag = (row_test.minCoeff() == 1) && (row_test.maxCoeff() == 1);

			if (!ColSumTest_flag || !RowSumTest_flag)
			{
				std::cout << " Error @ Scale_Counter = " << i_s << endl;
				for (int i_points = 0; i_points < num_points; i_points++)
				{
					std::cout << vectors[i_s].col(i_points).transpose() << " ,  " << values[i_s](i_points) << endl;
				}
				std::cout << endl;
			}

			vectors[i_s + 1] = (vectors[i_s + 1] * Identity_Close_Prod.transpose()).eval();
			values[i_s + 1] = (Identity_Close_Prod.cwiseAbs() * values[i_s + 1]).eval();

			Identity_Close_Prod = (Identity_Close_Prod * Identity_Close).eval();
			//
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << vectors[i_s + 1].col(i_points).transpose() << " ,  " << values[i_s + 1](i_points) << endl;
			}
			std::cout << endl;
			//

			// TEST!!
			MatrixXd Id_Test = (vectors[i_s + 2].transpose() * vectors[i_s + 1]);
			MatrixXd Id2_Test = Id_Test.cwiseProduct(Id_Test);
			MatrixXd Identity_Close_Test = MatrixXd::Constant(num_points, num_points, 0);
			for (int i = 0; i < num_points; i++)
				for (int j = 0; j < num_points; j++)
					if (Id2_Test(i, j) == Id2_Test.row(i).maxCoeff())
					{
						Identity_Close_Test(i, j) = Id_Test(i, j) / abs(Id_Test(i, j));
					}

			if (!Identity_Close_Test.isIdentity())
			{
				std::cout << Identity_Close_Test << endl;
				std::cout << endl;
			}
			//

		}

	}


}

void Compute::Smooth(Eigen::VectorXd* values, Eigen::MatrixXd* vectors, int scale_counter, int num_points)
{
	MatrixXd Identity_Close = MatrixXd::Constant(num_points, num_points, 0);
	bool ColSumTest_flag = false;
	bool RowSumTest_flag = false;

	MatrixXd Id = (vectors[scale_counter - 1].transpose() * vectors[scale_counter]);
	MatrixXd Id2 = Id.cwiseProduct(Id);
	for (int i = 0; i < num_points; i++)
		for (int j = 0; j < num_points; j++)
			if (Id2(i, j) == Id2.row(i).maxCoeff())
				Identity_Close(i, j) = Id(i, j) / abs(Id(i, j));

	ArrayXd col_test = Identity_Close.cwiseAbs().colwise().sum().array();
	ArrayXd row_test = Identity_Close.cwiseAbs().rowwise().sum().array();

	ArrayXb col_flag = (col_test == col_test.maxCoeff());
	ArrayXb row_flag = (row_test == row_test.maxCoeff());
	ColSumTest_flag = (col_test.minCoeff() == 1) && (col_test.maxCoeff() == 1);
	RowSumTest_flag = (row_test.minCoeff() == 1) && (row_test.maxCoeff() == 1);

	if (ColSumTest_flag && RowSumTest_flag)
	{
		vectors[scale_counter] = (vectors[scale_counter] * Identity_Close.transpose()).eval();
		values[scale_counter] = (Identity_Close.cwiseAbs() * values[scale_counter]).eval();
	}
	else
	{
		//Debug!
		std::cout << " INEEEEEEEEEEEEEE @Scale_Counter = " << scale_counter << endl;
		std::cout << (10 * Id).array().round() / 10 << endl;
		std::cout << endl;
		std::cout << (10 * Id2).array().round() / 10 << endl;
		std::cout << endl;
		std::cout << Identity_Close << endl;
		std::cout << endl;
		std::cout << endl;
		for (int i_points = 0; i_points < num_points; i_points++)
		{
			std::cout << Laplacian_Orthonormal_Transformation[scale_counter].col(i_points).transpose() << " ,  " << Laplacian_Energy[scale_counter](i_points) << endl;
		}
		std::cout << endl;
		// Fix
		if ((row_test.minCoeff() == 1) && (col_test.minCoeff() == 1))
		{
			int col_flag_counter = 0;
			for (int i = 0; i < num_points; i++)
			{
				if (col_flag(i))
				{
					col_flag_counter++;
					int row_flag_counter = 0;
					for (int j = 0; j < num_points; j++)
					{
						if (row_flag(j))
						{
							row_flag_counter++;

							if (row_flag_counter == col_flag_counter)
								Identity_Close(j, i) = 1;
							else
								Identity_Close(j, i) = 0;
						}
					}
				}
			}
		}
		else
		{ 
			for (int i_points = 0; i_points < num_points; i_points++)
			{
				std::cout << Laplacian_Orthonormal_Transformation[scale_counter - 1].col(i_points).transpose() << " ,  " << Laplacian_Energy[scale_counter - 1](i_points) << endl;
			}
			std::cout << endl;
		}
		std::cout << endl;
		std::cout << Identity_Close << endl;
		std::cout << endl;

	}
}

void Compute::Smooth(Eigen::VectorXd* values, Eigen::MatrixXcd* vectors, int scale_counter, int num_points)
{
	MatrixXd Identity_Close = MatrixXd::Constant(num_points, num_points, 0);
	bool ColSumTest_flag = false;
	bool RowSumTest_flag = false;

	//MatrixXcd Id = (vectors[scale_counter - 1].conjugate().transpose() * vectors[scale_counter]);
	//MatrixXcd Id2 = Id.cwiseProduct(Id.conjugate);
	MatrixXd Id = (vectors[scale_counter - 1].conjugate().transpose() * vectors[scale_counter]).real(); // This is not always right! The above should be corrected & replaced.
	MatrixXd Id2 = Id.cwiseProduct(Id);
	for (int i = 0; i < num_points; i++)
		for (int j = 0; j < num_points; j++)
			if (Id2(i, j) == Id2.row(i).maxCoeff())
				Identity_Close(i, j) = Id(i, j) / abs(Id(i, j));

	VectorXd col_test = Identity_Close.cwiseAbs().colwise().sum();
	VectorXd row_test = Identity_Close.cwiseAbs().rowwise().sum();
	ColSumTest_flag = (col_test.minCoeff() == 1) && (col_test.maxCoeff() == 1);
	RowSumTest_flag = (row_test.minCoeff() == 1) && (row_test.maxCoeff() == 1);

	if (ColSumTest_flag && RowSumTest_flag)
	{
		vectors[scale_counter] = (vectors[scale_counter] * Identity_Close.transpose()).eval();
		values[scale_counter] = (Identity_Close.cwiseAbs() * values[scale_counter]).eval();
	}
	else
	{
		////Debug!
		//std::cout << " INEEEEEEEEEEEEEE @Scale_Counter = " << scale_counter << endl;
		//std::cout << (10 * Id).array().round() / 10 << endl;
		//std::cout << endl;
		//std::cout << (10 * Id2).array().round() / 10 << endl;
		//std::cout << endl;
		//std::cout << Identity_Close << endl;
		//std::cout << endl;
		//std::cout << endl;
		//for (int i_points = 0; i_points < num_points; i_points++)
		//{
		//	std::cout << vectors[scale_counter].col(i_points).transpose() << " ,  " << values[scale_counter](i_points) << endl;
		//}
		//std::cout << endl;
	}
}

Compute::Compute(int num_points, long num_scale_bins, int num_tilelayers, int perturb_order)
{
	Peturb_Order = perturb_order;
	Number_Pairs = num_points * (num_points - 1) / 2;

	Laplacian_Energy = new VectorXd[num_scale_bins];
	Laplacian_Energy_Derivative = new VectorXd[num_scale_bins - 1];
	Laplacian_Energy_Derivative_smoothed = new VectorXd[num_scale_bins - 1];

	Commutator_Energy = new VectorXd[num_scale_bins];
	Mass_Vector = new VectorXd[num_scale_bins];
	Energy_Vector = new VectorXd[num_scale_bins];
	Laplacian_Orthonormal_Transformation = new MatrixXd[num_scale_bins];
	Energy_Orthonormal_Transformation = new MatrixXd[num_scale_bins];
	Commutator_Orthonormal_Transformation = new MatrixXcd[num_scale_bins];
	Commutator_Orthonormal_Transformation_Real = new MatrixXd[num_scale_bins];
	Commutator_Orthonormal_Transformation_Imag = new MatrixXd[num_scale_bins];
	//
	Correlation_Operator = MatrixXd::Constant(num_points, num_points, 0);
	Commutator = MatrixXcd::Constant(num_points, num_points, 0);
	//
	ClassicalEnergy_Hamiltonian = MatrixXd::Constant(num_points, num_points, 0);
	PotentialEnergy_Operator = MatrixXd::Constant(num_points, num_points, 0);


	Laplacian_Multiplicity = new Multiplicity[num_scale_bins];

	Vac = VectorXd::Constant(num_points, 1);
	Laplacian = new MatrixXd[num_scale_bins];

	Tile_Interface = new WaveParticle_Interface(num_tilelayers);
	Num_TileHubs = Tile_Interface->Tile->Num_Points;
	Tile_Dimension = Tile_Interface->Tile->dimension;
	Tile_Positions = Tile_Interface->Tile->positions;
}

Compute::Compute(int num_points, long num_scale_bins, int num_tilelayers) :Compute::Compute(num_points, num_scale_bins, num_tilelayers, 2) {}

Compute::Compute(int num_points, long num_scale_bins):Compute::Compute(num_points, num_scale_bins, 4) {}

Compute::~Compute()
{
	delete[] Laplacian;
	delete[] Laplacian_Energy;
	delete[] Laplacian_Energy_Derivative_smoothed;
	delete[] Laplacian_Energy_Derivative;
	delete[] Commutator_Energy;
	delete[] Mass_Vector;
	delete[] Energy_Vector;
	delete[] Laplacian_Orthonormal_Transformation;
	delete[] Energy_Orthonormal_Transformation;
	delete[] Commutator_Orthonormal_Transformation;
	delete[] Commutator_Orthonormal_Transformation_Real;
	delete[] Commutator_Orthonormal_Transformation_Imag;
	delete[] Laplacian_Multiplicity;
	delete[] KineticEnergy_Vector;
	delete[] ClassicalEnergy_Exchange;
	delete Hierarchical_Clusters;
	delete Tile_Interface;

	for (int i = 0; i < Num_TileHubs; i++)
		delete[] Tile_Positions[i];
	delete[] Tile_Positions;
}

int main()
{
	srand(time(NULL));

	bool error_calculation = false;

	const long Number_Points = 8;
	const long Dimension = 2;

	const long Number_Scale_Bins = Number_Points * Number_Points * 10;
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
	double*** Gravitational_Force_Operator = Force_Operator(new Gravitation(), initial_Positions, Number_Points, Dimension, Force_Vector);

	double* Power_Vector = Vector_Product(initial_Velocities, Force_Vector, Number_Points, Dimension);
	// Calculating Particle Dynamics


	// Total Energy & Momentum Calculation
	double* initial_Total_Momentum = Momentum(initial_Velocities, Masses, Number_Points, Dimension);

	//double kinetic_energy;
	//Kinetic_Energy(initial_Velocities, Masses, Number_Points, Dimension, kinetic_energy);
	//double initial_Total_Energy = Potential_Energy(new Gravitation(), initial_Positions, Number_Points, Dimension) + kinetic_energy;
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
	if (false)
	{
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
	}

	double** x = new double*[Number_Points];
	int Calculated_Dimension = Embedding_Dimension(ScaleInv_Distance_Matrix, Number_Points, x);
	if (false)
	{
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
	}


	Compute Q_Compute = Compute(Number_Points, Number_Scale_Bins, 10);
	bool eigenvectors_flag = true;
	bool smooth_flag = false;
	bool perturb_flag = false;

	long long computation_time1 = 0;
	auto start_time1 = high_resolution_clock::now();
	
	Q_Compute.Run(new Spring(), Positions, Velocities, Masses, Number_Points, Dimension, dt, Number_Scale_Bins, eigenvectors_flag, perturb_flag, smooth_flag);
	//Verlet(Spring::Force, Positions, Velocities, Masses, dt, Number_Points, Dimension);

	auto elapsed_time1 = high_resolution_clock::now() - start_time1;
	computation_time1 = duration_cast<microseconds>(elapsed_time1).count();

	int median_scale_index = (int)((Q_Compute.Hierarchical_Clusters->Min_Saturation_Scale + Q_Compute.Hierarchical_Clusters->Max_Vacuum_Scale) *  Number_Scale_Bins / Q_Compute.Max_Distance / 2);
	//int median_scale_index = Number_Scale_Bins;

	if (eigenvectors_flag)
	{
		MatrixXd midEnergy_EigenVectors = MatrixXd::Constant(Number_Points, Number_Points, 0);
		midEnergy_EigenVectors.col(0) = Q_Compute.Laplacian_Orthonormal_Transformation[Number_Scale_Bins - 1].col(0);

		double half_Energy = (double)Number_Points / 2;
		for (int i_p = 1; i_p < Number_Points; i_p++)
		{
			double residual = INFINITY;
			for (long i_s = 0; i_s < Number_Scale_Bins; i_s++)
			{
				if ((residual < abs(half_Energy - Q_Compute.Laplacian_Energy[i_s](i_p))) && (residual < .01))
				{
					midEnergy_EigenVectors.col(i_p) = Q_Compute.Laplacian_Orthonormal_Transformation[i_s - 1].col(i_p);

					if (false)
					{
						std::cout << "midEnergy EigenValue Test : " << endl;
						std::cout << half_Energy << ", " << Q_Compute.Laplacian_Energy[i_s - 1](i_p) << ", " << (i_s - 1) << endl;
					}
					break;
				}
				else
				{
					residual = abs(half_Energy - Q_Compute.Laplacian_Energy[i_s](i_p));
				}
			}
		}

		if (false)
		{
			std::cout << endl;
			std::cout << "midEnergy Eigen Vectors : " << endl;
			std::cout << midEnergy_EigenVectors << endl;
			std::cout << "midEnergy Orthonormal Check : " << endl;
			std::cout << midEnergy_EigenVectors.transpose() * midEnergy_EigenVectors << endl;


			std::cout << "Eigenvectors of median scale at t=0 : " << endl;
			std::cout << Q_Compute.Laplacian_Orthonormal_Transformation[median_scale_index] << endl;
			//std::cout << "Commutator Eigenvectors of median scale at t=0 : " << endl;
			//std::cout << Q_Compute.Commutator_Orthonormal_Transformation[median_scale_index] << endl;
			//for (int i_points = 0; i_points < Number_Points; i_points++)
			//{
				//std::cout << Q_Compute.Laplacian_Orthonormal_Transformation[median_scale_index].col(i_points).transpose() << " ,  " << Q_Compute.Laplacian_Energy[median_scale_index](i_points) << endl;
			//}
		}
	}
	if (false)
	{
		std::cout << endl;
		std::cout << "Laplacian Energy of median scale at t=0 : " << endl;
		std::cout << Q_Compute.Laplacian_Energy[median_scale_index] << endl;
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
	}

	if (false)
	{
		long long computation_time2 = 0;
		auto start_time2 = high_resolution_clock::now();

		Q_Compute.Run(new Gravitation(), Positions, Velocities, Masses, Number_Points, Dimension, dt, Number_Scale_Bins, eigenvectors_flag, perturb_flag, smooth_flag);
		Verlet(new Spring(), Positions, Velocities, Masses, dt, Number_Points, Dimension);

		auto elapsed_time2 = high_resolution_clock::now() - start_time2;
		computation_time2 = duration_cast<microseconds>(elapsed_time2).count();

		if (false)
		{
			std::cout << "Eigenvectors of median scale at t=1 : " << endl;
			for (int i_points = 0; i_points < Number_Points; i_points++)
			{
				std::cout << Q_Compute.Laplacian_Orthonormal_Transformation[median_scale_index].col(i_points).transpose() << " ,  " << Q_Compute.Laplacian_Energy[median_scale_index](i_points) << endl;
			}
			std::cout << endl;

			//for (int i_perturb = 0; i_perturb < Q_Compute.Peturb_Order; i_perturb++)
			//{
			//	std::cout << "Eigenvectors of median scale at t=1 (perturbed: " << i_perturb << " ): " << endl;
			//	for (int i_points = 0; i_points < Number_Points; i_points++)
			//	{
			//		double mm = min((Q_Compute.Laplacian_Orthonormal_Transformation_P[i_perturb][median_scale_index].col(i_points).transpose() - Q_Compute.Laplacian_Orthonormal_Transformation[median_scale_index].col(i_points).transpose()).norm(), (Q_Compute.Laplacian_Orthonormal_Transformation_P[i_perturb][median_scale_index].col(i_points).transpose() + Q_Compute.Laplacian_Orthonormal_Transformation[median_scale_index].col(i_points).transpose()).norm());
			//		std::cout << mm << " ,  " << abs(Q_Compute.Energy_Vector_P[i_perturb][median_scale_index](i_points)- Q_Compute.Energy_Vector[median_scale_index](i_points)) / Q_Compute.Energy_Vector[median_scale_index](i_points) << endl;
			//	}
			//	std::cout << endl;
			//}
		}
	}

	// Hexagonal Dome Wave States 
	double noise = 0.005;
	double Scale_Distance = 1.1;
	TriangleTile Tile = TriangleTile(4, noise, false);
	TriangleTile Tile_Noise = TriangleTile(4, noise, true);

	Wave_State Hexagon_Waves = Tile.Hexagon_Wave(Scale_Distance);
	Wave_State Hexagon_Waves_noise = Tile_Noise.Hexagon_Wave(Scale_Distance);
	Wave_State Hexagon_Waves0 = Tile.Hexagon_Wave();
	Wave_State Hexagon_Waves0_noise = Tile_Noise.Hexagon_Wave();
	//


	if (false)
	{
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
		std::cout << "Computation Time 1 : " << computation_time1 / 1e6 << endl;
		std::cout << endl;
		//std::cout << "Computation Time 2 : " << computation_time2 / 1e6 << endl;
	}

	std::cout << endl;
	std::cout << "Hexagon Waves w/o Noise Energies : " << Hexagon_Waves.Energy_Vector << endl;
	std::cout << endl;
	std::cout << "Hexagon Waves w Noise Energies : " << Hexagon_Waves_noise.Energy_Vector << endl;
	std::cout << endl;
	std::cout << "Hexagon Waves at MidScale w/o Noise Energies : " << Hexagon_Waves0.Energy_Vector << endl;
	std::cout << endl;
	std::cout << "Hexagon Waves at MidScale w Noise Energies : " << Hexagon_Waves0_noise.Energy_Vector << endl;

	return 0;
}

