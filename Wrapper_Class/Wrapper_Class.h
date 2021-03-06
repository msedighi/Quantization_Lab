#pragma once
#include "D:\Projects\Quantization_Engine\Quantization_Engine\TriangleTile.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\TriangleTile.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Particle_to_Wave.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Particle_to_Wave.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Clusters.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Clusters.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Particle_Dynamics.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Particle_Dynamics.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Data_Structures.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Data_Structures.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Quantization_Engine.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Quantization_Engine.cpp"

using namespace System;

namespace WrapperClass {
	public ref class Quantization
	{
	public:
		Quantization(int num_points, long num_scale_bins);
		void Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag, bool smooth_flag);

		double Min_Scale;
		double Min_NonVacScale;
		double Max_NonVacScale;
		double Max_Scale;

		bool** Dendogram_Original;
		double* Scale_Ladder_Original;
		bool** Dendogram_Dual;
		double* Scale_Ladder_Dual;

		double** Laplacian_Energy;
		double** Commutator_Energy;
		double** Mass_Vector;
		double** Energy_Vector;
		double*** Laplacian_Orthonormal_Transformation;
		double*** Energy_Orthonormal_Transformation;
		double*** Commutator_Orthonormal_Transformation_Real;
		double*** Commutator_Orthonormal_Transformation_Imag;

		// Classical Variables
		double ClassicalEnergy;
		double ClassicalKineticEnergy;
		double ClassicalPotentialEnergy;

		double* KineticEnergy_Vector;
		double* PotentialEnergy_Vector;
		double* ClassicalLocalEnergy_Vector;
		double* ClassicalEnergy_Exchange;

		double* ClassicalLaplacian_Energy;
		double** ClassicalLaplacian_EigenStates;

		double* ClassicalHamiltonian_Energy;
		double* ClassicalHamiltonian_wVacuum_Energy;
		double** ClassicalHamiltonian_EigenStates;
		double** ClassicalHamiltonian_wVacuum_EigenStates;
		double* ClassicalParticle_Entropy;
		//

		double** Laplacian_Energy_Derivative;
		double** Laplacian_Energy_Derivative_smoothed;

		// Triangle-Tile Waves
		int Num_TileHubs;
		int Tile_Dimension;
		double** Tile_Positions;
		double* Tile_ScalarField;
		double** Tile_VectorField;
		//

	private:
		Compute* Quantization_Engine_Pointer;
		int Num_TileLayers = 5;
	};
}
