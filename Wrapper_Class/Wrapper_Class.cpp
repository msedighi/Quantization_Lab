#include "stdafx.h"

#include "Wrapper_Class.h"

WrapperClass::Quantization::Quantization(int num_points, long num_scale_bins)
{
	Quantization_Engine_Pointer = new Compute(num_points, num_scale_bins, Num_TileLayers);

	Laplacian_Energy = new double*[num_scale_bins];
	Commutator_Energy = new double*[num_scale_bins];
	Energy_Vector = new double*[num_scale_bins];
	Mass_Vector = new double*[num_scale_bins];
	Laplacian_Orthonormal_Transformation = new double**[num_scale_bins];
	Energy_Orthonormal_Transformation = new double**[num_scale_bins];
	Commutator_Orthonormal_Transformation_Real = new double**[num_scale_bins];
	Commutator_Orthonormal_Transformation_Imag = new double**[num_scale_bins];

	Laplacian_Energy_Derivative = new double*[num_scale_bins - 1];
	Laplacian_Energy_Derivative_smoothed = new double*[num_scale_bins - 1];

	KineticEnergy_Vector = new double[num_points];
	PotentialEnergy_Vector = new double[num_points];
	ClassicalLocalEnergy_Vector = new double[num_points];
	ClassicalEnergy_Exchange = new double[num_points];

	ClassicalLaplacian_Energy = new double[num_points];
	ClassicalLaplacian_EigenStates = new double*[num_points];
	ClassicalHamiltonian_Energy = new double[num_points];
	ClassicalHamiltonian_wVacuum_Energy = new double[num_points];
	ClassicalHamiltonian_EigenStates = new double*[num_points];
	ClassicalHamiltonian_wVacuum_EigenStates = new double*[num_points];
	ClassicalParticle_Entropy = new double[num_points];;
	for (int i_p = 0; i_p < num_points; i_p++)
	{
		ClassicalLaplacian_EigenStates[i_p] = new double[num_points];
		ClassicalHamiltonian_EigenStates[i_p] = new double[num_points];
		ClassicalHamiltonian_wVacuum_EigenStates[i_p] = new double[num_points];
	}

	for (long i_s = 0; i_s < num_scale_bins; i_s++)
	{
		Laplacian_Energy[i_s] = new double[num_points];
		Commutator_Energy[i_s] = new double[num_points];
		Energy_Vector[i_s] = new double[num_points];
		Mass_Vector[i_s] = new double[num_points];
		Laplacian_Orthonormal_Transformation[i_s] = new double*[num_points];
		Energy_Orthonormal_Transformation[i_s] = new double*[num_points];
		Commutator_Orthonormal_Transformation_Real[i_s] = new double*[num_points];
		Commutator_Orthonormal_Transformation_Imag[i_s] = new double*[num_points];

		if (i_s < (num_scale_bins - 1))
		{
			Laplacian_Energy_Derivative[i_s] = new double[num_points];
			Laplacian_Energy_Derivative_smoothed[i_s] = new double[num_points];
		}

		for (int i_p = 0; i_p < num_points; i_p++)
		{
			Laplacian_Orthonormal_Transformation[i_s][i_p] = new double[num_points];
			Energy_Orthonormal_Transformation[i_s][i_p] = new double[num_points];
			Commutator_Orthonormal_Transformation_Real[i_s][i_p] = new double[num_points];
			Commutator_Orthonormal_Transformation_Imag[i_s][i_p] = new double[num_points];
		}
	}

	Num_TileHubs = Quantization_Engine_Pointer->Num_TileHubs;
	Tile_Dimension = Quantization_Engine_Pointer->Tile_Dimension;
	Tile_Positions = Quantization_Engine_Pointer->Tile_Positions;
	Tile_ScalarField = new double[Num_TileHubs];
	Tile_VectorField = new double*[Num_TileHubs];
	for (int i_p = 0; i_p < Num_TileHubs; i_p++)
	{
		Tile_VectorField[i_p] = new double[Tile_Dimension];
	}
	//Tile_Positions = new double*[Num_TileHubs];
	//for (int i_p = 0; i_p < Num_TileHubs; i_p++)
	//{
	//	Tile_Positions[i_p] = new double[Tile_Dimension];
	//	for (int i_d = 0; i_d < Tile_Dimension; i_d++)
	//		Tile_Positions[i_p][i_d] = Quantization_Engine_Pointer->Tile_Positions[i_p][i_d];
	//}

}

void WrapperClass::Quantization::Run(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag, bool smooth_flag)
{
	Quantization_Engine_Pointer->Run(positions, velocities, masses, num_points, dimension, dt, num_scale_bins, eigenvectors_flag, perturb_flag, smooth_flag);

	ClassicalEnergy = Quantization_Engine_Pointer->ClassicalEnergy;
	ClassicalKineticEnergy = Quantization_Engine_Pointer->ClassicalKineticEnergy;
	ClassicalPotentialEnergy = Quantization_Engine_Pointer->ClassicalPotentialEnergy;

	for (int i_p = 0; i_p < num_points; i_p++)
	{
		KineticEnergy_Vector[i_p] = (Quantization_Engine_Pointer->KineticEnergy_Vector)[i_p];
		PotentialEnergy_Vector[i_p] = (Quantization_Engine_Pointer->PotentialEnergy_Vector)[i_p];
		ClassicalLocalEnergy_Vector[i_p] = (Quantization_Engine_Pointer->ClassicalLocalEnergy_Vector)[i_p];
		ClassicalEnergy_Exchange[i_p] = (Quantization_Engine_Pointer->ClassicalEnergy_Exchange)[i_p];
		ClassicalHamiltonian_Energy[i_p] = (Quantization_Engine_Pointer->ClassicalHamiltonian_Energy)[i_p];
		ClassicalHamiltonian_wVacuum_Energy[i_p] = (Quantization_Engine_Pointer->ClassicalHamiltonian_wVacuum_Energy)[i_p];
		ClassicalLaplacian_Energy[i_p] = (Quantization_Engine_Pointer->ClassicalLaplacian_Energy)[i_p];
		ClassicalParticle_Entropy[i_p] = (Quantization_Engine_Pointer->ClassicalParticle_Entropy)[i_p];

		for (int j_p = 0; j_p < num_points; j_p++)
		{
			ClassicalHamiltonian_EigenStates[i_p][j_p] = (Quantization_Engine_Pointer->ClassicalHamiltonian_EigenStates)(i_p, j_p);
			ClassicalHamiltonian_wVacuum_EigenStates[i_p][j_p] = (Quantization_Engine_Pointer->ClassicalHamiltonian_wVacuum_EigenStates)(i_p, j_p);
			ClassicalLaplacian_EigenStates[i_p][j_p] = (Quantization_Engine_Pointer->ClassicalLaplacian_EigenStates)(i_p, j_p);
		}
	}
	for (long i_s = 0; i_s < num_scale_bins; i_s++)
	{
		for (int i_p = 0; i_p < num_points; i_p++)
		{
			Laplacian_Energy[i_s][i_p] = (Quantization_Engine_Pointer->Laplacian_Energy)[i_s][i_p];
			Commutator_Energy[i_s][i_p] = (Quantization_Engine_Pointer->Commutator_Energy)[i_s][i_p];
			Energy_Vector[i_s][i_p] = (Quantization_Engine_Pointer->Energy_Vector)[i_s][i_p];
			Mass_Vector[i_s][i_p] = (Quantization_Engine_Pointer->Mass_Vector)[i_s][i_p];

			if (smooth_flag)
			{
				if (i_s < (num_scale_bins - 1))
				{
					Laplacian_Energy_Derivative[i_s][i_p] = (Quantization_Engine_Pointer->Laplacian_Energy_Derivative)[i_s][i_p];
					Laplacian_Energy_Derivative_smoothed[i_s][i_p] = (Quantization_Engine_Pointer->Laplacian_Energy_Derivative_smoothed)[i_s][i_p];
				}
			}

			if (eigenvectors_flag)
			{
				for (int j_p = 0; j_p < num_points; j_p++)
				{
					Laplacian_Orthonormal_Transformation[i_s][i_p][j_p] = (Quantization_Engine_Pointer->Laplacian_Orthonormal_Transformation)[i_s](i_p, j_p);
					Energy_Orthonormal_Transformation[i_s][i_p][j_p] = (Quantization_Engine_Pointer->Energy_Orthonormal_Transformation)[i_s](i_p, j_p);
					Commutator_Orthonormal_Transformation_Real[i_s][i_p][j_p] = (Quantization_Engine_Pointer->Commutator_Orthonormal_Transformation_Real)[i_s](i_p, j_p);
					Commutator_Orthonormal_Transformation_Imag[i_s][i_p][j_p] = (Quantization_Engine_Pointer->Commutator_Orthonormal_Transformation_Imag)[i_s](i_p, j_p);
				}
			}
		}
	}

	Dendogram_Original = Quantization_Engine_Pointer->Hierarchical_Clusters->Dendogram->Structure;
	Scale_Ladder_Original = Quantization_Engine_Pointer->Hierarchical_Clusters->Dendogram->Scale;
	Dendogram_Dual = Quantization_Engine_Pointer->Hierarchical_Clusters->Dendogram_Dual->Structure;
	Scale_Ladder_Dual = Quantization_Engine_Pointer->Hierarchical_Clusters->Dendogram_Dual->Scale;

	Min_Scale = Quantization_Engine_Pointer->Min_Distance;
	Min_NonVacScale = Quantization_Engine_Pointer->Hierarchical_Clusters->Max_Vacuum_Scale;
	Max_NonVacScale = Quantization_Engine_Pointer->Hierarchical_Clusters->Min_Saturation_Scale;
	Max_Scale = Quantization_Engine_Pointer->Max_Distance;

	if (false) // tile_flag = false
	{
		for (int i_p = 0; i_p < Num_TileHubs; i_p++)
		{
			for (int i_d = 0; i_d < Tile_Dimension; i_d++)
			{
				Tile_ScalarField[i_p] = Quantization_Engine_Pointer->Tile_Fields->Scalar_Field[i_p];
				Tile_VectorField[i_p][i_d] = Quantization_Engine_Pointer->Tile_Fields->Vector_Field(i_p, i_d);
			}
		}
	}
}
