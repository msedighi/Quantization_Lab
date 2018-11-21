#include "stdafx.h"

#include "Wrapper_Class.h"

WrapperClass::Quantization::Quantization(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag)
{
	Quantization_Engine_Pointer = new Compute(num_points, num_scale_bins);
	Quantization_Engine_Pointer->Run(positions, velocities, masses, num_points, dimension, dt, num_scale_bins, eigenvectors_flag, perturb_flag);

	Laplacian_Energy = new double*[num_scale_bins + 2];
	Commutator_Energy = new double*[num_scale_bins + 2];
	Energy_Vector = new double*[num_scale_bins + 2];
	Mass_Vector = new double*[num_scale_bins + 2];
	Orthonormal_Transformation = new double**[num_scale_bins + 2];
	for (long i_s = 0; i_s < (num_scale_bins + 2); i_s++)
	{
		Laplacian_Energy[i_s] = new double[num_points];
		Commutator_Energy[i_s] = new double[num_points];
		Energy_Vector[i_s] = new double[num_points];
		Mass_Vector[i_s] = new double[num_points];
		Orthonormal_Transformation[i_s] = new double*[num_points];
		for (int i_p = 0; i_p < num_points; i_p++)
		{
			Laplacian_Energy[i_s][i_p] = (Quantization_Engine_Pointer->Laplacian_Energy)[i_s][i_p];
			Commutator_Energy[i_s][i_p] = (Quantization_Engine_Pointer->Commutator_Energy)[i_s][i_p];
			Energy_Vector[i_s][i_p] = (Quantization_Engine_Pointer->Energy_Vector)[i_s][i_p];
			Mass_Vector[i_s][i_p] = (Quantization_Engine_Pointer->Mass_Vector)[i_s][i_p];

			if (eigenvectors_flag)
			{
				Orthonormal_Transformation[i_s][i_p] = new double[num_points];
				for (int j_p = 0; j_p < num_points; j_p++)
				{
					Orthonormal_Transformation[i_s][i_p][j_p] = (Quantization_Engine_Pointer->Orthonormal_Transformation)[i_s](i_p, j_p);
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
}

