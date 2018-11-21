#pragma once
#include "D:\Projects\Quantization_Engine\Quantization_Engine\TriangleTile.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\TriangleTile.cpp"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Wave_States.h"
#include "D:\Projects\Quantization_Engine\Quantization_Engine\Wave_States.cpp"
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
		Quantization(double** positions, double** velocities, double* masses, int num_points, int dimension, double dt, long num_scale_bins, bool eigenvectors_flag, bool perturb_flag);
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
		double*** Orthonormal_Transformation;
	private:
		Compute* Quantization_Engine_Pointer;
	};
}
