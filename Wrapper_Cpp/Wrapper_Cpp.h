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

namespace Wrapper_Cpp
{
	class Wrapper_Class
	{
	public:
		Wrapper_Class(double** positions, double** velocities, double* masses, int num_points, int dimension, int num_steps);
		double* Min_Scale;
		double* Min_NonVac_Scale;
		double* Max_NonVac_Scale;
		double* Max_Scale;

		~Wrapper_Class();
	private:
		Compute* Quantization_Engine_Pointer;
	};
}


