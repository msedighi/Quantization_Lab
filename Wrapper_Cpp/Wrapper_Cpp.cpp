// Wrapper_Cpp.cpp : Defines the exported functions for the DLL application.
//
#include "stdafx.h"
#include "Wrapper_Cpp.h"

Wrapper_Cpp::Wrapper_Class::Wrapper_Class(double** positions, double** velocities, double* masses, int num_points, int dimension, int num_steps)
{
	Quantization_Engine_Pointer = new Compute(positions, velocities, masses, num_points, dimension, num_steps);
	Min_Scale = Quantization_Engine_Pointer->Min_Distance;
	Max_Scale = Quantization_Engine_Pointer->Max_Distance;
	Min_NonVac_Scale = Quantization_Engine_Pointer->Min_NonVacuum_Range;
	Max_NonVac_Scale = Quantization_Engine_Pointer->Max_NonVacuum_Range;
}

Wrapper_Cpp::Wrapper_Class::~Wrapper_Class()
{}



