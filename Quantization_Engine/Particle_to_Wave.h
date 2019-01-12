#pragma once
#ifndef PARTICLE_TO_WAVE_H
#define PARTICLE_TO_WAVE_H
#endif // !1

#include "Eigen/Dense"
#include "TriangleTile.h"
#include "Particle_Dynamics.h"

struct Tile_Fields
{
	Eigen::VectorXd Scalar_Field;
	Eigen::MatrixXd Vector_Field;
	Tile_Fields(int num_hubs, int dim)
	{
		Scalar_Field = Eigen::VectorXd::Constant(num_hubs, 0);
		Vector_Field = Eigen::MatrixXd::Constant(num_hubs, dim, 0);
	}
};

class WaveParticle_Interface
{	
public:
	TriangleTile* Tile;

	Tile_Fields* Run(Global_Variables* global_variables, int num_points, int dim);
	Eigen::VectorXd Run(Eigen::VectorXd global_points_state, int num_points);
	Eigen::VectorXd Run(Eigen::VectorXd points_state, Eigen::MatrixXd orthogonal_transformation, int num_points);
	Eigen::VectorXd Run(double* points_state, int num_points, Eigen::MatrixXd orthogonal_transformation);

	WaveParticle_Interface(int Num_Layers);
	~WaveParticle_Interface();
private:
	Wave_State Hexagon_Waves;

	Eigen::VectorXd acting_global_points_state;
	Eigen::MatrixXd acting_global_vector_state;
	int num_points;
};