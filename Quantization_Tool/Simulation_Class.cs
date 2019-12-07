using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using WrapperClass;
using System.Runtime.InteropServices;
using System.Xml;
using System.Xml.Serialization;

namespace Quantization_Tool
{
    public class Simulation
    {
        public double Time_Range;
        public double dt = 0.005; //Default
        //private int Num_TimeSteps = 1;

        public double[] Coordinate_Range;
        public double[] Speed_Range;

        WrapperClass.Quantization Q;
        Output_Variables Out;
        [XmlIgnore]
        public Tile_Variables Tile;

        GCHandle[] handle_Positions;
        GCHandle[] handle_Velocities;

        public Simulation(State_Variables state)
        {
            Coordinate_Range = new double[state.Dimension];
            Speed_Range = new double[state.Dimension];

            Out = new Output_Variables(state);
            // Fixing the pointer handles 
            handle_Positions = new GCHandle[state.Num_Points];
            handle_Velocities = new GCHandle[state.Num_Points];
            unsafe
            {
                // Un-managed Code!
                Q = new Quantization(state.Num_Points, state.Num_ScaleBins);
            }

            Tile = new Tile_Variables(Q.Num_TileHubs, Q.Tile_Dimension);
            unsafe
            {
                // Un-managed Code!
                for (uint i_p = 0; i_p < Tile.Num_TileHubs; i_p++)
                {
                    for (uint i_d = 0; i_d < Tile.Tile_Dimension; i_d++)
                    {
                        Tile.Tile_Positions[i_p][i_d] = Q.Tile_Positions[i_p][i_d];
                    }
                }
            }
        }

        public Simulation() { }

        public Output_Variables Evolve(State_Variables state, bool eigenvectors_flag, bool perturb_flag, bool smooth_flag)
        {
            for (uint i_p = 0; i_p < state.Num_Points; i_p++)
            {
                handle_Positions[i_p] = GCHandle.Alloc(state.Positions[i_p], GCHandleType.Pinned);
                handle_Velocities[i_p] = GCHandle.Alloc(state.Velocities[i_p], GCHandleType.Pinned);
            }
            //
            unsafe
            {
                double*[] _positions = new double*[state.Num_Points];
                double*[] _velocities = new double*[state.Num_Points];

                for (uint i_p = 0; i_p < state.Num_Points; i_p++)
                {
                    _positions[i_p] = (double*)handle_Positions[i_p].AddrOfPinnedObject().ToPointer();
                    _velocities[i_p] = (double*)handle_Velocities[i_p].AddrOfPinnedObject().ToPointer();
                }
                fixed (double** p = &_positions[0], v = &_velocities[0])
                {
                    fixed (double* m = &state.Mass_Ratios[0])
                    {
                        // Un-managed Code!
                        Q.Run(p, v, m, state.Num_Points, state.Dimension, dt, state.Num_ScaleBins, eigenvectors_flag, perturb_flag, smooth_flag);
                    }
                }

                for (uint i_s = 0; i_s < state.Num_ScaleBins; i_s++)
                {
                    for (uint i_p = 0; i_p < state.Num_Points; i_p++)
                    {
                        Out.Laplacian_Energy[i_s][i_p] = Q.Laplacian_Energy[i_s][i_p];
                        Out.Commutator_Energy[i_s][i_p] = Q.Commutator_Energy[i_s][i_p];
                        Out.Mass_Vector[i_s][i_p] = Q.Mass_Vector[i_s][i_p];
                        Out.Energy_Vector[i_s][i_p] = Q.Energy_Vector[i_s][i_p];

                        if (eigenvectors_flag)
                        {
                            if ((i_s < (state.Num_ScaleBins - 1)) && (smooth_flag))
                            {
                                Out.Laplacian_Energy_Derivative[i_s][i_p] = Q.Laplacian_Energy_Derivative[i_s][i_p];
                                Out.Laplacian_Energy_Derivative_smoothed[i_s][i_p] = Q.Laplacian_Energy_Derivative_smoothed[i_s][i_p];
                            }
                            for (uint j_p = 0; j_p < state.Num_Points; j_p++)
                            {
                                Out.Laplacian_Orthonormal_Transformation[i_s][i_p][j_p] = Q.Laplacian_Orthonormal_Transformation[i_s][i_p][j_p];
                                Out.Energy_Orthonormal_Transformation[i_s][i_p][j_p] = Q.Energy_Orthonormal_Transformation[i_s][i_p][j_p];
                                Out.Commutator_Orthonormal_Transformation_Real[i_s][i_p][j_p] = Q.Commutator_Orthonormal_Transformation_Real[i_s][i_p][j_p];
                                Out.Commutator_Orthonormal_Transformation_Imag[i_s][i_p][j_p] = Q.Commutator_Orthonormal_Transformation_Imag[i_s][i_p][j_p];
                            }
                        }
                    }
                }

                for (uint i_p = 0; i_p < state.Num_Points; i_p++)
                {
                    Out.Scale_Ladder_Original[i_p] = Q.Scale_Ladder_Original[i_p];
                    Out.Scale_Ladder_Dual[i_p] = Q.Scale_Ladder_Dual[i_p];

                    Out.ClassicalEnergy_Exchange[i_p] = Q.ClassicalEnergy_Exchange[i_p];
                    Out.KineticEnergy_Vector[i_p] = Q.KineticEnergy_Vector[i_p] / 2;
                    Out.PotentialEnergy_Vector[i_p] = Q.PotentialEnergy_Vector[i_p] / 2;
                    //Out.ClassicalEnergy_Vector[i_p] = Out.KineticEnergy_Vector[i_p] + Out.PotentialEnergy_Vector[i_p];
                    //Out.ClassicalEnergy_Vector[i_p] = Out.KineticEnergy_Vector[i_p];
                    Out.ClassicalEnergy_Vector[i_p] = Q.ClassicalLocalEnergy_Vector[i_p];

                    Out.ClassicalLaplacian_Energy[i_p] = Q.ClassicalLaplacian_Energy[i_p];
                    Out.ClassicalHamiltonian_Energy[i_p] = Q.ClassicalHamiltonian_Energy[i_p];
                    Out.ClassicalHamiltonian_wVacuum_Energy[i_p] = Q.ClassicalHamiltonian_wVacuum_Energy[i_p];
                    Out.ClassicalParticle_Entropy[i_p] = Q.ClassicalParticle_Entropy[i_p];
                    for (uint j_p = 0; j_p < state.Num_Points; j_p++)
                    {
                        Out.Dendogram_Original[i_p][j_p] = Q.Dendogram_Original[i_p][j_p];
                        Out.Dendogram_Dual[i_p][j_p] = Q.Dendogram_Dual[i_p][j_p];

                        Out.ClassicalHamiltonian_EigenStates[i_p][j_p] = Q.ClassicalHamiltonian_EigenStates[i_p][j_p];
                        Out.ClassicalHamiltonian_wVacuum_EigenStates[i_p][j_p] = Q.ClassicalHamiltonian_wVacuum_EigenStates[i_p][j_p];
                        Out.ClassicalLaplacian_EigenStates[i_p][j_p] = Q.ClassicalLaplacian_EigenStates[i_p][j_p];
                    }
                }

                for (uint i_p = 0; i_p < Tile.Num_TileHubs; i_p++)
                {
                    Tile.Tile_ScalarField[i_p] = Q.Tile_ScalarField[i_p];
                    for (uint i_d = 0; i_d < Tile.Tile_Dimension; i_d++)
                    {
                        Tile.Tile_VectorField[i_p][i_d] = Q.Tile_VectorField[i_p][i_d];
                    }
                }
            }

            // Releasing the pointer handles 
            for (uint i_p = 0; i_p < state.Num_Points; i_p++)
            {
                handle_Positions[i_p].Free();
                handle_Velocities[i_p].Free();
            }
            //
            Out.Min_Scale = Q.Min_Scale;
            Out.Max_Scale = Q.Max_Scale;
            Out.Min_NonVacScale = Q.Min_NonVacScale;
            Out.Max_NonVacScale = Q.Max_NonVacScale;

            Out.ClassicalEnergy = Q.ClassicalEnergy / state.Num_Points;
            Out.ClassicalKineticEnergy = Q.ClassicalKineticEnergy / state.Num_Points;
            Out.ClassicalPotentialEnergy = Q.ClassicalPotentialEnergy / state.Num_Points;

            return Out;
        }

    }
}
