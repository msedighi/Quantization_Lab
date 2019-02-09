using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Serialization;

namespace Quantization_Tool
{
    public class State_Variables
    {
        public int Num_Points;
        public int Dimension;
        public int Number_Pairs;
        public int Num_ScaleBins;
        public int Num_EnergyBins;

        public double Mass_Ratio_Max;
        public double[] Mass_Ratios;

        public double[][] Positions;
        public double[][] Velocities;

        public State_Variables(int num_points, int dimension, int num_scalebins) : this(num_points, dimension)
        {
            Num_ScaleBins = num_scalebins;
        }

        public State_Variables(int num_points, int dimension)
        {
            Num_Points = num_points;
            Dimension = dimension;
            Number_Pairs = num_points * (num_points - 1) / 2;
            Num_ScaleBins = 5 * num_points * num_points;
            Num_EnergyBins = num_points * 100;

            Mass_Ratios = new double[num_points];
            Positions = new double[num_points][];
            Velocities = new double[num_points][];
            for (uint i = 0; i < num_points; i++)
            {
                Positions[i] = new double[dimension];
                Velocities[i] = new double[dimension];
            }
        }

        public State_Variables() { }
    }

    public class Output_Variables
    {
        public double Min_Scale;
        public double Min_NonVacScale;
        public double Max_NonVacScale;
        public double Max_Scale;

        public double[] Scale_Ladder_Original;
        public double[] Scale_Ladder_Dual;
        public bool[][] Dendogram_Original;
        public bool[][] Dendogram_Dual;

        public double[][] Laplacian_Energy;
        public double[][] Commutator_Energy;
        public double[][] Mass_Vector;
        public double[][] Energy_Vector;
        public double[][][] Laplacian_Orthonormal_Transformation;
        public double[][][] Energy_Orthonormal_Transformation;
        public double[][][] Commutator_Orthonormal_Transformation_Real;
        public double[][][] Commutator_Orthonormal_Transformation_Imag;

        // Classical Variables
        public double ClassicalEnergy;
        public double ClassicalKineticEnergy;
        public double ClassicalPotentialEnergy;

        public double[] KineticEnergy_Vector;
        public double[] PotentialEnergy_Vector;
        public double[] ClassicalEnergy_Exchange;
        public double[] ClassicalEnergy_Vector;

        public double[] ClassicalLaplacian_Energy;
        public double[][] ClassicalLaplacian_EigenStates;
        public double[] ClassicalHamiltonian_Energy;
        public double[][] ClassicalHamiltonian_EigenStates;
        public double[] ClassicalHamiltonian_wVacuum_Energy;
        public double[][] ClassicalHamiltonian_wVacuum_EigenStates;
        //

        [XmlIgnore]
        public double[][] Laplacian_Energy_Derivative;
        [XmlIgnore]
        public double[][] Laplacian_Energy_Derivative_smoothed;

        public Output_Variables(State_Variables state)
        {
            Scale_Ladder_Original = new double[state.Num_Points];
            Scale_Ladder_Dual = new double[state.Num_Points];
            Dendogram_Original = new bool[state.Num_Points][];
            Dendogram_Dual = new bool[state.Num_Points][];

            KineticEnergy_Vector = new double[state.Num_Points];
            PotentialEnergy_Vector = new double[state.Num_Points];
            ClassicalEnergy_Exchange = new double[state.Num_Points];
            ClassicalEnergy_Vector = new double[state.Num_Points];

            ClassicalLaplacian_Energy = new double[state.Num_Points];
            ClassicalLaplacian_EigenStates = new double[state.Num_Points][];
            ClassicalHamiltonian_Energy = new double[state.Num_Points];
            ClassicalHamiltonian_wVacuum_Energy = new double[state.Num_Points];
            ClassicalHamiltonian_EigenStates = new double[state.Num_Points][];
            ClassicalHamiltonian_wVacuum_EigenStates = new double[state.Num_Points][];

            for (int i = 0; i < state.Num_Points; i++)
            {
                Dendogram_Original[i] = new bool[state.Num_Points];
                Dendogram_Dual[i] = new bool[state.Num_Points];

                ClassicalLaplacian_EigenStates[i] = new double[state.Num_Points];
                ClassicalHamiltonian_EigenStates[i] = new double[state.Num_Points];
                ClassicalHamiltonian_wVacuum_EigenStates[i] = new double[state.Num_Points];
            }
            Laplacian_Energy = new double[state.Num_ScaleBins][];
            Commutator_Energy = new double[state.Num_ScaleBins][];
            Mass_Vector = new double[state.Num_ScaleBins][];
            Energy_Vector = new double[state.Num_ScaleBins][];
            Laplacian_Orthonormal_Transformation = new double[state.Num_ScaleBins][][];
            Energy_Orthonormal_Transformation = new double[state.Num_ScaleBins][][];
            Commutator_Orthonormal_Transformation_Real = new double[state.Num_ScaleBins][][];
            Commutator_Orthonormal_Transformation_Imag = new double[state.Num_ScaleBins][][];

            Laplacian_Energy_Derivative = new double[state.Num_ScaleBins - 1][];
            Laplacian_Energy_Derivative_smoothed = new double[state.Num_ScaleBins - 1][];
            for (int i = 0; i < state.Num_ScaleBins; i++)
            {
                Laplacian_Energy[i] = new double[state.Num_Points];
                Commutator_Energy[i] = new double[state.Num_Points];
                Mass_Vector[i] = new double[state.Num_Points];
                Energy_Vector[i] = new double[state.Num_Points];
                Laplacian_Orthonormal_Transformation[i] = new double[state.Num_Points][];
                Energy_Orthonormal_Transformation[i] = new double[state.Num_Points][];
                Commutator_Orthonormal_Transformation_Real[i] = new double[state.Num_Points][];
                Commutator_Orthonormal_Transformation_Imag[i] = new double[state.Num_Points][];

                for (int j = 0; j < state.Num_Points; j++)
                {
                    Laplacian_Orthonormal_Transformation[i][j] = new double[state.Num_Points];
                    Energy_Orthonormal_Transformation[i][j] = new double[state.Num_Points];
                    Commutator_Orthonormal_Transformation_Real[i][j] = new double[state.Num_Points];
                    Commutator_Orthonormal_Transformation_Imag[i][j] = new double[state.Num_Points];
                }

                if (i < (state.Num_ScaleBins - 1))
                {
                    Laplacian_Energy_Derivative[i] = new double[state.Num_Points];
                    Laplacian_Energy_Derivative_smoothed[i] = new double[state.Num_Points];
                }
            }
        }

        public Output_Variables() { }
    }

    public class Tile_Variables
    {
        public int Num_TileHubs;
        public int Tile_Dimension;
        public double[][] Tile_Positions;
        public double[] Tile_ScalarField;
        public double[][] Tile_VectorField;

        public Tile_Variables(int num_hubs, int dimension)
        {
            Num_TileHubs = num_hubs;
            Tile_Dimension = dimension;
            Tile_ScalarField = new double[num_hubs];
            Tile_VectorField = new double[num_hubs][];
            Tile_Positions = new double[num_hubs][];
            for (uint i_p = 0; i_p < num_hubs; i_p++)
            {
                Tile_Positions[i_p] = new double[dimension];
                Tile_VectorField[i_p] = new double[dimension];
            }
        }
    }

}
