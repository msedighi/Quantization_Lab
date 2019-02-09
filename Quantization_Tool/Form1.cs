using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Runtime.InteropServices;
using System.IO;
using OxyPlot;

namespace Quantization_Tool
{
    public partial class Form1 : Form
    {
        private PlotModel Points_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] Points_Data;
        private PointsPlotWiring PointsPlot_Wiring;

        private PlotModel Tile_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] Tile_Data;
        private PointsPlotWiring TilePlot_Wiring;

        private PlotModel ClassicalHamiltonianTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries ClassicalEnergy_Time = new OxyPlot.Series.ScatterSeries();
        private OxyPlot.Series.ScatterSeries[] ClassicalHamiltonianTime_Data;
        private EnergyTimePlotWiring ClassicalHamiltonianTimePlot_Wiring;

        private PlotModel ClassicalHamiltonian_wVacuumTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] ClassicalHamiltonian_wVacuum_Time_Data;
        private DataPlotWiring ClassicalHamiltonian_wVacuumTimePlot_Wiring;

        private PlotModel ClassicalLaplacianTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] ClassicalLaplacianTime_Data;
        private DataPlotWiring ClassicalLaplacianTimePlot_Wiring;

        private PlotModel ClassicalEnergyExchangeTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] ClassicalEnergyExchangeTime_Data;
        private DataPlotWiring ClassicalEnergyExchangeTimePlot_Wiring;

        private PlotModel ClassicalPotentialEnergyTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries ClassicalTotalPotentialEnergy_Time = new OxyPlot.Series.ScatterSeries();
        private OxyPlot.Series.ScatterSeries[] ClassicalPotentialEnergyTime_Data;
        private EnergyTimePlotWiring ClassicalPotentialEnergyTimePlot_Wiring;

        private PlotModel ClassicalKineticEnergyTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries ClassicalTotalKineticEnergy_Time = new OxyPlot.Series.ScatterSeries();
        private OxyPlot.Series.ScatterSeries[] ClassicalKineticEnergyTime_Data;
        private EnergyTimePlotWiring ClassicalKineticEnergyTimePlot_Wiring;

        private PlotModel ClassicalEnergyTime_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries ClassicalTotalEnergy_Time = new OxyPlot.Series.ScatterSeries();
        private OxyPlot.Series.ScatterSeries[] ClassicalEnergyTime_Data;
        private EnergyTimePlotWiring ClassicalEnergyTimePlot_Wiring;

        private PlotModel ScaleTime_Plot = new PlotModel();
        private OxyPlot.Series.LineSeries ScaleTime_MinScale = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries ScaleTime_Min_NonVacScale = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries ScaleTime_Max_NonVacScale = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries ScaleTime_MaxScale = new OxyPlot.Series.LineSeries();

        private PlotModel MassScale_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] MassScale_Data;
        private DataPlotWiring MassScalePlot_Wiring;
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Max_Mass = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Min_Mass = new OxyPlot.Series.LineSeries();

        private PlotModel EnergyScale_Laplacian_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] EnergyScale_Laplacian_Data;
        private DataPlotWiring LaplacianScalePlot_Wiring;
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Max_Laplacian = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Min_Laplacian = new OxyPlot.Series.LineSeries();

        private PlotModel EnergyScale_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] EnergyScale_Data;
        private DataPlotWiring EnergyScalePlot_Wiring;
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Max = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Min = new OxyPlot.Series.LineSeries();

        private PlotModel EnergyScale_LaplacianDerivative_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] EnergyScale_LaplacianDerivative_Data;

        private PlotModel EnergyScale_LaplacianDerivative_Plot_smoothed = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] EnergyScale_LaplacianDerivative_Data_smoothed;

        private PlotModel EnergyScale_Commutator_Plot = new PlotModel();
        private OxyPlot.Series.ScatterSeries[] EnergyScale_Commutator_Data;
        private DataPlotWiring CommutatorScalePlot_Wiring;
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Max_Commutator = new OxyPlot.Series.LineSeries();
        private OxyPlot.Series.LineSeries NonVacScale_Boundary_Min_Commutator = new OxyPlot.Series.LineSeries();

        private PlotModel EnergyScale_Laplacian_HeatMap_Plot = new PlotModel();
        private OxyPlot.Series.HeatMapSeries EnergyScale_Laplacian_HeatMap = new OxyPlot.Series.HeatMapSeries();
        private double[,] EnergyScale_Laplacian_HeatMap_Data;
        private double energy_laplacian_variable;

        private PlotModel EnergyScale_Commutator_HeatMap_Plot = new PlotModel();
        private OxyPlot.Series.HeatMapSeries EnergyScale_Commutator_HeatMap = new OxyPlot.Series.HeatMapSeries();
        private double[,] EnergyScale_Commutator_HeatMap_Data;
        private double energy_commutator_variable;

        private PlotModel EnergyScale_HeatMap_Plot = new PlotModel();
        private OxyPlot.Series.HeatMapSeries EnergyScale_HeatMap = new OxyPlot.Series.HeatMapSeries();
        private double[,] EnergyScale_HeatMap_Data;
        private double energy_variable;

        private PlotModel MassScale_HeatMap_Plot = new PlotModel();
        private OxyPlot.Series.HeatMapSeries MassScale_HeatMap = new OxyPlot.Series.HeatMapSeries();
        private double[,] MassScale_HeatMap_Data;
        private double mass_variable;

        private PlotModel Selected_Plot;
        private OxyPlot.Series.ScatterSeries Selected_Series;
        private OxyPlot.Series.ScatterPoint Selected_Point;
        private int Selected_Index;

        private double Max_HeatMap_4plot;
        private double Max_HeatMap_4plot1;
        private double Max_Commutator_Energy;
        private double Min_Energy;

        private Simulation Sim_Variables;
        private State_Variables Initial_State;
        private State_Variables State;
        private Output_Variables output_Variables;
        private bool Eigenvectors_flag;
        private bool Perturb_flag;
        private bool Smooth_flag;
        private bool VaccumScaleBounray_flag;
        private bool PlotColoring_flag;

        private int TimeLimit_4plot;
        private int SimTime_Counter;
        private ManualResetEvent MRE;
        private Thread Sim_Thread;
        private Thread Graph_Thread;

        private string mydoc_path = System.Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);
        private string output_dir = "D:\\Projects\\Quantization_Log";
        private string output_path;
        private Data_Output DO = new Data_Output();

        public Form1()
        {
            InitializeComponent();

            // Initialize Values
            Random rand = new Random();
            int _Num_Points = 6;
            int _Dimension = 2;
            //int _Num_ScaleBins = 200;
            int _Num_ScaleBins = 20;
            Eigenvectors_flag = false;
            Perturb_flag = false;
            Smooth_flag = false;

            State = new State_Variables(_Num_Points, _Dimension, _Num_ScaleBins);
            Initial_State = new State_Variables(State.Num_Points, State.Dimension, State.Num_ScaleBins);

            Max_HeatMap_4plot = 160.0;
            Max_HeatMap_4plot1 = 20.0;
            Max_Commutator_Energy = State.Num_Points;
            Min_Energy = 0;
            VaccumScaleBounray_flag = false;
            PlotColoring_flag = false;

            TimeLimit_4plot = 200;
            SimTime_Counter = 0;
            RunButton_State = false;
            PauseButton_State = false;
            Sim_Variables = new Simulation(State);
            for (int i_d = 0; i_d < _Dimension; i_d++)
            {
                Sim_Variables.Coordinate_Range[i_d] = 10.0;
                Sim_Variables.Speed_Range[i_d] = 2.0;
            }

            Sim_Variables.Time_Range = 2.0;

            // This:
            Sim_Variables.dt = 0.002;
            //Sim_Variables.Num_TimeSteps = 100;
            // Or, 
            // Choose one and find the other : 
            // Num_Time_Bins = (uint)(Time_Range / dt);
            // dt = (double)(Time_Range / Num_Time_Bins);

            State.Mass_Ratio_Max = 10.0; // Default (should be > 1)            
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                // This:
                State.Mass_Ratios[i_p] = 1.0;
                // Or, 
                // Ask Mass_Ratio_Max (should be > 1) :
                //State.Mass_Ratios[i_p] = (State.Mass_Ratio_Max - 1.0) * rand.NextDouble() + 1.0;
                // Or,
                // Ask all the Mass Ratios!
            }

            double Hubble = 0;
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                for (uint i_d = 0; i_d < State.Dimension; i_d++)
                {
                    State.Positions[i_p][i_d] = Sim_Variables.Coordinate_Range[i_d] * (rand.NextDouble() - 0.5);
                    State.Velocities[i_p][i_d] = Sim_Variables.Speed_Range[i_d] * (rand.NextDouble() - 0.5);
                    //State.Velocities[i_p][i_d] = Hubble * State.Positions[i_p][i_d];

                    Initial_State.Positions[i_p][i_d] = State.Positions[i_p][i_d];
                    Initial_State.Velocities[i_p][i_d] = State.Velocities[i_p][i_d];
                }
            }

            output_Variables = new Output_Variables(State);

            // Writing to Files
            string[] dirs = Directory.GetDirectories(output_dir,"*", SearchOption.TopDirectoryOnly);
            output_path = output_dir + "\\" + (dirs.Length + 1).ToString();
            Directory.CreateDirectory(output_path);
            string State_FileName = output_path + "\\State_Input.txt";
            string Simulation_FileName = output_path + "\\Simulation_Input.txt";

            DO.Write_XML(State_FileName, State);
            DO.Write_XML(Simulation_FileName, Sim_Variables);

            // Choosing Force function ??
            // Probably have a choice between major force functions: Gravitation, Spring, Lenard-Jones, Hard (or Stiky) Balls, ..

            //
            // Points Plot
            //
            Initialize_PointsPlot(Points_Plot, "Points");
            Points_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];

            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                Points_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                Points_Data[i_p].MarkerType = MarkerType.Circle;
                Points_Data[i_p].MarkerFill = OxyColors.Black;
                Points_Data[i_p].MarkerSize = 3;
                if (State.Dimension == 1)
                    Points_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(State.Positions[i_p][0], 0));
                else
                    Points_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(State.Positions[i_p][0], State.Positions[i_p][1]));
            }

            PointsPlot_Wiring = new PointsPlotWiring(Points_Plot, Points_Data);
            plot_ToolStripContainer_BottomRight.Attach(Points_Plot);
            PointsPlot_Wiring.Activate();
            //

            // Tile Plot
            //
            Initialize_PointsPlot(Tile_Plot, "Triangular Tile");
            Tile_Data = new OxyPlot.Series.ScatterSeries[Sim_Variables.Tile.Num_TileHubs];

            for (uint i_p = 0; i_p < Sim_Variables.Tile.Num_TileHubs; i_p++)
            {
                Tile_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                Tile_Data[i_p].MarkerType = MarkerType.Circle; 
                Tile_Data[i_p].MarkerFill = OxyColors.Black;
                Tile_Data[i_p].MarkerSize = 10;
                //Tile_Data[i_p].MarkerSize = 2;
                if (Sim_Variables.Tile.Tile_Dimension == 1)
                    Tile_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(Sim_Variables.Tile.Tile_Positions[i_p][0], 0));
                else if (Sim_Variables.Tile.Tile_Dimension == 2)
                    Tile_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(Sim_Variables.Tile.Tile_Positions[i_p][0], Sim_Variables.Tile.Tile_Positions[i_p][1]));
                else
                    throw new Exception("Tile cannot be more than 2 dimensions");
            }

            TilePlot_Wiring = new PointsPlotWiring(Tile_Plot, Tile_Data);
            //plot_ToolStripContainer_TopRight.Attach(Tile_Plot);
            TilePlot_Wiring.Activate();
            //

            // Classical Energy Plots
            Initialize_ClassicalEnergyTimePlot(ClassicalHamiltonianTime_Plot, "Classical Hamitonian");
            Initialize_ClassicalEnergyTimePlot(ClassicalHamiltonian_wVacuumTime_Plot, "Classical Hamitonian w/ Vaccum");
            Initialize_ClassicalEnergyTimePlot(ClassicalLaplacianTime_Plot, "Classical Laplacian");
            Initialize_ClassicalEnergyTimePlot(ClassicalPotentialEnergyTime_Plot, "Potential Energy");
            Initialize_ClassicalEnergyTimePlot(ClassicalKineticEnergyTime_Plot, "Kinetic Energy");
            Initialize_ClassicalEnergyTimePlot(ClassicalEnergyExchangeTime_Plot, "Energy Exchange");
            Initialize_ClassicalEnergyTimePlot(ClassicalEnergyTime_Plot, "Classical Energy");

            ClassicalHamiltonianTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalHamiltonianTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalHamiltonianTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalHamiltonianTime_Data[i_p].MarkerSize = 3;
                //ClassicalHamiltonianTime_Data[i_p].MarkerFill = OxyColors.Black;
            }
            ClassicalEnergy_Time.MarkerType = MarkerType.Circle;
            ClassicalEnergy_Time.MarkerSize = 3;
            ClassicalEnergy_Time.MarkerFill = OxyColors.Black;

            ClassicalHamiltonian_wVacuum_Time_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalHamiltonian_wVacuum_Time_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalHamiltonian_wVacuum_Time_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalHamiltonian_wVacuum_Time_Data[i_p].MarkerSize = 3;
                ClassicalHamiltonian_wVacuum_Time_Data[i_p].MarkerFill = OxyColors.Black;
            }

            ClassicalLaplacianTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalLaplacianTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalLaplacianTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalLaplacianTime_Data[i_p].MarkerSize = 3;
                ClassicalLaplacianTime_Data[i_p].MarkerFill = OxyColors.Black;
            }

            ClassicalPotentialEnergyTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalPotentialEnergyTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalPotentialEnergyTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalPotentialEnergyTime_Data[i_p].MarkerSize = 3;
               // ClassicalPotentialEnergyTime_Data[i_p].MarkerFill = OxyColors.Black;
            }
            ClassicalTotalPotentialEnergy_Time.MarkerType = MarkerType.Circle;
            ClassicalTotalPotentialEnergy_Time.MarkerSize = 3;
            ClassicalTotalPotentialEnergy_Time.MarkerFill = OxyColors.Black;

            ClassicalKineticEnergyTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalKineticEnergyTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalKineticEnergyTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalKineticEnergyTime_Data[i_p].MarkerSize = 3;
                //ClassicalKineticEnergyTime_Data[i_p].MarkerFill = OxyColors.Black;
            }
            ClassicalTotalKineticEnergy_Time.MarkerType = MarkerType.Circle;
            ClassicalTotalKineticEnergy_Time.MarkerSize = 3;
            ClassicalTotalKineticEnergy_Time.MarkerFill = OxyColors.Black;

            ClassicalEnergyExchangeTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalEnergyExchangeTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalEnergyExchangeTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalEnergyExchangeTime_Data[i_p].MarkerSize = 3;
                //ClassicalEnergyExchangeTime_Data[i_p].MarkerFill = OxyColors.Black;
            }

            ClassicalEnergyTime_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            for (uint i_p = 0; i_p < State.Num_Points; i_p++)
            {
                ClassicalEnergyTime_Data[i_p] = new OxyPlot.Series.ScatterSeries();
                ClassicalEnergyTime_Data[i_p].MarkerType = MarkerType.Circle;
                ClassicalEnergyTime_Data[i_p].MarkerSize = 3;
                //ClassicalEnergyTime_Data[i_p].MarkerFill = OxyColors.Black;
            }
            ClassicalTotalEnergy_Time.MarkerType = MarkerType.Circle;
            ClassicalTotalEnergy_Time.MarkerSize = 3;
            ClassicalTotalEnergy_Time.MarkerFill = OxyColors.Black;

            plot_ToolStripContainer_TopMiddle.Attach(ClassicalHamiltonianTime_Plot);
            plot_ToolStripContainer_TopLeft.Attach(ClassicalLaplacianTime_Plot);
            plot_ToolStripContainer_TopRight.Attach(ClassicalHamiltonian_wVacuumTime_Plot);
            //plot_ToolStripContainer_TopRight.Attach(ClassicalPotentialEnergyTime_Plot);
            //plot_ToolStripContainer_BottomRight.Attach(ClassicalKineticEnergyTime_Plot);
            plot_ToolStripContainer_BottomMiddle.Attach(ClassicalEnergyExchangeTime_Plot);
            plot_ToolStripContainer_BottomLeft.Attach(ClassicalEnergyTime_Plot);

            ClassicalHamiltonianTimePlot_Wiring = new EnergyTimePlotWiring(ClassicalHamiltonianTime_Plot, ClassicalHamiltonianTime_Data, ClassicalEnergy_Time);
            ClassicalHamiltonian_wVacuumTimePlot_Wiring = new DataPlotWiring(ClassicalHamiltonian_wVacuumTime_Plot, ClassicalHamiltonian_wVacuum_Time_Data);
            ClassicalLaplacianTimePlot_Wiring = new DataPlotWiring(ClassicalLaplacianTime_Plot, ClassicalLaplacianTime_Data);
            ClassicalPotentialEnergyTimePlot_Wiring = new EnergyTimePlotWiring(ClassicalPotentialEnergyTime_Plot, ClassicalPotentialEnergyTime_Data, ClassicalTotalPotentialEnergy_Time);
            ClassicalKineticEnergyTimePlot_Wiring = new EnergyTimePlotWiring(ClassicalKineticEnergyTime_Plot, ClassicalKineticEnergyTime_Data, ClassicalTotalKineticEnergy_Time);
            ClassicalEnergyExchangeTimePlot_Wiring = new DataPlotWiring(ClassicalEnergyExchangeTime_Plot, ClassicalEnergyExchangeTime_Data);
            ClassicalEnergyTimePlot_Wiring = new EnergyTimePlotWiring(ClassicalEnergyTime_Plot, ClassicalEnergyTime_Data, ClassicalTotalEnergy_Time);

            ClassicalHamiltonianTimePlot_Wiring.Activate();
            ClassicalHamiltonian_wVacuumTimePlot_Wiring.Activate();
            ClassicalLaplacianTimePlot_Wiring.Activate();
            ClassicalPotentialEnergyTimePlot_Wiring.Activate();
            ClassicalKineticEnergyTimePlot_Wiring.Activate();
            ClassicalEnergyExchangeTimePlot_Wiring.Activate();
            ClassicalEnergyTimePlot_Wiring.Activate();

            ClassicalHamiltonianTimePlot_Wiring.Show_ClassicalEnergy();
            ClassicalPotentialEnergyTimePlot_Wiring.Show_ClassicalEnergy();
            ClassicalKineticEnergyTimePlot_Wiring.Show_ClassicalEnergy();
            ClassicalEnergyTimePlot_Wiring.Show_ClassicalEnergy();
            //

            Initialize_ScaleTimePlot(ScaleTime_Plot);
            ScaleTime_MinScale.Color = OxyColors.Red;
            ScaleTime_Min_NonVacScale.Color = OxyColors.Black;
            ScaleTime_Max_NonVacScale.Color = OxyColors.Black;
            ScaleTime_MaxScale.Color = OxyColors.Blue;

            ScaleTime_Plot.Series.Add(ScaleTime_MinScale);
            ScaleTime_Plot.Series.Add(ScaleTime_Min_NonVacScale);
            ScaleTime_Plot.Series.Add(ScaleTime_Max_NonVacScale);
            ScaleTime_Plot.Series.Add(ScaleTime_MaxScale);

            //plot_plotView.Model = EnergyScale_Laplacian_Plot;
            //plot_ToolStripContainer_TopRight.Attach(EnergyScale_Laplacian_Plot);
            EnergyScale_Laplacian_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_LaplacianEnergyScalePlot(EnergyScale_Laplacian_Plot, EnergyScale_Laplacian_Data, NonVacScale_Boundary_Max_Laplacian, NonVacScale_Boundary_Min_Laplacian, State, PlotColoring_flag, VaccumScaleBounray_flag);


            //plot_ToolStripContainer_BottomMiddle.Attach(EnergyScale_Commutator_Plot);
            EnergyScale_Commutator_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_CommutatorEnergyScalePlot(EnergyScale_Commutator_Plot, EnergyScale_Commutator_Data, NonVacScale_Boundary_Max_Commutator, NonVacScale_Boundary_Min_Commutator, State, PlotColoring_flag, VaccumScaleBounray_flag);
            //plot_ToolStripContainer_TopMiddle.Attach(MassScale_Plot);
            MassScale_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_MassScalePlot(MassScale_Plot, MassScale_Data, NonVacScale_Boundary_Max_Mass, NonVacScale_Boundary_Min_Mass, State, PlotColoring_flag, VaccumScaleBounray_flag);
            //plot_ToolStripContainer_BottomRight.Attach(EnergyScale_Plot);
            EnergyScale_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_EnergyScalePlot(EnergyScale_Plot, EnergyScale_Data, NonVacScale_Boundary_Max, NonVacScale_Boundary_Min, State, PlotColoring_flag, VaccumScaleBounray_flag);

            //bottomMiddle_Plot_SplitContainer.Attach(EnergyScale_LaplacianDerivative_Plot);
            EnergyScale_LaplacianDerivative_Data = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_LaplacianEnergyDerivativeScalePlot(EnergyScale_LaplacianDerivative_Plot, EnergyScale_LaplacianDerivative_Data, State);

            //bottomRight_Plot_SplitContainer.Attach(EnergyScale_LaplacianDerivative_Plot_smoothed);
            EnergyScale_LaplacianDerivative_Data_smoothed = new OxyPlot.Series.ScatterSeries[State.Num_Points];
            Initialize_LaplacianEnergyDerivativeScalePlot(EnergyScale_LaplacianDerivative_Plot_smoothed, EnergyScale_LaplacianDerivative_Data_smoothed, State);

            //plot_plotView.Model = EnergyScale_Commutator_HeatMap_Plot;
            //bottomMiddle_Plot_SplitContainer.Attach(EnergyScale_Commutator_HeatMap_Plot);
            Initialize_EnergyScaleHeatMap(EnergyScale_Commutator_HeatMap_Plot, EnergyScale_Commutator_HeatMap, State);
            EnergyScale_Commutator_HeatMap_Plot.Title = "Commutator Energy vs Scale";
            EnergyScale_Commutator_HeatMap.Y0 = -State.Num_EnergyBins;
            EnergyScale_Commutator_HeatMap.Y1 = State.Num_EnergyBins;

            plot_plotView.Model = EnergyScale_Laplacian_HeatMap_Plot;
            //topRight_plotView.Model = EnergyScale_Laplacian_HeatMap_Plot;
            Initialize_EnergyScaleHeatMap(EnergyScale_Laplacian_HeatMap_Plot, EnergyScale_Laplacian_HeatMap, State);
            EnergyScale_Laplacian_HeatMap_Plot.Title = "Laplacian Energy vs Scale";
            EnergyScale_Laplacian_HeatMap.Y0 = 0.0;
            EnergyScale_Laplacian_HeatMap.Y1 = State.Num_EnergyBins;

            //plot_plotView.Model = EnergyScale_HeatMap_Plot;
            //bottomLeft_plotView.Model = EnergyScale_HeatMap_Plot;
            Initialize_EnergyScaleHeatMap(EnergyScale_HeatMap_Plot, EnergyScale_HeatMap, State);
            EnergyScale_HeatMap_Plot.Title = "Energy vs Scale";
            EnergyScale_HeatMap.Y0 = 0.0;
            EnergyScale_HeatMap.Y1 = State.Num_EnergyBins;

            //topLeft_plotView.Model = MassScale_HeatMap_Plot;
            Initialize_EnergyScaleHeatMap(MassScale_HeatMap_Plot, MassScale_HeatMap, State);
            MassScale_HeatMap_Plot.Title = "Mass vs Scale";
            MassScale_HeatMap.Y0 = 0.0;
            MassScale_HeatMap.Y1 = State.Num_EnergyBins;

            // Threading
            MRE = new ManualResetEvent(false);

           
            //Sim_Thread = new Thread(new ThreadStart(Run_Simulation));
            //Sim_Thread.IsBackground = true;

            //Graph_Thread = new Thread(new ThreadStart(Plot_Graphs));
            //Graph_Thread.IsBackground = true;

            // Running the Simulation Thread
            //Sim_Thread.Start();
            // Running the Graph Thread
            //Graph_Thread.Start();

            //Plot_Graphs();

        }

        private void Plot_Graphs(int t)
        {
            string output_filename;
            Output_Variables current_output;
            //int t = 0;
            //while (t < SimTime_Counter)
            //{
            //    MRE.WaitOne();

            //    t++;
                output_filename = output_path + "\\Quantization_Output_" + t + ".txt";
                current_output = DO.Output_XML(output_filename);

                if (t > TimeLimit_4plot)
                {
                    ScaleTime_MinScale.Points.RemoveAt(0);
                    ScaleTime_Min_NonVacScale.Points.RemoveAt(0);
                    ScaleTime_Max_NonVacScale.Points.RemoveAt(0);
                    ScaleTime_MaxScale.Points.RemoveAt(0);
                }
                ScaleTime_MinScale.Points.Add(new DataPoint(t, current_output.Min_Scale));
                ScaleTime_Min_NonVacScale.Points.Add(new DataPoint(t, current_output.Min_NonVacScale));
                ScaleTime_Max_NonVacScale.Points.Add(new DataPoint(t, current_output.Max_NonVacScale));
                ScaleTime_MaxScale.Points.Add(new DataPoint(t, current_output.Max_Scale));
                // Classical energy plots
                ClassicalEnergy_Time.Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalEnergy));
                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalHamiltonianTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalHamiltonian_Energy[i_p]));
                }

                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalHamiltonian_wVacuum_Time_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalHamiltonian_wVacuum_Energy[i_p]));
                }

                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalLaplacianTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalLaplacian_Energy[i_p]));
                }

                ClassicalTotalPotentialEnergy_Time.Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalPotentialEnergy));
                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalPotentialEnergyTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.PotentialEnergy_Vector[i_p]));
                }

                ClassicalTotalKineticEnergy_Time.Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalKineticEnergy));
                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalKineticEnergyTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.KineticEnergy_Vector[i_p]));
                }

                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalEnergyExchangeTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalEnergy_Exchange[i_p]));
                }

                ClassicalTotalEnergy_Time.Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalEnergy));
                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    ClassicalEnergyTime_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(t, current_output.ClassicalEnergy_Vector[i_p]));
                }
                //
                NonVacScale_Boundary_Min_Laplacian.Points.Clear();
                NonVacScale_Boundary_Max_Laplacian.Points.Clear();
                NonVacScale_Boundary_Min_Laplacian.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, 0));
                NonVacScale_Boundary_Min_Laplacian.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));
                NonVacScale_Boundary_Max_Laplacian.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, 0));
                NonVacScale_Boundary_Max_Laplacian.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));

                NonVacScale_Boundary_Min_Commutator.Points.Clear();
                NonVacScale_Boundary_Max_Commutator.Points.Clear();
                NonVacScale_Boundary_Min_Commutator.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, -Max_Commutator_Energy));
                NonVacScale_Boundary_Min_Commutator.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, Max_Commutator_Energy));
                NonVacScale_Boundary_Max_Commutator.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, -Max_Commutator_Energy));
                NonVacScale_Boundary_Max_Commutator.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, Max_Commutator_Energy));

                NonVacScale_Boundary_Min.Points.Clear();
                NonVacScale_Boundary_Max.Points.Clear();
                NonVacScale_Boundary_Min.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, Min_Energy));
                NonVacScale_Boundary_Min.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));
                NonVacScale_Boundary_Max.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, Min_Energy));
                NonVacScale_Boundary_Max.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));

                NonVacScale_Boundary_Min_Mass.Points.Clear();
                NonVacScale_Boundary_Max_Mass.Points.Clear();
                NonVacScale_Boundary_Min_Mass.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, 0));
                NonVacScale_Boundary_Min_Mass.Points.Add(new OxyPlot.DataPoint((current_output.Min_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));
                NonVacScale_Boundary_Max_Mass.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, 0));
                NonVacScale_Boundary_Max_Mass.Points.Add(new OxyPlot.DataPoint((current_output.Max_NonVacScale / current_output.Max_Scale) * State.Num_ScaleBins, State.Num_Points));

                for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    Points_Data[i_p].Points.Clear();

                    if (State.Dimension == 1)
                        Points_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(State.Positions[i_p][0], 0));
                    else
                        Points_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(State.Positions[i_p][0], State.Positions[i_p][1]));


                    EnergyScale_Laplacian_Data[i_p].Points.Clear();
                    EnergyScale_Commutator_Data[i_p].Points.Clear();
                    EnergyScale_Data[i_p].Points.Clear();
                    MassScale_Data[i_p].Points.Clear();

                    EnergyScale_LaplacianDerivative_Data[i_p].Points.Clear();
                    EnergyScale_LaplacianDerivative_Data_smoothed[i_p].Points.Clear();
                    for (uint i_s = 0; i_s < State.Num_ScaleBins; i_s++)
                    {
                        //if (i_p > 0) //Excluding the Vacuum
                        EnergyScale_Laplacian_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, current_output.Laplacian_Energy[i_s][i_p]));

                        if ((i_s < (State.Num_ScaleBins - 1)) && (Smooth_flag))
                        {
                            EnergyScale_LaplacianDerivative_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, current_output.Laplacian_Energy_Derivative[i_s][i_p]));
                            EnergyScale_LaplacianDerivative_Data_smoothed[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, Math.Abs(current_output.Laplacian_Energy_Derivative_smoothed[i_s][i_p] - current_output.Laplacian_Energy_Derivative[i_s][i_p])));
                        }

                        EnergyScale_Commutator_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, current_output.Commutator_Energy[i_s][i_p]));
                        EnergyScale_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, current_output.Energy_Vector[i_s][i_p]));
                        MassScale_Data[i_p].Points.Add(new OxyPlot.Series.ScatterPoint(i_s, current_output.Mass_Vector[i_s][i_p]));

                        if (current_output.Commutator_Energy[i_s][i_p] > Max_Commutator_Energy)
                            Max_Commutator_Energy = current_output.Commutator_Energy[i_s][i_p];

                        if (current_output.Energy_Vector[i_s][i_p] < Min_Energy)
                            Min_Energy = current_output.Energy_Vector[i_s][i_p];
                    }
                }

                // Tile Plot Coloring
                TilePlot_Wiring.PointsPlot_Coloring(Sim_Variables.Tile.Tile_ScalarField);
                //TilePlot_Wiring.PointsPlot_Vectorizing(Sim_Variables.Tile.Tile_VectorField);

                //PointsPlot_Wiring.PointsPlot_Coloring(current_output.ClassicalEnergy_Exchange);
                //

                EnergyScale_Commutator_Plot.Axes[1].Maximum = Max_Commutator_Energy;
                EnergyScale_Commutator_Plot.Axes[1].Minimum = -Max_Commutator_Energy;
                EnergyScale_Plot.Axes[1].Minimum = Min_Energy;

                // Laplacian Energy Levels
                EnergyScale_Laplacian_HeatMap_Data = new double[State.Num_ScaleBins, State.Num_EnergyBins + 1];
                for (uint i_e = 0; i_e < (State.Num_EnergyBins + 1); i_e++)
                {
                    energy_laplacian_variable = (double)(i_e * State.Num_Points) / (double)State.Num_EnergyBins;
                    for (uint i_s = 0; i_s < State.Num_ScaleBins; i_s++)
                    {
                        EnergyScale_Laplacian_HeatMap_Data[i_s, i_e] = 0.0;
                        for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                        {
                            EnergyScale_Laplacian_HeatMap_Data[i_s, i_e] += -Math.Log(Math.Abs(energy_laplacian_variable - current_output.Laplacian_Energy[i_s][i_p]));
                        }

                        if (EnergyScale_Laplacian_HeatMap_Data[i_s, i_e] > Max_HeatMap_4plot)
                            EnergyScale_Laplacian_HeatMap_Data[i_s, i_e] = Max_HeatMap_4plot;
                    }
                }
                EnergyScale_Laplacian_HeatMap.Data = EnergyScale_Laplacian_HeatMap_Data;

                // Commutator Energy Levels
                EnergyScale_Commutator_HeatMap_Data = new double[State.Num_ScaleBins, State.Num_EnergyBins + 1];
                for (uint i_e = 0; i_e < (State.Num_EnergyBins + 1); i_e++)
                {
                    energy_commutator_variable = (double)(i_e * Max_Commutator_Energy * 2) / (double)State.Num_EnergyBins - Max_Commutator_Energy;
                    for (uint i_s = 0; i_s < State.Num_ScaleBins; i_s++)
                    {
                        EnergyScale_Commutator_HeatMap_Data[i_s, i_e] = 0.0;
                        for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                        {
                            EnergyScale_Commutator_HeatMap_Data[i_s, i_e] += -Math.Log(Math.Abs(energy_commutator_variable - current_output.Commutator_Energy[i_s][i_p]));
                        }

                        if (EnergyScale_Commutator_HeatMap_Data[i_s, i_e] > Max_HeatMap_4plot)
                            EnergyScale_Commutator_HeatMap_Data[i_s, i_e] = Max_HeatMap_4plot;
                    }
                }
                EnergyScale_Commutator_HeatMap.Data = EnergyScale_Commutator_HeatMap_Data;

                // Energy Levels
                EnergyScale_HeatMap_Data = new double[State.Num_ScaleBins, State.Num_EnergyBins + 1];
                for (uint i_e = 0; i_e < (State.Num_EnergyBins + 1); i_e++)
                {
                    energy_variable = (double)(i_e * (State.Num_Points - Min_Energy)) / (double)State.Num_EnergyBins + Min_Energy;
                    for (uint i_s = 0; i_s < State.Num_ScaleBins; i_s++)
                    {
                        EnergyScale_HeatMap_Data[i_s, i_e] = 0.0;
                        for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                        {
                            EnergyScale_HeatMap_Data[i_s, i_e] += -Math.Log(Math.Abs(energy_variable - current_output.Energy_Vector[i_s][i_p]));
                        }

                        if (EnergyScale_HeatMap_Data[i_s, i_e] > Max_HeatMap_4plot1)
                            EnergyScale_HeatMap_Data[i_s, i_e] = Max_HeatMap_4plot1;
                    }
                }
                EnergyScale_HeatMap.Data = EnergyScale_HeatMap_Data;

                // Mass Levels
                MassScale_HeatMap_Data = new double[State.Num_ScaleBins, State.Num_EnergyBins + 1];
                for (uint i_e = 0; i_e < (State.Num_EnergyBins + 1); i_e++)
                {
                    mass_variable = (double)(i_e * State.Num_Points) / (double)State.Num_EnergyBins;
                    for (uint i_s = 0; i_s < State.Num_ScaleBins; i_s++)
                    {
                        MassScale_HeatMap_Data[i_s, i_e] = 0.0;
                        for (uint i_p = 0; i_p < State.Num_Points; i_p++)
                        {
                            MassScale_HeatMap_Data[i_s, i_e] += -Math.Log(Math.Abs(mass_variable - current_output.Mass_Vector[i_s][i_p]));
                        }

                        //if (double.IsInfinity(MassScale_HeatMap_Data[i_s, i_e]))
                        //{
                        //    MassScale_HeatMap_Data[i_s, i_e] = double.NaN;
                        //}

                        if (MassScale_HeatMap_Data[i_s, i_e] > Max_HeatMap_4plot)
                            MassScale_HeatMap_Data[i_s, i_e] = Max_HeatMap_4plot;
                    }
                }
                MassScale_HeatMap.Data = MassScale_HeatMap_Data;

                ScaleTime_Plot.InvalidatePlot(true);
                Points_Plot.InvalidatePlot(true);
                Tile_Plot.InvalidatePlot(true);

                EnergyScale_Laplacian_Plot.InvalidatePlot(true);
                EnergyScale_Commutator_Plot.InvalidatePlot(true);
                EnergyScale_Plot.InvalidatePlot(true);
                MassScale_Plot.InvalidatePlot(true);
                ClassicalHamiltonianTime_Plot.InvalidatePlot(true);
                ClassicalHamiltonian_wVacuumTime_Plot.InvalidatePlot(true);
                ClassicalLaplacianTime_Plot.InvalidatePlot(true);
                ClassicalPotentialEnergyTime_Plot.InvalidatePlot(true);
                ClassicalKineticEnergyTime_Plot.InvalidatePlot(true);
                ClassicalEnergyExchangeTime_Plot.InvalidatePlot(true);
                ClassicalEnergyTime_Plot.InvalidatePlot(true);

                EnergyScale_LaplacianDerivative_Plot_smoothed.InvalidatePlot(true);
                EnergyScale_LaplacianDerivative_Plot.InvalidatePlot(true);

                EnergyScale_Laplacian_HeatMap_Plot.InvalidatePlot(true);
                EnergyScale_Commutator_HeatMap_Plot.InvalidatePlot(true);
                EnergyScale_HeatMap_Plot.InvalidatePlot(true);
                MassScale_HeatMap_Plot.InvalidatePlot(true);

                if (false)
                {
                    //int frequency = (int)Math.Round(1000 * Math.Sqrt(current_output.Laplacian_Energy[State.Num_ScaleBins / 2][1]));
                    int frequency = (int)Math.Round(700 * Math.Sqrt(current_output.Commutator_Energy[State.Num_ScaleBins / 2][State.Num_Points - 1]));
                    if (frequency >= 37)
                        Console.Beep(frequency, 50);
                }
                //watch.Stop();
                //var elapsed = watch.ElapsedMilliseconds/1000.0;
            //}
        }

        private void Run_Simulation()
        {
            while ((SimTime_Counter * Sim_Variables.dt) < Sim_Variables.Time_Range)
            {
                MRE.WaitOne();

                // Timing:
                //var watch = System.Diagnostics.Stopwatch.StartNew();

                if (double.IsNaN(output_Variables.Laplacian_Energy[0][0]))
                {
                    int x = 0;
                }
                else
                {
                    output_Variables = Sim_Variables.Evolve(State, Eigenvectors_flag, Perturb_flag, Smooth_flag);
                }

                // Test
                string Output_FileName = output_path + "\\Quantization_Output_" + SimTime_Counter + ".txt";
                DO.Write_XML(Output_FileName, output_Variables);


                Plot_Graphs(SimTime_Counter);


                SimTime_Counter++;
                backgroundWorker_Simulation.ReportProgress(SimTime_Counter);

            }
        }

        private void backgroundWorker_Simulation_DoWork(object sender, DoWorkEventArgs e)
        {
            Run_Simulation();
        }

        private void backgroundWorker_Simulation_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            simulation_ProgressBar.Value = (int)Math.Round(100 * SimTime_Counter * Sim_Variables.dt / Sim_Variables.Time_Range);
        }

    }

}
