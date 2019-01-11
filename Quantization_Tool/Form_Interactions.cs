using System;
using System.Windows.Forms;
using OxyPlot;

namespace Quantization_Tool
{
    public partial class Form1 : Form
    {

        private int EnergyScalePlot_MouseDown_Base(PlotModel energyscale_plot, OxyMouseEventArgs e)
        {
            int i_s = (int)energyscale_plot.Axes[0].InverseTransform(e.Position.X);
            double y = energyscale_plot.Axes[1].InverseTransform(e.Position.Y);

            int selected_index = 0;
            double closest_dist = double.PositiveInfinity;
            for (int i_p = 0; i_p < State.Num_Points; i_p++)
            {
                double dist = Math.Abs((energyscale_plot.Series[i_p] as OxyPlot.Series.ScatterSeries).Points[i_s].Y - y);
                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    selected_index = i_p;
                }
            }

            if (Selected_Point != null)
            {
                Selected_Point.Size = 2;
                Selected_Plot.InvalidatePlot(true);
            }

            Selected_Plot = energyscale_plot;
            Selected_Series = energyscale_plot.Series[selected_index] as OxyPlot.Series.ScatterSeries;
            Selected_Point = Selected_Series.Points[i_s];
            Selected_Index = selected_index;

            Selected_Point.Size = 4;
            label_X.Text = "X: " + i_s;
            label_Y.Text = ", Y: " + Selected_Point.Y;

            return i_s;
        }
        private int EnergyScalePlot_MouseMove_Base(PlotModel energyscale_plot, OxyMouseEventArgs e)
        {
            Selected_Point.Size = 2;

            int i_s = (int)energyscale_plot.Axes[0].InverseTransform(e.Position.X);
            if (i_s > (State.Num_ScaleBins - 1))
                i_s = (State.Num_ScaleBins - 1);
            else if (i_s < 0)
                i_s = 0;

            Selected_Point = Selected_Series.Points[i_s];

            Selected_Point.Size = 4;
            Selected_Plot.InvalidatePlot(true);

            label_X.Text = "X: " + i_s;
            label_Y.Text = ", Y: " + Selected_Point.Y;

            return i_s;
        }
        private void EnergyScalePlot_MouseLeave(object sender, OxyMouseEventArgs e)
        {
            if (Selected_Point != null)
            {
                Selected_Point.Size = 2;
                Selected_Plot.InvalidatePlot(true);
                Selected_Point = null;
                Selected_Series = null;

                Points_Plot.Annotations.Clear();
                for (int i_p = 0; i_p < State.Num_Points; i_p++)
                {
                    Points_Data[i_p].MarkerFill = OxyColors.Black;
                    Points_Data[i_p].LabelFormatString = "";
                }
                Points_Plot.InvalidatePlot(true);
            }

        }

        private void MassScalePlot_MouseDown(object sender, OxyMouseDownEventArgs e)
        {
            EnergyScalePlot_MouseDown_Base(MassScale_Plot, e);
            // Eigen Vectors
            PointsPlot_Wiring.PointsPlot_Coloring(Selected_Index);
        }
        private void MassScalePlot_MouseMove(object sender, OxyMouseEventArgs e)
        {
            if ((Selected_Series != null) && (Selected_Plot == MassScale_Plot))
            {
                int scale_index = EnergyScalePlot_MouseMove_Base(MassScale_Plot, e);
                // Eigen Vectors
                PointsPlot_Wiring.PointsPlot_Coloring(Selected_Index);
            }
        }

        private void LaplacianEnergyScalePlot_MouseDown(object sender, OxyMouseDownEventArgs e)
        {
            int scale_index = EnergyScalePlot_MouseDown_Base(EnergyScale_Laplacian_Plot, e);
            // Eigen Vectors
            PointsPlot_Wiring.PointsPlot_Coloring(output_Variables.Laplacian_Orthonormal_Transformation, scale_index, Selected_Index);
        }
        private void LaplacianEnergyScalePlot_MouseMove(object sender, OxyMouseEventArgs e)
        {
            if ((Selected_Series != null) && (Selected_Plot == EnergyScale_Laplacian_Plot))
            {
                int scale_index = EnergyScalePlot_MouseMove_Base(EnergyScale_Laplacian_Plot, e);
                // Eigen Vectors
                PointsPlot_Wiring.PointsPlot_Coloring(output_Variables.Laplacian_Orthonormal_Transformation, scale_index, Selected_Index);
            }
        }

        private void EnergyScalePlot_MouseDown(object sender, OxyMouseDownEventArgs e)
        {
            int scale_index = EnergyScalePlot_MouseDown_Base(EnergyScale_Plot, e);
            // Eigen Vectors
            PointsPlot_Wiring.PointsPlot_Coloring(output_Variables.Energy_Orthonormal_Transformation, scale_index, Selected_Index);
        }
        private void EnergyScalePlot_MouseMove(object sender, OxyMouseEventArgs e)
        {
            if ((Selected_Series != null) && (Selected_Plot == EnergyScale_Plot))
            {
                int scale_index = EnergyScalePlot_MouseMove_Base(EnergyScale_Plot, e);
                // Eigen Vectors
                PointsPlot_Wiring.PointsPlot_Coloring(output_Variables.Energy_Orthonormal_Transformation, scale_index, Selected_Index);
            }
        }

        private void CommutatorEnergyScalePlot_MouseDown(object sender, OxyMouseDownEventArgs e)
        {
            int scale_index = EnergyScalePlot_MouseDown_Base(EnergyScale_Commutator_Plot, e);
            // Eigen Vectors
            PointsPlot_Wiring.PointsPlot_Vectorizing(output_Variables.Commutator_Orthonormal_Transformation_Real, output_Variables.Commutator_Orthonormal_Transformation_Imag, scale_index, Selected_Index);
        }
        private void CommutatorEnergyScalePlot_MouseMove(object sender, OxyMouseEventArgs e)
        {
            if ((Selected_Series != null) && (Selected_Plot == EnergyScale_Commutator_Plot))
            {
                int scale_index = EnergyScalePlot_MouseMove_Base(EnergyScale_Commutator_Plot, e);
                // Eigen Vectors
                PointsPlot_Wiring.PointsPlot_Vectorizing(output_Variables.Commutator_Orthonormal_Transformation_Real, output_Variables.Commutator_Orthonormal_Transformation_Imag, scale_index, Selected_Index);
            }
        }

        private void button_Run_Click(object sender, EventArgs e)
        {
            RunButton_State = !RunButton_State;
            if (RunButton_State)
            {
                button_Run.Text = "Stop";
                MRE.Set();
            }
            else
            {
                button_Run.Text = "Run";
                //Thread.Sleep(2000);
                MRE.Reset();
            }

        }

        private void toolStripMenuItem_SimOptions_Click(object sender, EventArgs e)
        {
            if (toolStripMenuItem_SimOptions.Checked)
                this.tableLayoutPanel_Back.Controls.Add(this.tableLayoutPanel_bottom, 0, 1);
            else
                this.tableLayoutPanel_Back.Controls.Remove(this.tableLayoutPanel_bottom);
        }

        private void toolStripMenuItem_PlotOptions_Click(object sender, EventArgs e)
        {
            if (toolStripMenuItem_PlotOptions.Checked)
            {
                plot_ToolStripContainer_TopLeft.TopToolStripPanelVisible = true;
                plot_ToolStripContainer_TopMiddle.TopToolStripPanelVisible = true;
                plot_ToolStripContainer_TopRight.TopToolStripPanelVisible = true;
                plot_ToolStripContainer_BottomLeft.TopToolStripPanelVisible = true;
                plot_ToolStripContainer_BottomMiddle.TopToolStripPanelVisible = true;
                plot_ToolStripContainer_BottomRight.TopToolStripPanelVisible = true;
            }
            else
            {
                plot_ToolStripContainer_TopLeft.TopToolStripPanelVisible = false;
                plot_ToolStripContainer_TopMiddle.TopToolStripPanelVisible = false;
                plot_ToolStripContainer_TopRight.TopToolStripPanelVisible = false;
                plot_ToolStripContainer_BottomLeft.TopToolStripPanelVisible = false;
                plot_ToolStripContainer_BottomMiddle.TopToolStripPanelVisible = false;
                plot_ToolStripContainer_BottomRight.TopToolStripPanelVisible = false;
            }
        }

    }
}