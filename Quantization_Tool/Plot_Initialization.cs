using System;
using System.Windows.Forms;
using OxyPlot;

namespace Quantization_Tool
{
    public partial class Form1 : Form
    {
        private void Initialize_PointsPlot(PlotModel points_plot)
        {
            //points_plot.PlotAreaBorderThickness = new OxyThickness(2.0);
            points_plot.PlotMargins = new OxyThickness(0.0);
            points_plot.Title = "Points";
            points_plot.TitleFontWeight = 5;

            OxyPlot.Axes.LinearAxis x = new OxyPlot.Axes.LinearAxis();
            x.PositionAtZeroCrossing = true;
            x.Position = OxyPlot.Axes.AxisPosition.Bottom;
            x.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            x.TextColor = OxyColors.Transparent;
            points_plot.Axes.Add(x);

            OxyPlot.Axes.LinearAxis y = new OxyPlot.Axes.LinearAxis();
            y.PositionAtZeroCrossing = true;
            y.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            y.TextColor = OxyColors.Transparent;
            points_plot.Axes.Add(y);
        }
        private void Initialize_ScaleTimePlot(PlotModel scaletime_plot)
        {
            //scaletime_plot.PlotAreaBorderThickness = new OxyThickness(0.0);
            scaletime_plot.PlotMargins = new OxyThickness(8.0);
            scaletime_plot.Title = "Scale vs Time";
            scaletime_plot.TitleFontWeight = 5;

            OxyPlot.Axes.LinearAxis x = new OxyPlot.Axes.LinearAxis();
            //x.Maximum = Sim_Variables.Time_Range;
            //x.Minimum = 0.0;
            x.PositionAtZeroCrossing = false;
            x.Position = OxyPlot.Axes.AxisPosition.Bottom;
            //x.MajorGridlineStyle = LineStyle.Solid;
            //x.MinorGridlineStyle = LineStyle.Dot;
            x.AxislineThickness = 0.5;
            x.AxislineColor = OxyColors.Black;
            x.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            x.AxisTickToLabelDistance = 0.0;
            x.FontSize = 8.0;
            scaletime_plot.Axes.Add(x);

            OxyPlot.Axes.LinearAxis y = new OxyPlot.Axes.LinearAxis();
            y.Minimum = 0.0;
            y.PositionAtZeroCrossing = true;
            y.MajorGridlineStyle = LineStyle.Solid;
            y.MinorGridlineStyle = LineStyle.Dot;
            y.AxislineThickness = 0.5;
            y.AxislineColor = OxyColors.Black;
            y.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            y.AxisTickToLabelDistance = 0.0;
            y.FontSize = 8.0;
            scaletime_plot.Axes.Add(y);
        }

        private void Initialize_ClassicalEnergyTimePlot(PlotModel scaletime_plot, string energy_name)
        {
            scaletime_plot.Title = energy_name + " vs Time";
            scaletime_plot.TitleFontWeight = 5;

            Initialize_EnergyTimePlot_Base(scaletime_plot);
        }
        private void Initialize_EnergyTimePlot_Base(PlotModel scaletime_plot)
        {
            //scaletime_plot.PlotAreaBorderThickness = new OxyThickness(0.0);
            scaletime_plot.PlotMargins = new OxyThickness(8.0);

            OxyPlot.Axes.LinearAxis x = new OxyPlot.Axes.LinearAxis();
            x.PositionAtZeroCrossing = false;
            x.Position = OxyPlot.Axes.AxisPosition.Bottom;
            x.AxislineThickness = 0.5;
            x.AxislineColor = OxyColors.Black;
            x.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            x.AxisTickToLabelDistance = 0.0;
            x.FontSize = 8.0;
            scaletime_plot.Axes.Add(x);

            OxyPlot.Axes.LinearAxis y = new OxyPlot.Axes.LinearAxis();
            //y.Minimum = 0.0;
            y.PositionAtZeroCrossing = true;
            y.MajorGridlineStyle = LineStyle.Solid;
            y.MinorGridlineStyle = LineStyle.Dot;
            y.AxislineThickness = 0.5;
            y.AxislineColor = OxyColors.Black;
            y.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            y.AxisTickToLabelDistance = 0.0;
            y.FontSize = 8.0;
            scaletime_plot.Axes.Add(y);
        }

        private void Initialize_EnergyScalePlot_Base_wBoundary(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, OxyPlot.Series.LineSeries boundary_max, OxyPlot.Series.LineSeries boundary_min, State_Variables state, bool color_flag, bool scaleBoundary_flag)
        {
            Initialize_EnergyScalePlot_Base(energyscale_plot, energyscale_data, state, color_flag);

            boundary_max.Color = OxyColors.Red;
            boundary_min.Color = OxyColors.Red;
            boundary_max.StrokeThickness = 1;
            boundary_min.StrokeThickness = 1;
            if (scaleBoundary_flag)
            {
                energyscale_plot.Series.Add(boundary_max);
                energyscale_plot.Series.Add(boundary_min);
            }
        }
        private void Initialize_EnergyScalePlot_Base(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, State_Variables state, bool color_flag)
        {
            //scaletime_plot.PlotAreaBorderThickness = new OxyThickness(0.0);
            energyscale_plot.PlotMargins = new OxyThickness(8.0);
            energyscale_plot.TitleFontWeight = 5;

            OxyPlot.Axes.LinearAxis x = new OxyPlot.Axes.LinearAxis();
            x.Maximum = State.Num_ScaleBins;
            x.Minimum = 0.0;
            x.PositionAtZeroCrossing = true;
            x.Position = OxyPlot.Axes.AxisPosition.Bottom;
            x.AxislineThickness = 0.5;
            x.AxislineColor = OxyColors.Black;
            x.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            x.AxisTickToLabelDistance = 0.0;
            x.FontSize = 8.0;
            energyscale_plot.Axes.Add(x);

            OxyPlot.Axes.LinearAxis y = new OxyPlot.Axes.LinearAxis();
            y.PositionAtZeroCrossing = true;
            y.MajorGridlineStyle = LineStyle.Solid;
            y.MinorGridlineStyle = LineStyle.Dot;
            y.AxislineThickness = 0.5;
            y.AxislineColor = OxyColors.Black;
            y.TickStyle = OxyPlot.Axes.TickStyle.Crossing;
            y.AxisTickToLabelDistance = 0.0;
            y.FontSize = 8.0;
            energyscale_plot.Axes.Add(y);

            for (int i = 0; i < state.Num_Points; i++)
            {
                energyscale_data[i] = new OxyPlot.Series.ScatterSeries();
                energyscale_data[i].MarkerType = MarkerType.Circle;
                energyscale_data[i].MarkerSize = 2;
                if (!color_flag)
                    energyscale_data[i].MarkerFill = OxyColors.Black;

                energyscale_plot.Series.Add(energyscale_data[i]);
            }

            //energyscale_plot.MouseLeave += EnergyScalePlot_MouseLeave;
        }

        private void Initialize_MassScalePlot(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, OxyPlot.Series.LineSeries boundary_max, OxyPlot.Series.LineSeries boundary_min, State_Variables state, bool color_flag, bool scaleBoundary_flag)
        {
            Initialize_EnergyScalePlot_Base_wBoundary(energyscale_plot, energyscale_data, boundary_max, boundary_min, state, color_flag, scaleBoundary_flag);
            energyscale_plot.Title = "Mass vs Scale";
            energyscale_plot.Axes[1].Maximum = State.Num_Points;
            energyscale_plot.Axes[1].Minimum = 0.0;

            energyscale_plot.MouseDown += MassScalePlot_MouseDown;
            energyscale_plot.MouseMove += MassScalePlot_MouseMove;
        }
        private void Initialize_EnergyScalePlot(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, OxyPlot.Series.LineSeries boundary_max, OxyPlot.Series.LineSeries boundary_min, State_Variables state, bool color_flag, bool scaleBoundary_flag)
        {
            Initialize_EnergyScalePlot_Base_wBoundary(energyscale_plot, energyscale_data, boundary_max, boundary_min, state, color_flag, scaleBoundary_flag);
            energyscale_plot.Title = "Energy vs Scale";
            energyscale_plot.Axes[1].Maximum = State.Num_Points;
            energyscale_plot.Axes[1].Minimum = 0;

            energyscale_plot.MouseDown += EnergyScalePlot_MouseDown;
            energyscale_plot.MouseMove += EnergyScalePlot_MouseMove;
        }
        private void Initialize_LaplacianEnergyScalePlot(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, OxyPlot.Series.LineSeries boundary_max, OxyPlot.Series.LineSeries boundary_min, State_Variables state, bool color_flag, bool scaleBoundary_flag)
        {
            Initialize_EnergyScalePlot_Base_wBoundary(energyscale_plot, energyscale_data, boundary_max, boundary_min, state, color_flag, scaleBoundary_flag);
            energyscale_plot.Title = "Laplacian Energy vs Scale";
            energyscale_plot.Axes[1].Maximum = State.Num_Points;
            energyscale_plot.Axes[1].Minimum = 0.0;

            energyscale_plot.MouseDown += LaplacianEnergyScalePlot_MouseDown;
            energyscale_plot.MouseMove += LaplacianEnergyScalePlot_MouseMove;
        }
        private void Initialize_CommutatorEnergyScalePlot(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, OxyPlot.Series.LineSeries boundary_max, OxyPlot.Series.LineSeries boundary_min, State_Variables state, bool color_flag, bool scaleBoundary_flag)
        {
            Initialize_EnergyScalePlot_Base_wBoundary(energyscale_plot, energyscale_data, boundary_max, boundary_min, state, color_flag, scaleBoundary_flag);
            energyscale_plot.Title = "Commutator Energy vs Scale";
            energyscale_plot.Axes[1].Maximum = State.Num_Points;
            energyscale_plot.Axes[1].Minimum = -State.Num_Points;

            energyscale_plot.MouseDown += CommutatorEnergyScalePlot_MouseDown;
            energyscale_plot.MouseMove += CommutatorEnergyScalePlot_MouseMove;
        }

        private void Initialize_LaplacianEnergyDerivativeScalePlot(PlotModel energyscale_plot, OxyPlot.Series.ScatterSeries[] energyscale_data, State_Variables state)
        {
            Initialize_EnergyScalePlot_Base(energyscale_plot, energyscale_data, state, true);
            energyscale_plot.Title = "Energy Derivative vs Scale";
            energyscale_plot.Axes[1].Maximum = 0.01;
            energyscale_plot.Axes[1].Minimum = 0;

            energyscale_plot.MouseDown += LaplacianEnergyScalePlot_MouseDown;
            energyscale_plot.MouseMove += LaplacianEnergyScalePlot_MouseMove;
        }

        private void Initialize_EnergyScaleHeatMap(PlotModel energyscale_heatmap_plot, OxyPlot.Series.HeatMapSeries energyscale_heatmap, State_Variables state)
        {
            //scaletime_plot.PlotAreaBorderThickness = new OxyThickness(0.0);
            energyscale_heatmap_plot.PlotMargins = new OxyThickness(0.0);
            energyscale_heatmap_plot.TitleFontWeight = 5;

            OxyPlot.Axes.LinearColorAxis xy = new OxyPlot.Axes.LinearColorAxis();
            xy.Palette = OxyPalettes.Jet(500);
            //xy.HighColor = OxyColors.Gray;
            //xy.LowColor = OxyColors.Black;
            //xy.Position = OxyPlot.Axes.AxisPosition.Right;
            energyscale_heatmap_plot.Axes.Add(xy);

            energyscale_heatmap.Interpolate = true;
            energyscale_heatmap.RenderMethod = OxyPlot.Series.HeatMapRenderMethod.Bitmap;
            energyscale_heatmap_plot.Series.Add(energyscale_heatmap);

            energyscale_heatmap.X0 = 0.0;
            energyscale_heatmap.X1 = state.Num_ScaleBins;
        }


    }
}