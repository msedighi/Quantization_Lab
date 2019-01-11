using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OxyPlot;
using System.Windows.Forms;
using System.Drawing;

namespace Quantization_Tool
{
    public enum Plot_Choice { Points, Laplacian, Mass, Similarity, Commutator, ScaleTime };

    public class Plot_ToolStripDropDownButton : ToolStripDropDownButton
    {
        private ToolStripMenuItem pointsToolStripMenuItem;
        private ToolStripSeparator toolStripSeparator1;
        private ToolStripMenuItem laplacianEnergyToolStripMenuItem;
        private ToolStripMenuItem massToolStripMenuItem;
        private ToolStripMenuItem similarityOperatorToolStripMenuItem;
        private ToolStripMenuItem commutatorOperatorToolStripMenuItem;
        private ToolStripSeparator toolStripSeparator2;
        private ToolStripMenuItem scaleTimePlotToolStripMenuItem;

        public event EventHandler Plot_Chosen;

        public Plot_ToolStripDropDownButton()
        {
            pointsToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator1 = new ToolStripSeparator();
            laplacianEnergyToolStripMenuItem = new ToolStripMenuItem();
            massToolStripMenuItem = new ToolStripMenuItem();
            similarityOperatorToolStripMenuItem = new ToolStripMenuItem();
            commutatorOperatorToolStripMenuItem = new ToolStripMenuItem();
            toolStripSeparator2 = new ToolStripSeparator();
            scaleTimePlotToolStripMenuItem = new ToolStripMenuItem();
            //
            // Initiation
            //
            this.DisplayStyle = ToolStripItemDisplayStyle.Text;
            this.DropDownItems.AddRange(new ToolStripItem[] {
            this.pointsToolStripMenuItem,
            this.toolStripSeparator1,
            this.laplacianEnergyToolStripMenuItem,
            this.massToolStripMenuItem,
            this.similarityOperatorToolStripMenuItem,
            this.commutatorOperatorToolStripMenuItem,
            this.toolStripSeparator2,
            this.scaleTimePlotToolStripMenuItem});
            this.Size = new Size(90, 42);
            this.Text = "Plot";
            // 
            // pointsToolStripMenuItem
            // 
            this.pointsToolStripMenuItem.Size = new Size(389, 42);
            this.pointsToolStripMenuItem.Text = "Points";
            this.pointsToolStripMenuItem.Name = Plot_Choice.Points.ToString();
            this.pointsToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            // 
            // toolStripSeparator1
            // 
            this.toolStripSeparator1.Size = new Size(386, 6);
            //
            // ToolStrip DropDown Menu Items
            //
            this.laplacianEnergyToolStripMenuItem.Size = new Size(389, 42);
            this.laplacianEnergyToolStripMenuItem.Text = "Laplacian Operator";
            this.laplacianEnergyToolStripMenuItem.Name = Plot_Choice.Laplacian.ToString();
            this.laplacianEnergyToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            // 
            // massToolStripMenuItem
            // 
            this.massToolStripMenuItem.Size = new Size(389, 42);
            this.massToolStripMenuItem.Text = "Mass Operator";
            this.massToolStripMenuItem.Name = Plot_Choice.Mass.ToString();
            this.massToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            // 
            // similarityOperatorToolStripMenuItem
            // 
            this.similarityOperatorToolStripMenuItem.Size = new Size(389, 42);
            this.similarityOperatorToolStripMenuItem.Text = "Similarity Operator";
            this.similarityOperatorToolStripMenuItem.Name = Plot_Choice.Similarity.ToString();
            this.similarityOperatorToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            // 
            // commutatorOperatorToolStripMenuItem
            // 
            this.commutatorOperatorToolStripMenuItem.Size = new Size(389, 42);
            this.commutatorOperatorToolStripMenuItem.Text = "Commutator Operator";
            this.commutatorOperatorToolStripMenuItem.Name = Plot_Choice.Commutator.ToString();
            this.commutatorOperatorToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            // 
            // toolStripSeparator2
            // 
            this.toolStripSeparator2.Size = new Size(386, 6);
            // 
            // scaleTimePlotToolStripMenuItem
            // 
            this.scaleTimePlotToolStripMenuItem.Size = new Size(389, 42);
            this.scaleTimePlotToolStripMenuItem.Text = "Scale-Time Plot";
            this.scaleTimePlotToolStripMenuItem.Name = Plot_Choice.ScaleTime.ToString();
            this.scaleTimePlotToolStripMenuItem.Click += Plot_ToolStripMenuItem_Click;
            //

        }

        public void ChoosePlot(Plot_Choice p)
        {
            pointsToolStripMenuItem.CheckState = CheckState.Unchecked;
            laplacianEnergyToolStripMenuItem.CheckState = CheckState.Unchecked;
            massToolStripMenuItem.CheckState = CheckState.Unchecked;
            similarityOperatorToolStripMenuItem.CheckState = CheckState.Unchecked;
            commutatorOperatorToolStripMenuItem.CheckState = CheckState.Unchecked;
            scaleTimePlotToolStripMenuItem.CheckState = CheckState.Unchecked;

            switch (p)
            {
                case Plot_Choice.Points:
                    pointsToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
                case Plot_Choice.Laplacian:
                    laplacianEnergyToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
                case Plot_Choice.Mass:
                    massToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
                case Plot_Choice.Similarity:
                    similarityOperatorToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
                case Plot_Choice.Commutator:
                    commutatorOperatorToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
                case Plot_Choice.ScaleTime:
                    scaleTimePlotToolStripMenuItem.CheckState = CheckState.Checked;
                    break;
            }

            this.Plot_Chosen(p, new EventArgs());

        }

        private void Plot_ToolStripMenuItem_Click(object sender, EventArgs e)
        {
            ToolStripMenuItem Chosen_PlotMenuItem = sender as ToolStripMenuItem;
            switch (Chosen_PlotMenuItem.Name)
            {
                case "Points":
                    ChoosePlot(Plot_Choice.Points);
                    break;
                case "Laplacian":
                    ChoosePlot(Plot_Choice.Laplacian);
                    break;
                case "Mass":
                    ChoosePlot(Plot_Choice.Mass);
                    break;
                case "Similarity":
                    ChoosePlot(Plot_Choice.Similarity);
                    break;
                case "Commutator":
                    ChoosePlot(Plot_Choice.Commutator);
                    break;
                case "ScaleTime":
                    ChoosePlot(Plot_Choice.ScaleTime);
                    break;
            }
        }
    }

    public class View_ToolStripDropDownButton : ToolStripDropDownButton
    {
        internal ToolStripMenuItem optionsPlotToolStripMenuItem;

        public View_ToolStripDropDownButton()
        {
            this.optionsPlotToolStripMenuItem = new ToolStripMenuItem();

            this.optionsPlotToolStripMenuItem.Text = "Options";
            this.optionsPlotToolStripMenuItem.CheckOnClick = true;

            this.DisplayStyle = ToolStripItemDisplayStyle.Text;
            this.DropDownItems.AddRange(new ToolStripItem[] {
            this.optionsPlotToolStripMenuItem});
            this.Text = "View";

        }
    }

    public class Plot_ToolStrip : ToolStrip
    {
        private ToolStripSeparator toolStripSeparator1;

        public Plot_ToolStrip()
        {
            this.toolStripSeparator1 = new ToolStripSeparator();

            //this.SuspendLayout();

            //this.Dock = DockStyle.None;

            //this.ResumeLayout(true);
        }

        public void Attach(Plot_ToolStripDropDownButton plot_button, View_ToolStripDropDownButton view_button)
        {
            this.Items.AddRange(new ToolStripItem[]
            {
                plot_button,
                this.toolStripSeparator1,
                view_button
            });
        }
    }

    public class Plot_SplitContainer : SplitContainer
    {
        private TableLayoutPanel Options_Table;
        private OxyPlot.WindowsForms.PlotView Plot_View;

        public Plot_SplitContainer()
        {
            Options_Table = new TableLayoutPanel();
            Plot_View = new OxyPlot.WindowsForms.PlotView();

            //
            this.Options_Table.SuspendLayout();
            this.Panel1.SuspendLayout();
            this.Panel2.SuspendLayout();
            this.SuspendLayout();
            //
            // 
            this.BorderStyle = BorderStyle.FixedSingle;
            this.FixedPanel = FixedPanel.Panel1;
            this.IsSplitterFixed = true;
            this.SplitterDistance = 100;
            // 
            // Panel1
            // 
            this.Panel1.Controls.Add(Options_Table);
            this.Panel1MinSize = 100;
            this.Panel1Collapsed = true;
            // 
            // Panel2
            // 
            this.Panel2.Controls.Add(Plot_View);
            //
            // Options Table
            //
            Options_Table.ColumnCount = 1;
            Options_Table.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100F));
            Options_Table.Dock = DockStyle.Fill;
            Options_Table.Location = new Point(0, 0);
            Options_Table.RowCount = 2;
            Options_Table.RowStyles.Add(new RowStyle(SizeType.Absolute, 50F));
            Options_Table.RowStyles.Add(new RowStyle(SizeType.Percent, 100F));
            //
            // Plot View
            //
            Plot_View.Dock = DockStyle.Fill;
            Plot_View.Location = new Point(0, 0);
            Plot_View.PanCursor = Cursors.Hand;
            Plot_View.ZoomHorizontalCursor = Cursors.SizeWE;
            Plot_View.ZoomRectangleCursor = Cursors.SizeNWSE;
            Plot_View.ZoomVerticalCursor = Cursors.SizeNS;
            //
            // ToolStrip
            //
            this.Options_Table.ResumeLayout(true);
            this.Panel1.ResumeLayout(false);
            this.Panel2.ResumeLayout(false);
            this.ResumeLayout(false);
            //

        }

        public void Attach(PlotModel plotModel)
        {
            this.Plot_View.Model = plotModel;
        }

        public void Options_Panel(bool state)
        {
            this.Panel1Collapsed = state;
        }
    }

    public class Plot_ToolStripContainer : ToolStripContainer
    {
        private Plot_SplitContainer Plot_Container;
        private Plot_ToolStrip ToolStrip_Plot;
        private Plot_ToolStripDropDownButton ToolStripDropDownButton_Plot;
        private View_ToolStripDropDownButton ToolStripDropDownButton_View;

        public event EventHandler Plot_Changed;
        private Plot_Choice chosen_plot;
        public Plot_Choice Chosen_Plot
        {
            get
            { return chosen_plot; }
            set
            {
                if (chosen_plot != value)
                {
                    chosen_plot = value;
                    this.Plot_Changed(chosen_plot, new EventArgs());
                }
            }
        }

        public Plot_ToolStripContainer()
        {
            this.Plot_Container = new Plot_SplitContainer();
            this.ToolStrip_Plot = new Plot_ToolStrip();
            this.ToolStripDropDownButton_Plot = new Plot_ToolStripDropDownButton();
            this.ToolStripDropDownButton_View = new View_ToolStripDropDownButton();

            this.SuspendLayout();
            this.Plot_Container.SuspendLayout();
            this.ToolStrip_Plot.SuspendLayout();

            this.LeftToolStripPanelVisible = false;
            this.RightToolStripPanelVisible = false;
            this.BottomToolStripPanelVisible = false;

            this.Plot_Container.Dock = DockStyle.Fill;
            this.ContentPanel.Controls.Add(Plot_Container);

            this.ToolStrip_Plot.Attach(ToolStripDropDownButton_Plot, ToolStripDropDownButton_View);

            this.TopToolStripPanel.Controls.Add(ToolStrip_Plot);

            // Event
            this.ToolStripDropDownButton_Plot.Plot_Chosen += ToolStripDropDownButton_Plot_Chosen; ;
            //
            this.ToolStrip_Plot.ResumeLayout(true);
            this.Plot_Container.ResumeLayout(true);
            this.ResumeLayout(true);

            ToolStripDropDownButton_View.optionsPlotToolStripMenuItem.Click += OptionsPlotToolStripMenuItem_Click;
        }

        private void OptionsPlotToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (ToolStripDropDownButton_View.optionsPlotToolStripMenuItem.Checked)
                Plot_Container.Panel1Collapsed = false;
            else
                Plot_Container.Panel1Collapsed = true;
        }

        private void ToolStripDropDownButton_Plot_Chosen(object sender, EventArgs e)
        {
            Plot_Choice p = (Plot_Choice)sender;

            this.Chosen_Plot = p;
        }

        public void Attach(PlotModel plotModel)
        {
            this.Plot_Container.Attach(plotModel);
        }

        public void Options_Panel(bool state)
        {
            this.Plot_Container.Panel1Collapsed = state;
        }
    }

    internal class Plot_ControlView_Wiring
    {
        private PlotModel plotModel;
        private Plot_ToolStripContainer plotContainer;

        public Plot_ControlView_Wiring(PlotModel pm, Plot_ToolStripContainer pc)
        {
            plotModel = pm;
            plotContainer = pc;
        }

        public void Apply(Plot_Choice p)
        {
            switch (p)
            {
                case Plot_Choice.Points:
                    break;
                case Plot_Choice.Laplacian:
                    break;
                case Plot_Choice.Mass:
                    break;
                case Plot_Choice.Similarity:
                    break;
                case Plot_Choice.Commutator:
                    break;
                case Plot_Choice.ScaleTime:
                    break;
            }

        }
    }

    internal class DataPlotWiring
    {
        protected PlotModel plotModel;
        protected OxyPlot.Series.ScatterSeries[] plotData;
        protected int numPoints = 0;

        public bool Activated;

        public DataPlotWiring(PlotModel pm, OxyPlot.Series.ScatterSeries[] pd)
        {
            plotData = pd;
            plotModel = pm;
            Activated = false;

            numPoints = plotData.Count();
        }

        public void Activate()
        {
            for (uint i_p = 0; i_p < numPoints; i_p++)
            {
                plotModel.Series.Add(plotData[i_p]);
            }
            this.Activated = true;
        }

        public void Deactivate()
        {
            this.plotModel.Series.Clear();
            this.Activated = false;
        }
    }

    internal class EnergyTimePlotWiring : DataPlotWiring
    {
        private OxyPlot.Series.ScatterSeries Total_Energy;
        public bool ClassicalEnergy_Line { get; private set; }

        public EnergyTimePlotWiring(PlotModel pm, OxyPlot.Series.ScatterSeries[] energy, OxyPlot.Series.ScatterSeries total_energy) : base(pm, energy)
        {
            Total_Energy = total_energy;
            ClassicalEnergy_Line = false;
        }

        public void Show_ClassicalEnergy()
        {
            plotModel.Series.Add(Total_Energy);
            ClassicalEnergy_Line = true;
        }

        public void Hide_ClassicalEnergy()
        {
            plotModel.Series.Remove(Total_Energy);
            ClassicalEnergy_Line = false;
        }
    }

    internal class PointsPlotWiring : DataPlotWiring
    {
        public PointsPlotWiring(PlotModel pm, OxyPlot.Series.ScatterSeries[] pd) : base(pm, pd) { }

        public void PointsPlot_Coloring(double[][,] transformation, int i_s, int selected_index)
        {
            plotModel.Annotations.Clear();
            for (int i_p = 0; i_p < numPoints; i_p++)
            {
                plotData[i_p].MarkerFill = OxyPalettes.BlueWhiteRed(200).Colors[100 + (int)Math.Round(transformation[i_s][i_p, selected_index] * 99)];
                plotData[i_p].LabelFormatString = Math.Round(transformation[i_s][i_p, selected_index] * 100).ToString();
            }
            plotModel.InvalidatePlot(true);
        }

        public void PointsPlot_Coloring(int selected_index)
        {
            plotModel.Annotations.Clear();
            for (int i_p = 0; i_p < numPoints; i_p++)
            {
                if (i_p == selected_index)
                {
                    plotData[i_p].MarkerFill = OxyColors.Black;
                    plotData[i_p].LabelFormatString = "100";
                }
                else
                {
                    plotData[i_p].MarkerFill = OxyColors.WhiteSmoke;
                    plotData[i_p].LabelFormatString = "0";
                }
            }
            plotModel.InvalidatePlot(true);
        }

        public void PointsPlot_Vectorizing(double[][,] transformation_x, double[][,] transformation_y, int i_s, int selected_index)
        {
            OxyPlot.Annotations.ArrowAnnotation[] vector_field = new OxyPlot.Annotations.ArrowAnnotation[numPoints];
            plotModel.Annotations.Clear();
            for (int i_p = 0; i_p < numPoints; i_p++)
            {
                vector_field[i_p] = new OxyPlot.Annotations.ArrowAnnotation();
                vector_field[i_p].StartPoint = new DataPoint(plotData[i_p].Points[0].X, plotData[i_p].Points[0].Y);
                vector_field[i_p].EndPoint = new DataPoint(plotData[i_p].Points[0].X + transformation_x[i_s][i_p, selected_index], plotData[i_p].Points[0].Y + transformation_y[i_s][i_p, selected_index]);
                vector_field[i_p].HeadLength = 4;
                vector_field[i_p].HeadWidth = 1;

                double vector_length = Math.Sqrt(transformation_x[i_s][i_p, selected_index] * transformation_x[i_s][i_p, selected_index] + transformation_y[i_s][i_p, selected_index] * transformation_y[i_s][i_p, selected_index]);
                vector_field[i_p].Color = OxyPalettes.Gray(100).Colors[99 - (int)Math.Round(vector_length * 99)];
                plotData[i_p].MarkerFill = OxyPalettes.Gray(100).Colors[99 - (int)Math.Round(vector_length * 99)];

                plotModel.Annotations.Add(vector_field[i_p]);
            }
            plotModel.InvalidatePlot(true);
        }

    }



}
