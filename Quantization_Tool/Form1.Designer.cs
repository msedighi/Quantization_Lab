namespace Quantization_Tool
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(Form1));
            this.button_Run = new System.Windows.Forms.Button();
            this.tableLayoutPanel_Plots = new System.Windows.Forms.TableLayoutPanel();
            this.plot_ToolStripContainer_TopMiddle = new Quantization_Tool.Plot_ToolStripContainer();
            this.topLeft_panel = new System.Windows.Forms.Panel();
            this.plot_ToolStripContainer_TopLeft = new Quantization_Tool.Plot_ToolStripContainer();
            this.plot_ToolStripContainer_TopRight = new Quantization_Tool.Plot_ToolStripContainer();
            this.plot_ToolStripContainer_BottomLeft = new Quantization_Tool.Plot_ToolStripContainer();
            this.plot_ToolStripContainer_BottomMiddle = new Quantization_Tool.Plot_ToolStripContainer();
            this.plot_ToolStripContainer_BottomRight = new Quantization_Tool.Plot_ToolStripContainer();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.tabControl1 = new System.Windows.Forms.TabControl();
            this.tabPage2 = new System.Windows.Forms.TabPage();
            this.panel_Plot = new System.Windows.Forms.Panel();
            this.plot_plotView = new OxyPlot.WindowsForms.PlotView();
            this.tableLayoutPanel_Back = new System.Windows.Forms.TableLayoutPanel();
            this.tableLayoutPanel_bottom = new System.Windows.Forms.TableLayoutPanel();
            this.flowLayoutPanel1 = new System.Windows.Forms.FlowLayoutPanel();
            this.label_X = new System.Windows.Forms.Label();
            this.label_Y = new System.Windows.Forms.Label();
            this.simulation_ProgressBar = new System.Windows.Forms.ProgressBar();
            this.button_Pause = new System.Windows.Forms.Button();
            this.BottomToolStripPanel = new System.Windows.Forms.ToolStripPanel();
            this.TopToolStripPanel = new System.Windows.Forms.ToolStripPanel();
            this.RightToolStripPanel = new System.Windows.Forms.ToolStripPanel();
            this.LeftToolStripPanel = new System.Windows.Forms.ToolStripPanel();
            this.ContentPanel = new System.Windows.Forms.ToolStripContentPanel();
            this.toolStrip_Main = new System.Windows.Forms.ToolStrip();
            this.toolStripDropDownButton_File = new System.Windows.Forms.ToolStripDropDownButton();
            this.toolStripMenuItem_NewFile = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem_OpenFile = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem_SaveFile = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem_SaveAsFile = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripSeparator1 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripMenuItem5 = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripSeparator2 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripDropDownButton_View = new System.Windows.Forms.ToolStripDropDownButton();
            this.toolStripMenuItem_PlotOptions = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem_SimOptions = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem_DataOptions = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripSeparator3 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripDropDownButton_Tools = new System.Windows.Forms.ToolStripDropDownButton();
            this.toolStripMenuItem_DBSetup = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripSeparator4 = new System.Windows.Forms.ToolStripSeparator();
            this.toolStripDropDownButton_Help = new System.Windows.Forms.ToolStripDropDownButton();
            this.toolStripMenuItem_About = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripContainer_Back = new System.Windows.Forms.ToolStripContainer();
            this.plot_SplitContainer1 = new Quantization_Tool.Plot_SplitContainer();
            this.backgroundWorker_Simulation = new System.ComponentModel.BackgroundWorker();
            this.tableLayoutPanel_Plots.SuspendLayout();
            this.plot_ToolStripContainer_TopMiddle.SuspendLayout();
            this.topLeft_panel.SuspendLayout();
            this.plot_ToolStripContainer_TopLeft.SuspendLayout();
            this.plot_ToolStripContainer_TopRight.SuspendLayout();
            this.plot_ToolStripContainer_BottomLeft.SuspendLayout();
            this.plot_ToolStripContainer_BottomMiddle.SuspendLayout();
            this.plot_ToolStripContainer_BottomRight.SuspendLayout();
            this.tabPage1.SuspendLayout();
            this.tabControl1.SuspendLayout();
            this.tabPage2.SuspendLayout();
            this.panel_Plot.SuspendLayout();
            this.tableLayoutPanel_Back.SuspendLayout();
            this.tableLayoutPanel_bottom.SuspendLayout();
            this.flowLayoutPanel1.SuspendLayout();
            this.toolStrip_Main.SuspendLayout();
            this.toolStripContainer_Back.ContentPanel.SuspendLayout();
            this.toolStripContainer_Back.TopToolStripPanel.SuspendLayout();
            this.toolStripContainer_Back.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.plot_SplitContainer1)).BeginInit();
            this.plot_SplitContainer1.SuspendLayout();
            this.SuspendLayout();
            // 
            // button_Run
            // 
            this.button_Run.Anchor = System.Windows.Forms.AnchorStyles.None;
            this.button_Run.Location = new System.Drawing.Point(1624, 18);
            this.button_Run.Name = "button_Run";
            this.button_Run.Size = new System.Drawing.Size(126, 58);
            this.button_Run.TabIndex = 0;
            this.button_Run.Text = "Run";
            this.button_Run.UseVisualStyleBackColor = true;
            this.button_Run.Click += new System.EventHandler(this.button_Run_Click);
            // 
            // tableLayoutPanel_Plots
            // 
            this.tableLayoutPanel_Plots.CellBorderStyle = System.Windows.Forms.TableLayoutPanelCellBorderStyle.Single;
            this.tableLayoutPanel_Plots.ColumnCount = 3;
            this.tableLayoutPanel_Plots.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 33.33333F));
            this.tableLayoutPanel_Plots.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 33.33333F));
            this.tableLayoutPanel_Plots.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 33.33333F));
            this.tableLayoutPanel_Plots.Controls.Add(this.plot_ToolStripContainer_TopMiddle, 1, 0);
            this.tableLayoutPanel_Plots.Controls.Add(this.topLeft_panel, 0, 0);
            this.tableLayoutPanel_Plots.Controls.Add(this.plot_ToolStripContainer_TopRight, 2, 0);
            this.tableLayoutPanel_Plots.Controls.Add(this.plot_ToolStripContainer_BottomLeft, 0, 1);
            this.tableLayoutPanel_Plots.Controls.Add(this.plot_ToolStripContainer_BottomMiddle, 1, 1);
            this.tableLayoutPanel_Plots.Controls.Add(this.plot_ToolStripContainer_BottomRight, 2, 1);
            this.tableLayoutPanel_Plots.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tableLayoutPanel_Plots.Location = new System.Drawing.Point(3, 3);
            this.tableLayoutPanel_Plots.Name = "tableLayoutPanel_Plots";
            this.tableLayoutPanel_Plots.RowCount = 2;
            this.tableLayoutPanel_Plots.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 50F));
            this.tableLayoutPanel_Plots.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 50F));
            this.tableLayoutPanel_Plots.Size = new System.Drawing.Size(3529, 1551);
            this.tableLayoutPanel_Plots.TabIndex = 7;
            // 
            // plot_ToolStripContainer_TopMiddle
            // 
            this.plot_ToolStripContainer_TopMiddle.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_TopMiddle.ContentPanel
            // 
            this.plot_ToolStripContainer_TopMiddle.ContentPanel.Size = new System.Drawing.Size(1169, 699);
            this.plot_ToolStripContainer_TopMiddle.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_TopMiddle.Location = new System.Drawing.Point(1180, 4);
            this.plot_ToolStripContainer_TopMiddle.Name = "plot_ToolStripContainer_TopMiddle";
            this.plot_ToolStripContainer_TopMiddle.Size = new System.Drawing.Size(1169, 768);
            this.plot_ToolStripContainer_TopMiddle.TabIndex = 10;
            this.plot_ToolStripContainer_TopMiddle.Text = "plot_ToolStripContainer1";
            // 
            // topLeft_panel
            // 
            this.topLeft_panel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.topLeft_panel.Controls.Add(this.plot_ToolStripContainer_TopLeft);
            this.topLeft_panel.Dock = System.Windows.Forms.DockStyle.Fill;
            this.topLeft_panel.Location = new System.Drawing.Point(4, 4);
            this.topLeft_panel.Name = "topLeft_panel";
            this.topLeft_panel.Size = new System.Drawing.Size(1169, 768);
            this.topLeft_panel.TabIndex = 7;
            // 
            // plot_ToolStripContainer_TopLeft
            // 
            this.plot_ToolStripContainer_TopLeft.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_TopLeft.ContentPanel
            // 
            this.plot_ToolStripContainer_TopLeft.ContentPanel.Size = new System.Drawing.Size(1167, 697);
            this.plot_ToolStripContainer_TopLeft.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_TopLeft.Location = new System.Drawing.Point(0, 0);
            this.plot_ToolStripContainer_TopLeft.Name = "plot_ToolStripContainer_TopLeft";
            this.plot_ToolStripContainer_TopLeft.Size = new System.Drawing.Size(1167, 766);
            this.plot_ToolStripContainer_TopLeft.TabIndex = 0;
            this.plot_ToolStripContainer_TopLeft.Text = "plot_ToolStripContainer1";
            // 
            // plot_ToolStripContainer_TopRight
            // 
            this.plot_ToolStripContainer_TopRight.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_TopRight.ContentPanel
            // 
            this.plot_ToolStripContainer_TopRight.ContentPanel.Size = new System.Drawing.Size(1169, 699);
            this.plot_ToolStripContainer_TopRight.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_TopRight.Location = new System.Drawing.Point(2356, 4);
            this.plot_ToolStripContainer_TopRight.Name = "plot_ToolStripContainer_TopRight";
            this.plot_ToolStripContainer_TopRight.Size = new System.Drawing.Size(1169, 768);
            this.plot_ToolStripContainer_TopRight.TabIndex = 11;
            this.plot_ToolStripContainer_TopRight.Text = "plot_ToolStripContainer1";
            // 
            // plot_ToolStripContainer_BottomLeft
            // 
            this.plot_ToolStripContainer_BottomLeft.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_BottomLeft.ContentPanel
            // 
            this.plot_ToolStripContainer_BottomLeft.ContentPanel.Size = new System.Drawing.Size(1169, 699);
            this.plot_ToolStripContainer_BottomLeft.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_BottomLeft.Location = new System.Drawing.Point(4, 779);
            this.plot_ToolStripContainer_BottomLeft.Name = "plot_ToolStripContainer_BottomLeft";
            this.plot_ToolStripContainer_BottomLeft.Size = new System.Drawing.Size(1169, 768);
            this.plot_ToolStripContainer_BottomLeft.TabIndex = 12;
            this.plot_ToolStripContainer_BottomLeft.Text = "plot_ToolStripContainer1";
            // 
            // plot_ToolStripContainer_BottomMiddle
            // 
            this.plot_ToolStripContainer_BottomMiddle.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_BottomMiddle.ContentPanel
            // 
            this.plot_ToolStripContainer_BottomMiddle.ContentPanel.Size = new System.Drawing.Size(1169, 699);
            this.plot_ToolStripContainer_BottomMiddle.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_BottomMiddle.Location = new System.Drawing.Point(1180, 779);
            this.plot_ToolStripContainer_BottomMiddle.Name = "plot_ToolStripContainer_BottomMiddle";
            this.plot_ToolStripContainer_BottomMiddle.Size = new System.Drawing.Size(1169, 768);
            this.plot_ToolStripContainer_BottomMiddle.TabIndex = 13;
            this.plot_ToolStripContainer_BottomMiddle.Text = "plot_ToolStripContainer1";
            // 
            // plot_ToolStripContainer_BottomRight
            // 
            this.plot_ToolStripContainer_BottomRight.Chosen_Plot = Quantization_Tool.Plot_Choice.Points;
            // 
            // plot_ToolStripContainer_BottomRight.ContentPanel
            // 
            this.plot_ToolStripContainer_BottomRight.ContentPanel.Size = new System.Drawing.Size(1169, 699);
            this.plot_ToolStripContainer_BottomRight.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_ToolStripContainer_BottomRight.Location = new System.Drawing.Point(2356, 779);
            this.plot_ToolStripContainer_BottomRight.Name = "plot_ToolStripContainer_BottomRight";
            this.plot_ToolStripContainer_BottomRight.Size = new System.Drawing.Size(1169, 768);
            this.plot_ToolStripContainer_BottomRight.TabIndex = 14;
            this.plot_ToolStripContainer_BottomRight.Text = "plot_ToolStripContainer1";
            // 
            // tabPage1
            // 
            this.tabPage1.Controls.Add(this.tableLayoutPanel_Plots);
            this.tabPage1.Location = new System.Drawing.Point(10, 47);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage1.Size = new System.Drawing.Size(3535, 1557);
            this.tabPage1.TabIndex = 0;
            this.tabPage1.Text = "tabPage1";
            this.tabPage1.UseVisualStyleBackColor = true;
            // 
            // tabControl1
            // 
            this.tabControl1.Controls.Add(this.tabPage1);
            this.tabControl1.Controls.Add(this.tabPage2);
            this.tabControl1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tabControl1.Location = new System.Drawing.Point(3, 3);
            this.tabControl1.Name = "tabControl1";
            this.tabControl1.SelectedIndex = 0;
            this.tabControl1.Size = new System.Drawing.Size(3555, 1614);
            this.tabControl1.SizeMode = System.Windows.Forms.TabSizeMode.FillToRight;
            this.tabControl1.TabIndex = 8;
            // 
            // tabPage2
            // 
            this.tabPage2.Controls.Add(this.panel_Plot);
            this.tabPage2.Location = new System.Drawing.Point(10, 47);
            this.tabPage2.Name = "tabPage2";
            this.tabPage2.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage2.Size = new System.Drawing.Size(3535, 1557);
            this.tabPage2.TabIndex = 1;
            this.tabPage2.Text = "tabPage2";
            this.tabPage2.UseVisualStyleBackColor = true;
            // 
            // panel_Plot
            // 
            this.panel_Plot.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panel_Plot.Controls.Add(this.plot_plotView);
            this.panel_Plot.Dock = System.Windows.Forms.DockStyle.Fill;
            this.panel_Plot.Location = new System.Drawing.Point(3, 3);
            this.panel_Plot.Name = "panel_Plot";
            this.panel_Plot.Size = new System.Drawing.Size(3529, 1551);
            this.panel_Plot.TabIndex = 0;
            // 
            // plot_plotView
            // 
            this.plot_plotView.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_plotView.Location = new System.Drawing.Point(0, 0);
            this.plot_plotView.Name = "plot_plotView";
            this.plot_plotView.PanCursor = System.Windows.Forms.Cursors.Hand;
            this.plot_plotView.Size = new System.Drawing.Size(3527, 1549);
            this.plot_plotView.TabIndex = 0;
            this.plot_plotView.Text = "plotView1";
            this.plot_plotView.ZoomHorizontalCursor = System.Windows.Forms.Cursors.SizeWE;
            this.plot_plotView.ZoomRectangleCursor = System.Windows.Forms.Cursors.SizeNWSE;
            this.plot_plotView.ZoomVerticalCursor = System.Windows.Forms.Cursors.SizeNS;
            // 
            // tableLayoutPanel_Back
            // 
            this.tableLayoutPanel_Back.ColumnCount = 1;
            this.tableLayoutPanel_Back.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel_Back.Controls.Add(this.tabControl1, 0, 0);
            this.tableLayoutPanel_Back.Controls.Add(this.tableLayoutPanel_bottom, 0, 1);
            this.tableLayoutPanel_Back.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tableLayoutPanel_Back.Location = new System.Drawing.Point(0, 0);
            this.tableLayoutPanel_Back.Name = "tableLayoutPanel_Back";
            this.tableLayoutPanel_Back.RowCount = 2;
            this.tableLayoutPanel_Back.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel_Back.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Absolute, 100F));
            this.tableLayoutPanel_Back.Size = new System.Drawing.Size(3561, 1720);
            this.tableLayoutPanel_Back.TabIndex = 9;
            // 
            // tableLayoutPanel_bottom
            // 
            this.tableLayoutPanel_bottom.ColumnCount = 4;
            this.tableLayoutPanel_bottom.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 45F));
            this.tableLayoutPanel_bottom.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 5F));
            this.tableLayoutPanel_bottom.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 5F));
            this.tableLayoutPanel_bottom.ColumnStyles.Add(new System.Windows.Forms.ColumnStyle(System.Windows.Forms.SizeType.Percent, 45F));
            this.tableLayoutPanel_bottom.Controls.Add(this.flowLayoutPanel1, 3, 0);
            this.tableLayoutPanel_bottom.Controls.Add(this.simulation_ProgressBar, 0, 0);
            this.tableLayoutPanel_bottom.Controls.Add(this.button_Run, 1, 0);
            this.tableLayoutPanel_bottom.Controls.Add(this.button_Pause, 2, 0);
            this.tableLayoutPanel_bottom.Dock = System.Windows.Forms.DockStyle.Fill;
            this.tableLayoutPanel_bottom.Location = new System.Drawing.Point(3, 1623);
            this.tableLayoutPanel_bottom.Name = "tableLayoutPanel_bottom";
            this.tableLayoutPanel_bottom.RowCount = 1;
            this.tableLayoutPanel_bottom.RowStyles.Add(new System.Windows.Forms.RowStyle(System.Windows.Forms.SizeType.Percent, 100F));
            this.tableLayoutPanel_bottom.Size = new System.Drawing.Size(3555, 94);
            this.tableLayoutPanel_bottom.TabIndex = 9;
            // 
            // flowLayoutPanel1
            // 
            this.flowLayoutPanel1.Controls.Add(this.label_X);
            this.flowLayoutPanel1.Controls.Add(this.label_Y);
            this.flowLayoutPanel1.Location = new System.Drawing.Point(1956, 3);
            this.flowLayoutPanel1.Name = "flowLayoutPanel1";
            this.flowLayoutPanel1.Size = new System.Drawing.Size(485, 88);
            this.flowLayoutPanel1.TabIndex = 1;
            // 
            // label_X
            // 
            this.label_X.AutoSize = true;
            this.label_X.Location = new System.Drawing.Point(3, 0);
            this.label_X.Name = "label_X";
            this.label_X.Size = new System.Drawing.Size(48, 29);
            this.label_X.TabIndex = 0;
            this.label_X.Text = "X : ";
            // 
            // label_Y
            // 
            this.label_Y.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Bottom | System.Windows.Forms.AnchorStyles.Left)));
            this.label_Y.AutoSize = true;
            this.label_Y.Location = new System.Drawing.Point(57, 0);
            this.label_Y.Name = "label_Y";
            this.label_Y.Size = new System.Drawing.Size(59, 29);
            this.label_Y.TabIndex = 1;
            this.label_Y.Text = ", Y : ";
            // 
            // simulation_ProgressBar
            // 
            this.simulation_ProgressBar.Anchor = System.Windows.Forms.AnchorStyles.Right;
            this.simulation_ProgressBar.Location = new System.Drawing.Point(1148, 25);
            this.simulation_ProgressBar.Name = "simulation_ProgressBar";
            this.simulation_ProgressBar.Size = new System.Drawing.Size(448, 44);
            this.simulation_ProgressBar.TabIndex = 2;
            // 
            // button_Pause
            // 
            this.button_Pause.Anchor = System.Windows.Forms.AnchorStyles.None;
            this.button_Pause.Location = new System.Drawing.Point(1801, 18);
            this.button_Pause.Name = "button_Pause";
            this.button_Pause.Size = new System.Drawing.Size(126, 58);
            this.button_Pause.TabIndex = 3;
            this.button_Pause.Text = "Pause";
            this.button_Pause.UseVisualStyleBackColor = true;
            this.button_Pause.Click += new System.EventHandler(this.button_Pause_Click);
            // 
            // BottomToolStripPanel
            // 
            this.BottomToolStripPanel.Location = new System.Drawing.Point(0, 0);
            this.BottomToolStripPanel.Name = "BottomToolStripPanel";
            this.BottomToolStripPanel.Orientation = System.Windows.Forms.Orientation.Horizontal;
            this.BottomToolStripPanel.RowMargin = new System.Windows.Forms.Padding(3, 0, 0, 0);
            this.BottomToolStripPanel.Size = new System.Drawing.Size(0, 0);
            // 
            // TopToolStripPanel
            // 
            this.TopToolStripPanel.Location = new System.Drawing.Point(0, 0);
            this.TopToolStripPanel.Name = "TopToolStripPanel";
            this.TopToolStripPanel.Orientation = System.Windows.Forms.Orientation.Horizontal;
            this.TopToolStripPanel.RowMargin = new System.Windows.Forms.Padding(3, 0, 0, 0);
            this.TopToolStripPanel.Size = new System.Drawing.Size(0, 0);
            // 
            // RightToolStripPanel
            // 
            this.RightToolStripPanel.Location = new System.Drawing.Point(0, 0);
            this.RightToolStripPanel.Name = "RightToolStripPanel";
            this.RightToolStripPanel.Orientation = System.Windows.Forms.Orientation.Horizontal;
            this.RightToolStripPanel.RowMargin = new System.Windows.Forms.Padding(3, 0, 0, 0);
            this.RightToolStripPanel.Size = new System.Drawing.Size(0, 0);
            // 
            // LeftToolStripPanel
            // 
            this.LeftToolStripPanel.Location = new System.Drawing.Point(0, 0);
            this.LeftToolStripPanel.MinimumSize = new System.Drawing.Size(200, 0);
            this.LeftToolStripPanel.Name = "LeftToolStripPanel";
            this.LeftToolStripPanel.Orientation = System.Windows.Forms.Orientation.Horizontal;
            this.LeftToolStripPanel.RowMargin = new System.Windows.Forms.Padding(3, 0, 0, 0);
            this.LeftToolStripPanel.Size = new System.Drawing.Size(200, 0);
            // 
            // ContentPanel
            // 
            this.ContentPanel.Size = new System.Drawing.Size(1169, 765);
            // 
            // toolStrip_Main
            // 
            this.toolStrip_Main.Dock = System.Windows.Forms.DockStyle.None;
            this.toolStrip_Main.ImageScalingSize = new System.Drawing.Size(36, 36);
            this.toolStrip_Main.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripDropDownButton_File,
            this.toolStripSeparator2,
            this.toolStripDropDownButton_View,
            this.toolStripSeparator3,
            this.toolStripDropDownButton_Tools,
            this.toolStripSeparator4,
            this.toolStripDropDownButton_Help});
            this.toolStrip_Main.Location = new System.Drawing.Point(3, 0);
            this.toolStrip_Main.Name = "toolStrip_Main";
            this.toolStrip_Main.Size = new System.Drawing.Size(406, 44);
            this.toolStrip_Main.TabIndex = 1;
            // 
            // toolStripDropDownButton_File
            // 
            this.toolStripDropDownButton_File.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.toolStripDropDownButton_File.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripMenuItem_NewFile,
            this.toolStripMenuItem_OpenFile,
            this.toolStripMenuItem_SaveFile,
            this.toolStripMenuItem_SaveAsFile,
            this.toolStripSeparator1,
            this.toolStripMenuItem5});
            this.toolStripDropDownButton_File.Image = ((System.Drawing.Image)(resources.GetObject("toolStripDropDownButton_File.Image")));
            this.toolStripDropDownButton_File.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripDropDownButton_File.Name = "toolStripDropDownButton_File";
            this.toolStripDropDownButton_File.Size = new System.Drawing.Size(81, 41);
            this.toolStripDropDownButton_File.Text = "File";
            // 
            // toolStripMenuItem_NewFile
            // 
            this.toolStripMenuItem_NewFile.Name = "toolStripMenuItem_NewFile";
            this.toolStripMenuItem_NewFile.Size = new System.Drawing.Size(362, 42);
            this.toolStripMenuItem_NewFile.Text = "New";
            // 
            // toolStripMenuItem_OpenFile
            // 
            this.toolStripMenuItem_OpenFile.Name = "toolStripMenuItem_OpenFile";
            this.toolStripMenuItem_OpenFile.Size = new System.Drawing.Size(362, 42);
            this.toolStripMenuItem_OpenFile.Text = "Open";
            // 
            // toolStripMenuItem_SaveFile
            // 
            this.toolStripMenuItem_SaveFile.Name = "toolStripMenuItem_SaveFile";
            this.toolStripMenuItem_SaveFile.Size = new System.Drawing.Size(362, 42);
            this.toolStripMenuItem_SaveFile.Text = "Save";
            // 
            // toolStripMenuItem_SaveAsFile
            // 
            this.toolStripMenuItem_SaveAsFile.Name = "toolStripMenuItem_SaveAsFile";
            this.toolStripMenuItem_SaveAsFile.Size = new System.Drawing.Size(362, 42);
            this.toolStripMenuItem_SaveAsFile.Text = "Save As";
            // 
            // toolStripSeparator1
            // 
            this.toolStripSeparator1.Name = "toolStripSeparator1";
            this.toolStripSeparator1.Size = new System.Drawing.Size(359, 6);
            // 
            // toolStripMenuItem5
            // 
            this.toolStripMenuItem5.Name = "toolStripMenuItem5";
            this.toolStripMenuItem5.Size = new System.Drawing.Size(362, 42);
            this.toolStripMenuItem5.Text = "toolStripMenuItem5";
            // 
            // toolStripSeparator2
            // 
            this.toolStripSeparator2.Name = "toolStripSeparator2";
            this.toolStripSeparator2.Size = new System.Drawing.Size(6, 44);
            // 
            // toolStripDropDownButton_View
            // 
            this.toolStripDropDownButton_View.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.toolStripDropDownButton_View.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripMenuItem_PlotOptions,
            this.toolStripMenuItem_SimOptions,
            this.toolStripMenuItem_DataOptions});
            this.toolStripDropDownButton_View.Image = ((System.Drawing.Image)(resources.GetObject("toolStripDropDownButton_View.Image")));
            this.toolStripDropDownButton_View.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripDropDownButton_View.Name = "toolStripDropDownButton_View";
            this.toolStripDropDownButton_View.Size = new System.Drawing.Size(98, 41);
            this.toolStripDropDownButton_View.Text = "View";
            // 
            // toolStripMenuItem_PlotOptions
            // 
            this.toolStripMenuItem_PlotOptions.Checked = true;
            this.toolStripMenuItem_PlotOptions.CheckOnClick = true;
            this.toolStripMenuItem_PlotOptions.CheckState = System.Windows.Forms.CheckState.Checked;
            this.toolStripMenuItem_PlotOptions.Name = "toolStripMenuItem_PlotOptions";
            this.toolStripMenuItem_PlotOptions.Size = new System.Drawing.Size(352, 42);
            this.toolStripMenuItem_PlotOptions.Text = "Plot Options";
            this.toolStripMenuItem_PlotOptions.Click += new System.EventHandler(this.toolStripMenuItem_PlotOptions_Click);
            // 
            // toolStripMenuItem_SimOptions
            // 
            this.toolStripMenuItem_SimOptions.Checked = true;
            this.toolStripMenuItem_SimOptions.CheckOnClick = true;
            this.toolStripMenuItem_SimOptions.CheckState = System.Windows.Forms.CheckState.Checked;
            this.toolStripMenuItem_SimOptions.Name = "toolStripMenuItem_SimOptions";
            this.toolStripMenuItem_SimOptions.Size = new System.Drawing.Size(352, 42);
            this.toolStripMenuItem_SimOptions.Text = "Simulation Options";
            this.toolStripMenuItem_SimOptions.Click += new System.EventHandler(this.toolStripMenuItem_SimOptions_Click);
            // 
            // toolStripMenuItem_DataOptions
            // 
            this.toolStripMenuItem_DataOptions.Name = "toolStripMenuItem_DataOptions";
            this.toolStripMenuItem_DataOptions.Size = new System.Drawing.Size(352, 42);
            this.toolStripMenuItem_DataOptions.Text = "Data Options";
            // 
            // toolStripSeparator3
            // 
            this.toolStripSeparator3.Name = "toolStripSeparator3";
            this.toolStripSeparator3.Size = new System.Drawing.Size(6, 44);
            // 
            // toolStripDropDownButton_Tools
            // 
            this.toolStripDropDownButton_Tools.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.toolStripDropDownButton_Tools.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripMenuItem_DBSetup});
            this.toolStripDropDownButton_Tools.Image = ((System.Drawing.Image)(resources.GetObject("toolStripDropDownButton_Tools.Image")));
            this.toolStripDropDownButton_Tools.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripDropDownButton_Tools.Name = "toolStripDropDownButton_Tools";
            this.toolStripDropDownButton_Tools.Size = new System.Drawing.Size(101, 41);
            this.toolStripDropDownButton_Tools.Text = "Tools";
            // 
            // toolStripMenuItem_DBSetup
            // 
            this.toolStripMenuItem_DBSetup.Name = "toolStripMenuItem_DBSetup";
            this.toolStripMenuItem_DBSetup.Size = new System.Drawing.Size(310, 42);
            this.toolStripMenuItem_DBSetup.Text = "Database Setup";
            // 
            // toolStripSeparator4
            // 
            this.toolStripSeparator4.Name = "toolStripSeparator4";
            this.toolStripSeparator4.Size = new System.Drawing.Size(6, 44);
            // 
            // toolStripDropDownButton_Help
            // 
            this.toolStripDropDownButton_Help.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.toolStripDropDownButton_Help.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripMenuItem_About});
            this.toolStripDropDownButton_Help.Image = ((System.Drawing.Image)(resources.GetObject("toolStripDropDownButton_Help.Image")));
            this.toolStripDropDownButton_Help.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripDropDownButton_Help.Name = "toolStripDropDownButton_Help";
            this.toolStripDropDownButton_Help.Size = new System.Drawing.Size(96, 41);
            this.toolStripDropDownButton_Help.Text = "Help";
            // 
            // toolStripMenuItem_About
            // 
            this.toolStripMenuItem_About.Name = "toolStripMenuItem_About";
            this.toolStripMenuItem_About.Size = new System.Drawing.Size(197, 42);
            this.toolStripMenuItem_About.Text = "About";
            // 
            // toolStripContainer_Back
            // 
            this.toolStripContainer_Back.BottomToolStripPanelVisible = false;
            // 
            // toolStripContainer_Back.ContentPanel
            // 
            this.toolStripContainer_Back.ContentPanel.Controls.Add(this.tableLayoutPanel_Back);
            this.toolStripContainer_Back.ContentPanel.Controls.Add(this.plot_SplitContainer1);
            this.toolStripContainer_Back.ContentPanel.Size = new System.Drawing.Size(3561, 1720);
            this.toolStripContainer_Back.Dock = System.Windows.Forms.DockStyle.Fill;
            this.toolStripContainer_Back.LeftToolStripPanelVisible = false;
            this.toolStripContainer_Back.Location = new System.Drawing.Point(0, 0);
            this.toolStripContainer_Back.Name = "toolStripContainer_Back";
            this.toolStripContainer_Back.RightToolStripPanelVisible = false;
            this.toolStripContainer_Back.Size = new System.Drawing.Size(3561, 1764);
            this.toolStripContainer_Back.TabIndex = 0;
            this.toolStripContainer_Back.Text = "toolStripContainer1";
            // 
            // toolStripContainer_Back.TopToolStripPanel
            // 
            this.toolStripContainer_Back.TopToolStripPanel.Controls.Add(this.toolStrip_Main);
            // 
            // plot_SplitContainer1
            // 
            this.plot_SplitContainer1.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.plot_SplitContainer1.Dock = System.Windows.Forms.DockStyle.Fill;
            this.plot_SplitContainer1.FixedPanel = System.Windows.Forms.FixedPanel.Panel1;
            this.plot_SplitContainer1.IsSplitterFixed = true;
            this.plot_SplitContainer1.Location = new System.Drawing.Point(0, 0);
            this.plot_SplitContainer1.Name = "plot_SplitContainer1";
            this.plot_SplitContainer1.Panel1Collapsed = true;
            this.plot_SplitContainer1.Panel1MinSize = 100;
            this.plot_SplitContainer1.Size = new System.Drawing.Size(3561, 1720);
            this.plot_SplitContainer1.SplitterDistance = 100;
            this.plot_SplitContainer1.TabIndex = 0;
            // 
            // backgroundWorker_Simulation
            // 
            this.backgroundWorker_Simulation.WorkerReportsProgress = true;
            this.backgroundWorker_Simulation.WorkerSupportsCancellation = true;
            this.backgroundWorker_Simulation.DoWork += new System.ComponentModel.DoWorkEventHandler(this.backgroundWorker_Simulation_DoWork);
            this.backgroundWorker_Simulation.ProgressChanged += new System.ComponentModel.ProgressChangedEventHandler(this.backgroundWorker_Simulation_ProgressChanged);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(14F, 29F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(3561, 1764);
            this.Controls.Add(this.toolStripContainer_Back);
            this.Name = "Form1";
            this.Text = "Form1";
            this.tableLayoutPanel_Plots.ResumeLayout(false);
            this.plot_ToolStripContainer_TopMiddle.ResumeLayout(false);
            this.plot_ToolStripContainer_TopMiddle.PerformLayout();
            this.topLeft_panel.ResumeLayout(false);
            this.plot_ToolStripContainer_TopLeft.ResumeLayout(false);
            this.plot_ToolStripContainer_TopLeft.PerformLayout();
            this.plot_ToolStripContainer_TopRight.ResumeLayout(false);
            this.plot_ToolStripContainer_TopRight.PerformLayout();
            this.plot_ToolStripContainer_BottomLeft.ResumeLayout(false);
            this.plot_ToolStripContainer_BottomLeft.PerformLayout();
            this.plot_ToolStripContainer_BottomMiddle.ResumeLayout(false);
            this.plot_ToolStripContainer_BottomMiddle.PerformLayout();
            this.plot_ToolStripContainer_BottomRight.ResumeLayout(false);
            this.plot_ToolStripContainer_BottomRight.PerformLayout();
            this.tabPage1.ResumeLayout(false);
            this.tabControl1.ResumeLayout(false);
            this.tabPage2.ResumeLayout(false);
            this.panel_Plot.ResumeLayout(false);
            this.tableLayoutPanel_Back.ResumeLayout(false);
            this.tableLayoutPanel_bottom.ResumeLayout(false);
            this.flowLayoutPanel1.ResumeLayout(false);
            this.flowLayoutPanel1.PerformLayout();
            this.toolStrip_Main.ResumeLayout(false);
            this.toolStrip_Main.PerformLayout();
            this.toolStripContainer_Back.ContentPanel.ResumeLayout(false);
            this.toolStripContainer_Back.TopToolStripPanel.ResumeLayout(false);
            this.toolStripContainer_Back.TopToolStripPanel.PerformLayout();
            this.toolStripContainer_Back.ResumeLayout(false);
            this.toolStripContainer_Back.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.plot_SplitContainer1)).EndInit();
            this.plot_SplitContainer1.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button button_Run;
        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel_Plots;
        private System.Windows.Forms.TabControl tabControl1;
        private System.Windows.Forms.TabPage tabPage1;
        private System.Windows.Forms.TabPage tabPage2;
        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel_Back;
        private System.Windows.Forms.Panel topLeft_panel;
        private System.Windows.Forms.TableLayoutPanel tableLayoutPanel_bottom;
        private System.Windows.Forms.Panel panel_Plot;
        private OxyPlot.WindowsForms.PlotView plot_plotView;
        private System.Windows.Forms.FlowLayoutPanel flowLayoutPanel1;
        private System.Windows.Forms.Label label_X;
        private System.Windows.Forms.Label label_Y;
        private System.Windows.Forms.ToolStrip toolStrip_Main;
        private System.Windows.Forms.ToolStripDropDownButton toolStripDropDownButton_File;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_NewFile;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_OpenFile;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_SaveFile;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_SaveAsFile;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator1;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem5;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator2;
        private System.Windows.Forms.ToolStripDropDownButton toolStripDropDownButton_View;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator3;
        private System.Windows.Forms.ToolStripDropDownButton toolStripDropDownButton_Tools;
        private System.Windows.Forms.ToolStripSeparator toolStripSeparator4;
        private System.Windows.Forms.ToolStripDropDownButton toolStripDropDownButton_Help;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_PlotOptions;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_SimOptions;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_DataOptions;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_DBSetup;
        private System.Windows.Forms.ToolStripMenuItem toolStripMenuItem_About;
        private System.Windows.Forms.ToolStripContainer toolStripContainer_Back;
        private Plot_SplitContainer plot_SplitContainer1;
        private Plot_ToolStripContainer plot_ToolStripContainer_TopMiddle;
        private Plot_ToolStripContainer plot_ToolStripContainer_TopLeft;
        private Plot_ToolStripContainer plot_ToolStripContainer_TopRight;
        private Plot_ToolStripContainer plot_ToolStripContainer_BottomLeft;
        private Plot_ToolStripContainer plot_ToolStripContainer_BottomMiddle;
        private Plot_ToolStripContainer plot_ToolStripContainer_BottomRight;
        private System.Windows.Forms.ToolStripPanel BottomToolStripPanel;
        private System.Windows.Forms.ToolStripPanel TopToolStripPanel;
        private System.Windows.Forms.ToolStripPanel RightToolStripPanel;
        private System.Windows.Forms.ToolStripPanel LeftToolStripPanel;
        private System.Windows.Forms.ToolStripContentPanel ContentPanel;
        private System.Windows.Forms.ProgressBar simulation_ProgressBar;
        private System.ComponentModel.BackgroundWorker backgroundWorker_Simulation;
        private System.Windows.Forms.Button button_Pause;
    }
}

