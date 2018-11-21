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
            this.button_Run = new System.Windows.Forms.Button();
            this.plotView1 = new OxyPlot.WindowsForms.PlotView();
            this.plotView2 = new OxyPlot.WindowsForms.PlotView();
            this.left_Plotpanel = new System.Windows.Forms.Panel();
            this.right_Plotpanel = new System.Windows.Forms.Panel();
            this.leftbottom_Plotpanel = new System.Windows.Forms.Panel();
            this.plotView3 = new OxyPlot.WindowsForms.PlotView();
            this.rightbottom_Plotpanel = new System.Windows.Forms.Panel();
            this.plotView4 = new OxyPlot.WindowsForms.PlotView();
            this.left_Plotpanel.SuspendLayout();
            this.right_Plotpanel.SuspendLayout();
            this.leftbottom_Plotpanel.SuspendLayout();
            this.rightbottom_Plotpanel.SuspendLayout();
            this.SuspendLayout();
            // 
            // button_Run
            // 
            this.button_Run.Location = new System.Drawing.Point(2360, 831);
            this.button_Run.Name = "button_Run";
            this.button_Run.Size = new System.Drawing.Size(153, 81);
            this.button_Run.TabIndex = 0;
            this.button_Run.Text = "Run";
            this.button_Run.UseVisualStyleBackColor = true;
            this.button_Run.Click += new System.EventHandler(this.button_Run_Click);
            // 
            // plotView1
            // 
            this.plotView1.Location = new System.Drawing.Point(3, 1);
            this.plotView1.Name = "plotView1";
            this.plotView1.PanCursor = System.Windows.Forms.Cursors.Hand;
            this.plotView1.Size = new System.Drawing.Size(1047, 791);
            this.plotView1.TabIndex = 1;
            this.plotView1.Text = "plotView1";
            this.plotView1.ZoomHorizontalCursor = System.Windows.Forms.Cursors.SizeWE;
            this.plotView1.ZoomRectangleCursor = System.Windows.Forms.Cursors.SizeNWSE;
            this.plotView1.ZoomVerticalCursor = System.Windows.Forms.Cursors.SizeNS;
            // 
            // plotView2
            // 
            this.plotView2.Location = new System.Drawing.Point(3, 4);
            this.plotView2.Name = "plotView2";
            this.plotView2.PanCursor = System.Windows.Forms.Cursors.Hand;
            this.plotView2.Size = new System.Drawing.Size(1047, 791);
            this.plotView2.TabIndex = 2;
            this.plotView2.Text = "plotView2";
            this.plotView2.ZoomHorizontalCursor = System.Windows.Forms.Cursors.SizeWE;
            this.plotView2.ZoomRectangleCursor = System.Windows.Forms.Cursors.SizeNWSE;
            this.plotView2.ZoomVerticalCursor = System.Windows.Forms.Cursors.SizeNS;
            // 
            // left_Plotpanel
            // 
            this.left_Plotpanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.left_Plotpanel.Controls.Add(this.plotView3);
            this.left_Plotpanel.Location = new System.Drawing.Point(43, 42);
            this.left_Plotpanel.Name = "left_Plotpanel";
            this.left_Plotpanel.Size = new System.Drawing.Size(1055, 800);
            this.left_Plotpanel.TabIndex = 3;
            // 
            // right_Plotpanel
            // 
            this.right_Plotpanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.right_Plotpanel.Controls.Add(this.plotView1);
            this.right_Plotpanel.Location = new System.Drawing.Point(1162, 45);
            this.right_Plotpanel.Name = "right_Plotpanel";
            this.right_Plotpanel.Size = new System.Drawing.Size(1055, 800);
            this.right_Plotpanel.TabIndex = 4;
            // 
            // leftbottom_Plotpanel
            // 
            this.leftbottom_Plotpanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.leftbottom_Plotpanel.Controls.Add(this.plotView2);
            this.leftbottom_Plotpanel.Location = new System.Drawing.Point(43, 866);
            this.leftbottom_Plotpanel.Name = "leftbottom_Plotpanel";
            this.leftbottom_Plotpanel.Size = new System.Drawing.Size(1055, 800);
            this.leftbottom_Plotpanel.TabIndex = 5;
            // 
            // plotView3
            // 
            this.plotView3.Location = new System.Drawing.Point(-1, 2);
            this.plotView3.Name = "plotView3";
            this.plotView3.PanCursor = System.Windows.Forms.Cursors.Hand;
            this.plotView3.Size = new System.Drawing.Size(1047, 791);
            this.plotView3.TabIndex = 0;
            this.plotView3.Text = "plotView3";
            this.plotView3.ZoomHorizontalCursor = System.Windows.Forms.Cursors.SizeWE;
            this.plotView3.ZoomRectangleCursor = System.Windows.Forms.Cursors.SizeNWSE;
            this.plotView3.ZoomVerticalCursor = System.Windows.Forms.Cursors.SizeNS;
            // 
            // rightbottom_Plotpanel
            // 
            this.rightbottom_Plotpanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.rightbottom_Plotpanel.Controls.Add(this.plotView4);
            this.rightbottom_Plotpanel.Location = new System.Drawing.Point(1162, 866);
            this.rightbottom_Plotpanel.Name = "rightbottom_Plotpanel";
            this.rightbottom_Plotpanel.Size = new System.Drawing.Size(1055, 800);
            this.rightbottom_Plotpanel.TabIndex = 6;
            // 
            // plotView4
            // 
            this.plotView4.Location = new System.Drawing.Point(3, 3);
            this.plotView4.Name = "plotView4";
            this.plotView4.PanCursor = System.Windows.Forms.Cursors.Hand;
            this.plotView4.Size = new System.Drawing.Size(1047, 791);
            this.plotView4.TabIndex = 0;
            this.plotView4.Text = "plotView4";
            this.plotView4.ZoomHorizontalCursor = System.Windows.Forms.Cursors.SizeWE;
            this.plotView4.ZoomRectangleCursor = System.Windows.Forms.Cursors.SizeNWSE;
            this.plotView4.ZoomVerticalCursor = System.Windows.Forms.Cursors.SizeNS;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(14F, 29F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(2916, 1723);
            this.Controls.Add(this.rightbottom_Plotpanel);
            this.Controls.Add(this.leftbottom_Plotpanel);
            this.Controls.Add(this.right_Plotpanel);
            this.Controls.Add(this.left_Plotpanel);
            this.Controls.Add(this.button_Run);
            this.Name = "Form1";
            this.Text = "Form1";
            this.left_Plotpanel.ResumeLayout(false);
            this.right_Plotpanel.ResumeLayout(false);
            this.leftbottom_Plotpanel.ResumeLayout(false);
            this.rightbottom_Plotpanel.ResumeLayout(false);
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button button_Run;
        private OxyPlot.WindowsForms.PlotView plotView1;
        private OxyPlot.WindowsForms.PlotView plotView2;
        private System.Windows.Forms.Panel left_Plotpanel;
        private System.Windows.Forms.Panel right_Plotpanel;
        private System.Windows.Forms.Panel leftbottom_Plotpanel;
        private OxyPlot.WindowsForms.PlotView plotView3;
        private System.Windows.Forms.Panel rightbottom_Plotpanel;
        private OxyPlot.WindowsForms.PlotView plotView4;
    }
}

