using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO;
using System.Web;
using System.Threading;
using System.Collections;
using MSFileReaderLib;


namespace DIADataExtracter
{
    public partial class DIAExtracter : Form
    {
        string[] filedirectory_array;
        string filedirectory;
        string extracttype;
        string ms1rms_checkbox;
        string ms2rms_checkbox;
        string deisotope_checkbox;
        string deconvolution_checkbox;
        public DIAExtracter()
        {
            InitializeComponent();
        }

        private void helpHToolStripMenuItem_Click(object sender, EventArgs e)
        {
            MessageBox.Show("RCPA_DIADataExtracter version 1.0.0\n"
                            + "E-mail: wssdandan2009@outlook.com\n"
                            + "Thank you for using RCPA_DIADataExtracter.");
        }

        private void selectfile_Click(object sender, EventArgs e)
        {
            //Select a DIA raw file
            DialogResult pathdialog = path.ShowDialog();
            if (pathdialog == DialogResult.OK)
            {
                //Add all files you select to the ListBox
                filedirectory_array = path.FileNames;
                foreach (string filename in this.filedirectory_array)
                {
                    this.pathway.Items.Add(filename);
                }
                
            }
        }

        private void MS1RMS_CheckedChanged(object sender, EventArgs e)
        {
            //RMS: Root Mean Square
            if (MS1RMS.Checked == true) ms1rms_checkbox = "ms1rms";
            else ms1rms_checkbox = "";
        }

        private void MS2RMS_CheckedChanged(object sender, EventArgs e)
        {
            if (MS2RMS.Checked == true) ms2rms_checkbox = "ms2rms";
            else ms2rms_checkbox = "";
        }

        private void deisotope_CheckedChanged(object sender, EventArgs e)
        {
            if (deisotope.Checked == true) deisotope_checkbox = "deisotope";
            else deisotope_checkbox = "";
        }

        private void deconvolution_CheckedChanged(object sender, EventArgs e)
        {
            if (deconvolution.Checked == true) deconvolution_checkbox = "deconvolution";
            else deconvolution_checkbox = "";
        }

        private void comboBox_Etype_SelectedIndexChanged(object sender, EventArgs e)
        {
            /* Chose the extract type:
               1. MS1MS2Data: obtain the ms1 and ms2 data with 8 columns information;
                 (ScanNum-scan number;
                  MSOrder-the scan serial number;
                  PeakNumber-the total number of peaks in one scan;
                  PeakNumberAfter-the number of peaks after pre-process in one scan;
                  RT-rentention time;
                  m/z-mass over charge;
                  Intensity-the intensity of peak,
                  Noise;
                  Baseline)
               2. LocalData(SMPR): obtain the ms1 and ms2 data with the first 4 columns information;
                  (ScanNum, MSOrder, PeakNumber, RT)*/
            switch (comboBox_Etype.SelectedIndex)
            {
                case 0: extracttype = "MS1MS2Data"; break;
                case 1: extracttype = "LocalData(SMPR)"; break;
            }

            if (extracttype == "LocalData(SMPR)")
            {
                ms12filter.Visible = false; //hide the filter window
            }
            else
            {
                ms12filter.Visible = true;
            }
        }

        private void DIAExtracter_Load_1(object sender, EventArgs e)
        {
            comboBox_Etype.Text = "MS1MS2Data";
        }

        //Show the progress that how many raw files left;
        private delegate void SetPos(int ipos_rawfile_num);
        private void SetTextMessage_rawfile_num(int ipos_rawfile_num)
        {
            if (this.InvokeRequired)
            {
                SetPos setpos = new SetPos(SetTextMessage_rawfile_num);
                this.Invoke(setpos, new object[] { ipos_rawfile_num });
            }
            else
            {
                this.rawfilenum.Text = "A:" + ipos_rawfile_num.ToString() + "/" + pathway.Items.Count.ToString();
            }
        }

        //Show the progress of dealing with one raw file;
        //private delegate void SetPos(int ipos);
        private void SetTextMessage(int ipos)
        {
            if (this.InvokeRequired)
            {
                SetPos setpos = new SetPos(SetTextMessage);
                this.Invoke(setpos, new object[] { ipos });
            }
            else
            {
                this.processlabel.Text ="B:"+ ipos.ToString() + "%/100%";
                this.progressBar1.Value = Convert.ToInt32(ipos);
            }
        }

        private void extracter_Click(object sender, EventArgs e)
        {
            //Open a new thread;
            Thread fThread = new Thread(new ThreadStart(SleepT));
            fThread.Start();
            Progress.Visible = true;
        }

        private void SleepT()
        {
            if (filedirectory_array.Length != 0)
            {
          #region
                for (int i_filedirectory = 0; i_filedirectory < pathway.Items.Count; i_filedirectory++)
                {
                    SetTextMessage_rawfile_num(i_filedirectory+1);

                    //Utilize the MSFileReader to obtain the information from the DIA data;
                    filedirectory = (string)pathway.Items[i_filedirectory];
                    MSFileReader_XRawfile rawfile = new MSFileReader_XRawfile();
                    IXRawfile5 rawfile5 = (IXRawfile5)new MSFileReader_XRawfile();
                    string rawpath = filedirectory;
                    rawfile.Open(rawpath);
                    rawfile.SetCurrentController(0, 1);
                    string rawpath5 = filedirectory;
                    rawfile5.Open(rawpath5);
                    rawfile5.SetCurrentController(0, 1);

                    string filename = Path.GetFileName(filedirectory);
                    string newfilepath = Path.GetDirectoryName(filedirectory);
                    #region
                    if (extracttype == "MS1MS2Data")
                    {
                        //creat the output file names
                        string outfilems1 = newfilepath + "\\" + "MS1Data_" + filename + ".txt";
                        StreamWriter swms1 = File.AppendText(outfilems1);
                        string outfilems2 = newfilepath + "\\" + "MS2Data_" + filename + ".txt";
                        StreamWriter swms2 = File.AppendText(outfilems2);

                        //write the 9 columns information in the .txt file
                        swms1.Write("ScanNum\tMSOrder\tPeakNumber\tPeakNumberAfter\tRT\tm/z\tIntensity\tNoise\tBaseline\n");
                        swms2.Write("ScanNum\tMSOrder\tPeakNumber\tPeakNumberAfter\tRT\tm/z\tIntensity\tNoise\tBaseline\n");
                        
                        int specnum = 0;
                        rawfile5.GetNumSpectra(ref specnum);
                        for (int i = 1; i <= specnum; i++)
                        {
                            System.Threading.Thread.Sleep(0);
                            SetTextMessage(100 * i / specnum);

                            int nOrder = 0;
                            rawfile5.GetMSOrderForScanNum(i, ref nOrder);
                            double scanRT = 0;
                            rawfile5.RTFromScanNum(i, ref scanRT);
                            object varlabels = null;
                            object varflags = null;
                            rawfile5.GetLabelData(ref varlabels, ref varflags, ref i);
                            var labels = (double[,])varlabels;
                            int dim = labels.GetLength(1);

                            double massNeutron = 1.00335; // massC13 - massC12
                            double mztol = Convert.ToDouble(ionmtoltextbox.Text);
                            double mzdecontol = Convert.ToDouble(ionmtoltextbox.Text);
                            double ms1sonthreshold = Convert.ToDouble(ms1sontextbox.Text);
                            double ms2sonthreshold = Convert.ToDouble(ms2sontextbox.Text);
                            double ms1resolutiontol = Convert.ToDouble(ms1unitoltextbox.Text);
                            double ms2resolutiontol = Convert.ToDouble(ms2unitoltextbox.Text);

                            if (nOrder == 1)
                            {
                                swms1.Write("{0}\t", i);
                                swms1.Write("{0}\t", nOrder);
                                swms1.Write("{0}\t", dim);

                                List<double> ms1mzlist = new List<double>();
                                List<double> ms1intensitylist = new List<double>();
                                List<double> ms1chargelist = new List<double>();
                                List<double> ms1noiselist = new List<double>();
                                List<double> ms1baselinelist = new List<double>();

                                for (int inx = 0; inx < dim; inx++)
                                {
                                    double dMassjms1 = labels[0, inx];
                                    double dInms1 = labels[1, inx];
                                    double dBasems1 = labels[3, inx];
                                    double dNoisems1 = labels[4, inx];
                                    double dChams1 = labels[5, inx];
                                    double sonms1 = 0;
                                    if (ms1rms_checkbox == "ms1rms")
                                    {
                                        //if checked, the S/N is changed to RMS S/N(root mean square signal over noise);
                                        sonms1 = (dInms1 - dBasems1) / (dNoisems1 - dBasems1);
                                    }
                                    else
                                    {
                                        sonms1 = dInms1 / dNoisems1;
                                    }

                                    if (sonms1 >= ms1sonthreshold)
                                    {
                                        ms1chargelist.Add(dChams1);
                                        ms1mzlist.Add(dMassjms1);
                                        ms1intensitylist.Add(dInms1);
                                        ms1noiselist.Add(dNoisems1);
                                        ms1baselinelist.Add(dBasems1);
                                    }
                                    else
                                    {
                                        //if the S/N is too low, set their values to zero
                                        ms1chargelist.Add(0);
                                        ms1mzlist.Add(0);
                                        ms1intensitylist.Add(0);
                                        ms1noiselist.Add(0);
                                        ms1baselinelist.Add(0);
                                    }
                                }

                                //remove the peaks, whose values of (m/z, intensity, charge, noise, baseline) are equal to zero;
                                List<int> ms1mzzeroindex = new List<int>();
                                for (int izero = 0; izero < ms1mzlist.Count; izero++)
                                {
                                    if (ms1mzlist[izero] == 0)
                                    {
                                        ms1mzzeroindex.Add(izero);
                                    }
                                }
                                for (int izero1 = 0; izero1 < ms1mzzeroindex.Count(); izero1++)
                                {
                                    ms1mzlist.RemoveAt(ms1mzzeroindex[izero1] - izero1);
                                    ms1intensitylist.RemoveAt(ms1mzzeroindex[izero1] - izero1);
                                    ms1baselinelist.RemoveAt(ms1mzzeroindex[izero1] - izero1);
                                    ms1noiselist.RemoveAt(ms1mzzeroindex[izero1] - izero1);
                                    ms1chargelist.RemoveAt(ms1mzzeroindex[izero1] - izero1);
                                }

                                //keep only one peak in a tolerance;
                                for (int ims1reso = 0; ims1reso < ms1mzlist.Count; ims1reso++)
                                {
                                    if (ms1mzlist[ims1reso] != 0)
                                    {
                                        double ms1mzreso = ms1mzlist[ims1reso] * ms1intensitylist[ims1reso];
                                        double ms1intensityreso = ms1intensitylist[ims1reso];
                                        double ms1baselinereso = ms1baselinelist[ims1reso];
                                        for (int jms1reso = ims1reso + 1; jms1reso < ms1mzlist.Count; jms1reso++)
                                        {
                                            double ms1resodiff = ms1mzlist[jms1reso] - ms1mzlist[jms1reso - 1];
                                            if (ms1resodiff > ms1resolutiontol) break;
                                            if (ms1resodiff <= ms1resolutiontol)
                                            {
                                                ms1mzreso = ms1mzreso + ms1mzlist[jms1reso] * ms1intensitylist[jms1reso];
                                                ms1intensityreso = ms1intensityreso + ms1intensitylist[jms1reso];
                                                ms1intensitylist[jms1reso] = 0;
                                            }
                                        }
                                        ms1mzlist[ims1reso] = ms1mzreso / ms1intensityreso;
                                        ms1intensitylist[ims1reso] = ms1intensityreso;
                                    }
                                }

                                List<int> ms1intensityzeroafterresoindex = new List<int>();
                                for (int izero = 0; izero < ms1mzlist.Count; izero++)
                                {
                                    if (ms1intensitylist[izero] == 0)
                                    {
                                        ms1intensityzeroafterresoindex.Add(izero);
                                    }
                                }
                                for (int izero1 = 0; izero1 < ms1intensityzeroafterresoindex.Count(); izero1++)
                                {
                                    ms1mzlist.RemoveAt(ms1intensityzeroafterresoindex[izero1] - izero1);
                                    ms1intensitylist.RemoveAt(ms1intensityzeroafterresoindex[izero1] - izero1);
                                    ms1baselinelist.RemoveAt(ms1intensityzeroafterresoindex[izero1] - izero1);
                                    ms1noiselist.RemoveAt(ms1intensityzeroafterresoindex[izero1] - izero1);
                                    ms1chargelist.RemoveAt(ms1intensityzeroafterresoindex[izero1] - izero1);
                                }

                                //if checked, the peaks are deisotoped, which could remove the isotope peaks;
                                if (deisotope_checkbox == "deisotope")
                                {
                                    List<int> ms1deisotopindexbefore = new List<int>();

                                    for (int ims1 = 0; ims1 < ms1mzlist.Count; ims1++)
                                    {
                                        if (ms1chargelist[ims1] != 0)
                                        {
                                            for (int jms1 = ims1 + 1; jms1 < ms1mzlist.Count; jms1++)
                                            {
                                                double ms1mzdiff = ms1mzlist[jms1] - ms1mzlist[ims1];
                                                double ms1tolvalue = Math.Abs(massNeutron / ms1chargelist[ims1] - ms1mzdiff);
                                                if (ms1tolvalue < mztol && ms1intensitylist[jms1] < ms1intensitylist[ims1])
                                                {
                                                    ms1deisotopindexbefore.Add(jms1);
                                                }
                                            }
                                        }
                                    }

                                    List<int> ms1deisotopindex = new List<int>();
                                    foreach (int ms1deiso in ms1deisotopindexbefore)
                                    {
                                        if (!ms1deisotopindex.Contains(ms1deiso))
                                        {
                                            ms1deisotopindex.Add(ms1deiso);
                                        }
                                    }
                                    ms1deisotopindex.Sort();
                                    for (int ide1 = 0; ide1 < ms1deisotopindex.Count; ide1++)
                                    {
                                        ms1mzlist.RemoveAt(ms1deisotopindex[ide1] - ide1);
                                        ms1intensitylist.RemoveAt(ms1deisotopindex[ide1] - ide1);
                                        ms1baselinelist.RemoveAt(ms1deisotopindex[ide1] - ide1);
                                        ms1noiselist.RemoveAt(ms1deisotopindex[ide1] - ide1);
                                        ms1chargelist.RemoveAt(ms1deisotopindex[ide1] - ide1);
                                    }
                                }

                                //if checked, the peaks are deconvoluted to one charge;
                                if (deconvolution_checkbox == "deconvolution")
                                {
                                    for (int idems1 = 0; idems1 < ms1mzlist.Count; idems1++)
                                    {
                                        if (ms1chargelist[idems1] != 0)
                                        {
                                            ms1mzlist[idems1] = ms1mzlist[idems1] * ms1chargelist[idems1] - ms1chargelist[idems1] * 1.0079;
                                        }
                                    }
                                    List<int> ms1deconvolutionindexbefore = new List<int>();
                                    for (int ideconms1 = 0; ideconms1 < ms1mzlist.Count; ideconms1++)
                                    {
                                        if (ms1mzlist[ideconms1] != 0)
                                        {
                                            double ms1mzdecon = ms1mzlist[ideconms1];
                                            double ms1intensitydecon = ms1intensitylist[ideconms1];
                                            double ms1baselinedecon = ms1baselinelist[ideconms1];
                                            int ims1deconcount = 1;
                                            for (int jdeconms1 = ideconms1 + 1; jdeconms1 < ms1mzlist.Count; jdeconms1++)
                                            {
                                                double ms1mzdecondiff = Math.Abs(ms1mzlist[jdeconms1] - ms1mzlist[ideconms1]);
                                                if (ms1mzdecondiff <= mzdecontol)
                                                {
                                                    ms1deconvolutionindexbefore.Add(jdeconms1);
                                                    ms1mzdecon = ms1mzdecon + ms1mzlist[jdeconms1];
                                                    ms1intensitydecon = ms1intensitydecon + ms1intensitylist[jdeconms1];
                                                    ims1deconcount++;
                                                }
                                            }
                                            ms1mzlist[ideconms1] = ms1mzdecon / ims1deconcount;
                                            ms1intensitylist[ideconms1] = ms1intensitydecon;
                                        }
                                    }

                                    List<int> ms1deconvolutionindex = new List<int>();
                                    foreach (int ms1deiso in ms1deconvolutionindexbefore)
                                    {
                                        if (!ms1deconvolutionindex.Contains(ms1deiso))
                                        {
                                            ms1deconvolutionindex.Add(ms1deiso);
                                        }
                                    }
                                    ms1deconvolutionindex.Sort();
                                    for (int ide1 = 0; ide1 < ms1deconvolutionindex.Count; ide1++)
                                    {
                                        ms1mzlist.RemoveAt(ms1deconvolutionindex[ide1] - ide1);
                                        ms1intensitylist.RemoveAt(ms1deconvolutionindex[ide1] - ide1);
                                        ms1baselinelist.RemoveAt(ms1deconvolutionindex[ide1] - ide1);
                                        ms1noiselist.RemoveAt(ms1deconvolutionindex[ide1] - ide1);
                                        ms1chargelist.RemoveAt(ms1deconvolutionindex[ide1] - ide1);
                                    }
                                }

                                //write the data into the .txt file;
                                swms1.Write("{0}\t", ms1mzlist.Count);
                                swms1.Write("{0}\t", scanRT);
                                for (int ims1mz = 0; ims1mz < ms1mzlist.Count; ims1mz++)
                                {
                                    swms1.Write("{0};", ms1mzlist[ims1mz]);
                                }
                                swms1.Write("\t");
                                for (int ims1in = 0; ims1in < ms1intensitylist.Count; ims1in++)
                                {
                                    swms1.Write("{0};", ms1intensitylist[ims1in]);
                                }
                                swms1.Write("\t");
                                for (int ims1in = 0; ims1in < ms1noiselist.Count; ims1in++)
                                {
                                    swms1.Write("{0};", ms1noiselist[ims1in]);
                                }
                                swms1.Write("\t");
                                for (int ims1base = 0; ims1base < ms1baselinelist.Count; ims1base++)
                                {
                                    swms1.Write("{0};", ms1baselinelist[ims1base]);
                                }
                                swms1.Write("\n");
                            }
                            else //extract the ms/ms data, the process is similar with above;
                            {
                                swms2.Write("{0}\t", i);
                                swms2.Write("{0}\t", nOrder);
                                swms2.Write("{0}\t", dim);

                                List<double> ms2mzlist = new List<double>();
                                List<double> ms2intensitylist = new List<double>();
                                List<double> ms2noiselist = new List<double>();
                                List<double> ms2chargelist = new List<double>();
                                List<double> ms2baselinelist = new List<double>();


                                for (int inx = 0; inx < dim; inx++)
                                {
                                    double dMassjms2 = labels[0, inx];
                                    double dInms2 = labels[1, inx];
                                    double dBasems2 = labels[3, inx];
                                    double dNoisems2 = labels[4, inx];
                                    double dChams2 = labels[5, inx];
                                    double sonms2 = 0;
                                    if (ms2rms_checkbox == "ms2rms")
                                    {
                                        sonms2 = (dInms2 - dBasems2) / (dNoisems2 - dBasems2);
                                    }
                                    else
                                    {
                                        sonms2 = dInms2 / dNoisems2;
                                    }

                                    if (sonms2 >= ms2sonthreshold)
                                    {
                                        ms2chargelist.Add(dChams2);
                                        ms2mzlist.Add(dMassjms2);
                                        ms2intensitylist.Add(dInms2);
                                        ms2noiselist.Add(dNoisems2);
                                        ms2baselinelist.Add(dBasems2);
                                    }
                                    else
                                    {
                                        ms2chargelist.Add(0);
                                        ms2mzlist.Add(0);
                                        ms2intensitylist.Add(0);
                                        ms2noiselist.Add(0);
                                        ms2baselinelist.Add(0);
                                    }
                                }

                                List<int> ms2mzzeroindex = new List<int>();
                                for (int izero = 0; izero < ms2mzlist.Count; izero++)
                                {
                                    if (ms2mzlist[izero] == 0)
                                    {
                                        ms2mzzeroindex.Add(izero);
                                    }
                                }
                                for (int izero1 = 0; izero1 < ms2mzzeroindex.Count(); izero1++)
                                {
                                    ms2mzlist.RemoveAt(ms2mzzeroindex[izero1] - izero1);
                                    ms2intensitylist.RemoveAt(ms2mzzeroindex[izero1] - izero1);
                                    ms2baselinelist.RemoveAt(ms2mzzeroindex[izero1] - izero1);
                                    ms2noiselist.RemoveAt(ms2mzzeroindex[izero1] - izero1);
                                    ms2chargelist.RemoveAt(ms2mzzeroindex[izero1] - izero1);
                                }

                                for (int ims2reso = 0; ims2reso < ms2mzlist.Count; ims2reso++)
                                {
                                    if (ms2mzlist[ims2reso] != 0)
                                    {
                                        double ms2mzreso = ms2mzlist[ims2reso] * ms2intensitylist[ims2reso];
                                        double ms2intensityreso = ms2intensitylist[ims2reso];
                                        double ms2baselinereso = ms2baselinelist[ims2reso];
                                        int iresocount = 1;
                                        for (int jms2reso = ims2reso + 1; jms2reso < ms2mzlist.Count; jms2reso++)
                                        {
                                            double ms2resodiff = ms2mzlist[jms2reso] - ms2mzlist[jms2reso - 1];
                                            if (ms2resodiff > ms2resolutiontol) break;
                                            if (ms2resodiff <= ms2resolutiontol)
                                            {
                                                ms2mzreso = ms2mzreso + ms2mzlist[jms2reso] * ms2intensitylist[jms2reso];
                                                ms2intensityreso = ms2intensityreso + ms2intensitylist[jms2reso];
                                                ms2intensitylist[jms2reso] = 0;
                                                iresocount++;
                                            }
                                        }
                                        ms2mzlist[ims2reso] = ms2mzreso / ms2intensityreso;
                                        ms2intensitylist[ims2reso] = ms2intensityreso;
                                    }
                                }

                                List<int> ms2intensityzeroafterresoindex = new List<int>();
                                for (int izero = 0; izero < ms2mzlist.Count; izero++)
                                {
                                    if (ms2intensitylist[izero] == 0)
                                    {
                                        ms2intensityzeroafterresoindex.Add(izero);
                                    }
                                }
                                for (int izero1 = 0; izero1 < ms2intensityzeroafterresoindex.Count(); izero1++)
                                {
                                    ms2mzlist.RemoveAt(ms2intensityzeroafterresoindex[izero1] - izero1);
                                    ms2intensitylist.RemoveAt(ms2intensityzeroafterresoindex[izero1] - izero1);
                                    ms2baselinelist.RemoveAt(ms2intensityzeroafterresoindex[izero1] - izero1);
                                    ms2noiselist.RemoveAt(ms2intensityzeroafterresoindex[izero1] - izero1);
                                    ms2chargelist.RemoveAt(ms2intensityzeroafterresoindex[izero1] - izero1);
                                }

                                if (deisotope_checkbox == "deisotope")
                                {
                                    List<int> ms2deisotopindexbefore = new List<int>();
                                    for (int ims2 = 0; ims2 < ms2mzlist.Count(); ims2++)
                                    {
                                        if (ms2chargelist[ims2] != 0)
                                        {
                                            for (int jms2 = ims2 + 1; jms2 < ms2mzlist.Count(); jms2++)
                                            {
                                                double ms2mzdiff = ms2mzlist[jms2] - ms2mzlist[ims2];
                                                double ms2tolvalue = Math.Abs(massNeutron / ms2chargelist[ims2] - ms2mzdiff);
                                                if (ms2tolvalue <= mztol && ms2intensitylist[jms2] < ms2intensitylist[ims2])
                                                {
                                                    ms2deisotopindexbefore.Add(jms2);
                                                }
                                            }
                                        }
                                    }

                                    List<int> ms2deisotopindex = new List<int>();
                                    foreach (int ms2deiso in ms2deisotopindexbefore)
                                    {
                                        if (!ms2deisotopindex.Contains(ms2deiso))
                                        {
                                            ms2deisotopindex.Add(ms2deiso);
                                        }
                                    }
                                    ms2deisotopindex.Sort();
                                    for (int ide2 = 0; ide2 < ms2deisotopindex.Count(); ide2++)
                                    {
                                        ms2mzlist.RemoveAt(ms2deisotopindex[ide2] - ide2);
                                        ms2intensitylist.RemoveAt(ms2deisotopindex[ide2] - ide2);
                                        ms2noiselist.RemoveAt(ms2deisotopindex[ide2] - ide2);
                                        ms2baselinelist.RemoveAt(ms2deisotopindex[ide2] - ide2);
                                        ms2chargelist.RemoveAt(ms2deisotopindex[ide2] - ide2);
                                    }
                                }

                                if (deconvolution_checkbox == "deconvolution")
                                {
                                    for (int idems2 = 0; idems2 < ms2mzlist.Count; idems2++)
                                    {
                                        if (ms2chargelist[idems2] != 0)
                                        {
                                            ms2mzlist[idems2] = ms2mzlist[idems2] * ms2chargelist[idems2] - (ms2chargelist[idems2] - 1) * 1.0079;
                                        }
                                    }
                                    List<int> ms2deconvolutionindexbefore = new List<int>();
                                    for (int ideconms2 = 0; ideconms2 < ms2mzlist.Count; ideconms2++)
                                    {
                                        if (ms2mzlist[ideconms2] != 0)
                                        {
                                            double ms2mzdecon = ms2mzlist[ideconms2];
                                            double ms2intensitydecon = ms2intensitylist[ideconms2];
                                            int ims2deconcount = 1;
                                            for (int jdeconms2 = ideconms2 + 1; jdeconms2 < ms2mzlist.Count; jdeconms2++)
                                            {
                                                double ms2mzdecondiff = Math.Abs(ms2mzlist[jdeconms2] - ms2mzlist[ideconms2]);
                                                if (ms2mzdecondiff <= mzdecontol)
                                                {
                                                    ms2deconvolutionindexbefore.Add(jdeconms2);
                                                    ms2mzdecon = ms2mzdecon + ms2mzlist[jdeconms2];
                                                    ms2intensitydecon = ms2intensitydecon + ms2intensitylist[jdeconms2];
                                                    ims2deconcount++;
                                                }
                                            }
                                            ms2mzlist[ideconms2] = ms2mzdecon / ims2deconcount;
                                            ms2intensitylist[ideconms2] = ms2intensitydecon;
                                        }
                                    }
                                    List<int> ms2deconvolutionindex = new List<int>();
                                    foreach (int ms2deiso in ms2deconvolutionindexbefore)
                                    {
                                        if (!ms2deconvolutionindex.Contains(ms2deiso))
                                        {
                                            ms2deconvolutionindex.Add(ms2deiso);
                                        }
                                    }
                                    ms2deconvolutionindex.Sort();
                                    for (int ide2 = 0; ide2 < ms2deconvolutionindex.Count; ide2++)
                                    {
                                        ms2mzlist.RemoveAt(ms2deconvolutionindex[ide2] - ide2);
                                        ms2intensitylist.RemoveAt(ms2deconvolutionindex[ide2] - ide2);
                                        ms2baselinelist.RemoveAt(ms2deconvolutionindex[ide2] - ide2);
                                        ms2noiselist.RemoveAt(ms2deconvolutionindex[ide2] - ide2);
                                        ms2chargelist.RemoveAt(ms2deconvolutionindex[ide2] - ide2);
                                    }
                                }

                                swms2.Write("{0}\t", ms2mzlist.Count);
                                swms2.Write("{0}\t", scanRT);

                                for (int ims2mz = 0; ims2mz < ms2mzlist.Count; ims2mz++)
                                {
                                    swms2.Write("{0};", ms2mzlist[ims2mz]);
                                }
                                swms2.Write("\t");
                                for (int ims2in = 0; ims2in < ms2intensitylist.Count; ims2in++)
                                {
                                    swms2.Write("{0};", ms2intensitylist[ims2in]);
                                }
                                swms2.Write("\t");
                                for (int ims2noise = 0; ims2noise < ms2noiselist.Count; ims2noise++)
                                {
                                    swms2.Write("{0};", ms2noiselist[ims2noise]);
                                }
                                swms2.Write("\t");
                                for (int ims2base = 0; ims2base < ms2baselinelist.Count; ims2base++)
                                {
                                    swms2.Write("{0};", ms2baselinelist[ims2base]);
                                }
                                swms2.Write("\n");
                            }
                        }
                        swms1.Flush();
                        swms1.Close();
                        swms2.Flush();
                        swms2.Close();
                    }
                    #endregion
                    #region
                    else if (extracttype == "LocalData(SMPR)")
                    {
                        // creat a new file name for "LocalData(SMPR)";
                        string outfile = newfilepath + "\\" + "LocalData(SMPR)_" + filename + ".txt";
                        StreamWriter sw = File.AppendText(outfile);

                        sw.Write("ScanNum\tMSOrder\tPeakNumber\tRT\n");

                        int specnum = 0;
                        rawfile.GetNumSpectra(ref specnum);
                        for (int i = 1; i <= specnum; i++)
                        {
                            System.Threading.Thread.Sleep(0);
                            SetTextMessage(100 * i / specnum);

                            int nOrder = 0;
                            rawfile5.GetMSOrderForScanNum(i, ref nOrder);
                            double scanRT = 0;
                            rawfile.RTFromScanNum(i, ref scanRT);

                            sw.Write("{0}\t", i);
                            sw.Write("{0}\t", nOrder);

                            object varlabels = null;
                            object varflags = null;
                            rawfile5.GetLabelData(ref varlabels, ref varflags, ref i);
                            var labels = (double[,])varlabels;
                            int dim = labels.GetLength(1);
                            sw.Write("{0}\t", dim);
                            sw.Write("{0}\n", scanRT);


                        }
                        sw.Flush();
                        sw.Close();
                    }
                    //
                    else
                    {
                        MessageBox.Show("Please Select a type below!");
                    }
                    #endregion
                }
                MessageBox.Show("Extract Done!");
            }
            else
            {
                MessageBox.Show("Please Select a File!");
            }
           #endregion
        }

        //remove the file you choose;
        private void remove_Click(object sender, EventArgs e)
        {
            object[] selected_objs=new object[pathway.SelectedItems.Count];
            pathway.SelectedItems.CopyTo(selected_objs,0);
            foreach (object item in selected_objs)
            {
                pathway.Items.Remove(item);
            }
            if (selected_objs.Length == 0)
            {
                MessageBox.Show("Please select a file!");
            }
        }

        //clear all the files in the pathway listBox;
        private void clear_Click(object sender, EventArgs e)
        {
            pathway.Items.Clear();
        }

        private void close_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void pathway_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void Progress_Enter(object sender, EventArgs e)
        {

        }

    }
}
