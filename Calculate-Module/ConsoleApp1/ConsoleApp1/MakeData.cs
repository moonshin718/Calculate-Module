using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ConsoleApp1
{
    class MakeData
    {
        string mainpath = Directory.GetParent(Environment.CurrentDirectory).Parent.FullName;
        public List<SiO2nm_on_Si_new> make_SiO2nm_on_Si_new(List<SiO2nm_on_Si> data)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\SiO2nm_on_Si_new.dat");
            List<SiO2nm_on_Si_new> L_siO2Nm_On_Si_News = new List<SiO2nm_on_Si_new>();
            string Header = "wavelength" + "\t" + "aoi" + "\t" + "alpha" + "\t" + "beta" + "\r\n";
            outputFile.Write(Header);
            Calculate cal = new Calculate();
            foreach (var item in data)
            {
                SiO2nm_on_Si_new SiO2nm_on_Si_new_data = new SiO2nm_on_Si_new();
                if (item.nm > 350 && item.nm < 1000)
                {
                    outputFile.Write(item.nm + "\t");
                    outputFile.Write(item.aoi + "\t");
                    double alpha = cal.calculateAlpha_exp(item.Psi, item.Delta);
                    double beta = cal.calculateAlpha_exp(item.Psi, item.Delta);
                    outputFile.Write(alpha + "\t");
                    outputFile.Write(beta + "\r\n");
                    SiO2nm_on_Si_new_data.nm = item.nm;
                    SiO2nm_on_Si_new_data.aoi = item.aoi;
                    SiO2nm_on_Si_new_data.alpha = alpha;
                    SiO2nm_on_Si_new_data.beta = beta;
                    L_siO2Nm_On_Si_News.Add(SiO2nm_on_Si_new_data);
                }


            }
            outputFile.Close();

            return L_siO2Nm_On_Si_News;
        }

        public List<Si_nm> make_nm_Si(List<Si> data)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\Si_nm.txt");
            List<Si_nm> L_Si_nm = new List<Si_nm>();
            string Header = "wavelength" + "\t" + "n" + "\t" + "k" + "\r\n";
            outputFile.Write(Header);
            foreach (var item in data)
            {
                Si_nm Si_nm_data = new Si_nm();
                double n = Math.Sqrt(0.5 * (Math.Sqrt((Math.Pow(item.e1, 2) + Math.Pow(item.e2, 2))) + item.e1));
                double k = Math.Sqrt(0.5 * (Math.Sqrt((Math.Pow(item.e1, 2) + Math.Pow(item.e2, 2))) - item.e1));
                double nm = 1240 / item.eV;

                outputFile.Write(nm + "\t");
                outputFile.Write(n + "\t");
                outputFile.Write(k + "\r\n");

                Si_nm_data.nm = nm;
                Si_nm_data.n = n;
                Si_nm_data.k = k;
                L_Si_nm.Add(Si_nm_data);
            }
            outputFile.Close();

            return L_Si_nm;
        }
        public List<SiO2_nm> make_nm_SiO2(List<SiO2> data)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\SIO2_nm.txt");
            List<SiO2_nm> L_SiO2_nm = new List<SiO2_nm>();

            string Header = "wavelength" + "\t" + "n" + "\t" + "k" + "\r\n";
            outputFile.Write(Header);


            foreach (var item in data)
            {
                SiO2_nm SiO2_nm_data = new SiO2_nm();
                double nm = item.angstroms * 0.1;
                outputFile.Write(nm + "\t");
                outputFile.Write(item.n + "\t");
                outputFile.Write(item.k + "\r\n");
                SiO2_nm_data.nm = nm;
                SiO2_nm_data.n = item.n;
                SiO2_nm_data.k = item.k;
                L_SiO2_nm.Add(SiO2_nm_data);

            }

            outputFile.Close();
            return L_SiO2_nm;
        }

        public List<Si_new> make_new_Si(List<Si_nm> data, List<SiO2nm_on_Si_new> wave)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\Si_new.txt");
            var n_itp = new Dictionary<double, double>();
            var k_itp = new Dictionary<double, double>();

            string Header = "wavelength" + "\t" + "n" + "\t" + "k" + "\r\n";
            outputFile.Write(Header);

            List<Si_new> L_si_new = new List<Si_new>();

            foreach (var item in data)
            {
                n_itp.Add(item.nm, item.n);
                k_itp.Add(item.nm, item.k);
            }

            var scaler_n = new SplineInterpolator(n_itp);
            var scaler_k = new SplineInterpolator(k_itp);

            double y_n, y_k;
            foreach (var item in wave)
            {

                Si_new inters = new Si_new();
                y_n = scaler_n.GetValue(item.nm);
                y_k = scaler_k.GetValue(item.nm);

                inters.nm = item.nm;
                inters.n = y_n;
                inters.k = y_k;
                L_si_new.Add(inters);
            }


            Si_new inter = new Si_new();
            inter.nm = 700.000;
            inter.n = scaler_n.GetValue(inter.nm);
            inter.k = scaler_k.GetValue(inter.nm);
            L_si_new.Add(inter);
            var sortList = L_si_new.OrderBy(x => x.nm).ToList();

            foreach (var item in sortList)
            {
                outputFile.Write(item.nm + "\t");
                outputFile.Write(item.n + "\t");
                outputFile.Write(item.k + "\r\n");
            }

            outputFile.Close();

            return L_si_new;
        }




        public List<SiO2_new> make_new_SiO2(List<SiO2_nm> data, List<SiO2nm_on_Si_new> wave)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\SiO2_new.txt");
            var n_itp = new Dictionary<double, double>();
            var k_itp = new Dictionary<double, double>();

            string Header = "wavelength" + "\t" + "n" + "\t" + "k" + "\r\n";
            outputFile.Write(Header);

            List<SiO2_new> L_sio2_new = new List<SiO2_new>();


            foreach (var item in data)
            {
                n_itp.Add(item.nm, item.n);
                k_itp.Add(item.nm, item.k);
            }

            var scaler_n = new SplineInterpolator(n_itp);
            var scaler_k = new SplineInterpolator(k_itp);

            double y_n, y_k;
            foreach (var item in wave)
            {
                SiO2_new inters = new SiO2_new();
                y_n = scaler_n.GetValue(item.nm);
                y_k = scaler_k.GetValue(item.nm);
                inters.nm = item.nm;
                inters.n = y_n;
                inters.k = y_k;
                L_sio2_new.Add(inters);
            }
            SiO2_new inter = new SiO2_new();
            inter.nm = 700.000;
            inter.n = scaler_n.GetValue(inter.nm);
            inter.k = scaler_k.GetValue(inter.nm);
            L_sio2_new.Add(inter);
            var sortList = L_sio2_new.OrderBy(x => x.nm).ToList();

            foreach (var item in sortList)
            {
                outputFile.Write(item.nm + "\t");
                outputFile.Write(item.n + "\t");
                outputFile.Write(item.k + "\r\n");
            }

            outputFile.Close();

            return L_sio2_new;
        }

        public List<SiN_new> make_new_SiN(List<SiN> data, List<SiO2nm_on_Si_new> wave)
        {
            StreamWriter outputFile = new StreamWriter(mainpath+@"\SiN_new.txt");
            var n_itp = new Dictionary<double, double>();
            var k_itp = new Dictionary<double, double>();

            string Header = "wavelength" + "\t" + "n" + "\t" + "k" + "\r\n";
            outputFile.Write(Header);

            List<SiN_new> L_SiN_new = new List<SiN_new>();

            foreach (var item in data)
            {

                n_itp.Add(item.nm, item.n);
                k_itp.Add(item.nm, item.k);
            }

            var scaler_n = new SplineInterpolator(n_itp);
            var scaler_k = new SplineInterpolator(k_itp);
            double y_n, y_k;
            foreach (var item in wave)
            {
                SiN_new inters = new SiN_new();
                y_n = scaler_n.GetValue(item.nm);
                y_k = scaler_k.GetValue(item.nm);
                inters.nm = item.nm;
                inters.n = y_n;
                inters.k = y_k;
                L_SiN_new.Add(inters);
            }
            SiN_new inter = new SiN_new();
            inter.nm = 700.000;
            inter.n = scaler_n.GetValue(inter.nm);
            inter.k = scaler_k.GetValue(inter.nm);
            L_SiN_new.Add(inter);
            var sortList = L_SiN_new.OrderBy(x => x.nm).ToList();
            foreach (var item in sortList)
            {
                outputFile.Write(item.nm + "\t");
                outputFile.Write(item.n + "\t");
                outputFile.Write(item.k + "\r\n");
            }

            outputFile.Close();

            return L_SiN_new;
        }
    }
}
