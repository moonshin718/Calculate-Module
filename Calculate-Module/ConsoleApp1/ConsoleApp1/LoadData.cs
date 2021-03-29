using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ConsoleApp1
{
    class LoadData
    {
        public List<Si> Load_Si()
        {
            List<Si> L_si = new List<Si>();

            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\Si.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                Si data_si = new Si();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_si.eV = tmp;
                tmp = double.Parse(Data[1]);
                data_si.e1 = tmp;
                tmp = double.Parse(Data[2]);
                data_si.e2 = tmp;
                L_si.Add(data_si);

            }
            return L_si;
        }

        public List<SiO2> Load_SIO2()
        {
            List<SiO2> L_sio2 = new List<SiO2>();


            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\SIO2.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiO2 data_sio2 = new SiO2();
                string[] Data = token[i].Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_sio2.angstroms = tmp;
                tmp = double.Parse(Data[1]);
                data_sio2.n = tmp;
                tmp = double.Parse(Data[2]);
                data_sio2.k = tmp;
                L_sio2.Add(data_sio2);
            }

            return L_sio2;
        }

        public List<SiN> Load_SiN()
        {
            List<SiN> L_SiN = new List<SiN>();

            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\SiN.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiN data_sin = new SiN();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_sin.nm = tmp;
                tmp = double.Parse(Data[1]);
                data_sin.n = tmp;
                tmp = double.Parse(Data[2]);
                data_sin.k = tmp;
                L_SiN.Add(data_sin);

            }
            return L_SiN;
        }

        public List<SiO2nm_on_Si> Load_SiO2nm_on_Si()
        {
            List<SiO2nm_on_Si> L_SiO2_nm_on_Si = new List<SiO2nm_on_Si>();
            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\SiO2 2nm_on_Si.dat";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiO2nm_on_Si data_SiO2_nm_on_Si = new SiO2nm_on_Si();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_SiO2_nm_on_Si.nm = tmp;
                tmp = double.Parse(Data[1]);
                data_SiO2_nm_on_Si.aoi = tmp;
                tmp = double.Parse(Data[2]);
                data_SiO2_nm_on_Si.Psi = tmp;
                tmp = double.Parse(Data[3]);
                data_SiO2_nm_on_Si.Delta = tmp;
                L_SiO2_nm_on_Si.Add(data_SiO2_nm_on_Si);

            }
            return L_SiO2_nm_on_Si;
        }

        public List<Si_new> Load_Si_new()
        {
            List<Si_new> L_si_new = new List<Si_new>();

            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\Si_new.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                Si_new data_si_new = new Si_new();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_si_new.nm = tmp;
                tmp = double.Parse(Data[1]);
                data_si_new.n = tmp;
                tmp = double.Parse(Data[2]);
                data_si_new.k = tmp;
                L_si_new.Add(data_si_new);

            }

            return L_si_new;
        }

        public List<SiO2_new> Load_SIO2_new()
        {
            List<SiO2_new> L_sio2_new = new List<SiO2_new>();


            String filename = @"C:\Users\bit\source\repos\ConsoleApp1\ConsoleApp1\SiO2_2nm_on_Si_new.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiO2_new data_sio2_new = new SiO2_new();
                string[] Data = token[i].Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                data_sio2_new.nm = tmp;
                tmp = double.Parse(Data[1]);
                data_sio2_new.n = tmp;
                tmp = double.Parse(Data[2]);
                data_sio2_new.k = tmp;
                L_sio2_new.Add(data_sio2_new);
            }

            return L_sio2_new;
        }
    }
}
