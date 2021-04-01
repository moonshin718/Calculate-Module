using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class FileRead
    {
        public List<SiO2_1000nm_on_Si> R_SiO2_1000nm_on_Si(string path)
        {
            List<SiO2_1000nm_on_Si> L_SiO2_1000nm_on_Si = new List<SiO2_1000nm_on_Si>();

            String filename = path + @"\SiO2_1000nm_on_Si.dat";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiO2_1000nm_on_Si SiO2_1000nm_on_Si_data = new SiO2_1000nm_on_Si();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                if(tmp>350&&tmp<1000)
                {
                    SiO2_1000nm_on_Si_data.nm = tmp;
                    tmp = double.Parse(Data[1]);
                    SiO2_1000nm_on_Si_data.aoi = tmp;
                    tmp = double.Parse(Data[2]);
                    SiO2_1000nm_on_Si_data.Psi = tmp;
                    tmp = double.Parse(Data[3]);
                    SiO2_1000nm_on_Si_data.Delta = tmp;
                    L_SiO2_1000nm_on_Si.Add(SiO2_1000nm_on_Si_data);
                }
            }
            return L_SiO2_1000nm_on_Si;
        }
        public List<Si_new> Read_Si_new(string path)
        {
            List<Si_new> L_Si_new = new List<Si_new>();

            String filename = path + @"\Si_new.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                Si_new si_data = new Si_new();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                si_data.nm = tmp;
                tmp = double.Parse(Data[1]);
                si_data.n = tmp;
                tmp = double.Parse(Data[2]);
                si_data.k = tmp;
                L_Si_new.Add(si_data);

            }
            return L_Si_new;
        }
        public List<SiO2_new> Read_SiO2_new(string path)
        {
            List<SiO2_new> L_SiO2_new = new List<SiO2_new>();

            String filename = path + @"\SiO2_new.txt";
            string textvalue = System.IO.File.ReadAllText(filename);
            string paserkey = "\r\n";
            string[] token = textvalue.Split(paserkey, StringSplitOptions.RemoveEmptyEntries);
            for (int i = 1; i < token.Length; i++)
            {
                SiO2_new sio2_data = new SiO2_new();
                string tab = "\t";
                string[] Data = token[i].Split(tab, StringSplitOptions.RemoveEmptyEntries);
                double tmp = double.Parse(Data[0]);
                sio2_data.nm = tmp;
                tmp = double.Parse(Data[1]);
                sio2_data.n = tmp;
                tmp = double.Parse(Data[2]);
                sio2_data.k = tmp;
                L_SiO2_new.Add(sio2_data);

            }
            return L_SiO2_new;
        }
    }
}
