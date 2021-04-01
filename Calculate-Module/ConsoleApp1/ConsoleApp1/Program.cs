using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;
using System.Reflection;

namespace ConsoleApp1
{

    class Program
    {

        static void Main(string[] args)
        {
            string mainpath = Directory.GetParent(Environment.CurrentDirectory).Parent.FullName;
            List<SiO2nm_on_Si> L_SiO2nm_on_Si = new List<SiO2nm_on_Si>(); List<SiO2nm_on_Si_new> L_SiO2nm_on_Si_new = new List<SiO2nm_on_Si_new>();
            List<Si> L_Si = new List<Si>(); List<SiO2> L_SiO2 = new List<SiO2>(); List<SiN> L_SiN = new List<SiN>();
            List<Si_nm> L_Si_nm = new List<Si_nm>(); List<SiO2_nm> L_SiO2_nm = new List<SiO2_nm>();
            List<Si_new> L_Si_new = new List<Si_new>(); List<SiO2_new> L_SiO2_new = new List<SiO2_new>(); List<SiN_new> L_SiN_new = new List<SiN_new>();
            LoadData loding = new LoadData();    MakeData making = new MakeData();    Calculate cal = new Calculate();
            double aoi_p = cal.DegreeToRadian(45);
            double theta65 = cal.DegreeToRadian(65);
            //데이터 로딩
            L_SiO2nm_on_Si = loding.Load_SiO2nm_on_Si();
            L_Si = loding.Load_Si();
            L_SiO2 = loding.Load_SIO2();
            L_SiN = loding.Load_SiN();

            //1-1 Psi,Delta->Alpha,Beta
            L_SiO2nm_on_Si_new = making.make_SiO2nm_on_Si_new(L_SiO2nm_on_Si);

            //물성값 ->nm,n,k 변환
            L_Si_nm = making.make_nm_Si(L_Si);
            L_SiO2_nm = making.make_nm_SiO2(L_SiO2);

            //1-2 데이터 보간 파장 일치
            L_Si_new = making.make_new_Si(L_Si_nm, L_SiO2nm_on_Si_new);
            L_SiO2_new = making.make_new_SiO2(L_SiO2_nm, L_SiO2nm_on_Si_new);
            L_SiN_new = making.make_new_SiN(L_SiN, L_SiO2nm_on_Si_new);


            //1-3
            //700nm에서 rs rp 
            int k = 0;
            for (int i = 0; i < L_Si_new.Count; i++)
            {
                if (L_Si_new[i].nm == 700)
                {
                    k = i;
                }
            }
            StreamWriter outputFile_700 = new StreamWriter(mainpath + @"\SiO2nm_on_Si_700.txt");
            string Header700 = "aoi" + "\t" + "Rs" + "\t" + "Rp" + "\r\n";
            outputFile_700.Write(Header700);
            for (int i = 40; i < 85; i += 2)
            {
                double radianAOI = cal.DegreeToRadian(i);
                Complex N2 = new Complex(L_Si_new[k].n, L_Si_new[k].k);
                Complex theta2 = cal.CalSnell(cal.DegreeToRadian(i), 1, N2);

                Complex rp = cal.Calrp(radianAOI, theta2, 1, N2);
                Complex rs = cal.Calrs(radianAOI, theta2, 1, N2);

                double R_r = cal.CalR_rate(rs);
                double R_p = cal.CalR_rate(rp);
                double alpha = cal.calculateAlpha_cal(radianAOI, rs, rp);
                double beta = cal.calculateBeta_cal(radianAOI, rs, rp);

                outputFile_700.Write(i + "\t");
                outputFile_700.Write(R_r + "\t");
                outputFile_700.Write(R_p + "\r\n");
            }
            outputFile_700.Close();
            //350~ 1000사이의 파장 rs rp
            foreach (var item in L_SiN_new)
            {
                Complex N2 = new Complex(item.n, item.k);
                Complex theta2 = cal.CalSnell(theta65, 1, N2);

                Complex rp = cal.Calrp(theta65, theta2, 1, N2);
                Complex rs = cal.Calrs(theta65, theta2, 1, N2);
                Complex p = rp / rs;
                double alpha = cal.calculateAlpha_cal(aoi_p, rs, rp);
                double beta = cal.calculateBeta_cal(aoi_p, rs, rp);
                //Console.WriteLine(alpha + "  " + beta);

            }



           

            //mse 계산
            double mse_65 = 0;
            List<CalSpectrum[]> L_CalSpectrum = new List<CalSpectrum[]>();
            
            for (int i = 40; i < 85; i += 5)
            {
                StreamWriter outputFile = new StreamWriter(mainpath + String.Format(@"\MSE\SiO2nm_on_Si_aoi{0}.txt", i));
                string Header = "wavelength" + "\t" + "aoi" + "\t" + "alpha" + "\t" + "beta" + "\r\n";
                outputFile.Write(Header);
                CalSpectrum[] Si_Cal = new CalSpectrum[L_SiO2nm_on_Si_new.Count];
                for (int j = 0; j < L_SiO2nm_on_Si_new.Count; j++)
                {

                    double radianAOI = cal.DegreeToRadian(i);
                    Complex N2 = new Complex(L_Si_new[j].n, L_Si_new[j].k);
                    Complex theta2 = cal.CalSnell(radianAOI, 1, N2);

                    Complex rp = cal.Calrp(radianAOI, theta2, 1, N2);
                    Complex rs = cal.Calrs(radianAOI, theta2, 1, N2);
                    double Rs_r = cal.CalR_rate(rs);
                    double Rs_p = cal.CalR_rate(rp);
                    double alpha = cal.calculateAlpha_cal(aoi_p, rs, rp);
                    double beta = cal.calculateBeta_cal(aoi_p, rs, rp);
                    Si_Cal[j] = new CalSpectrum(L_Si_new[j].nm, i, alpha, beta);

                    outputFile.Write(L_SiO2nm_on_Si_new[j].nm + "\t");
                    outputFile.Write(i + "\t");
                    outputFile.Write(alpha + "\t");
                    outputFile.Write(beta + "\r\n");

                }
                L_CalSpectrum.Add(Si_Cal);
                outputFile.Close();
            }
            for (int i = 0; i < L_CalSpectrum.Count; i++)
            {
                if (L_CalSpectrum[i][0].aoi == 65)
                {
                    mse_65 = cal.calculateMSE(L_CalSpectrum[i], L_SiO2nm_on_Si_new);
                }
            }

            Console.WriteLine(mse_65);
        }
    }
}

