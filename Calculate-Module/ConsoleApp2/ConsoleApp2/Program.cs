using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;

namespace ConsoleApp2
{
    class Program
    {
        static void Main(string[] args)
        {
            List<SiO2_1000nm_on_Si> L_SiO2_1000nm_on_Si = new List<SiO2_1000nm_on_Si>(); List<Si_new> L_Si_new = new List<Si_new>(); List<SiO2_new> L_SiO2_new = new List<SiO2_new>();
            string program_path = Directory.GetParent(Environment.CurrentDirectory).Parent.FullName;
            FileRead reader = new FileRead();
            Calculate cal = new Calculate();
            double aoi_p = cal.DegreeToRadian(45);

            //Si 기판(substrate)과 SiO2(thin film)의 물성 값을 읽기
            L_Si_new = reader.Read_Si_new(program_path);
            L_SiO2_new = reader.Read_SiO2_new(program_path);
            

            //2-1

            //“SiO2 1000nm_on_Si.dat”에서 스펙트럼 읽기
            L_SiO2_1000nm_on_Si = reader.R_SiO2_1000nm_on_Si(program_path);
            List<SiO2_1000nm_on_Si_new> L_SiO2_1000nm_on_Si_new = new List<SiO2_1000nm_on_Si_new>();
            List<CalSpectrum> L_CalSpectrum = new List<CalSpectrum>();
            //exp alphs beta 얻기
            foreach (var item in L_SiO2_1000nm_on_Si)
            {
               
              
            }

            //alpha beta 구하기 모델 세우기
            for (int i=0;i< L_SiO2_1000nm_on_Si.Count;i++)
            {
                double aoi_radian = cal.DegreeToRadian(L_SiO2_1000nm_on_Si[i].aoi);
                Complex N1 = new Complex(L_SiO2_new[i].n, -L_SiO2_new[i].k);
                Complex theta2_01 = cal.Cal_Snell(aoi_radian, 1, N1);
                Complex rp_01 = cal.Cal_rp(aoi_radian, theta2_01, 1, N1);
                Complex rs_01 = cal.Cal_rs(aoi_radian, theta2_01, 1, N1);

                Complex N2 = new Complex(L_Si_new[i].n, -L_Si_new[i].k);
                Complex theta2_12 = cal.Cal_Snell(theta2_01, N1, N2);
                Complex rp_12 = cal.Cal_rp(theta2_01, theta2_12, N1, N2);
                Complex rs_12 = cal.Cal_rs(theta2_01, theta2_12, N1, N2);

                
                Complex beta_thickness = cal.DegreeToRadian(360.0f)*(1000 / L_SiO2_1000nm_on_Si[i].nm) * N1 * Complex.Cos(theta2_01);
                Complex minus_i = new Complex(0, -1);
                Complex e_minus_i = Complex.Exp(minus_i);
                Complex e_minus_2bi = Complex.Pow(e_minus_i, 2 * beta_thickness);

                //무한 등비급수의 일반항으로 얻은 rs,rp
                Complex rp = (rp_01 + rp_12 * e_minus_2bi) / (1 + rp_01 * rp_12 * e_minus_2bi);
                Complex rs = (rs_01 + rs_12 * e_minus_2bi) / (1 + rs_01 * rs_12 * e_minus_2bi);

                //수렴식에 대해 alpha, beta 스펙트럼 -> exp alpha,beta
                if (L_SiO2_1000nm_on_Si[i].nm > 350 && L_SiO2_1000nm_on_Si[i].nm < 1000)
                {
                    SiO2_1000nm_on_Si_new item_data = new SiO2_1000nm_on_Si_new();
                    item_data.nm = L_SiO2_1000nm_on_Si[i].nm;
                    item_data.aoi = L_SiO2_1000nm_on_Si[i].aoi;
                    item_data.alpha = cal.calculateAlpha_cal(aoi_p, rs, rp);
                    item_data.beta = cal.calculateBeta_cal(aoi_p, rs, rp);

                    L_SiO2_1000nm_on_Si_new.Add(item_data);
                }
               

                //항의 개수가 정해저 있을때 rs,rp
                Complex rp_n = cal.Cal_rn(2000, rp_01, rp_12, e_minus_2bi);
                Complex rs_n = cal.Cal_rn(2000, rs_01, rs_12, e_minus_2bi);

                //무한급수의 항의 개수에 따른 alpha, beta
                
                double alpha_cal = cal.calculateAlpha_cal(aoi_p, rs_n, rp_n);
                double beta_cal = cal.calculateBeta_cal(aoi_p, rs_n, rp_n);
                CalSpectrum cal_data = new CalSpectrum(L_SiO2_1000nm_on_Si[i].nm, L_SiO2_1000nm_on_Si[i].aoi, alpha_cal, beta_cal);
                L_CalSpectrum.Add(cal_data);


                Console.WriteLine(L_CalSpectrum[i].alpha+"   "+ L_SiO2_1000nm_on_Si_new[i].alpha);
            }

            //mse 계산
            
            double mse = cal.calculateMSE(L_CalSpectrum, L_SiO2_1000nm_on_Si_new);

            Console.WriteLine(mse);

        }
    }
}
