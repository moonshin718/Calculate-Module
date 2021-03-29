﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Numerics;

namespace ConsoleApp1
{

    class Program
    {

        static void Main(string[] args)
        {
            List<SiO2nm_on_Si> L_SiO2nm_on_Si = new List<SiO2nm_on_Si>();         List<SiO2nm_on_Si_new> L_SiO2nm_on_Si_new = new List<SiO2nm_on_Si_new>();
            List<Si> L_Si = new List<Si>();       List<SiO2> L_SiO2 = new List<SiO2>();        List<SiN> L_SiN = new List<SiN>();
            List<Si_nm> L_Si_nm = new List<Si_nm>();         List<SiO2_nm> L_SiO2_nm = new List<SiO2_nm>();
            List<Si_new> L_Si_new = new List<Si_new>();      List<SiO2_new> L_SiO2_new = new List<SiO2_new>();     List<SiN_new> L_SiN_new = new List<SiN_new>();
            LoadData loding = new LoadData(); MakeData making = new MakeData();
            Calculate cal=new Calculate();
            L_SiO2nm_on_Si = loding.Load_SiO2nm_on_Si();
            L_Si = loding.Load_Si();
            L_SiO2 = loding.Load_SIO2();
            L_SiN = loding.Load_SiN();
            L_SiO2nm_on_Si_new = making.make_SiO2nm_on_Si_new(L_SiO2nm_on_Si);
            L_Si_nm = making.make_nm_Si(L_Si);
            L_SiO2_nm=making.make_nm_SiO2(L_SiO2);

            L_Si_new=making.make_new_Si(L_Si_nm, L_SiO2nm_on_Si_new);
            L_SiO2_new = making.make_new_SiO2(L_SiO2_nm, L_SiO2nm_on_Si_new);
            L_SiN_new = making.make_new_SiN(L_SiN, L_SiO2nm_on_Si_new);

            double theta45 = cal.DegreeToRadian(45);
            double theta65 = cal.DegreeToRadian(65);

            //foreach (var item in L_SiN_new)
            //{
            //    Complex N2 = new Complex(item.n, item.k);
            //    Complex theta2 = cal.CalSnell(theta65, 1, N2);

            //    Complex rp = cal.Calrp(theta65, theta2, 1, N2);
            //    Complex p=rp/rs;
            //    double alpha = cal.calculateAlpha_cal(theta45, rs, rp);
            //    double beta = cal.calculateBeta_cal(theta45, rs, rp);
            //    Console.WriteLine(alpha+"  "+ beta);

            //}
            int k=0;
            for(int i=0;i< L_Si_new.Count;i++)
            { 
                if(L_Si_new[i].nm==700)
                {
                    k = i;
                }
            }

            for (int i = 40; i < 85; i+=2)
            {
                Complex N2 = new Complex(L_Si_new[k].n, L_Si_new[k].k);
                Complex theta2 = cal.CalSnell(cal.DegreeToRadian(i), 1, N2);

                Complex rp = cal.Calrp(cal.DegreeToRadian(i), theta2, 1, N2);
                Complex rs = cal.Calrs(cal.DegreeToRadian(i), theta2, 1, N2);

                double Rs_r = cal.CalRs_rate(rs);
                double Rs_p = cal.CalRp_rate(rp);
                double alpha = cal.calculateAlpha_cal(cal.DegreeToRadian(i), rs, rp);
                double beta = cal.calculateBeta_cal(cal.DegreeToRadian(i), rs, rp);
                //Console.WriteLine(i+"\t"+Rs_r + "\t" + Rs_p);
                Console.WriteLine(Rs_r + " " + Rs_p);
            }


        }
    }
}