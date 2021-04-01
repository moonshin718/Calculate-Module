using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace ConsoleApp1
{
    class Calculate
    {
        public double DegreeToRadian(double angle)
        {
            return (Math.PI * angle) / 180.0f;
        }
        public double calculateAlpha_exp(double psi, double delta)
        {
            double result = 0;
            double tan45 = Math.Tan(DegreeToRadian(45));
            double TanPsi = Math.Tan(psi);
            double TanDelta = Math.Tan(delta);
            result = (Math.Pow(TanPsi, 2) - Math.Pow(tan45, 2)) / (Math.Pow(TanPsi, 2) + Math.Pow(tan45, 2));
            return result;
        }
        public double calculateBeta_exp(double psi, double delta) 
        {
            double result = 0;
            double tan45 = Math.Tan(DegreeToRadian(45));
            double TanPsi = Math.Tan(psi);
            double TanDelta = Math.Tan(delta);
            result = (2 * Math.Cos(delta) * TanPsi * tan45) / (Math.Pow(TanPsi, 2) + Math.Pow(tan45, 2));
            return result;
        }

        // aoi= radian
        public Complex CalSnell(double aoi, Complex N1, Complex N2)
        {

            Complex aoi_k = Complex.Asin((N1 / N2) * Complex.Sin(aoi));
            return aoi_k;
        }
        //theta1= radian
        public Complex Calrp(double theta1, Complex theta2, Complex N1, Complex N2)
        {
            Complex rp = (N2 * Complex.Cos(theta1) - N1 * Complex.Cos(theta2)) / (N2 * Complex.Cos(theta1) + N1 * Complex.Cos(theta2));
            return rp;
        }
        //theta1= radian
        public Complex Calrs(double theta1, Complex theta2, Complex N1, Complex N2)
        {
            Complex rs = (N1 * Complex.Cos(theta1) - N2 * Complex.Cos(theta2)) / (N1 * Complex.Cos(theta1) + N2 * Complex.Cos(theta2));
            return rs;
        }

        public double CalR_rate(Complex r)
        {
            double R_r = Math.Pow(r.Magnitude, 2);

            return R_r;
        }

        //p= radian
        public double calculateAlpha_cal(double aoi_p, Complex rs, Complex rp)
        {
            double Rs_r = CalR_rate(rs);
            double Rp_r = CalR_rate(rp);
            double alpha = (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) - Rs_r * Math.Pow(Math.Sin(aoi_p), 2)) / (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) + Rs_r * Math.Pow(Math.Sin(aoi_p), 2));

            return alpha;
        }

        //p= radian
        public double calculateBeta_cal(double aoi_p, Complex rs, Complex rp)
        {
            double Rs_r = CalR_rate(rs);
            double Rp_r = CalR_rate(rp);
            double beta = (2 * (rp * Complex.Conjugate(rs)).Real * Math.Cos(aoi_p) * Math.Sin(aoi_p)) / (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) + Rs_r * Math.Pow(Math.Sin(aoi_p), 2));

            return beta;
        }

        public double calculateMSE(CalSpectrum[] calSpectrums, List<SiO2nm_on_Si_new> data)
        {
            double mse = 0;
            for (int i = 0; i < calSpectrums.Length; i++)
            {
                mse += Math.Pow(data[i].alpha - calSpectrums[i].alpha, 2) + Math.Pow(data[i].beta - calSpectrums[i].beta, 2);
            }

            return (1 / (double)calSpectrums.Length) * mse;
        }
    }
}
