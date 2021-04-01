using System;
using System.Collections.Generic;
using System.Numerics;
using System.Text;

namespace ConsoleApp2
{
    public class Calculate
    {
        public double calculateAlpha_exp(double psi, double delta)
        {
            double Alpha_exp = 0;
            double tan45 = Math.Tan(DegreeToRadian(45));
            double TanPsi = Math.Tan(psi);
            double TanDelta = Math.Tan(delta);
            Alpha_exp = (Math.Pow(TanPsi, 2) - Math.Pow(tan45, 2)) / (Math.Pow(TanPsi, 2) + Math.Pow(tan45, 2));
            return Alpha_exp;
        }
        public double calculateBeta_exp(double psi, double delta)
        {
            double Beta_exp = 0;
            double tan45 = Math.Tan(DegreeToRadian(45));
            double TanPsi = Math.Tan(psi);
            double TanDelta = Math.Tan(delta);
            Beta_exp = (2 * Math.Cos(delta) * TanPsi * tan45) / (Math.Pow(TanPsi, 2) + Math.Pow(tan45, 2));
            return Beta_exp;
        }
        public double DegreeToRadian(double angle)
        {
            return (Math.PI * angle) / 180.0f;
        }
        public Complex Cal_Snell(Complex theta0, Complex N1, Complex N2)
        {

            Complex theta1 = Complex.Asin((N1 / N2) * Complex.Sin(theta0));
            return theta1;
        }
        public Complex Cal_rp(Complex theta1, Complex theta2, Complex N0, Complex N1)
        {
            Complex rp = (N1 * Complex.Cos(theta1) - N0 * Complex.Cos(theta2)) / (N1 * Complex.Cos(theta1) + N0 * Complex.Cos(theta2));
            return rp;
        }
        public Complex Cal_rs(Complex theta1, Complex theta2, Complex N0, Complex N1)
        {
            Complex rs = (N0 * Complex.Cos(theta1) - N1 * Complex.Cos(theta2)) / (N0 * Complex.Cos(theta1) + N1 * Complex.Cos(theta2));
            return rs;
        }
        public Complex Cal_rn(double n,Complex r01,Complex r12,Complex e_minus_2bi)
        {
            Complex r_n = new Complex();
            Complex C_ratio = new Complex();

            for (int i=0;i<n;i++)
            {
                C_ratio +=  Complex.Pow((-r01) * r12 * e_minus_2bi, i);
            }
            r_n = r01+ ((1 - Complex.Pow(r01, 2)) * r12 * e_minus_2bi) * C_ratio;
            return r_n;
        }

        public double CalR_rate(Complex r)
        {
            double R_r = Math.Pow(r.Magnitude, 2);

            return R_r;
        }
        public double calculateAlpha_cal(double aoi_p, Complex rs, Complex rp)
        {
            double Rs_r = CalR_rate(rs);
            double Rp_r = CalR_rate(rp);
            double alpha = (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) - Rs_r * Math.Pow(Math.Sin(aoi_p), 2)) / (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) + Rs_r * Math.Pow(Math.Sin(aoi_p), 2));

            return alpha;
        }
        public double calculateBeta_cal(double aoi_p, Complex rs, Complex rp)
        {
            double Rs_r = CalR_rate(rs);
            double Rp_r = CalR_rate(rp);
            double beta = (2 * (rp * Complex.Conjugate(rs)).Real * Math.Cos(aoi_p) * Math.Sin(aoi_p)) / (Rp_r * Math.Pow(Math.Cos(aoi_p), 2) + Rs_r * Math.Pow(Math.Sin(aoi_p), 2));

            return beta;
        }
        public double calculateMSE(List<CalSpectrum> L_CalSpectrum, List<SiO2_1000nm_on_Si_new> data)
        {
            double mse = 0;
            for (int i = 0; i < L_CalSpectrum.Count; i++)
            {
                mse += Math.Pow(data[i].alpha - L_CalSpectrum[i].alpha, 2) + Math.Pow(data[i].beta - L_CalSpectrum[i].beta, 2);
            }

            return (1 / (double)L_CalSpectrum.Count) * mse;
        }

    }
}
