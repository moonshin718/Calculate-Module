using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp2_2
{
    class MSpectrum
    {
        private double wavelength;
        private double aoi;
        private double alpha;
        private double beta;
        public double Wavelength
        {
            set
            {
                wavelength = value;
            }
            get
            {
                return wavelength;
            }
        }
        public double AOI
        {
            set
            {
                aoi = value;
            }
            get
            {
                return aoi;
            }
        }
        public double Alpha
        {
            set
            {
                alpha = value;
            }
            get
            {
                return alpha;
            }
        }
        public double Beta
        {
            set
            {
                beta = value;
            }
            get
            {
                return beta;
            }
        }
        public MSpectrum(double wavelength, double aoi, double alpha, double beta)
        {
            this.wavelength = wavelength;
            this.aoi = aoi;
            this.alpha = alpha;
            this.beta = beta;
        }
        public override string ToString()
        {
            return this.wavelength + "\t" + this.aoi + "\t" + this.alpha + "\t" + this.beta;
        }
    }
    class Material
    {
        private double wavelength;
        private double refractIdx;
        private double ExtinctionCoeff;
        public double Wavelength
        {
            set
            {
                wavelength = value;
            }
            get
            {
                return wavelength;
            }
        }
        public double RefractIdx
        {
            set
            {
                refractIdx = value;
            }
            get
            {
                return refractIdx;
            }
        }
        public double ExtinctCoeff
        {
            set
            {
                ExtinctionCoeff = value;
            }
            get
            {
                return ExtinctionCoeff;
            }
        }
        public Material(double wavelength, double refractIdx, double extinctionCoeff)
        {
            this.wavelength = wavelength;
            this.refractIdx = refractIdx;
            this.ExtinctionCoeff = extinctionCoeff;
        }
        public override string ToString()
        {
            return this.wavelength + "\t" + this.refractIdx + "\t" + this.ExtinctionCoeff;
        }
    }
    class Reflectance
    {
        double Lambda;
        double AOI;
        Complex r_p;
        Complex r_s;
        Complex rho;
        double Rp;
        double Rs;
        public double lambda
        {
            get
            {
                return Lambda;
            }
            set
            {
                Lambda = value;
            }
        }
        public double aoi
        {
            get
            {
                return AOI;
            }
            set
            {
                AOI = value;
            }
        }
        public Complex Rcoeff_p
        {
            get
            {
                return r_p;
            }
            set
            {
                r_p = value;
            }
        }
        public Complex Rcoeff_s
        {
            get
            {
                return r_s;
            }
            set
            {
                r_s = value;
            }
        }
        public Complex Rho
        {
            get
            {
                return rho;
            }
            set
            {
                rho = value;
            }
        }
        public double R_P
        {
            get
            {
                return Rp;
            }
            set
            {
                Rp = value;
            }
        }
        public double R_S
        {
            get
            {
                return Rs;
            }
            set
            {
                Rs = value;
            }
        }
        public Reflectance(double lambda, double aOI, Complex p_r, Complex s_r, Complex Rho, double rp, double rs)
        {
            Lambda = lambda;
            AOI = aOI;
            r_p = p_r;
            r_s = s_r;
            rho = Rho;
            Rp = rp;
            Rs = rs;
        }
        public override string ToString()
        {
            return Lambda + "\t" + AOI + "\t" + r_p + "\t" + r_s + "\t" + rho + "\t" + Rp + "\t" + Rs;
        }
    }
    class Error
    {
        double wavelength;
        double diffAlpha;
        double diffBeta;
        double mse;
        public Error(double wavelength, double diffAlpha, double diffBeta, double mse)
        {
            this.wavelength = wavelength;
            this.diffAlpha = diffAlpha;
            this.diffBeta = diffBeta;
            this.mse = mse;
        }
        public override string ToString()
        {
            return this.wavelength + "\t" + this.diffAlpha + "\t" + this.diffBeta + "\t" + this.mse;
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            string path = "C:/Users/mhshi/Documents/Data/SiO2_1000nm_on_Si.dat";
            string path2 = "C:/Users/mhshi/Documents/Data/Si_new.txt";
            string path3 = "C:/Users/mhshi/Documents/Data/SiO2_new.txt";
            string[] textData = System.IO.File.ReadAllLines(path);
            string[] textData2 = System.IO.File.ReadAllLines(path2);
            string[] textData3 = System.IO.File.ReadAllLines(path3);
            string[] replace = { "", ",", "\t", "\n" };
            List<MSpectrum> SiO2_1000nm_on_Si_List = new List<MSpectrum>();
            for (int i = 1; i < textData.Length; i++)
            {
                double lambda = double.Parse(textData[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double aoi = double.Parse(textData[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double psi = double.Parse(textData[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                double delta = double.Parse(textData[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[3]);
                double psiRadian = degToRad(psi);
                double deltaRadian = degToRad(delta);
                double alpha = calcAlpha(psiRadian, degToRad(45));
                double beta = calcBeta(deltaRadian, psiRadian, degToRad(45));
                MSpectrum SiO2_1000nm_on_Si = new MSpectrum(lambda, aoi, alpha, beta);
                SiO2_1000nm_on_Si_List.Add(SiO2_1000nm_on_Si);
            }

            List<Material> Si_Data_List = new List<Material>();
            for (int i = 1; i < textData2.Length; i++) {
                double lambda = double.Parse(textData2[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double n = double.Parse(textData2[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double k = double.Parse(textData2[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                Material Si_data = new Material(lambda, n, k);
                Si_Data_List.Add(Si_data);
            }

            List<Material> SiO2_Data_List = new List<Material>();
            for (int i = 1; i < textData2.Length; i++)
            {
                double lambda = double.Parse(textData3[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double n = double.Parse(textData3[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double k = double.Parse(textData3[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                Material SiO2_data = new Material(lambda, n, k);
                SiO2_Data_List.Add(SiO2_data);
            }
            double mse = 0;

            string header = "Thickness(nm)"+ "\t" + "mse";
            List<string> mse_list = new List<string>();
            mse_list.Add(header);
            for (double i = 700; i <= 1300; i+=5)
            {
                
                mse = 0;
                for (int j = 0; j < SiO2_1000nm_on_Si_List.Count; j ++)
                {
                    double lambda = SiO2_1000nm_on_Si_List[j].Wavelength;
                    double aoi = SiO2_1000nm_on_Si_List[j].AOI;
                    double aoiRadian = degToRad(aoi);
                    double alphaExp = SiO2_1000nm_on_Si_List[j].Alpha;
                    double betaExp = SiO2_1000nm_on_Si_List[j].Beta;

                    double Si_n = Si_Data_List[j].RefractIdx;
                    double Si_k = Si_Data_List[j].ExtinctCoeff;
                    double SiO2_n = SiO2_Data_List[j].RefractIdx;
                    double SiO2_k = SiO2_Data_List[j].ExtinctCoeff;

                    Complex Si_N = new Complex(Si_n, -Si_k);
                    Complex SiO2_N = new Complex(SiO2_n, -SiO2_k);
                    Complex Air_N = new Complex(1, -0);

                    Complex theta_0 = new Complex(aoiRadian, 0);
                    Complex theta_1 = Complex.Asin((Air_N * Complex.Sin(theta_0)) / SiO2_N);
                    Complex theta_2 = Complex.Asin((SiO2_N * Complex.Sin(theta_1)) / Si_N);

                    Complex r01_p = refCoeff_p(Air_N, SiO2_N, theta_0, theta_1);
                    Complex r01_s = refCoeff_s(Air_N, SiO2_N, theta_0, theta_1);
                    Complex r12_p = refCoeff_p(SiO2_N, Si_N, theta_1, theta_2);
                    Complex r12_s = refCoeff_s(SiO2_N, Si_N, theta_1, theta_2);

                    Complex phaseBeta = PhaseThick(i, lambda, SiO2_N, theta_1);
                    Complex doubleBeta = 2 * phaseBeta;
                    Complex expBase = new Complex(0, -1);
                    Complex exp = Complex.Exp(expBase);
                    Complex expBeta = Complex.Pow(exp, doubleBeta);

                    Complex r_p = (r01_p + (r12_p * expBeta)) / (1 + r01_p * (r12_p * expBeta));
                    Complex r_s = (r01_s + (r12_s * expBeta)) / (1 + r01_s * (r12_s * expBeta));

                    double alphaCal = calcAlpha(r_s, r_p, 45.0);
                    double betaCal = calcBeta(r_s, r_p, 45.0);

                    mse += MSE(alphaExp, alphaCal, betaExp, betaCal);
                    //if (i == 1000) {
                    //    Console.WriteLine(lambda + " " + i + " " + (alphaExp - alphaCal) + " " + (betaExp - betaCal) + " " + (1 / (double)SiO2_1000nm_on_Si_List.Count) * mse);
                    //}
                    
                    
                }
                mse_list.Add(i + "\t" +  (1 / (double)SiO2_1000nm_on_Si_List.Count) * mse);
            }

            System.IO.File.WriteAllLines("C:/Users/mhshi/Documents/Data/2_2_mse.txt", mse_list);
        }
        static double degToRad(double degree)
        {
            return (Math.PI * (degree / 180.0f));
        }
        // Calculates the value of alpha using the value of Psi and Polarizer Angle(Polar);
        static double calcAlpha(double Psi, double Polar)
        {
            double tanPsi = Math.Tan(Psi);
            double tanAOI = Math.Tan(Polar);
            double num = (Math.Pow(tanPsi, 2) - (Math.Pow(tanAOI, 2)));
            double deno = (Math.Pow(tanPsi, 2) + (Math.Pow(tanAOI, 2)));
            return num / deno;
        }
        static double calcAlpha(Complex r_s, Complex r_p, double P)
        {
            double num = ((Math.Pow(r_p.Magnitude, 2) * Math.Pow(Math.Cos(degToRad(P)), 2))) - ((Math.Pow(r_s.Magnitude, 2) * Math.Pow(Math.Sin(degToRad(P)), 2)));
            double deno = ((Math.Pow(r_p.Magnitude, 2) * Math.Pow(Math.Cos(degToRad(P)), 2))) + ((Math.Pow(r_s.Magnitude, 2) * Math.Pow(Math.Sin(degToRad(P)), 2)));
            double Alpha = num / deno;
            return Alpha;
        }
        static double calcBeta(Complex r_s, Complex r_p, double P)
        {
            double num = ((2 * (r_p * Complex.Conjugate(r_s)).Real) * Math.Cos(degToRad(P))) * Math.Sin(degToRad(P));
            double deno = (Math.Pow(r_p.Magnitude, 2) * Math.Pow(Math.Cos(degToRad(P)), 2)) + (Math.Pow(r_s.Magnitude, 2) * Math.Pow(Math.Sin(degToRad(P)), 2));
            double Beta = num / deno;
            return Beta;
        }
        // Calculates the value of beta using the value of Psi, Delta and Polarizer Angle(Polar)
        static double calcBeta(double Delta, double Psi, double Polar)
        {
            double tanPsi = Math.Tan(Psi);
            double tanPolar = Math.Tan(Polar);
            double deno = (2 * Math.Cos(Delta) * tanPsi * tanPolar);
            double num = ((Math.Pow(tanPsi, 2) + (Math.Pow(tanPolar, 2))));
            return deno / num;
        }
        // Convert the value of Photon Energy to wavelength in nm
        static double cvtElVtoLambda(double PhotonEn)
        {
            return (1.99 * Math.Pow(10, -25)) / (PhotonEn * (1.6 * Math.Pow(10, -19)));
        }
        // Convert the value of Angstrom to wavelength in nm
        static double cvtANGtoLambda(double Angstrom)
        {
            return 0.1 * Angstrom;
        }
        // Convert epsilons to Refractive Index
        static double cvtRefractIdx(double dielectric_1, double dielectric_2)
        {
            return Math.Sqrt(0.5 * (Math.Sqrt(Math.Pow(dielectric_1, 2) + Math.Pow(dielectric_2, 2)) + dielectric_1));
        }
        // Convert epsilons to Extinction Coefficient
        static double cvtExtinctCoeff(double dielectric_1, double dielectric_2)
        {
            return Math.Sqrt(0.5 * (Math.Sqrt(Math.Pow(dielectric_1, 2) + Math.Pow(dielectric_2, 2)) - dielectric_1));
        }
        static double MSE(double alphaExp, double alphaCal, double betaExp, double betaCal)
        {
            double mse = Math.Pow(alphaExp - alphaCal, 2) + Math.Pow(betaExp - betaCal, 2);
            //Console.WriteLine("alpha 오차 값 " + (alphaExp - alphaCal) + " " + "beta 오차 " +(betaExp-betaCal));
            return mse;
        }
        static Complex refCoeff_p(Complex N0, Complex N1, Complex theta0, Complex theta1)
        {
            Complex r01p = (N1 * Complex.Cos(theta0) - N0 * Complex.Cos(theta1)) / (N1 * Complex.Cos(theta0) + N0 * Complex.Cos(theta1));
            return r01p;
        }
        static Complex thetaRef(Complex N0, Complex N1, Complex theta0)
        {
            Complex theta1 = Complex.Asin((N0 * Complex.Sin(theta0)) / N1);
            return theta1;
        }
        static Complex refCoeff_s(Complex N0, Complex N1, Complex theta0, Complex theta1)
        {
            Complex r01s = (N0 * Complex.Cos(theta0) - N1 * Complex.Cos(theta1)) / (N0 * Complex.Cos(theta0) + N1 * Complex.Cos(theta1));
            return r01s;
        }
        static Complex transCoeff_p(Complex N0, Complex N1, Complex theta0, Complex theta1)
        {
            Complex t01p = (2 * N0 * Complex.Cos(theta0)) / (N1 * Complex.Cos(theta0) + N0 * Complex.Cos(theta1));
            return t01p;
        }
        static Complex transCoeff_s(Complex N0, Complex N1, Complex theta0, Complex theta1)
        {
            Complex t01s = (2 * N0 * Complex.Cos(theta0)) / (N0 * Complex.Cos(theta0) + N1 * Complex.Cos(theta1));
            return t01s;
        }
        static Complex PhaseThick(double Thickness, double Wavelength, Complex N, Complex Angle)
        {
            Complex Beta = (2 * Math.PI) * (Thickness / Wavelength) * N * Complex.Cos(Angle);
            return Beta;
        }
        static Complex refCoeff_p(Complex r01p, Complex r12p, Complex Beta)
        {
            return 0;
        }
        static Complex refCoeff_s(Complex r01s, Complex r12s, Complex Beta)
        {
            return 0;
        }
        static Complex refCoeff_p(Complex t01p, Complex t10p, Complex r01p, Complex r12p, Complex Beta, int n)
        {
            return 0;
        }
        static Complex refCoeff_s(Complex t01s, Complex t10s, Complex r01s, Complex r12s, Complex Beta, int n)
        {
            return 0;
        }
        static Complex[,] ComplexMatrixMultiply(Complex[,] A, Complex[,] B)
        {
            int rA = A.GetLength(0);
            int cA = A.GetLength(1);
            int rB = B.GetLength(0);
            int cB = B.GetLength(1);
            Complex temp = 0;
            Complex[,] result = new Complex[rA, cB];
            if (cA != rB)
            {
                Console.WriteLine("Matrices cannot be multiplied");
            }
            else
            {
                for (int i = 0; i < rA; i++)
                {
                    for (int j = 0; j < cB; j++)
                    {
                        temp = 0;
                        for (int k = 0; k < cA; k++)
                        {
                            temp += A[i, k] * B[k, j];
                        }
                        result[i, j] = temp;
                    }
                }
            }
            return result;
        }
        static Complex[,] Interface_ij(Complex r_ij, Complex t_ij)
        {
            Complex[,] Interface_01 = { { 1 / t_ij, r_ij / t_ij }, { r_ij / t_ij, 1 / t_ij } };
            return Interface_01;
        }
        static Complex[,] Layer_j(Complex Angle, double Thickness, double Wavelength, Complex N)
        {
            Complex Beta_1 = 2 * Math.PI * (Thickness / Wavelength) * N * Complex.Cos(Angle);
            Complex expBase_11 = new Complex(0, 1);
            Complex expBase_22 = new Complex(0, -1);
            Complex exp_11 = Complex.Exp(expBase_11);
            Complex exp_22 = Complex.Exp(expBase_22);
            Complex expBeta_11 = Complex.Pow(exp_11, Beta_1);
            Complex expBeta_22 = Complex.Pow(exp_22, Beta_1);
            Complex[,] Layer_1 = { { expBeta_11, 0 }, { 0, expBeta_22 } };
            return Layer_1;
        }
        static Complex[,] sMat_p(Complex N_1, Complex N_2, Complex theta_1, Complex theta_2, Complex theta_3, double d1, double d2, double lambda1, double lambda2, Complex[,] s_p, int i, int n)
        {
            if (i != n)
            {
                if (i % 2 == 0)
                {
                    Complex r_ij_p = refCoeff_p(N_1, N_2, theta_1, theta_2);
                    Complex t_ij_p = transCoeff_p(N_1, N_2, theta_1, theta_2);
                    Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
                    Complex[,] layer_j = Layer_j(theta_2, d1, lambda1, N_2);
                    s_p = ComplexMatrixMultiply(s_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
                }
                else
                {
                    Complex r_ij_p = refCoeff_p(N_2, N_1, theta_2, theta_3);
                    Complex t_ij_p = transCoeff_p(N_2, N_1, theta_2, theta_3);
                    Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
                    Complex[,] layer_j = Layer_j(theta_3, d2, lambda2, N_1);
                    s_p = ComplexMatrixMultiply(s_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
                }
                return ComplexMatrixMultiply(s_p, sMat_p(N_1, N_2, theta_1, theta_2, theta_3, d1, d2, lambda1, lambda2, s_p, i + 1, n));
            }
            else
            {
                return s_p;
            }
        }
        static Complex[,] sMat_s(Complex N_1, Complex N_2, Complex theta_1, Complex theta_2, Complex theta_3, double d1, double d2, double lambda1, double lambda2, Complex[,] s_s, int i, int n)
        {
            if (i != n)
            {
                if (i % 2 == 0)
                {
                    Complex r_ij_s = refCoeff_s(N_1, N_2, theta_1, theta_2);
                    Complex t_ij_s = transCoeff_s(N_1, N_2, theta_1, theta_2);
                    Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
                    Complex[,] layer_j = Layer_j(theta_2, d1, lambda1, N_2);
                    s_s = ComplexMatrixMultiply(s_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
                }
                else
                {
                    Complex r_ij_s = refCoeff_s(N_2, N_1, theta_2, theta_3);
                    Complex t_ij_s = transCoeff_s(N_2, N_1, theta_2, theta_3);
                    Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
                    Complex[,] layer_j = Layer_j(theta_3, d2, lambda2, N_1);
                    s_s = ComplexMatrixMultiply(s_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
                }
                return ComplexMatrixMultiply(s_s, sMat_p(N_1, N_2, theta_1, theta_2, theta_3, d1, d2, lambda1, lambda2, s_s, i + 1, n));
            }
            else
            {
                return s_s;
            }
        }
    }
}
