using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp3
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
            string[] Si_DataPath = System.IO.File.ReadAllLines(@"C:/Users/mhshi/OneDrive/문서/Data/Si_new.txt");
            string[] SiN_DataPath = System.IO.File.ReadAllLines(@"C:/Users/mhshi/OneDrive/문서/Data/SiN_new.txt");
            string[] SiO2_DataPath = System.IO.File.ReadAllLines(@"C:/Users/mhshi/OneDrive/문서/Data/SiO2_new.txt");

            string[] replace = { "", ",", "\t", "\n" };

            List<Material> Si_Data = new List<Material>();
            List<Material> SiN_Data = new List<Material>();
            List<Material> SiO2_Data = new List<Material>();

            for (int i = 1; i < SiO2_DataPath.Length; i++) {
                double lambda_Si = double.Parse(Si_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double n_Si = double.Parse(Si_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double k_Si = double.Parse(Si_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                Material Si_Material = new Material(lambda_Si, n_Si, k_Si);
                Si_Data.Add(Si_Material);

                double lambda_SiN = double.Parse(SiN_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double n_SiN = double.Parse(SiN_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double k_SiN = double.Parse(SiN_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                Material SiN_Material = new Material(lambda_SiN, n_SiN, k_SiN);
                SiN_Data.Add(SiN_Material);
                
                double lambda_SiO2 = double.Parse(SiO2_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[0]);
                double n_SiO2 = double.Parse(SiO2_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[1]);
                double k_SiO2 = double.Parse(SiO2_DataPath[i].Split(replace, StringSplitOptions.RemoveEmptyEntries)[2]);
                Material SiO2_Material = new Material(lambda_SiO2, n_SiO2, k_SiO2);
                SiO2_Data.Add(SiO2_Material);
            }
            #region 이중 For 문 + 재귀함수 + 병렬 처리
            // Double ForLoop + Recursive Function + Parallel

            //int from = 0;
            //int to = 100;
            //int taskCount = 10;
            //Complex N_0 = new Complex(1, -0);
            //Complex theta_0 = new Complex(degToRad(65), 0);
            //Complex Sintheta_0 = Complex.Sin(theta_0);

            //Func<object, List<MSpectrum>> TFthicknessFunc = (objRange) =>
            //{
            //    int[] range = (int[])objRange;
            //    List<MSpectrum> multiCal = new List<MSpectrum>();
            //    for (int i = range[0]; i < range[1]; i++)
            //    {
            //        DateTime start = DateTime.Now;
            //        for (int n = 0; n < SiN_Data.Count; n++)
            //        {
            //            double n_SiO2_350 = SiO2_Data[n].RefractIdx;
            //            double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
            //            double n_SiN_350 = SiN_Data[n].RefractIdx;
            //            double k_SiN_350 = SiN_Data[n].ExtinctCoeff;

            //            Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
            //            Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

            //            Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
            //            Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

            //            Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
            //            Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
            //            Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
            //            Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
            //            Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

            //            Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
            //            Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
            //            Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);

            //            SMatrix_p = sMat_p(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_p, 1, 1000);
            //            SMatrix_s = sMat_s(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_s, 1, 1000);

            //            Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
            //            Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

            //            Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

            //            Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
            //            Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

            //            SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
            //            SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);


            //            double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            Console.WriteLine(i + " " + SiN_Data[n].Wavelength + "," + alpha + "," + beta);
            //            MSpectrum calSpectrum = new MSpectrum(SiN_Data[n].Wavelength, 65, alpha, beta);
            //            multiCal.Add(calSpectrum);
            //        }
            //        DateTime end = DateTime.Now;
            //        TimeSpan duration = end - start;
            //        Console.WriteLine(duration);
            //    }
            //    return multiCal;
            //};

            //Task<List<MSpectrum>>[] tasks = new Task<List<MSpectrum>>[taskCount];
            //int currentFrom = from;
            //int currentTo = to / tasks.Length;
            //for (int i = 0; i < tasks.Length; i++)
            //{
            //    Console.WriteLine("Task[{0}] : {1} ~ {2}", i, currentFrom, currentTo);

            //    tasks[i] = new Task<List<MSpectrum>>(TFthicknessFunc, new int[] { currentFrom, currentTo });
            //    currentFrom = currentTo + 1;
            //    if (i == tasks.Length - 2)
            //    {
            //        currentTo = to;
            //    }
            //    else
            //    {
            //        currentTo = currentTo + (to / tasks.Length);
            //    }
            //}
            //Console.ReadLine();
            //Console.WriteLine("Start!");
            //DateTime startTime = DateTime.Now;
            //Console.WriteLine(startTime);
            //foreach (Task<List<MSpectrum>> task in tasks)
            //    task.Start();

            //List<MSpectrum> total = new List<MSpectrum>();
            //foreach (Task<List<MSpectrum>> task in tasks)
            //{

            //    task.Wait();
            //    total.AddRange(task.Result.ToArray());

            //}

            //DateTime endTime = DateTime.Now;
            //TimeSpan elapsed = endTime - startTime;
            //Console.WriteLine();
            //Console.WriteLine("{0} {1}", elapsed, total.Count);
            //Console.ReadLine();
            //Console.Clear();
            //for (int i = 0; i < total.Count; i++)
            //{
            //    Console.WriteLine(i + " " + total[i].ToString());
            //}
            #endregion

            #region  삼중 For 문 + 병렬 처리
            // Triple ForLoop + Parallel
            //int from = 0;
            //int to = 100;
            //int taskCount = 10;
            //Complex N_0 = new Complex(1, -0);
            //Complex theta_0 = new Complex(degToRad(65), 0);
            //Complex Sintheta_0 = Complex.Sin(theta_0);

            //Func<object, List<MSpectrum>> TFthicknessFunc = (objRange) =>
            //{
            //    int[] range = (int[])objRange;
            //    List<MSpectrum> multiCal = new List<MSpectrum>();
            //    for (int i = range[0]; i < range[1]; i++)
            //    {
            //        DateTime start = DateTime.Now;
            //        for (int n = 0; n < SiN_Data.Count; n++)
            //        {
            //            double n_SiO2_350 = SiO2_Data[n].RefractIdx;
            //            double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
            //            double n_SiN_350 = SiN_Data[n].RefractIdx;
            //            double k_SiN_350 = SiN_Data[n].ExtinctCoeff;



            //            Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
            //            Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

            //            Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
            //            Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

            //            Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
            //            Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
            //            Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
            //            Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
            //            Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

            //            Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
            //            Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
            //            Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);



            //            for (int j = 1; j <= 1000; j++)
            //            {
            //                if (j % 2 == 0)
            //                {
            //                    Complex r_ij_p = refCoeff_p(N_1, N_2, theta_1, theta_2);
            //                    Complex r_ij_s = refCoeff_s(N_1, N_2, theta_1, theta_2);
            //                    Complex t_ij_p = transCoeff_p(N_1, N_2, theta_1, theta_2);
            //                    Complex t_ij_s = transCoeff_s(N_1, N_2, theta_1, theta_2);
            //                    Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
            //                    Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
            //                    Complex[,] layer_j = Layer_j(theta_2, 20, SiN_Data[n].Wavelength, N_2);
            //                    SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
            //                    SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
            //                }
            //                else
            //                {
            //                    Complex r_ij_p = refCoeff_p(N_2, N_1, theta_2, theta_3);
            //                    Complex r_ij_s = refCoeff_s(N_2, N_1, theta_2, theta_3);
            //                    Complex t_ij_p = transCoeff_p(N_2, N_1, theta_2, theta_3);
            //                    Complex t_ij_s = transCoeff_s(N_2, N_1, theta_2, theta_3);
            //                    Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
            //                    Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
            //                    Complex[,] layer_j = Layer_j(theta_3, 30, SiO2_Data[n].Wavelength, N_1);
            //                    SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
            //                    SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
            //                }
            //            }

            //            Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
            //            Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

            //            Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

            //            Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
            //            Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

            //            SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
            //            SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);


            //            double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            //Console.WriteLine(i + " " + SiN_Data[n].Wavelength + "," + alpha + "," + beta);
            //            MSpectrum calSpectrum = new MSpectrum(SiN_Data[n].Wavelength, 65, alpha, beta);
            //            multiCal.Add(calSpectrum);
            //        }
            //        DateTime end = DateTime.Now;
            //        TimeSpan duration = end - start;
            //        Console.WriteLine(duration);
            //    }
            //    return multiCal;
            //};

            //Task<List<MSpectrum>>[] tasks = new Task<List<MSpectrum>>[taskCount];
            //int currentFrom = from;
            //int currentTo = to / tasks.Length;
            //for (int i = 0; i < tasks.Length; i++)
            //{
            //    Console.WriteLine("Task[{0}] : {1} ~ {2}", i, currentFrom, currentTo);

            //    tasks[i] = new Task<List<MSpectrum>>(TFthicknessFunc, new int[] { currentFrom, currentTo });
            //    currentFrom = currentTo + 1;
            //    if (i == tasks.Length)
            //    {
            //        currentTo = to;
            //    }
            //    else
            //    {
            //        currentTo = currentTo + (to / tasks.Length);
            //    }
            //}
            //Console.ReadLine();
            //Console.WriteLine("Start!");
            //DateTime startTime = DateTime.Now;
            //Console.WriteLine(startTime);
            //foreach (Task<List<MSpectrum>> task in tasks)
            //    task.Start();

            //List<MSpectrum> total = new List<MSpectrum>();
            //foreach (Task<List<MSpectrum>> task in tasks)
            //{

            //    task.Wait();
            //    total.AddRange(task.Result.ToArray());

            //}

            //DateTime endTime = DateTime.Now;
            //TimeSpan elapsed = endTime - startTime;
            //Console.WriteLine();
            //Console.WriteLine("{0} {1}", elapsed, total.Count);
            //Console.ReadLine();
            //Console.Clear();
            //for (int i = 0; i < total.Count; i++)
            //{
            //    Console.WriteLine(i + " " + total[i].ToString());
            //}
            #endregion

            #region 삼중 For 문
            /// Triple ForLoop
            //DateTime startTime = DateTime.Now;
            //for (int cnt = 0; cnt < 100; cnt++)
            //{

            //    for (int n = 0; n < SiN_Data.Count; n++)
            //    {
            //        double n_SiO2_350 = SiO2_Data[n].RefractIdx;
            //        double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
            //        double n_SiN_350 = SiN_Data[n].RefractIdx;
            //        double k_SiN_350 = SiN_Data[n].ExtinctCoeff;

            //        Complex N_0 = new Complex(1, -0);
            //        Complex theta_0 = new Complex(degToRad(65), 0);
            //        Complex Sintheta_0 = Complex.Sin(theta_0);

            //        Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
            //        Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

            //        Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
            //        Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
            //        Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
            //        Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

            //        Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
            //        Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
            //        Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
            //        Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
            //        Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

            //        Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
            //        Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
            //        Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);

            //        Complex theta_4 = Complex.Asin((N_1 * Complex.Sin(theta_3)) / N_2);
            //        Complex theta_5 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);
            //        //Console.WriteLine(theta_0 + " " + theta_1 + " " + theta_2 + " " + theta_3 + " " + theta_4 + " " + theta_5);

            //        for (int i = 1; i <= 1000; i++)
            //        {
            //            if (i % 2 == 0)
            //            {
            //                Complex r_ij_p = refCoeff_p(N_1, N_2, theta_1, theta_2);
            //                Complex r_ij_s = refCoeff_s(N_1, N_2, theta_1, theta_2);
            //                Complex t_ij_p = transCoeff_p(N_1, N_2, theta_1, theta_2);
            //                Complex t_ij_s = transCoeff_s(N_1, N_2, theta_1, theta_2);
            //                Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
            //                Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
            //                Complex[,] layer_j = Layer_j(theta_2, 20, SiN_Data[n].Wavelength, N_2);
            //                SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
            //                SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
            //            }
            //            else
            //            {
            //                Complex r_ij_p = refCoeff_p(N_2, N_1, theta_2, theta_3);
            //                Complex r_ij_s = refCoeff_s(N_2, N_1, theta_2, theta_3);
            //                Complex t_ij_p = transCoeff_p(N_2, N_1, theta_2, theta_3);
            //                Complex t_ij_s = transCoeff_s(N_2, N_1, theta_2, theta_3);
            //                Complex[,] interface_ij_p = Interface_ij(r_ij_p, t_ij_p);
            //                Complex[,] interface_ij_s = Interface_ij(r_ij_s, t_ij_s);
            //                Complex[,] layer_j = Layer_j(theta_3, 30, SiO2_Data[n].Wavelength, N_1);
            //                SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p, layer_j));
            //                SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s, layer_j));
            //            }
            //        }

            //        Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
            //        Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

            //        Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
            //        Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
            //        Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
            //        Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

            //        Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
            //        Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

            //        SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
            //        SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);

            //        //for (int i = 0; i < 2; i++)
            //        //{
            //        //    Console.Write("|");
            //        //    for (int j = 0; j < 2; j++)
            //        //    {
            //        //        Console.Write(" " + SMatrix_p[i, j] + " ");
            //        //    }
            //        //    Console.WriteLine("|");
            //        //}
            //        //Console.WriteLine();
            //        double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //        double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //        Console.WriteLine(cnt + "," + SiN_Data[n].Wavelength + "," + alpha + "," + beta);

            //    }

            //}
            //DateTime endTime = DateTime.Now;
            //TimeSpan elapsed = endTime - startTime;
            //Console.WriteLine(elapsed);
            #endregion

            #region 이중 For 문 + 재귀 함수
            // DoubleLoop + Recursive Function
            //DateTime start = DateTime.Now;
            //for (int cnt = 0; cnt < 100; cnt++)
            //{
            //    DateTime startTime = DateTime.Now;
            //    for (int n = 0; n < SiN_Data.Count; n++)
            //    {
            //        double n_SiO2_350 = SiO2_Data[n].RefractIdx;
            //        double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
            //        double n_SiN_350 = SiN_Data[n].RefractIdx;
            //        double k_SiN_350 = SiN_Data[n].ExtinctCoeff;

            //        Complex N_0 = new Complex(1, -0);
            //        Complex theta_0 = new Complex(degToRad(65), 0);
            //        Complex Sintheta_0 = Complex.Sin(theta_0);

            //        Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
            //        Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

            //        Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
            //        Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
            //        Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
            //        Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

            //        Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
            //        Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
            //        Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
            //        Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
            //        Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

            //        Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
            //        Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
            //        Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);

            //        Complex theta_4 = Complex.Asin((N_1 * Complex.Sin(theta_3)) / N_2);
            //        Complex theta_5 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);
            //        //Console.WriteLine(theta_0 + " " + theta_1 + " " + theta_2 + " " + theta_3 + " " + theta_4 + " " + theta_5);

            //        SMatrix_p = sMat_p(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_p, 1, 1000);
            //        SMatrix_s = sMat_s(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_s, 1, 1000);

            //        Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
            //        Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

            //        Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
            //        Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
            //        Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
            //        Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

            //        Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
            //        Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

            //        SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
            //        SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);

            //        //for (int i = 0; i < 2; i++)
            //        //{
            //        //    Console.Write("|");
            //        //    for (int j = 0; j < 2; j++)
            //        //    {
            //        //        Console.Write(" " + SMatrix_p[i, j] + " ");
            //        //    }
            //        //    Console.WriteLine("|");
            //        //}
            //        //Console.WriteLine();
            //        double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //        double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //        //Console.WriteLine(cnt + "," + SiN_Data[n].Wavelength + "," + alpha + "," + beta);

            //    }
            //    DateTime endTime = DateTime.Now;
            //    TimeSpan elapsed = endTime - startTime;
            //    Console.WriteLine(elapsed);
            //}
            //DateTime end = DateTime.Now;
            //TimeSpan duration = end - start;
            //Console.WriteLine(duration);
            #endregion

            #region 삼중 For 문 + 병렬 처리(반복 횟수 감축)
            //int from = 0;
            //int to = 100;
            //int taskCount = 5;
            //Complex N_0 = new Complex(1, -0);
            //Complex theta_0 = new Complex(degToRad(65), 0);
            //Complex Sintheta_0 = Complex.Sin(theta_0);

            //Func<object, List<MSpectrum>> TFthicknessFunc = (objRange) =>
            //{
            //    int[] range = (int[])objRange;
            //    List<MSpectrum> multiCal = new List<MSpectrum>();
            //    for (int i = range[0]; i < range[1]; i++)
            //    {
            //        DateTime start = DateTime.Now;
            //        for (int n = 0; n < SiN_Data.Count; n++)
            //        {
            //            double n_SiO2_350 = SiO2_Data[n].RefractIdx;
            //            double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
            //            double n_SiN_350 = SiN_Data[n].RefractIdx;
            //            double k_SiN_350 = SiN_Data[n].ExtinctCoeff;



            //            Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
            //            Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

            //            Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
            //            Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
            //            Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

            //            Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
            //            Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
            //            Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
            //            Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
            //            Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

            //            Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
            //            Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
            //            Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);



            //            for (int j = 1; j <= 500; j++)
            //            {

            //                Complex r_ij_p_SIN = refCoeff_p(N_1, N_2, theta_1, theta_2);
            //                Complex r_ij_s_SIN = refCoeff_s(N_1, N_2, theta_1, theta_2);
            //                Complex t_ij_p_SIN = transCoeff_p(N_1, N_2, theta_1, theta_2);
            //                Complex t_ij_s_SIN = transCoeff_s(N_1, N_2, theta_1, theta_2);
            //                Complex[,] interface_ij_p_SIN = Interface_ij(r_ij_p_SIN, t_ij_p_SIN);
            //                Complex[,] interface_ij_s_SIN = Interface_ij(r_ij_s_SIN, t_ij_s_SIN);
            //                Complex[,] layer_j_SIN = Layer_j(theta_2, 20, SiN_Data[n].Wavelength, N_2);
            //                SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p_SIN, layer_j_SIN));
            //                SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s_SIN, layer_j_SIN));


            //                Complex r_ij_p_SIO2 = refCoeff_p(N_2, N_1, theta_2, theta_3);
            //                Complex r_ij_s_SIO2 = refCoeff_s(N_2, N_1, theta_2, theta_3);
            //                Complex t_ij_p_SIO2 = transCoeff_p(N_2, N_1, theta_2, theta_3);
            //                Complex t_ij_s_SIO2 = transCoeff_s(N_2, N_1, theta_2, theta_3);
            //                Complex[,] interface_ij_p_SI02 = Interface_ij(r_ij_p_SIO2, t_ij_p_SIO2);
            //                Complex[,] interface_ij_s_SIO2 = Interface_ij(r_ij_s_SIO2, t_ij_s_SIO2);
            //                Complex[,] layer_j_SIO2 = Layer_j(theta_3, 30, SiO2_Data[n].Wavelength, N_1);
            //                SMatrix_p = ComplexMatrixMultiply(SMatrix_p, ComplexMatrixMultiply(interface_ij_p_SI02, layer_j_SIO2));
            //                SMatrix_s = ComplexMatrixMultiply(SMatrix_s, ComplexMatrixMultiply(interface_ij_s_SIO2, layer_j_SIO2));

            //            }

            //            Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
            //            Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

            //            Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
            //            Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

            //            Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
            //            Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

            //            SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
            //            SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);


            //            double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
            //            //Console.WriteLine(i + " " + SiN_Data[n].Wavelength + "," + alpha + "," + beta);
            //            MSpectrum calSpectrum = new MSpectrum(SiN_Data[n].Wavelength, 65, alpha, beta);
            //            multiCal.Add(calSpectrum);
            //        }
            //        DateTime end = DateTime.Now;
            //        TimeSpan duration = end - start;
            //        Console.WriteLine(duration);
            //    }
            //    return multiCal;
            //};

            //Task<List<MSpectrum>>[] tasks = new Task<List<MSpectrum>>[taskCount];
            //int currentFrom = from;
            //int currentTo = to / tasks.Length;
            //for (int i = 0; i < tasks.Length; i++)
            //{
            //    Console.WriteLine("Task[{0}] : {1} ~ {2}", i, currentFrom, currentTo);

            //    tasks[i] = new Task<List<MSpectrum>>(TFthicknessFunc, new int[] { currentFrom, currentTo });
            //    currentFrom = currentTo + 1;
            //    if (i == tasks.Length)
            //    {
            //        currentTo = to;
            //    }
            //    else
            //    {
            //        currentTo = currentTo + (to / tasks.Length);
            //    }
            //}
            //Console.ReadLine();
            //Console.WriteLine("Start!");
            //DateTime startTime = DateTime.Now;
            //Console.WriteLine(startTime);
            //foreach (Task<List<MSpectrum>> task in tasks)
            //    task.Start();

            //List<MSpectrum> total = new List<MSpectrum>();
            //foreach (Task<List<MSpectrum>> task in tasks)
            //{

            //    task.Wait();
            //    total.AddRange(task.Result.ToArray());

            //}

            //DateTime endTime = DateTime.Now;
            //TimeSpan elapsed = endTime - startTime;
            //Console.WriteLine();
            //Console.WriteLine("{0} {1}", elapsed, total.Count);
            //Console.ReadLine();
            //Console.Clear();
            //for (int i = 0; i < total.Count; i++)
            //{
            //    Console.WriteLine(i + " " + total[i].ToString());
            //}
            #endregion

            #region 이중 For 문 + 재귀 함수 + 병렬 처리 (반복 횟수 감축)
            int from = 0;
            int to = 100;
            int taskCount = 10;
            Complex N_0 = new Complex(1, -0);
            Complex theta_0 = new Complex(degToRad(65), 0);
            Complex Sintheta_0 = Complex.Sin(theta_0);

            Func<object, List<MSpectrum>> TFthicknessFunc = (objRange) =>
            {
                int[] range = (int[])objRange;
                List<MSpectrum> multiCal = new List<MSpectrum>();
                for (int i = range[0]; i < range[1]; i++)
                {
                    DateTime start = DateTime.Now;
                    for (int n = 0; n < SiN_Data.Count; n++)
                    {
                        double n_SiO2_350 = SiO2_Data[n].RefractIdx;
                        double k_SiO2_350 = SiO2_Data[n].ExtinctCoeff;
                        double n_SiN_350 = SiN_Data[n].RefractIdx;
                        double k_SiN_350 = SiN_Data[n].ExtinctCoeff;

                        Complex N_1 = new Complex(n_SiO2_350, -k_SiO2_350);
                        Complex theta_1 = Complex.Asin((N_0 * Sintheta_0) / N_1);

                        Complex r01_p = refCoeff_p(N_0, N_1, theta_0, theta_1);
                        Complex t01_p = transCoeff_p(N_0, N_1, theta_0, theta_1);
                        Complex r01_s = refCoeff_s(N_0, N_1, theta_0, theta_1);
                        Complex t01_s = transCoeff_s(N_0, N_1, theta_0, theta_1);

                        Complex[,] Interface_01_p = Interface_ij(r01_p, t01_p);
                        Complex[,] Layer_1 = Layer_j(theta_1, 30, SiO2_Data[n].Wavelength, N_1);
                        Complex[,] SMatrix_p = ComplexMatrixMultiply(Interface_01_p, Layer_1);
                        Complex[,] Interface_01_s = Interface_ij(r01_s, t01_s);
                        Complex[,] SMatrix_s = ComplexMatrixMultiply(Interface_01_s, Layer_1);

                        Complex N_2 = new Complex(n_SiN_350, -k_SiN_350);
                        Complex theta_2 = Complex.Asin((N_1 * Complex.Sin(theta_1)) / N_2);
                        Complex theta_3 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_1);

                        SMatrix_p = sMat_p_short(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_p, 1, 500);
                        SMatrix_s = sMat_s_short(N_1, N_2, theta_1, theta_2, theta_3, 20, 30, SiN_Data[n].Wavelength, SiO2_Data[n].Wavelength, SMatrix_s, 1, 500);

                        Complex N_201 = new Complex(Si_Data[n].RefractIdx, -Si_Data[n].ExtinctCoeff);
                        Complex theta201 = Complex.Asin((N_2 * Complex.Sin(theta_2)) / N_201);

                        Complex r_200_201_p = refCoeff_p(N_2, N_201, theta_2, theta201);
                        Complex r_200_201_s = refCoeff_s(N_2, N_201, theta_2, theta201);
                        Complex t_200_201_p = transCoeff_p(N_2, N_201, theta_2, theta201);
                        Complex t_200_201_s = transCoeff_s(N_2, N_201, theta_2, theta201);

                        Complex[,] Interface_200_201_p = Interface_ij(r_200_201_p, t_200_201_p);
                        Complex[,] Interface_200_201_s = Interface_ij(r_200_201_s, t_200_201_s);

                        SMatrix_p = ComplexMatrixMultiply(SMatrix_p, Interface_200_201_p);
                        SMatrix_s = ComplexMatrixMultiply(SMatrix_s, Interface_200_201_s);


                        double alpha = calcAlpha(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
                        double beta = calcBeta(SMatrix_s[1, 0] / SMatrix_s[0, 0], SMatrix_p[1, 0] / SMatrix_p[0, 0], 45);
                        //Console.WriteLine(i + " " + SiN_Data[n].Wavelength + "," + alpha + "," + beta);
                        MSpectrum calSpectrum = new MSpectrum(SiN_Data[n].Wavelength, 65, alpha, beta);
                        multiCal.Add(calSpectrum);
                    }
                    DateTime end = DateTime.Now;
                    TimeSpan duration = end - start;
                    Console.WriteLine(duration);
                }
                return multiCal;
            };

            Task<List<MSpectrum>>[] tasks = new Task<List<MSpectrum>>[taskCount];
            int currentFrom = from;
            int currentTo = to / tasks.Length;
            for (int i = 0; i < tasks.Length; i++)
            {
                Console.WriteLine("Task[{0}] : {1} ~ {2}", i, currentFrom, currentTo);

                tasks[i] = new Task<List<MSpectrum>>(TFthicknessFunc, new int[] { currentFrom, currentTo });
                currentFrom = currentTo + 1;
                if (i == tasks.Length - 2)
                {
                    currentTo = to;
                }
                else
                {
                    currentTo = currentTo + (to / tasks.Length);
                }
            }
            Console.ReadLine();
            Console.WriteLine("Start!");
            DateTime startTime = DateTime.Now;
            Console.WriteLine(startTime);
            foreach (Task<List<MSpectrum>> task in tasks)
                task.Start();

            List<MSpectrum> total = new List<MSpectrum>();
            foreach (Task<List<MSpectrum>> task in tasks)
            {

                task.Wait();
                total.AddRange(task.Result.ToArray());

            }

            DateTime endTime = DateTime.Now;
            TimeSpan elapsed = endTime - startTime;
            Console.WriteLine();
            Console.WriteLine("{0} {1}", elapsed, total.Count);
            Console.ReadLine();
            Console.Clear();
            for (int i = 0; i < total.Count; i++)
            {
                Console.WriteLine(i + " " + total[i].ToString());
            }
            #endregion

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
        static Complex[,] Interface_ij(Complex r_ij, Complex t_ij) {
            Complex[,] Interface_01 = { { 1 / t_ij, r_ij / t_ij }, { r_ij / t_ij, 1 / t_ij } };
            return Interface_01;
        }
        static Complex[,] Layer_j(Complex Angle, double Thickness, double Wavelength, Complex N) {
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
        static Complex[,] sMat_p(Complex N_1, Complex N_2,  Complex theta_1, Complex theta_2, Complex theta_3, double d1, double d2, double lambda1, double lambda2, Complex[,] s_p,  int i, int n)
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
        static Complex[,] sMat_p_short(Complex N_1, Complex N_2, Complex theta_1, Complex theta_2, Complex theta_3, double d1, double d2, double lambda1, double lambda2, Complex[,] s_p, int i, int n)
        {
            if (i != n)
            {
                
                Complex r_ij_p_SiN = refCoeff_p(N_1, N_2, theta_1, theta_2);
                Complex t_ij_p_SiN = transCoeff_p(N_1, N_2, theta_1, theta_2);
                Complex[,] interface_ij_p_SiN = Interface_ij(r_ij_p_SiN, t_ij_p_SiN);
                Complex[,] layer_j_SiN = Layer_j(theta_2, d1, lambda1, N_2);
                s_p = ComplexMatrixMultiply(s_p, ComplexMatrixMultiply(interface_ij_p_SiN, layer_j_SiN));
                
                Complex r_ij_p_SiO2 = refCoeff_p(N_2, N_1, theta_2, theta_3);
                Complex t_ij_p_SiO2 = transCoeff_p(N_2, N_1, theta_2, theta_3);
                Complex[,] interface_ij_p_SiO2 = Interface_ij(r_ij_p_SiO2, t_ij_p_SiO2);
                Complex[,] layer_j_SiO2 = Layer_j(theta_3, d2, lambda2, N_1);
                s_p = ComplexMatrixMultiply(s_p, ComplexMatrixMultiply(interface_ij_p_SiO2, layer_j_SiO2));
                
                return ComplexMatrixMultiply(s_p, sMat_p(N_1, N_2, theta_1, theta_2, theta_3, d1, d2, lambda1, lambda2, s_p, i + 1, n));
            }
            else
            {
                return s_p;
            }
        }
        static Complex[,] sMat_s_short(Complex N_1, Complex N_2, Complex theta_1, Complex theta_2, Complex theta_3, double d1, double d2, double lambda1, double lambda2, Complex[,] s_s, int i, int n)
        {
            if (i != n)
            {
                
                Complex r_ij_s_SiN = refCoeff_s(N_1, N_2, theta_1, theta_2);
                Complex t_ij_s_SiN = transCoeff_s(N_1, N_2, theta_1, theta_2);
                Complex[,] interface_ij_s_SiN = Interface_ij(r_ij_s_SiN, t_ij_s_SiN);
                Complex[,] layer_j_SiN = Layer_j(theta_2, d1, lambda1, N_2);
                s_s = ComplexMatrixMultiply(s_s, ComplexMatrixMultiply(interface_ij_s_SiN, layer_j_SiN));
                                
                Complex r_ij_s_SiO2 = refCoeff_s(N_2, N_1, theta_2, theta_3);
                Complex t_ij_s_SiO2 = transCoeff_s(N_2, N_1, theta_2, theta_3);
                Complex[,] interface_ij_s_SiO2 = Interface_ij(r_ij_s_SiO2, t_ij_s_SiO2);
                Complex[,] layer_j_SiO2 = Layer_j(theta_3, d2, lambda2, N_1);
                s_s = ComplexMatrixMultiply(s_s, ComplexMatrixMultiply(interface_ij_s_SiO2, layer_j_SiO2));
                
                return ComplexMatrixMultiply(s_s, sMat_p(N_1, N_2, theta_1, theta_2, theta_3, d1, d2, lambda1, lambda2, s_s, i + 1, n));
            }
            else
            {
                return s_s;
            }
        }
    }
}
