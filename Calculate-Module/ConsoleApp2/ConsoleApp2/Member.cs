using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp2
{
    class Member
    {
    }
    public class Si_new
    {
        public double nm = 0, n = 0, k = 0;
    }

    public class SiO2_new
    {
        public double nm = 0, n = 0, k = 0;
    }
    public class SiO2_1000nm_on_Si
    {
        public double nm = 0, aoi = 0, Psi = 0, Delta = 0;
    }
    public class SiO2_1000nm_on_Si_new
    {
        public double nm=0, aoi=0, alpha=0, beta=0;
    }
    public class CalSpectrum
    {
        public double wavelength;
        public double aoi;
        public double alpha;
        public double beta;

        public CalSpectrum(double wavelength, double aoi, double alpha, double beta)
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


}
