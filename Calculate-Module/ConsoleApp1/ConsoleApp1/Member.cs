using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleApp1
{
    class Member
    {

    }

    class Si
    {
        public double eV, e1, e2;
    }

    class SiO2
    {
        public double angstroms, n, k;
    }

    class SiN
    {
        public double nm, n, k;
    }

    class SiO2nm_on_Si
    {
        public double nm, aoi, Psi, Delta;
    }

    class CalSpectrum
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
