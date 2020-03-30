﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Structure.Elements;
using MathNet.Numerics.LinearAlgebra;
using BH.oM.Structure.MaterialFragments;
using BH.oM.Physical;

namespace BH.Engine.FEM.Compute
{
    public static partial class Compute
    {
        public static double GreenBarElementStress(Bar aBar, Vector<double> ed)
        {
            double E = (aBar.SectionProperty.Material as IIsotropic).YoungsModulus;
            double A = aBar.SectionProperty.Area;

            // initial coordinates
            double[,] _ec = { { aBar.StartNode.Position.X, aBar.EndNode.Position.X }, { aBar.StartNode.Position.Y, aBar.EndNode.Position.Y }, { aBar.StartNode.Position.Z, aBar.EndNode.Position.Z } };
            Matrix<double> ec = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfArray(_ec);

            // initial length
            Vector<double> xo = ec.Column(1) - ec.Column(0);

            // bar displacement
            Vector<double> u = ed.SubVector(3, 3) - ed.SubVector(0, 3);
            Vector<double> x = xo + u;

            // current length
            Vector<double> lo = (x.ToRowMatrix() * x).PointwiseSqrt();
            double l0 = lo.ElementAt(0);

            // element stress (green)
            Vector<double> ee = (1 / Math.Pow(l0, 2)) * (xo.ToColumnMatrix().Transpose() * u + u.ToColumnMatrix().Transpose() * u / 2);
            Vector<double> es = E * A * ee;

            return es.At(0);

        }
    }
}


    
        