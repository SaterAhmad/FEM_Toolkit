using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Structure.Elements;
using BH.oM.Structure.MaterialFragments;
using BH.oM.Physical;
using BH.Engine.Structure;
using MathNet.Numerics.LinearAlgebra;

namespace BH.Engine.FEM
{
    public static partial class Compute
    {
        public static Matrix<double> GreenBarStiffnessMatrix(Bar aBar, double N, Vector<double> ed, bool isActive)
        {
            double Ao;
            double E;

            if (isActive)
            {
                Ao = aBar.SectionProperty.Area;
                E = (aBar.SectionProperty.Material as IIsotropic).YoungsModulus;
            }
            else
            {
                Ao = 1;
                E = 1;
                N = 1;
            }
            

            // intitial bar coordinates
            double[,] _ec = { { aBar.StartNode.Position.X, aBar.EndNode.Position.X }, { aBar.StartNode.Position.Y, aBar.EndNode.Position.Y }, { aBar.StartNode.Position.Z, aBar.EndNode.Position.Z } };
            Matrix<double> ec = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.OfArray(_ec);

            // initial length
            Vector<double> xo = ec.Column(1) - ec.Column(0);

            // bar displacement 
            Vector<double> u = ed.SubVector(3, 3) - ed.SubVector(0, 3);
            Vector<double> x = xo + u;

            // initial length
            double lo = Structure.Query.Length(aBar);

            // direction ...
            Matrix<double> k = x.ToColumnMatrix() * x.ToRowMatrix();
            Matrix<double> k_ = k.Stack(-1 * k).Transpose();
            Matrix<double> k_final = k_.Stack(-1 * k_);

            // elemental stiffness
            Matrix<double> K1 = (E * Ao / Math.Pow(lo, 3)) * k_final;

            // geometric stiffness
            Matrix<double> _I_ = MathNet.Numerics.LinearAlgebra.Double.DiagonalMatrix.CreateIdentity(3);
            Matrix<double> _I = _I_.Stack(-1 * _I_).Transpose();
            Matrix<double> I = _I.Stack(-1 * _I);
            Matrix<double> K2 = (N / lo) * I;

            // return element stiffness matrix
            Matrix<double> Ke = K1 + K2;
            return Ke;
        }

    }
}
