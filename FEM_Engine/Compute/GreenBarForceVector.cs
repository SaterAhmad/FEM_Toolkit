using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Structure.Elements;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BH.Engine.FEM
{
    public static partial class Compute
    {
        public static Matrix<double> GreenBarForceVector(Bar aBar, double N, Vector<double> ed)
        {
            double[,] ecArray = { { aBar.StartNode.Position.X, aBar.EndNode.Position.X }, { aBar.StartNode.Position.Y, aBar.EndNode.Position.Y }, { aBar.StartNode.Position.Z, aBar.EndNode.Position.Z } };
            Matrix<double> ec = DenseMatrix.OfArray(ecArray);

            Vector<double> xo = ec.Column(1) - ec.Column(0);

            Vector<double> edStart = ed.SubVector(0, 3);
            Vector<double> edEnd = ed.SubVector(3, 3);

            Vector<double> u = edEnd - edStart;
            Vector<double> x = xo + u;

            double deltaX = aBar.EndNode.Position.X - aBar.StartNode.Position.X;
            double deltaY = aBar.EndNode.Position.Y - aBar.StartNode.Position.Y;
            double deltaZ = aBar.EndNode.Position.Z - aBar.StartNode.Position.Z;

            double lo = Math.Sqrt(Math.Pow(deltaX, 2) + Math.Pow(deltaY, 2) + Math.Pow(deltaZ, 2));

            Matrix<double> ef = (N / lo) * (-1 * x).ToColumnMatrix().Stack(x.ToColumnMatrix());
            return ef;
        }
    }
}
