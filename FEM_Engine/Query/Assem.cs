using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace BH.Engine.FEM
{
    public static partial class Query
    {
        public static void Assem(ref Matrix<double> K, Matrix<double> Ke, ref Vector<double> f_int, Vector<double> fe, double[] edofRow)
        {
            for (int i = 0; i < edofRow.Length; i++)
            {
                for (int j = 0; j < edofRow.Length; j++)
                {
                    // Stiffness matrix assembly
                    K[(int)edofRow[i]-1, (int)edofRow[j]-1] = K[(int)edofRow[i]-1, (int)edofRow[j]-1] + Ke[i, j];
                }
                // Force matrix assembly
                f_int[(int)edofRow[i]-1] = f_int[(int)edofRow[i]-1] + fe[i];
            }

        }
    }
}
