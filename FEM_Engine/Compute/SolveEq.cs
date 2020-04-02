using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.FEM;
using BH.oM.Structure.Elements;
using BH.oM.Structure.Loads;
using BH.oM.Structure.Results;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BH.Engine.FEM
{
    public static partial class Compute
    {
        public static void SolveEq(Matrix<double> K, Vector<double> f_int, Vector<double> f_ext_current, ref Vector<double> u, List<double> freedof, out double conv)
        {
            // Partitioning
            Matrix<double> K_part;
            Vector<double> f_part;
            Vector<double> f_ext_part;

            Query.Partioning(K, f_int, f_ext_current, freedof, out K_part, out f_part, out f_ext_part);


            conv = (Math.Pow(f_part.Norm(2), 2) / (1 + Math.Pow(f_ext_part.Norm(2), 2)));

            Vector<double> delta_u = K_part.Solve(f_part);

            // a(freedof) = a(freedof) + delta_a;
            for (int j = 0; j < freedof.Count; j++)
            {
                int a = (int)freedof[j] - 1;
                u[a] = u[a] + delta_u[j];
            }
        }
    }
}
