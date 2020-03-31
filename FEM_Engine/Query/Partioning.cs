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
        public static void Partioning(Matrix<double> K, Vector<double> f_int, Vector<double> f_ext, List<double> freedof, out Matrix<double> K_part, out Vector<double> f_part, out Vector<double> f_ext_part)
        {

            Matrix<double> K_part_temp = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(freedof.Count, freedof.Count, 0);
            Vector<double> f_part_temp = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(freedof.Count, 0);
            Vector<double> f_ext_part_temp = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(freedof.Count, 0);

            for (int j = 0; j < freedof.Count; j++)
            {
                int a = (int)freedof[j] - 1;
                for (int k = 0; k < freedof.Count; k++)
                {
                    int b = (int)freedof[k] - 1;
                    K_part_temp[j, k] = K[a, b];
                }
                f_part_temp[j] = f_ext[a] - f_int[a];
                f_ext_part_temp[j] = f_ext[a];
            }

            K_part = K_part_temp;
            f_part = f_part_temp;
            f_ext_part = f_ext_part_temp;
        }
    }
}
