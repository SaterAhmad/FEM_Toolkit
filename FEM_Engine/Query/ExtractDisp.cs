using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace BH.Engine.FEM
{
    public static partial class Query
    {
        public static Matrix<double> ExtractDisp(Matrix<double> dof, Vector<double> u)
        {
            Matrix<double> ed = DenseMatrix.Create(dof.RowCount, 6, 0);

            for ( int i = 0; i < dof.RowCount; i++)
            {
                
                for ( int j = 0; j < dof.ColumnCount; j++)
                {
                    int ind = (int)dof.Row(i).At(j)-1;
                    
                    ed[i, j] = u[ind];
                }
            }
            return ed;
        }
    }
}
