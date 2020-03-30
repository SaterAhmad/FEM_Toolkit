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
        public static void Assem(Matrix<double> Edof, Matrix<double> K, Matrix<double> Ke, Vector<double> f, Vector<double> fe, int elementNr)
        {
            int n = Edof.RowCount;
            int nie = Edof.ColumnCount;

            Matrix<double> t = Edof.RemoveColumn(0);

            for (int i = 0; i<nie; i++)
            {
                int a = (int)t[elementNr-1,i];
                
                for (int j = 0; j < n; j++)
                {
                    int b = (int)t[elementNr-1,j];

                    K[a,b] = K[a,b] + Ke[i, j];
                }

                f[a] = f[a] + fe[i];
            }
            
        }
    }