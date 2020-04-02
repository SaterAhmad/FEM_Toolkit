﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.FEM;
using BH.oM.Structure.Constraints;
using BH.oM.Structure.Elements;
using BH.oM.Structure.Loads;
using BH.Engine.FEM;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace BH.Engine.FEM
{
    public static partial class Compute
    {
        public static List<FEMResult> Run(List<Bar> bars, List<Node> constraints, List<PointLoad> loads, List<BarPrestressLoad> prestress, Settings settings)
        {
            int loadsteps = settings.loadsteps;
            int maxIter = settings.maxIter;
            double tol = settings.tol;

            int nEL = bars.Count;

            // Create Edof and get unique Dofs
            Matrix<double> edof;
            Vector<double> bc;
            Vector<double> f_ext;
            Vector<double> pres;

            Query.Edof_BC_fext(out edof, out bc, out f_ext, out pres, bars, constraints, loads, prestress);

            double[,] edofArray = edof.ToArray();

            List<double> dofList = edofArray.Cast<double>().ToList();
            List<double> uniqueDofs = dofList.Distinct().ToList();

            Matrix<double> dofs = edof.RemoveColumn(0);

            // Extract freedofs
            List<double> freedof = uniqueDofs.Except(bc).ToList();

            // Create Initial element displacements and forces
            Matrix<double> ed = DenseMatrix.Create(nEL, 6, 0);
            Vector<double> es = DenseVector.Create(nEL, 0);

            // Allocate space for global diplacement vector a
            Vector<double> u = DenseVector.Create(uniqueDofs.Count, 0);

            int count = 0;

            for (int step = 1; step <= loadsteps; step++)
            {
                Vector<double> f_ext_current = step / loadsteps * f_ext;
                Vector<double> pres_current = step / loadsteps * pres;
                double conv = 1;

                while (conv>tol)
                {
                    count = count + 1;

                    // add Initial prestress
                    es += pres_current;

                    // Allocate space for global stiffness matrix
                    Matrix<double> K = DenseMatrix.Create(uniqueDofs.Count, uniqueDofs.Count, 0);

                    // Allocate space for global force vetor f_int
                    Vector<double> f_int = DenseVector.Create(uniqueDofs.Count, 0);

                    for (int i = 0; i < nEL; i++)
                    {
                        Matrix<double> Ke = GreenBarStiffnessMatrix(bars[i], es.At(i), ed.Row(i));
                        Vector<double> fe = GreenBarForceVector(bars[i], es.At(i), ed.Row(i));

                        double[] edofRow = dofs.Row(i).ToArray();

                        // Assebmly Stiffness matrix and internal force vector
                        Query.Assem(ref K, Ke, ref f_int, fe, edofRow);

                    }

                    SolveEq(K, f_int, f_ext_current, ref u, freedof, out conv);

                    // ed = extract(edof,a)
                    for (int i = 0; i < nEL; i++)
                    {
                        int a = (int)edof.Row(i).At(0) - 1;
                        for (int j = 0; j < dofs.Row(i).Count; j++)
                        {
                            int b = (int)dofs.Row(i).At(j) - 1;
                            ed[a, j] = u[b];
                        }
                    }

                    // Extract Element Forces
                    for (int i = 0; i < nEL; i++)
                    {
                        Bar aBar = bars[i];
                        es[i] = GreenBarElementStress(aBar, ed.Row(i));
                    }

                    if (count > maxIter)
                    {
                        break;
                    }
                }
            }

            List<FEMResult> outResult = FEM.Query.Results(bars, loads, es, ed);
            return outResult;
        }

    }
}