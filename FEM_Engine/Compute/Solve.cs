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
        public static List<FEMResult> Solve(List<Bar> bars, List<Node> constraints, List<PointLoad> loads, int maxIter, double tol)
        {
            int nEL = bars.Count;

            // Create Edof and get unique Dofs
            Matrix<double> edof;
            Vector<double> bc;
            Vector<double> f_ext;

            Query.Edof_BC_fext(out edof, out bc, out f_ext, bars, constraints, loads);

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

            double conv = 1;
            //double tol = 0.000000001;

            while (conv > tol)
            {
                count = count + 1;

                // Allocate space for global stiffness matrix
                Matrix<double> K = DenseMatrix.Create(uniqueDofs.Count, uniqueDofs.Count, 0);

                // Allocate space for global force vetor f_int
                Vector<double> f_int = DenseVector.Create(uniqueDofs.Count,0);

                for (int i = 0; i < nEL; i++)
                {
                    Matrix<double> Ke = GreenBarStiffnessMatrix(bars[i], es.At(i), ed.Row(i));
                    Vector<double> fe = GreenBarForceVector(bars[i], es.At(i), ed.Row(i));

                    double[] edofRow = dofs.Row(i).ToArray();

                    // Assebmly Stiffness matrix and internal force vector
                    Query.Assem(ref K, Ke, ref f_int, fe, edofRow);

                }

                // Partitioning
                Matrix<double> K_part;
                Vector<double> f_part;
                Vector<double> f_ext_part;

                Query.Partioning(K, f_int, f_ext, freedof, out K_part, out f_part, out f_ext_part);


                conv = (Math.Pow(f_part.Norm(2), 2) / (1 + Math.Pow(f_ext_part.Norm(2), 2)));

                // delta_a = inv(K(freedof,freedof))*f(freedof)
                Vector<double> delta_u = K_part.Solve(f_part);
                //Vector<double> delta_u = K_part.Inverse() * f_part;

                // a(freedof) = a(freedof) + delta_a;
                for (int j = 0; j < freedof.Count; j++)
                {
                    int a = (int)freedof[j] - 1;
                    u[a] = u[a] + delta_u[j];
                }

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
                    BH.oM.Structure.Elements.Bar aBar = (BH.oM.Structure.Elements.Bar)bars[i];
                    es[i] = GreenBarElementStress(aBar, ed.Row(i));
                }

                if (count > maxIter)
                {
                    break;
                }
            }

            List<FEMResult> outResults = new List<FEMResult>();

            for (int i = 0; i < nEL; i++)
            {
                FEMResult outResult1 = new FEMResult();

                BarDeformation barDef1 = new BarDeformation();
                BarDisplacement barDisp1 = new BarDisplacement();
                BarForce barForce1 = new BarForce();
                BarStrain barStrain1 = new BarStrain();
                BarStress barStress1 = new BarStress();

                barDisp1.UX = ed[i, 0];
                barDisp1.UY = ed[i, 1];
                barDisp1.UZ = ed[i, 2];
                barDisp1.ObjectId = bars[i].Name;
                barDisp1.ResultCase = loads[0].Loadcase.Number;
                barDisp1.Position = 0;
                barDisp1.Divisions = 2;

                barForce1.FX = es[i];
                barForce1.ObjectId = bars[i].Name;
                barForce1.ResultCase = loads[0].Loadcase.Number;
                barForce1.Position = 0;
                barForce1.Divisions = 2;

                outResult1.barDisplacement = barDisp1;
                outResult1.barForce = barForce1;

                FEMResult outResult2 = new FEMResult();

                BarDeformation barDef2 = new BarDeformation();
                BarDisplacement barDisp2 = new BarDisplacement();
                BarForce barForce2 = new BarForce();
                BarStrain barStrain2 = new BarStrain();
                BarStress barStress2 = new BarStress();

                barDisp2.UX = ed[i, 3];
                barDisp2.UY = ed[i, 4];
                barDisp2.UZ = ed[i, 5];
                barDisp2.ObjectId = bars[i].Name;
                barDisp2.ResultCase = loads[0].Loadcase.Number;
                barDisp2.Position = 1;
                barDisp2.Divisions = 2;

                barForce2.FX = es[i];
                barForce2.ObjectId = bars[i].Name;
                barForce2.ResultCase = loads[0].Loadcase.Number;
                barForce2.Position = 1;
                barForce2.Divisions = 2;

                outResult2.barDisplacement = barDisp2;
                outResult2.barForce = barForce2;
               
                outResults.Add(outResult1);
                outResults.Add(outResult2);
            }

            //return u;
            return outResults;          
        }
    }
}
