using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Base;
using BH.oM.Geometry;
using BH.oM.Structure.Elements;
using BH.oM.Structure.Loads;
using BH.Engine.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BH.Engine.FEM
{
    public static partial class Query
    {
        public static void Edof_BC_fext(out Matrix<double> Edof, out Vector<double> BC, out Vector<double> Fext, out Vector<double> pres, out Matrix<double> bcSprings, List<Bar> bars, List<Node> constraints, List<PointLoad> loads, List<BarPrestressLoad> ploads)
        {

            int nEL = bars.Count;
            int nNode = constraints.Count;

            List<int> dofs = new List<int>();
            List<Point> elPoints = new List<Point>();

            for (int i = 0; i < nEL; i++)
            {
                Bar aBar = bars[i];
                elPoints.Add(aBar.StartNode.Position);
                elPoints.Add(aBar.EndNode.Position);
            }

            double tol = 1e-5;

            List<Point> cullDup = Geometry.Compute.CullDuplicates(elPoints, tol);

            Matrix<double> nodeDofs = DenseMatrix.Create(cullDup.Count, 3, 0);

            for (int i = 0; i < cullDup.Count; i++)
            {
                Vector<double> nodeDof = DenseVector.Create(3, 0);
                nodeDof[0] = 3 * i + 1;
                nodeDof[1] = 3 * i + 2;
                nodeDof[2] = 3 * i + 3;

                nodeDofs.SetRow(i, nodeDof);
            }


            // create edof
            Vector<double> e = DenseVector.Create(nEL, 0);
            Matrix<double> edof = DenseMatrix.Create(nEL, 6, 0);

            for (int i = 0; i < nEL; i++)
            {
                Bar aBar = bars[i];
                
                int indStart = cullDup.IndexOf(Geometry.Query.ClosestPoint(aBar.StartNode.Position, cullDup));
                int indEnd = cullDup.IndexOf(Geometry.Query.ClosestPoint(aBar.EndNode.Position, cullDup));
                dofs.Add(indStart);

                e[i] = i + 1;

                edof.SetRow(i, nodeDofs.Row(indStart).ToColumnMatrix().Stack(nodeDofs.Row(indEnd).ToColumnMatrix()).Column(0));

                //test of indexing
                string id = "FEMID: " + i.ToString();
                bars[i].Tags.Add(id);
                
            }

            Vector<double> nDof = DenseVector.Create(1, 0);
            Vector<double> temp = DenseVector.Create(2, 0);
            Matrix<double> bc = DenseMatrix.Create(1, 1, 0);
            Matrix<double> bc_stiffness = DenseMatrix.Create(1, 2, 0);

            for (int i = 0; i < constraints.Count; i++)
            {
                Node aNode = constraints[i];

                int nodeIndex = cullDup.IndexOf(Geometry.Query.ClosestPoint(aNode.Position, cullDup));

                // Fixed
                if (aNode.Support.TranslationX.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = nodeDofs.Row(nodeIndex).At(0);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }
                if (aNode.Support.TranslationY.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = nodeDofs.Row(nodeIndex).At(1);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }
                if (aNode.Support.TranslationZ.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = nodeDofs.Row(nodeIndex).At(2);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }

                // Spring
                if (aNode.Support.TranslationX.Equals(BH.oM.Structure.Constraints.DOFType.Spring))
                {
                    temp[0] = nodeDofs.Row(nodeIndex).At(0);
                    temp[1] = aNode.Support.TranslationalStiffnessX;
                    int index = bc_stiffness.Column(0).Count;
                    bc_stiffness = bc_stiffness.InsertRow(index, temp);
                }
                if (aNode.Support.TranslationY.Equals(BH.oM.Structure.Constraints.DOFType.Spring))
                {
                    temp[0] = nodeDofs.Row(nodeIndex).At(1);
                    temp[1] = aNode.Support.TranslationalStiffnessY;
                    int index = bc_stiffness.Column(0).Count;
                    bc_stiffness = bc_stiffness.InsertRow(index, temp);

                }
                if (aNode.Support.TranslationZ.Equals(BH.oM.Structure.Constraints.DOFType.Spring))
                {

                    temp[0] = nodeDofs.Row(nodeIndex).At(2);
                    temp[1] = aNode.Support.TranslationalStiffnessZ;
                    int index = bc_stiffness.Column(0).Count;
                    bc_stiffness = bc_stiffness.InsertRow(index, temp);

                }
            }

            Vector<double> fext = DenseVector.Create(cullDup.Count * 3, 0);

            for (int i = 0; i < loads.Count; i++)
            {

                PointLoad ptLoad = loads[i];
                List<Node> loadNodes = ptLoad.Objects.Elements;
                int a = loadNodes.Count;

                for (int j = 0; j < a; j++)
                {
                    BH.oM.Structure.Elements.Node aNode = loadNodes[j];
                    int nodeIndex = cullDup.IndexOf(Geometry.Query.ClosestPoint(aNode.Position, cullDup));

                    fext[(int)nodeDofs.Row(nodeIndex).At(0) - 1] = ptLoad.Force.X;
                    fext[(int)nodeDofs.Row(nodeIndex).At(1) - 1] = ptLoad.Force.Y;
                    fext[(int)nodeDofs.Row(nodeIndex).At(2) - 1] = ptLoad.Force.Z;
                }
            }



            pres = DenseVector.Create(edof.RowCount, 0);

            for (int i = 0; i < ploads.Count; i++)
            {

                double Load = ploads[i].Prestress;
                List <Bar> loadElements = ploads[i].Objects.Elements;

               for (int j = 0; j < loadElements.Count; j++)
                {
                    int barIndex = bars.IndexOf(loadElements[j]);
                    pres[barIndex] = Load;
                }

            }


            edof = edof.InsertColumn(0, e);
            Edof = edof;
            BC = bc.RemoveRow(0).Column(0);
            Fext = fext;
           

            if (bc_stiffness.RowCount > 1) 
            {
                bcSprings = bc_stiffness.RemoveRow(0); 
            }
            else
            {
                bcSprings = null;
            }
            
            
        }
    }
}
