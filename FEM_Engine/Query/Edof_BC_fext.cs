using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Base;
using BH.oM.Geometry;
using BH.oM.Structure.Elements;
using BH.Engine.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace BH.Engine.FEM
{
    public static partial class Query
    {
        public static void EDOF_BC_fext(out Matrix<double> Edof, out Vector<double> Bc, out Vector<double> Fext, List<Bar> bars, List<Node> constraints, List<Node> loads)
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

            double tol = 1e-6;

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
            }

            Vector<double> nDof = DenseVector.Create(1, 0);
            Matrix<double> bc = DenseMatrix.Create(1, 1, 0);

            for (int i = 0; i < constraints.Count; i++)
            {
                Node aNode = constraints[i];

                int nodeIndex = cullDup.IndexOf(Geometry.Query.ClosestPoint(aNode.Position, cullDup));

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
            }

            Vector<double> f_ext = DenseVector.Create(cullDup.Count * 3, 0);

            for (int i = 0; i < loads.Count; i++)
            {

                BH.oM.Structure.Loads.PointLoad ptLoad = loads[i];
                List<BH.oM.Structure.Elements.Node> loadNodes = ptLoad.Objects.Elements;
                int a = loadNodes.Count;

                for (int j = 0; j < a; j++)
                {
                    BH.oM.Structure.Elements.Node aNode = loadNodes[j];
                    int nodeIndex = cull_dup.IndexOf(BH.Engine.Geometry.Query.ClosestPoint(aNode.Position, cull_dup));

                    f_ext[(int)node_dofs.Row(nodeIndex).At(0) - 1] = ptLoad.Force.X;
                    f_ext[(int)node_dofs.Row(nodeIndex).At(1) - 1] = ptLoad.Force.Y;
                    f_ext[(int)node_dofs.Row(nodeIndex).At(2) - 1] = ptLoad.Force.Z;
                }
            }

            Edof = edof;
        }
        static void getEDOFandBC(out Matrix<double> EDOF, out Vector<double> BC, out Vector<double> fext, List<Object> barList, List<Object> constraintNodes, List<object> loads)
        {
            int nel = barList.Count;
            int nNode = constraintNodes.Count;
            List<int> dofs = new List<int>();

            List<BH.oM.Geometry.Point> el_nodes = new List<BH.oM.Geometry.Point>();
            Matrix<double> edof = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(nel, 6, 0);
            Matrix<double> bc = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(1, 1, 0);


            for (int i = 0; i < nel; i++)
            {
                BH.oM.Structure.Elements.Bar aBar = (BH.oM.Structure.Elements.Bar)barList[i];

                el_nodes.Add(aBar.StartNode.Position);
                el_nodes.Add(aBar.EndNode.Position);
            }

            double tol = 1e-6;
            List<BH.oM.Geometry.Point> cull_dup = BH.Engine.Geometry.Compute.CullDuplicates(el_nodes, tol);

            MathNet.Numerics.LinearAlgebra.Matrix<double> node_dofs = MathNet.Numerics.LinearAlgebra.Double.DenseMatrix.Create(cull_dup.Count, 3, 0);

            for (int i = 0; i < cull_dup.Count; i++)
            {
                Vector<double> node_dof = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(3, 0);
                node_dof[0] = 3 * i + 1;
                node_dof[1] = 3 * i + 2;
                node_dof[2] = 3 * i + 3;

                node_dofs.SetRow(i, node_dof);
            }

            Vector<double> e = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(nel, 0);

            for (int i = 0; i < nel; i++)
            {
                BH.oM.Structure.Elements.Bar aBar = (BH.oM.Structure.Elements.Bar)barList[i];

                int ind_Start = cull_dup.IndexOf(BH.Engine.Geometry.Query.ClosestPoint(aBar.StartNode.Position, cull_dup));
                int ind_End = cull_dup.IndexOf(BH.Engine.Geometry.Query.ClosestPoint(aBar.EndNode.Position, cull_dup));
                dofs.Add(ind_Start);

                e[i] = i + 1;

                edof.SetRow(i, node_dofs.Row(ind_Start).ToColumnMatrix().Stack(node_dofs.Row(ind_End).ToColumnMatrix()).Column(0));
            }

            Vector<double> nDof = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(1, 0);

            for (int i = 0; i < constraintNodes.Count; i++)
            {
                BH.oM.Structure.Elements.Node aNode = (BH.oM.Structure.Elements.Node)constraintNodes[i];

                int nodeIndex = cull_dup.IndexOf(BH.Engine.Geometry.Query.ClosestPoint(aNode.Position, cull_dup));

                if (aNode.Support.TranslationX.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = node_dofs.Row(nodeIndex).At(0);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }
                if (aNode.Support.TranslationY.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = node_dofs.Row(nodeIndex).At(1);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }
                if (aNode.Support.TranslationZ.Equals(BH.oM.Structure.Constraints.DOFType.Fixed))
                {
                    nDof[0] = node_dofs.Row(nodeIndex).At(2);
                    int index = bc.Column(0).Count;
                    bc = bc.InsertRow(index, nDof);
                }
            }

            Vector<double> f_ext = MathNet.Numerics.LinearAlgebra.Double.DenseVector.Create(cull_dup.Count * 3, 0);

            for (int i = 0; i < loads.Count; i++)
            {

                BH.oM.Structure.Loads.PointLoad ptLoad = (BH.oM.Structure.Loads.PointLoad)loads[i];
                List<BH.oM.Structure.Elements.Node> loadNodes = ptLoad.Objects.Elements;
                int a = loadNodes.Count;

                for (int j = 0; j < a; j++)
                {
                    BH.oM.Structure.Elements.Node aNode = loadNodes[j];
                    int nodeIndex = cull_dup.IndexOf(BH.Engine.Geometry.Query.ClosestPoint(aNode.Position, cull_dup));

                    f_ext[(int)node_dofs.Row(nodeIndex).At(0) - 1] = ptLoad.Force.X;
                    f_ext[(int)node_dofs.Row(nodeIndex).At(1) - 1] = ptLoad.Force.Y;
                    f_ext[(int)node_dofs.Row(nodeIndex).At(2) - 1] = ptLoad.Force.Z;
                }
            }

            edof = edof.InsertColumn(0, e);
            EDOF = edof;
            BC = bc.RemoveRow(0).Column(0);
            fext = f_ext;
        }
    }
}
