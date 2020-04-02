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

namespace BH.Engine.FEM
{
    public static partial class Query
    {
        public static List<FEMResult> Results(List<Bar> bars, List<PointLoad> loads, Vector<double> es, Matrix<double> ed)
        {
            int nEL = bars.Count;
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
            return outResults;
        }
    }
}
