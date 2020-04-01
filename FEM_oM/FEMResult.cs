using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Base;
using BH.oM.Structure.Results;

namespace BH.oM.FEM
{
    public class FEMResult : BHoMObject
    {
        public BarDeformation barDeformation { get; set; }
        public BarDisplacement barDisplacement { get; set; }
        public BarForce barForce { get; set; }
        public BarStrain barStrain { get; set; }
        public BarStress barStress { get; set; }      
    }
}
