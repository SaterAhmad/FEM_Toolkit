using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BH.oM.Base;

namespace BH.oM.FEM
{
    public class Settings : BHoMObject
    {
        public int loadsteps { get; set; }
        public int maxIter { get; set; }
        public double tol { get; set; }        
    }
}
