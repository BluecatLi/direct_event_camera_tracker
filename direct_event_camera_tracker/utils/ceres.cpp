#include "ceres.h"

using namespace ceres;
using namespace Eigen;

////////////////////////////////////////////////////////////////////////////////

namespace CeresCaster
{
    // double

    double toDouble(const double& val) { return val; }
    float  toFloat(const double& val) { return (float) val; }
};

////////////////////////////////////////////////////////////////////////////////
