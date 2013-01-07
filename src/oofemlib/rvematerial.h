/*
 * RVEMaterial.h
 *
 *  Created on: Mar 26, 2010
 *      Author: carl
 */

#ifndef rvematerial_h
#define rvematerial_h

#include <cstdio>
#include <cstdlib>

#include "material.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "oofem_limits.h"

namespace oofem {

/**
 * @author Carl Sandström
 */
class RVEMaterial //: virtual public Material
{
private:
    int stdoutFID;
    fpos_t stdoutPos;

protected:
    /// Name of .in file containing the RVE
    std::string rveFilename;
    std::string rveLogFilename;
    /// Type of boundary condition.
    int BCType;

public:

    EngngModel *rve;

    // Constructor
    RVEMaterial(int n, Domain *d) { };// : Material(n, d) { };

    // Destructor
    virtual ~RVEMaterial() { free (rve); };

    int SupressRVEoutput;

    IRResultType initializeFrom(InputRecord *ir);

    void suppressStdout();
    void enableStdout();

    const char *giveClassName() const { return "RVEMaterial"; };
    classType giveClassID() const { return MaterialClass; };
};

}

#endif // rvematerial_h
