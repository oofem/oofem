/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef rvematerial_h
#define rvematerial_h

#include <cstdio>
#include <cstdlib>

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "oofem_limits.h"

///@name Input fields for RVEMaterial
//@{
#define _IFT_RVEMaterial_bctype "bctype"
#define _IFT_RVEMaterial_supressoutput "supressoutput"
#define _IFT_RVEMaterial_fileName "file"
//@}

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
    ~RVEMaterial() { free (rve); };

    int SupressRVEoutput;

    IRResultType initializeFrom(InputRecord *ir);

    void suppressStdout();
    void enableStdout();

    const char *giveClassName() const { return "RVEMaterial"; };
    classType giveClassID() const { return RVEMaterialClass; };
};

}

#endif // rvematerial_h
