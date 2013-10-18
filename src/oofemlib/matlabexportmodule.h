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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef matlabexportmodule_h_
#define matlabexportmodule_h_

#include <vector>

#include "exportmodule.h"

///@name Input fields for MatlabExportModule
//@{
#define _IFT_MatlabExportModule_Name "matlab"
#define _IFT_MatlabExportModule_mesh "mesh"
#define _IFT_MatlabExportModule_data "data"
#define _IFT_MatlabExportModule_area "area"
#define _IFT_MatlabExportModule_specials "specials"
#define _IFT_MatlabExportModule_ReactionForces "reactionforces"
#define _IFT_MatlabExportModule_DofManList "dofmanlist"
//@}

namespace oofem {
/**
 * (Under development) The Matlab export module enables oofem to export the results to a textfile containing the description of the mesh used
 * along with the pertinent results.
 *
 * @author Carl Sandström
 */
class OOFEM_EXPORT MatlabExportModule : public ExportModule
{
protected:
    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    /// list of primary unknowns to export
    IntArray primaryVarsToExport;
    std :: string functionname;

    FILE *giveOutputStream(TimeStep *);
    std :: vector <double> smax;
    std :: vector <double> smin;

    double Area, Volume;
    int ndim;

    bool exportMesh;
    bool exportData;
    bool exportArea;
    bool exportSpecials;
    bool exportReactionForces;

private:
    void computeArea();

    // Export reaction forces
    IntArray reactionForcesDofManList; // Holds which dof managers reaction forces should be exported from.

public:
    MatlabExportModule(int n, EngngModel *e);
    virtual ~MatlabExportModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep, bool forcedOutput=false);
    virtual void initialize();
    virtual void terminate();

    void doOutputMesh(TimeStep *tStep, FILE *FID);
    void doOutputData(TimeStep *tStep, FILE *FID);
    void doOutputSpecials(TimeStep *tStep, FILE *FID);
    void doOutputReactionForces(TimeStep *tStep, FILE *FID);

    virtual const char *giveClassName() const { return "MatlabExportModule"; };
    virtual const char *giveInputRecordName() const { return _IFT_MatlabExportModule_Name; }


};
} // end namespace oofem
#endif // matlabexportmodule_h_
