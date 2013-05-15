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

#ifndef matlabexportmodule_h_
#define matlabexportmodule_h_

#include "exportmodule.h"

///@name Input fields for MatlabExportModule
//@{
#define _IFT_MatlabExportModule_Name "matlab"
#define _IFT_MatlabExportModule_mesh "mesh"
#define _IFT_MatlabExportModule_data "data"
#define _IFT_MatlabExportModule_area "area"
#define _IFT_MatlabExportModule_specials "specials"
//@}

namespace oofem {
/**
 * (Under development) The Matlab export module enables oofem to export the results to a textfile containing the description of the mesh used
 * along with the pertinent results.
 *
 * @author Carl Sandström
 */
class MatlabExportModule : public ExportModule
{
protected:
    /// list of InternalStateType values, identifying the selected vars for export
    IntArray internalVarsToExport;
    /// list of primary unknowns to export
    IntArray primaryVarsToExport;
    std :: string functionname;

    FILE *giveOutputStream(TimeStep *);
    double xmax, xmin, ymax, ymin;
    double Area;

    bool exportMesh;
    bool exportData;
    bool exportArea;
    bool exportSpecials;

private:
    void computeArea();

public:
    MatlabExportModule(int n, EngngModel *e);
    virtual ~MatlabExportModule();
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void doOutput(TimeStep *tStep);
    virtual void initialize();
    virtual void terminate();

    void doOutputMesh(TimeStep *tStep,  FILE *FID);
    void doOutputData(TimeStep *tStep,  FILE *FID);
    void doOutputSpecials(TimeStep *tStep,      FILE *FID);

    virtual const char *giveClassName() const { return "MatlabExportModule"; };
};
} // end namespace oofem
#endif // matlabexportmodule_h_
