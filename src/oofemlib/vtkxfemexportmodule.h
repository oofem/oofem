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

#ifndef vtkxfemexportmodule_h
#define vtkxfemexportmodule_h

#include "vtkxmlexportmodule.h"
#include "intarray.h"
#include "nodalrecoverymodel.h"
#include "interface.h"
#include "internalstatevaluetype.h"
#include "integrationrule.h"
#include "xfem/xfemmanager.h"
#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef _WIN32
 #define NULL_DEVICE "NUL:"
#else
 #define NULL_DEVICE "/dev/null"
#endif

#include <string>
#include <list>

using namespace std;
namespace oofem {
class Node;

///@name Input fields for VTK XML export module
//@{
#define _IFT_VTKXMLXFemExportModule_Name "vtkxmlxfem"
//@}

/**
 * Represents VTK (Visualization Toolkit) export module for Xfem. 
 * It uses VTK (.vtu) file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible
 * nonsmooth character of some internal variables at region boundaries.
 * The exported variables are dermined by XFemManager (FemMan->vtkExportFields).
 * Each region and each enrichment item is exported as a single piece. 
 * When region contains composite cells, these are assumed to be
 * exported in individual subsequent pieces after the default one for the particular region.
 */
class OOFEM_EXPORT VTKXMLXFemExportModule : public VTKXMLExportModule
{
protected:

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKXMLXFemExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKXMLXFemExportModule();

    void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    void terminate() override {}
    const char *giveClassName() const override { return "VTKXMLXFemExportModule"; }

protected:

    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

    bool writeXFEMVars(ExportRegion &vtkPiece, int field, int enrItIndex);
    void getNodalVariableFromXFEMST(FloatArray &answer, Node *node, TimeStep *tStep, XFEMStateType xfemstype, Set &region, EnrichmentItem *ei);
    void exportIntVars(ExportRegion &vtkPiece, Set& region, int field, int enrItIndex,  IntArray& internalVarsToExport, NodalRecoveryModel& smoother, TimeStep *tStep);
    void giveDataHeaders(std::string &pointHeader, std::string &cellHeader) override;     // returns the headers
};


} // end namespace oofem
#endif // vtkxfemexportmodule_h
