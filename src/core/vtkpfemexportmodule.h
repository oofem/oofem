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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef vtkpfemexportmodule_h
#define vtkpfemexportmodule_h

#include "vtkbaseexportmodule.h"
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

///@name Input fields for VTK XML export module

using namespace std;
namespace oofem {
class Node;

#define _IFT_VTKPFEMXMLExportModule_Name "vtkpfem"          

class OOFEM_EXPORT VTKPFEMXMLExportModule : public VTKBaseExportModule
{
protected:
    std::ofstream fileStream;
public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    VTKPFEMXMLExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~VTKPFEMXMLExportModule();

    void initializeFrom(InputRecord &ir) override;

    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
protected:
    /// Returns the filename for the given time step.
    std::string giveOutputFileName(TimeStep *tStep);

    /// Returns the output stream for given solution step.
    std::ofstream giveOutputStream(TimeStep *tStep);

};

} // end namespace oofem
#endif // vtkxmlexportmodule_h
