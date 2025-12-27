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

#include "vtkmemoryexportmodule.h"
#include "element.h"
#include "timestep.h"
#include "engngm.h"
#include "node.h"
#include "nodalaveragingrecoverymodel.h"
#include "zznodalrecoverymodel.h"
#include "classfactory.h"

namespace oofem {
REGISTER_ExportModule(VTKMemoryExportModule)

VTKMemoryExportModule::VTKMemoryExportModule(int n, EngngModel *e) : VTKBaseExportModule(n,e) {}
VTKMemoryExportModule::~VTKMemoryExportModule() {}

void
VTKMemoryExportModule::initializeFrom(InputRecord &ir)
{
  VTKBaseExportModule::initializeFrom(ir);
  
  IR_GIVE_OPTIONAL_FIELD(ir, cellVarsToExport, _IFT_VTKXMLExportModule_cellvars); // Macro - see internalstatetype.h
  IR_GIVE_OPTIONAL_FIELD(ir, internalVarsToExport, _IFT_VTKXMLExportModule_vars); // Macro - see internalstatetype.h
  IR_GIVE_OPTIONAL_FIELD(ir, primaryVarsToExport, _IFT_VTKXMLExportModule_primvars); // Macro - see unknowntype.h
  IR_GIVE_OPTIONAL_FIELD(ir, externalForcesToExport, _IFT_VTKXMLExportModule_externalForces); // Macro - see unknowntype.h
  IR_GIVE_OPTIONAL_FIELD(ir, ipInternalVarsToExport, _IFT_VTKXMLExportModule_ipvars); // Macro - see internalstatetype.h
}

void
VTKMemoryExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
  if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
    return;
  }

  int nPiecesToExport = this->giveNumberOfRegions(); //old name: region, meaning: sets
  ZZNodalRecoveryModel smoother(emodel->giveDomain(1));
  NodalAveragingRecoveryModel primVarSmoother(emodel->giveDomain(1));

  this->vtkPieces.resize(nPiecesToExport);
  
  // loop over regular pieces only (no support for composite elements at present)
  for ( int pieceNum = 1; pieceNum <= nPiecesToExport; pieceNum++ ) {
    ExportRegion& p = this->vtkPieces[pieceNum-1];
    p.clear();
    // Fills a data struct (VTKPiece) with all the necessary data.
    Set* region = this->giveRegionSet(pieceNum);
    this->setupVTKPiece(p, tStep, *region);
    // Export primary, internal and XFEM variables as nodal quantities
    this->exportPrimaryVars(p, *region, primaryVarsToExport, primVarSmoother, tStep);
    this->exportIntVars(p, *region, internalVarsToExport, smoother, tStep);
    this->exportExternalForces(p, *region, externalForcesToExport, tStep);
    this->exportCellVars(p, *region, cellVarsToExport, tStep);
    
  }
}

std::vector< ExportRegion>& 
VTKMemoryExportModule::getExportRegions() {
  return this->vtkPieces;
}


} // end namespace oofem
