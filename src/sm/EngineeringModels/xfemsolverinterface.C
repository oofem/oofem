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

#include "sm/EngineeringModels/xfemsolverinterface.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "sm/xfem/xfemstructuralelementinterface.h"
#include "sm/mappers/primvarmapper.h"
#include "timestep.h"
#include "structengngmodel.h"
#include "staticstructural.h"
#include "primaryfield.h"
#include "domain.h"
#include "xfem/xfemmanager.h"
#include "element.h"
#include "matstatmapperint.h"
#include "nummet.h"
#include "floatarray.h"
#include "exportmodulemanager.h"
#include "vtkxmlexportmodule.h"

namespace oofem {

XfemSolverInterface::XfemSolverInterface() :
    mNeedsVariableMapping(false)
{ }

void XfemSolverInterface::propagateXfemInterfaces(TimeStep *tStep, StructuralEngngModel &ioEngngModel, bool iRecomputeStepAfterCrackProp)
{
    int domainInd = 1;
    Domain *domain = ioEngngModel.giveDomain(domainInd);

    if ( domain->hasXfemManager() ) {
        XfemManager *xMan = domain->giveXfemManager();
        bool frontsHavePropagated = false;
        if ( xMan->hasInitiationCriteria() ) {
            // TODO: generalise this?
            // Intitiate delaminations (only implemented for listbasedEI/delamination. Treated the same way as propagation)
            xMan->initiateFronts(frontsHavePropagated,tStep);
        }

        if ( xMan->hasPropagatingFronts() ) {
            // Propagate crack tips
            xMan->propagateFronts(frontsHavePropagated);
        }

        bool eiWereNucleated = false;
        if ( xMan->hasNucleationCriteria() ) {
        	xMan->nucleateEnrichmentItems(eiWereNucleated);
        }

        for ( auto &elem : domain->giveElements() ) {
            ////////////////////////////////////////////////////////
            // Map state variables for enriched elements
            XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( elem.get() );

            if ( xfemElInt ) {
                xfemElInt->XfemElementInterface_updateIntegrationRule();
            }
        }

        if ( frontsHavePropagated || eiWereNucleated ) {
            mNeedsVariableMapping = false;

            ioEngngModel.giveDomain(1)->postInitialize();
            ioEngngModel.forceEquationNumbering();

            if ( iRecomputeStepAfterCrackProp ) {
                OOFEM_LOG_RELEVANT("Recomputing time step.\n");
                ioEngngModel.forceEquationNumbering();
                ioEngngModel.solveYourselfAt(tStep);
                ioEngngModel.updateYourself( tStep );
                ioEngngModel.terminate( tStep );
            }
        }
    }
}


} /* namespace oofem */
