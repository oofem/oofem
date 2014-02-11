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

#include "mmashapefunctprojection.h"
#include "dofmanager.h"
#include "gausspoint.h"
#include "element.h"
#include "domain.h"
#include "spatiallocalizer.h"
#include "timestep.h"
#include "nodalaveragingrecoverymodel.h"

#include <cstdlib>

namespace oofem {
MMAShapeFunctProjection :: MMAShapeFunctProjection() : MaterialMappingAlgorithm()
{
    stateCounter = 0;
    //smootherList(0);
    domain = NULL;
}

MMAShapeFunctProjection :: ~MMAShapeFunctProjection()
{ }

void
MMAShapeFunctProjection :: __init(Domain *dold, IntArray &varTypes, FloatArray &coords, int region, TimeStep *tStep, bool iCohesiveZoneGP)
//(Domain* dold, IntArray& varTypes, GaussPoint* gp, TimeStep* tStep)
{
    int nvar = varTypes.giveSize();
    // check time stemp
    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return;
    }


    // Project Gauss point components to nodes on old mesh
    if ( this->smootherList.giveSize() != nvar ) {
        this->smootherList.clear();
        this->smootherList.growTo(nvar);
        for ( int ivar = 1; ivar <= nvar; ivar++ ) {
            this->smootherList.put( ivar, new NodalAveragingRecoveryModel(dold) );
        }
    }

    this->intVarTypes = varTypes;
    for ( int ivar = 1; ivar <= nvar; ivar++ ) {
        this->smootherList.at(ivar)->recoverValues( ( InternalStateType ) varTypes.at(ivar), tStep );
    }

    // remember time stemp
    stateCounter = tStep->giveSolutionStateCounter();
    this->domain = dold;
}


void
MMAShapeFunctProjection :: finish(TimeStep *tStep)
{
    this->smootherList.clear();
    stateCounter = -1;
}

int
MMAShapeFunctProjection :: mapVariable(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    Element *elem = gp->giveElement();
    int nnodes = elem->giveNumberOfDofManagers();
    MMAShapeFunctProjectionInterface :: nodalValContainerType container(nnodes);
    MMAShapeFunctProjectionInterface *interface;
    const FloatArray *nvec;

    if ( ( interface = static_cast< MMAShapeFunctProjectionInterface * >
                       ( elem->giveInterface(MMAShapeFunctProjectionInterfaceType) ) ) == NULL ) {
        abort();
    }

    int indx = this->intVarTypes.findFirstIndexOf( ( int ) type );
    if ( indx ) {
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            container.put(inode, new FloatArray);
            this->smootherList.at(indx)->giveNodalVector( nvec, elem->giveDofManager(inode)->giveNumber(),
                                                          elem->giveRegionNumber() );
            * ( container.at(inode) ) = * nvec;
        }

        interface->MMAShapeFunctProjectionInterface_interpolateIntVarAt(answer, ( * gp->giveCoordinates() ),
                                                                        MMAShapeFunctProjectionInterface :: coordType_local,
                                                                        container, type, tStep);
    } else {
        OOFEM_ERROR("MMAShapeFunctProjection::mapVariable: var not initialized");
    }

    return 1;
}


int
MMAShapeFunctProjection :: __mapVariable(FloatArray &answer, FloatArray &coords,
                                         InternalStateType type, TimeStep *tStep)
{
    Element *elem = domain->giveSpatialLocalizer()->giveElementContainingPoint(coords);
    if ( !elem ) {
        OOFEM_ERROR("MMAShapeFunctProjection::__mapVariable: no suitable source found");
    }

    int nnodes = elem->giveNumberOfDofManagers();
    MMAShapeFunctProjectionInterface :: nodalValContainerType container(nnodes);
    MMAShapeFunctProjectionInterface *interface;
    const FloatArray *nvec;

    if ( ( interface = static_cast< MMAShapeFunctProjectionInterface * >
                       ( elem->giveInterface(MMAShapeFunctProjectionInterfaceType) ) )  == NULL ) {
        abort();
    }

    int indx = this->intVarTypes.findFirstIndexOf( ( int ) type );
    if ( indx ) {
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            container.put(inode, new FloatArray);
            this->smootherList.at(indx)->giveNodalVector( nvec, elem->giveDofManager(inode)->giveNumber(),
                                                          elem->giveRegionNumber() );
            * ( container.at(inode) ) = * nvec;
        }

        interface->MMAShapeFunctProjectionInterface_interpolateIntVarAt(answer, coords,
                                                                        MMAShapeFunctProjectionInterface :: coordType_global,
                                                                        container, type, tStep);
    } else {
        OOFEM_ERROR("MMAShapeFunctProjection::__mapVariable: var not initialized");
    }

    return 1;
}

int
MMAShapeFunctProjection :: mapStatus(MaterialStatus &oStatus) const
{
    OOFEM_ERROR("ERROR: MMAShapeFunctProjection :: mapStatus() is not implemented yet.")

    return 0;
}
} // end namespace oofem
