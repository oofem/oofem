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
#include "feinterpol.h"
#include "nodalaveragingrecoverymodel.h"

#include <cstdlib>

namespace oofem {
MMAShapeFunctProjection :: MMAShapeFunctProjection() : MaterialMappingAlgorithm()
{
    stateCounter = 0;
    domain = NULL;
}

MMAShapeFunctProjection :: ~MMAShapeFunctProjection()
{ }

void
MMAShapeFunctProjection :: __init(Domain *dold, IntArray &varTypes, FloatArray &coords, Set &elemSet, TimeStep *tStep, bool iCohesiveZoneGP)
//(Domain* dold, IntArray& varTypes, GaussPoint* gp, TimeStep* tStep)
{
    int nvar = varTypes.giveSize();
    // check time stemp
    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return;
    }


    // Project Gauss point components to nodes on old mesh
    if ( (int)this->smootherList.size() != nvar ) {
        this->smootherList.clear();
        this->smootherList.reserve(nvar);
        for ( int ivar = 1; ivar <= nvar; ivar++ ) {
            this->smootherList.emplace_back( new NodalAveragingRecoveryModel(dold) );
        }
    }

    this->intVarTypes = varTypes;
    for ( int ivar = 1; ivar <= nvar; ivar++ ) {
        this->smootherList[ivar-1]->recoverValues(elemSet, ( InternalStateType ) varTypes.at(ivar), tStep);
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
    std::vector< FloatArray > container;
    const FloatArray *nvec;

    int indx = this->intVarTypes.findFirstIndexOf( ( int ) type );
    if ( indx ) {
        container.reserve(nnodes);
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            this->smootherList[indx-1]->giveNodalVector( nvec, elem->giveDofManager(inode)->giveNumber() );
            container.emplace_back(*nvec);
        }

        this->interpolateIntVarAt(answer, elem, gp->giveNaturalCoordinates(),
                                  container, type, tStep);
    } else {
        OOFEM_ERROR("var not initialized");
    }

    return 1;
}


int
MMAShapeFunctProjection :: __mapVariable(FloatArray &answer, FloatArray &coords,
                                         InternalStateType type, TimeStep *tStep)
{
    FloatArray lcoords, closest;
    Element *elem = domain->giveSpatialLocalizer()->giveElementClosestToPoint(coords, closest, lcoords);
    if ( !elem ) {
        OOFEM_ERROR("no suitable source found");
    }

    int nnodes = elem->giveNumberOfDofManagers();
    std::vector< FloatArray > container;
    const FloatArray *nvec;

    int indx = this->intVarTypes.findFirstIndexOf( ( int ) type );
    if ( indx ) {
        container.reserve(nnodes);
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            this->smootherList.at(indx)->giveNodalVector( nvec, elem->giveDofManager(inode)->giveNumber() );
            container.emplace_back(*nvec);
        }

        this->interpolateIntVarAt(answer, elem, lcoords, container, type, tStep);
    } else {
        OOFEM_ERROR("var not initialized");
    }

    return 1;
}


int
MMAShapeFunctProjection :: mapStatus(MaterialStatus &oStatus) const
{
    OOFEM_ERROR("not implemented yet.")

    return 0;
}


void
MMAShapeFunctProjection :: interpolateIntVarAt(FloatArray &answer, Element *elem, const FloatArray &lcoords,
                                               std :: vector< FloatArray > &list, InternalStateType type, TimeStep *tStep) const
{
    FloatArray n;

    elem->giveInterpolation()->evalN( n, lcoords, FEIElementGeometryWrapper(elem) );

    answer.resize(0);
    for ( int i = 0; i < n.giveSize(); ++i ) {
        answer.add(n[i], list[i]);
    }
}

} // end namespace oofem
