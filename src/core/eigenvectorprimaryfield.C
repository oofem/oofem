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

#include "eigenvectorprimaryfield.h"
#include "timestep.h"
#include "domain.h"
#include "dofmanager.h"
#include "dof.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "engngm.h"
#include "set.h"
#include "boundarycondition.h"
#include "initialcondition.h"
#include "element.h"
#include "activebc.h"


namespace oofem {
EigenVectorPrimaryField :: EigenVectorPrimaryField(EngngModel *a, int idomain, FieldType ft, int nHist):
    DofDistributedPrimaryField(a, idomain, ft, nHist)
{ }


EigenVectorPrimaryField :: ~EigenVectorPrimaryField()
{ }


double
EigenVectorPrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    if ( mode != VM_Total ) OOFEM_ERROR("Only VM_Total is applicable to eigen vector fields");
    return dof->giveUnknownsDictionaryValue(tStep, mode);
}


void
EigenVectorPrimaryField :: updateAll(const FloatMatrix &eigenVectors, const UnknownNumberingScheme &s)
{
    Domain *d = emodel->giveDomain(domainIndx);

    for (int i = 1; i <= eigenVectors.giveNumberOfColumns(); ++i ) {
        TimeStep step(i, nullptr, 0, 1.0, 0.0, 0);
        this->applyBoundaryCondition(&step);
    }

    auto set_values = [&eigenVectors, &s](DofManager &dman) {
        for ( Dof *dof: dman ) {
            if ( !dof->isPrimaryDof() ) continue;
            int eqNum = dof->giveEquationNumber(s);
            if ( eqNum > 0 ) {
                for (int i = 1; i <= eigenVectors.giveNumberOfColumns(); ++i ) {
                    // updateUnknownsDict should be changed to let the field have direct input rather than calling such a method.
                    // Until then, we "fake" a time step, relying on the fact that only the time-step number matters here.
                    TimeStep step(i, nullptr, 0, 1.0, 0.0, 0);
                    dof->updateUnknownsDictionary(&step, VM_Total, eigenVectors.at(eqNum, i));
                }
            }
        }
    };

    for ( auto &node : d->giveDofManagers() ) {
        set_values(*node);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            set_values(*elem->giveInternalDofManager(i));
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            set_values(*bc->giveInternalDofManager(i));
        }
    }
}


void
EigenVectorPrimaryField :: applyDefaultInitialCondition()
{ }


void
EigenVectorPrimaryField :: advanceSolution(TimeStep *tStep)
{ }

} // end namespace oofem
