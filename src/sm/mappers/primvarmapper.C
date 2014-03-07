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

/*
 * primvarmapper.C
 *
 * @author: Erik Svenning
 */

#include "primvarmapper.h"
#include "domain.h"
#include "dofmanager.h"
#include "linsystsolvertype.h"
#include "structuralelement.h"
#include "engngm.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"

#include <fstream>

namespace oofem {
PrimaryVariableMapper :: PrimaryVariableMapper() { }

PrimaryVariableMapper :: ~PrimaryVariableMapper() { }

LSPrimaryVariableMapper :: LSPrimaryVariableMapper()
{ }

LSPrimaryVariableMapper :: ~LSPrimaryVariableMapper()
{ }

void LSPrimaryVariableMapper :: mapPrimaryVariables(FloatArray &oU, Domain &iOldDom, Domain &iNewDom, ValueModeType iMode, TimeStep &iTStep)
{
    EngngModel *engngMod = iNewDom.giveEngngModel();


    const int dim = iNewDom.giveNumberOfSpatialDimensions();

    int numElNew = iNewDom.giveNumberOfElements();

    // Count dofs
    int numDofsNew = engngMod->giveNumberOfDomainEquations( 1, engngMod->giveUnknownNumberingScheme(EID_MomentumBalance) );


    oU.resize(numDofsNew);
    oU.zero();

    FloatArray du(numDofsNew);
    du.zero();

    FloatArray res(numDofsNew);
    SparseMtrx *K = classFactory.createSparseMtrx(SMT_Skyline);
    SparseLinearSystemNM *solver = classFactory.createSparseLinSolver(ST_Direct, & iOldDom, engngMod);


    K->buildInternalStructure( engngMod, 1, EID_MomentumBalance, engngMod->giveUnknownNumberingScheme(EID_MomentumBalance) );

    int maxIter = 1;

    for ( int iter = 0; iter < maxIter; iter++ ) {
        K->zero();
        res.zero();

        for ( int elIndex = 1; elIndex <= numElNew; elIndex++ ) {
            StructuralElement *elNew = dynamic_cast< StructuralElement * >( iNewDom.giveElement(elIndex) );
            if ( elNew == NULL ) {
                OOFEM_SIMPLE_ERROR("In LSPrimaryVariableMapper::mapPrimaryVariables(): Failed to cast Element new to StructuralElement.\n");
            }

            ///////////////////////////////////
            // Compute residual

            // Count element dofs
            int numElNodes = elNew->giveNumberOfDofManagers();
            int numElDofs = 0;
            for ( int i = 1; i <= numElNodes; i++ ) {
                numElDofs += elNew->giveDofManager(i)->giveNumberOfDofs();
            }

            FloatArray elRes(numElDofs);
            elRes.zero();

            IntArray elDofsGlob;
            elNew->giveLocationArray( elDofsGlob, EID_MomentumBalance, engngMod->giveUnknownNumberingScheme(EID_MomentumBalance) );


            // Loop over Gauss points
            for ( int intRuleInd = 0; intRuleInd < elNew->giveNumberOfIntegrationRules(); intRuleInd++ ) {
                IntegrationRule *iRule = elNew->giveIntegrationRule(intRuleInd);
                int numGP = iRule->giveNumberOfIntegrationPoints();

                for ( int gpInd = 0; gpInd < numGP; gpInd++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(gpInd);

                    // New N-matrix
                    FloatMatrix NNew;
                    elNew->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), NNew);


                    //////////////
                    // Global coordinates of GP
                    const int nDofMan = elNew->giveNumberOfDofManagers();

                    FloatArray Nc;
                    FEInterpolation *interp = elNew->giveInterpolation();
                    const FloatArray &localCoord = * ( gp->giveCoordinates() );
                    interp->evalN( Nc, localCoord, FEIElementGeometryWrapper(elNew) );

                    const IntArray &elNodes = elNew->giveDofManArray();

                    FloatArray globalCoord(dim);
                    globalCoord.zero();

                    for ( int i = 1; i <= nDofMan; i++ ) {
                        DofManager *dMan = elNew->giveDofManager(i);

                        for ( int j = 1; j <= dim; j++ ) {
                            globalCoord.at(j) += Nc.at(i) * dMan->giveCoordinate(j);
                        }
                    }
                    //////////////


                    // Localize element and point in the old domain
                    FloatArray localCoordOld(dim), pointCoordOld(dim);
                    StructuralElement *elOld = dynamic_cast< StructuralElement * >( iOldDom.giveSpatialLocalizer()->giveElementClosestToPoint(localCoordOld, pointCoordOld, globalCoord, 0) );
                    if ( elOld == NULL ) {
                        OOFEM_SIMPLE_ERROR("In LSPrimaryVariableMapper::mapPrimaryVariables(): Failed to cast Element old to StructuralElement.\n");
                    }


                    // Compute N-Matrix for the old element
                    FloatMatrix NOld;
                    elOld->computeNmatrixAt(localCoordOld, NOld);

                    // Fetch nodal displacements for the new element
                    FloatArray nodeDispNew( elDofsGlob.giveSize() );


                    int dofsPassed = 1;
                    for ( int i = 1; i <= elNodes.giveSize(); i++ ) {
                        DofManager *dMan = elNew->giveDofManager(i);

                        for ( int j = 1; j <= dMan->giveNumberOfDofs(); j++ ) {
                            if ( elDofsGlob.at(dofsPassed) != 0 ) {
                                nodeDispNew.at(dofsPassed) = oU.at( elDofsGlob.at(dofsPassed) );
                            } else {
                                Dof *dof = dMan->giveDof(j);

                                if ( dof->hasBc(& iTStep) ) {
                                    nodeDispNew.at(dofsPassed) = dof->giveBcValue(iMode, & iTStep);
                                }
                            }

                            dofsPassed++;
                        }
                    }


                    FloatArray newDisp;
                    newDisp.beProductOf(NNew, nodeDispNew);


                    // Fetch nodal displacements for the old element
                    FloatArray nodeDisp;
                    elOld->computeVectorOf(EID_MomentumBalance, iMode, & iTStep, nodeDisp);

                    FloatArray oldDisp;
                    oldDisp.beProductOf(NOld, nodeDisp);

                    FloatArray temp, du;
                    du.beDifferenceOf(oldDisp, newDisp);
                    temp.beTProductOf(NNew, du);
                    double dV = elNew->computeVolumeAround(gp);
                    elRes.add(dV, temp);
                }
            }


            ///////////////////////////////////
            // Compute matrix
            FloatMatrix me;

            double mass = 0.0;
            double density = 1.0;
            elNew->computeConsistentMassMatrix(me, & iTStep, mass, & density);

            K->assemble(elDofsGlob, me);
            res.assemble(elRes, elDofsGlob);
        }

        //		printf("iter: %d res norm: %e\n", iter, res.computeNorm() );


        // Solve
        solver->solve(K, & res, & du);
        oU.add(du);
    }

    delete solver;
    delete K;
}
} /* namespace oofem */
