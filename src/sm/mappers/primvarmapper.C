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

#include "primvarmapper.h"
#include "domain.h"
#include "dofmanager.h"
#include "linsystsolvertype.h"
#include "../sm/Elements/structuralelement.h"
#include "engngm.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "sparsemtrx.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"
#include "timestep.h"
#include "activebc.h"
#include "prescribedgradientbcweak.h"
#include "prescribedgradientbcneumann.h"
#include "unknownnumberingscheme.h"

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
    EModelDefaultEquationNumbering num;


    const int dim = iNewDom.giveNumberOfSpatialDimensions();

    int numElNew = iNewDom.giveNumberOfElements();

    // Count dofs
    int numDofsNew = engngMod->giveNumberOfDomainEquations( 1, num );


    oU.resize(numDofsNew);
    oU.zero();

    FloatArray du(numDofsNew);
    du.zero();

    FloatArray res(numDofsNew);

    std :: unique_ptr< SparseMtrx > K;
    std :: unique_ptr< SparseLinearSystemNM > solver;

    solver.reset( classFactory.createSparseLinSolver(ST_Petsc, & iOldDom, engngMod) );
    if (!solver) {
        solver.reset( classFactory.createSparseLinSolver(ST_Direct, & iOldDom, engngMod) );
    }
    K.reset( classFactory.createSparseMtrx(solver->giveRecommendedMatrix(true)) );

    K->buildInternalStructure( engngMod, iNewDom.giveNumber(), num );

    int maxIter = 1;

    for ( int iter = 0; iter < maxIter; iter++ ) {
        K->zero();
        res.zero();


        // Contribution from elements
        for ( int elIndex = 1; elIndex <= numElNew; elIndex++ ) {
            StructuralElement *elNew = dynamic_cast< StructuralElement * >( iNewDom.giveElement(elIndex) );
            if ( elNew == NULL ) {
                OOFEM_ERROR("Failed to cast Element new to StructuralElement.");
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
            elNew->giveLocationArray( elDofsGlob, num );


            // Loop over Gauss points
            for ( int intRuleInd = 0; intRuleInd < elNew->giveNumberOfIntegrationRules(); intRuleInd++ ) {

                for ( GaussPoint *gp: *elNew->giveIntegrationRule(intRuleInd) ) {

                    // New N-matrix
                    FloatMatrix NNew;
                    elNew->computeNmatrixAt(gp->giveNaturalCoordinates(), NNew);


                    //////////////
                    // Global coordinates of GP
                    const FloatArray &localCoord = gp->giveNaturalCoordinates();
                    FloatArray globalCoord;
                    elNew->computeGlobalCoordinates(globalCoord, localCoord);
                    //////////////


                    // Localize element and point in the old domain
                    FloatArray localCoordOld(dim), pointCoordOld(dim);
                    StructuralElement *elOld = dynamic_cast< StructuralElement * >( iOldDom.giveSpatialLocalizer()->giveElementClosestToPoint(localCoordOld, pointCoordOld, globalCoord, 0) );
                    if ( elOld == NULL ) {
                        OOFEM_ERROR("Failed to cast Element old to StructuralElement.");
                    }


                    // Compute N-Matrix for the old element
                    FloatMatrix NOld;
                    elOld->computeNmatrixAt(localCoordOld, NOld);

                    // Fetch nodal displacements for the new element
                    FloatArray nodeDispNew( elDofsGlob.giveSize() );


                    int dofsPassed = 1;
                    for ( int elNode: elNew->giveDofManArray() ) {
//                        DofManager *dMan = iNewDom.giveNode(elNode);
                        DofManager *dMan = iNewDom.giveDofManager(elNode);

                        for ( Dof *dof: *dMan ) {
                            if ( elDofsGlob.at(dofsPassed) != 0 ) {
                                nodeDispNew.at(dofsPassed) = oU.at( elDofsGlob.at(dofsPassed) );
                            } else {
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
                    FloatArray nodeDispOld;
                    dofsPassed = 1;
                    IntArray elDofsGlobOld;
                    elOld->giveLocationArray( elDofsGlobOld, num );

//                    elOld->computeVectorOf(iMode, &(iTStep), nodeDisp);
                    int numElNodesOld = elOld->giveNumberOfDofManagers();
                    for(int nodeIndOld = 1; nodeIndOld <= numElNodesOld; nodeIndOld++) {
                        DofManager *dManOld = elOld->giveDofManager(nodeIndOld);

                        for ( Dof *dof: *dManOld ) {
                            if ( elDofsGlobOld.at(dofsPassed) != 0 ) {
                                FloatArray dofUnknowns;

                                if(dof->giveEqn() > 0) {
                                    dof->giveUnknowns(dofUnknowns, iMode, &iTStep);

#ifdef DEBUG
                                    if(!dofUnknowns.isFinite()) {
                                        OOFEM_ERROR("!dofUnknowns.isFinite()")
                                    }

                                    if(dofUnknowns.giveSize() < 1) {
                                        OOFEM_ERROR("dofUnknowns.giveSize() < 1")
                                    }
#endif
                                    nodeDispOld.push_back(dofUnknowns.at(1));
                                }
                                else {
                                    // TODO: Why does this case occur?
                                    nodeDispOld.push_back(0.0);
                                }
                            } else {
                                if ( dof->hasBc(& iTStep) ) {
//                                    printf("hasBC.\n");
#ifdef DEBUG
                                    if(!std::isfinite(dof->giveBcValue(iMode, & iTStep))) {
                                        OOFEM_ERROR("!std::isfinite(dof->giveBcValue(iMode, & iTStep))")
                                    }
#endif
                                    nodeDispOld.push_back( dof->giveBcValue(iMode, & iTStep) );
                                }
                                else {
//                                    printf("Unhandled case in LSPrimaryVariableMapper :: mapPrimaryVariables().\n");
                                    nodeDispOld.push_back( 0.0 );
                                }
                            }

                            dofsPassed++;
                        }

                    }


                    FloatArray oldDisp;
                    oldDisp.beProductOf(NOld, nodeDispOld);

                    FloatArray temp, duu;

#ifdef DEBUG
                    if(!oldDisp.isFinite()) {
                        OOFEM_ERROR("!oldDisp.isFinite()")
                    }

                    if(!newDisp.isFinite()) {
                        OOFEM_ERROR("!newDisp.isFinite()")
                    }
#endif

                    duu.beDifferenceOf(oldDisp, newDisp);
                    temp.beTProductOf(NNew, duu);
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

#ifdef DEBUG
            if(!me.isFinite()) {
                OOFEM_ERROR("!me.isFinite()")
            }
#endif

            K->assemble(elDofsGlob, me);
            res.assemble(elRes, elDofsGlob);
        }


        // Contribution from active boundary conditions
        int numBC = iNewDom.giveNumberOfBoundaryConditions();
        for(int bcInd = 1; bcInd <= numBC; bcInd++) {

            PrescribedGradientBCWeak *activeBC = dynamic_cast<PrescribedGradientBCWeak*> ( iNewDom.giveBc(bcInd) );

            if(activeBC != NULL) {
                IntArray tractionRows;
                activeBC->giveTractionLocationArray(tractionRows, num);

                // TODO: Compute correct mass matrix and add residual contribution.

                FloatMatrix mNode(1,1);
                mNode.beUnitMatrix();

                for(int tracDofInd : tractionRows) {
                    const IntArray tracDofArray = {tracDofInd};
                    K->assemble(tracDofArray, tracDofArray, mNode);
                }

            }


            PrescribedGradientBCNeumann *neumannBC = dynamic_cast<PrescribedGradientBCNeumann*> ( iNewDom.giveBc(bcInd) );

            if(neumannBC != NULL) {
                IntArray stressRows;
                neumannBC->giveStressLocationArray(stressRows, num);

                FloatMatrix massMtrxBc( stressRows.giveSize(), stressRows.giveSize() );
                massMtrxBc.beUnitMatrix(); // TODO: Compute correct mass matrix and add residual contribution.
                K->assemble(stressRows, massMtrxBc);
            }

        }

#ifdef DEBUG
        if(!res.isFinite()) {
            OOFEM_ERROR("!res.isFinite()")
        }
#endif

//        printf("iter: %d res norm: %e\n", iter, res.computeNorm() );
        K->assembleBegin();
        K->assembleEnd();
//        K->writeToFile("Kmapping.txt");

        // Solve
        solver->solve(*K, res, du);

#ifdef DEBUG
        if(!du.isFinite()) {
            OOFEM_ERROR("!du.isFinite()")
        }
#endif

        oU.add(du);
    }
}
} /* namespace oofem */
