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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "planstrssxfem.h"
#include "engngm.h"
#include "dof.h"
#include "masterdof.h"
#include "planstrss.h"
#include "structuralmaterial.h"
#include "patchintegrationrule.h"
#include "interfacetype.h"
#include "xfemelementinterface.h"
#include "structuralcrosssection.h"
#include "enrichmentitem.h"

namespace oofem {
Interface *
PlaneStress2dXfem :: giveInterface(InterfaceType interface)
{
    if ( interface != XfemElementInterfaceType ) {
        return PlaneStress2d :: giveInterface(interface);
    } else if ( interface == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else {
        return NULL;
    }
}


void PlaneStress2dXfem :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    int i;
    // evaluation of N,dNdx
    FloatMatrix dNdx;
    FloatArray N;
    interpolation.evaldNdx(dNdx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    interpolation.evalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    FloatMatrix *simple = new FloatMatrix(3, 8);
    simple->zero();

    // assemble standard FEM part of strain-displacement matrix
    for ( i = 1; i <= 4; i++ ) {
        simple->at(1, 2 * i - 1) = dNdx.at(i, 1);
        simple->at(2, 2 * i - 0) = dNdx.at(i, 2);
    }

    for ( i = 1; i <= 4; i++ ) {
        simple->at(3, 2 * i - 1) = dNdx.at(i, 2);
        simple->at(3, 2 * i - 0) = dNdx.at(i, 1);
    }


    // assemble xfem part of strain-displacement matrix
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    int counter = 0;
    AList< FloatMatrix >additionals;
    additionals.put(1, simple);
    counter += simple->giveNumberOfColumns();
    // loop over enrichment items
    for ( i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *er = xf->giveEnrichmentItem(i);
        int erndofs = er->giveNumberOfDofs();
        // enrichment function at the gauss point
        double efgp = er->giveEnrichmentFunction()->evaluateFunctionAt(gp, er);
        // derivative of enrichment function at the gauss point
        FloatArray efgpD;
        er->giveEnrichmentFunction()->evaluateDerivativeAt(efgpD, gp, er);
        // adds up the number of the dofs from an enrichment item
        // for each node
        for ( int j = 1; j <= this->giveNumberOfDofManagers(); j++ ) {
            if ( er->isDofManEnriched( dofManArray.at(j) ) ) {
                FloatArray *nodecoords = domain->giveDofManager( dofManArray.at(j) )->giveCoordinates();
                // ef is a FloatArray containing the value of EnrichmentFunction in a specific for all enriched dofs
                double efnode = er->giveEnrichmentFunction()->evaluateFunctionAt(nodecoords, er);
                // matrix to be added anytime a node is enriched
                FloatMatrix *toAdd = new FloatMatrix(3, erndofs);
                toAdd->zero();
                FloatArray help;
                help.resize(2);
                for ( int p = 1; p <= 2; p++ ) {
                    help.at(p) = dNdx.at(j, p) * ( efgp - efnode ) + N.at(j) * efgpD.at(p);
                }

                for ( int k = 1; k <= erndofs; k++ ) {
                    toAdd->at(k, k) = help.at(k);
                    if ( k == 1 ) {
                        toAdd->at(3, k) = help.at(2);
                    }

                    if ( k == 2 ) {
                        toAdd->at(3, k) = help.at(1);
                    }
                }

                int sz = additionals.giveSize();
                additionals.put(sz + 1, toAdd);
                counter += toAdd->giveNumberOfColumns();
            }
        }
    }

    answer.resize(3, counter);
    answer.zero();
    int columns = 1;
    for ( int i = 1; i <= additionals.giveSize(); i++ ) {
        for ( int j = 1; j <= additionals.at(i)->giveNumberOfColumns(); j++ ) {
            for ( int k = 1; k <= 3; k++ ) {
                answer.at(k, columns) = additionals.at(i)->at(k, j);
            }

            columns++;
        }
    }
}

void PlaneStress2dXfem :: giveLocationArray(IntArray &locationArray, EquationID, const UnknownNumberingScheme &s) const
{
    IntArray interactedEI;
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    xf->getInteractedEI( interactedEI, const_cast< PlaneStress2dXfem * >(this) );

    // for all enrichment items which element interacts
    for ( int i = 1; i <= ( const_cast< PlaneStress2dXfem * >(this) )->giveNumberOfDofManagers(); i++ ) {
        DofManager *dm = this->giveDomain()->giveDofManager( dofManArray.at(i) );
        for ( int j = 1; j <= xf->giveNumberOfEnrichmentItems(); j++ ) {
            EnrichmentItem *er = xf->giveEnrichmentItem(j);
            if ( er->isDofManEnriched( dofManArray.at(i) ) ) {
                IntArray *dofIdAr = er->getDofIdArray();
                for ( int k = 1; k <= dofIdAr->giveSize(); k++ ) {
                    if ( dm->hasDofID( (DofIDItem)dofIdAr->at(k) ) == false ) {
                        int sz = dm->giveNumberOfDofs();
                        Dof *df = new MasterDof( sz + 1, dm, 0, 0, (DofIDItem)dofIdAr->at(k) );
                        int eqN = xf->giveFictPosition( dofManArray.at(i) )->at(k);
                        df->setEquationNumber(eqN);
                        dm->appendDof(df);
                    }
                }
            }
        }
    }

    locationArray.resize(0);
    IntArray enriched;
    enriched.resize(0);
    for ( int i = 1; i <= ( const_cast< PlaneStress2dXfem * >(this) )->giveNumberOfDofManagers(); i++ ) {
        DofManager *dm = ( const_cast< PlaneStress2dXfem * >(this) )->giveDomain()->giveDofManager( dofManArray.at(i) );
        for ( int j = 1; j <= dm->giveNumberOfDofs(); j++ ) {
            int eqN = dm->giveDof(j)->giveEquationNumber(s);
            if ( j <= 2 ) {
                locationArray.followedBy(eqN);
            } else {
                enriched.followedBy(eqN);
            }
        }
    }

    locationArray.followedBy(enriched);
}


int PlaneStress2dXfem :: computeNumberOfDofs(EquationID ut)
{
    int ret = 0;
    for ( int i = 1; i <= giveNumberOfDofManagers(); i++ ) {
        DofManager *dm = this->giveDomain()->giveDofManager( dofManArray.at(i) );
        for ( int j = 1; j <= dm->giveNumberOfDofs(); j++ ) {
            ret++;
        }
    }

    return ret;
}


void
PlaneStress2dXfem :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    PlaneStress2d :: giveDofManDofIDMask(inode, EID_MomentumBalance, answer);
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    for ( int i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *er = xf->giveEnrichmentItem(i);
        if ( er->isDofManEnriched( dofManArray.at(inode) ) ) {
            IntArray *dofIdAr = er->getDofIdArray();
            for ( int j = 1; j <= dofIdAr->giveSize(); j++ ) {
                answer.followedBy( dofIdAr->at(j) );
            }
        }
    }
}

void PlaneStress2dXfem :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( xf->isInteracted(this) ) {
        PatchIntegrationRule *pir = ( PatchIntegrationRule * ) gp->giveIntegrationRule();
        StructuralMaterial *sm = ( StructuralMaterial * ) this->giveDomain()->giveMaterial( pir->giveMaterial() );
        sm->giveCharacteristicMatrix(answer, ReducedForm, rMode, gp, tStep);
    } else {
        PlaneStress2d :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
    }
}


void
PlaneStress2dXfem :: computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs in the element local c.s.
{
    int i, j, k, nDofs, size;
    IntArray elementNodeMask;
    FloatArray vec;
    answer.resize( size = this->computeNumberOfGlobalDofs(type) );
    k = 0;
    int m = 0;
    FloatArray p1(size - 8);

    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, type, elementNodeMask);
        this->giveDofManager(i)->giveUnknownVector(vec, elementNodeMask, type, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            if ( j <= 2 ) {
                answer.at(++k) = vec.at(j);
            } else {
                p1.at(++m) = vec.at(j);
            }
        }
    }

    for ( int i = 1; i <= m; i++ ) {
        answer.at(k + i) = p1.at(i);
    }
    // Rotate it as well? but this element doesn't support local coordinate systems anyway.
}


void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatArray Epsilon;
    this->computeStrainVector(Epsilon, gp, stepN);
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( xf->isInteracted(this) ) {
        PatchIntegrationRule *pir = ( PatchIntegrationRule * ) gp->giveIntegrationRule();
        StructuralMaterial *sm = ( StructuralMaterial * ) this->giveDomain()->giveMaterial( pir->giveMaterial() );
        sm->giveRealStressVector(answer, ReducedForm, gp, Epsilon, stepN);
    } else {
        StructuralCrossSection *cs = ( StructuralCrossSection * ) this->giveCrossSection();
        cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
    }
}

void PlaneStress2dXfem :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    this->computeStiffnessMatrix_withIRulesAsSubcells(answer, rMode, tStep);
}

void
PlaneStress2dXfem :: giveInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, int useUpdatedGpRecord) {
    this->giveInternalForcesVector_withIRulesAsSubcells(answer, tStep, useUpdatedGpRecord);
}


double PlaneStress2dXfem :: computeArea() const
{
    FloatArray *node1 = this->giveDofManager(1)->giveCoordinates();
    FloatArray *node2 = this->giveDofManager(2)->giveCoordinates();
    FloatArray *node3 = this->giveDofManager(3)->giveCoordinates();
    double a = node1->distance(node2);
    double b = node2->distance(node3);
    return a * b;
}

#ifdef __OOFEG
void PlaneStress2dXfem :: drawRawGeometry(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( !xf->isInteracted(this) ) {
        PlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( numberOfIntegrationRules > 1 ) {
            int i;
            PatchIntegrationRule *iRule;
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >(integrationRulesArray [ i ]);
                if ( iRule ) {
                    iRule->givePatch()->draw(context);
                }
            }
        } else {
            PlaneStress2d :: drawRawGeometry(context);
        }
    }
}

void PlaneStress2dXfem :: drawScalar(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( !xf->isInteracted(this) ) {
        PlaneStress2d :: drawScalar(context);
    } else {
        if ( context.giveIntVarMode() == ISM_local ) {
            int i, j, indx, result = 1;
            double val;
            FloatArray s(3), v;
            IntArray map;

            this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
            if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
                return;
            }


            TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
            PatchIntegrationRule *iRule;
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >(integrationRulesArray [ i ]);

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(0);
                    result += giveIPValue(v, gp, context.giveIntVarType(), tStep);
                    val += v.at(indx);
                }

                val /= iRule->getNumberOfIntegrationPoints();
 #endif
                s.at(1) = s.at(2) = s.at(3) = val;
                iRule->givePatch()->drawWD(context, s);
            }
        } else {
            PlaneStress2d :: drawScalar(context);
        }
    }
}
#endif
} // end namespace oofem
