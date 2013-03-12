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
#include "xfemmanager.h"

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

#if 0

    // Standard (continuous) part of the B-matrix
    FloatMatrix Bc;
    PlaneStress2d :: computeBmatrixAt(gp, Bc, li, ui);

    // Enriched (discontinuous) part of the B-matrix
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);

    // loop over enrichment items
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        
        // Evaluate enrichment function and its gradient at a given gauss point     
        EnrichmentFunction *ef = ei->giveEnrichmentFunction(1);
        double valEF = ef->evaluateFunctionAt(gp, ei);
        FloatArray gradEF;
        ef->evaluateDerivativeAt(gradEF, gp, ei);

        // nabla*ef
        FloatMatrix nablaEF(3,2);
        nablaEF.zero();
        nablaEF.at(1,1) = gradEF.at(1);
        nablaEF.at(2,2) = gradEF.at(2);
        nablaEF.at(3,1) = gradEF.at(2);
        nablaEF.at(3,2) = gradEF.at(1);

        FloatMatrix Nc;
        PlaneStress2d ::computeNmatrixAt(gp, Nc);

        FloatMatrix Bd, B(2,16);
        Bd.beProductOf(nablaEF,Nc);
        Bd.add(valEF,Bc);

        //answer.resize(3,16);
        //answer.setSubMatrix(Bc,1,1);
        //answer.setSubMatrix(Bd,1,9);

        B.resize(3,16);
        B.setSubMatrix(Bc,1,1);
        B.setSubMatrix(Bd,1,9);
        IntArray maskCol(4), maskRow(3);
        maskCol.setValues(8,  1, 2, 3, 4, 5, 6, 7, 8);
        maskRow.setValues(3,  1, 2, 3);
        // Remove part which is not enriched
        for ( int j = 1; j <= 4; j++ ) {
            DofManager *dMan = this->giveDofManager(j);
            if ( ei->isDofManEnriched(dMan) ) {
                maskCol.followedBy(8+(j-1)*2 + 1);
                maskCol.followedBy(8+(j-1)*2 + 2);
            }
        }

        //answer.resize(maskRow.giveSize(), maskCol.giveSize());
        answer.beSubMatrixOf(B, maskRow, maskCol);    
        
    }
}

#else

    FloatMatrix dNdx;
    FloatArray N;
    interpolation.evaldNdx( dNdx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    interpolation.evalN(     N  , * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    FloatMatrix *simple = new FloatMatrix(3, 8);
    simple->zero();

    // assemble standard FEM part of strain-displacement matrix
    for ( int i = 1; i <= 4; i++ ) {
        simple->at(1, 2 * i - 1) = dNdx.at(i, 1);
        simple->at(2, 2 * i - 0) = dNdx.at(i, 2);
    }

    for ( int i = 1; i <= 4; i++ ) {
        simple->at(3, 2 * i - 1) = dNdx.at(i, 2);
        simple->at(3, 2 * i - 0) = dNdx.at(i, 1);
    }


    // assemble xfem part of strain-displacement matrix
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    int counter = 0;
    AList< FloatMatrix >additionals;
    additionals.put(1, simple);
    counter += simple->giveNumberOfColumns();
    // loop over enrichment items
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        int erndofs = 2 ; //ei->giveNumberOfDofs();
        
        // enrichment function at the gauss point     
        EnrichmentFunction *ef = ei->giveEnrichmentFunction(1);
        double efgp = ef->evaluateFunctionAt(gp, ei);

        // derivative of enrichment function at the gauss point
        FloatArray efgpD;
        ef->evaluateDerivativeAt(efgpD, gp, ei);

        // adds up the number of the dofs from an enrichment item
        // this part is for the construction of a shifted enrichment
        for ( int j = 1; j <= this->giveNumberOfDofManagers(); j++ ) {
            DofManager *dMan = this->giveDofManager(j);
            if ( ei->isDofManEnriched( dMan ) ) {
                FloatArray *nodecoords = dMan->giveCoordinates();
                
                // should ask after specific EF in a loop
                double efnode = ef->evaluateFunctionAt(nodecoords, ei);

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
#endif


int PlaneStress2dXfem :: computeNumberOfDofs(EquationID ut)
{
    int ndofs = 0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        ndofs += this->giveDofManager(i)->giveNumberOfDofs();
    }
    return ndofs;
}


void
PlaneStress2dXfem :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    // Returns the total id mask of the dof manager - regular id's + enriched id's
    this->giveDofManager(inode)->giveCompleteMasterDofIDArray(answer);
}



void PlaneStress2dXfem :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    if ( xMan->isElementEnriched(this) ) {
        Inclusion *ei = static_cast< Inclusion * > ( xMan->giveEnrichmentItem(1) );
        StructuralMaterial *sm = static_cast< StructuralMaterial * >( ei->giveMaterial() );
        sm->giveCharacteristicMatrix(answer, ReducedForm, rMode, gp, tStep);
    } else {
        PlaneStress2d :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
    }
}

void
PlaneStress2dXfem :: computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown. Reorders them such that the regular (continuous) dofs
// comes first and then the enriched dofs.
{
    int k, m, nDofs, size;
    IntArray elementNodeMask;
    FloatArray vec;
    answer.resize( size = this->computeNumberOfGlobalDofs(type) );
    FloatArray p1(size - 8);

    m = k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, type, elementNodeMask);
        this->giveDofManager(i)->giveUnknownVector(vec, elementNodeMask, type, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
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
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    if ( xMan->isElementEnriched(this) ) {
        Inclusion *ei = static_cast< Inclusion * > ( xMan->giveEnrichmentItem(1) );
        StructuralMaterial *sm = static_cast< StructuralMaterial * >( ei->giveMaterial() );
        sm->giveRealStressVector(answer, ReducedForm, gp, Epsilon, stepN);
    } else {
        StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( this->giveCrossSection() );
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

Element_Geometry_Type 
PlaneStress2dXfem :: giveGeometryType() const 
{ 
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    if ( xMan->isElementEnriched(this) ) {
        return EGT_Composite;
    } else {
        return EGT_quad_1; 
    }
    
}


#ifdef __OOFEG
void PlaneStress2dXfem :: drawRawGeometry(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    if ( !xf->isElementEnriched(this) ) {
        PlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( numberOfIntegrationRules > 1 ) {
            int i;
            PatchIntegrationRule *iRule;
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );
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
    if ( !xf->isElementEnriched(this) ) {
        PlaneStress2d :: drawScalar(context);
    } else {
        if ( context.giveIntVarMode() == ISM_local ) {
            int indx, ans, result = 1;
            double val;
            FloatArray s(3), v;
            IntArray map;

            ans = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
            if ( ( !ans ) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
                return;
            }


            TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
            PatchIntegrationRule *iRule;
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( int j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
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
