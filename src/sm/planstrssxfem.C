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
#include "structuralmaterial.h"
#include "xfemelementinterface.h"
#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "enrichmentdomain.h"
#include "structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#ifdef __OOFEG
#include "patchintegrationrule.h"
#endif

namespace oofem {

Interface *
PlaneStress2dXfem :: giveInterface(InterfaceType it)
{
    if ( it == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else if ( it == VTKXMLExportModuleElementInterfaceType ) {
        return ( VTKXMLExportModuleElementInterface * ) this;
    } else {
        return PlaneStress2d :: giveInterface(it);
    }
}


void
PlaneStress2dXfem :: computeGaussPoints()
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    EnrichmentItem *ei = xMan->giveEnrichmentItem(1);
    EnrichmentDomain *ed = ei->giveEnrichmentDomain(1);
    Inclusion *iei = dynamic_cast< Inclusion * > (ei); 
    EDBGCircle *edc = dynamic_cast< EDBGCircle *> (ed);
    // update integration rule only if element is cut by the ed
    if ( iei && edc && edc->computeNumberOfIntersectionPoints(this) > 0) { 
        this->XfemElementInterface_updateIntegrationRule();
    } else {
        PlaneStress2d ::computeGaussPoints();
    }
}

///@todo: Computation of N and B can be moved to the xfem interface and thus handle all continuum elements in the same way
void PlaneStress2dXfem :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    FloatMatrix dNdx;
    FloatArray N;
    interpolation.evaldNdx( dNdx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    interpolation.evalN(     N  , * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
 
    FloatMatrix Bc[4];
    // Assemble standard FEM part of strain-displacement matrix
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        FloatMatrix &BNode = Bc[i-1];
        BNode.resize(3, 2);
        BNode.zero();
        BNode.at(1, 1) = dNdx.at(i, 1);
        BNode.at(2, 2) = dNdx.at(i, 2);
        BNode.at(3, 1) = dNdx.at(i, 2);
        BNode.at(3, 2) = dNdx.at(i, 1);
    }

    // Assemble xfem part of strain-displacement matrix
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    FloatMatrix Bd[4];

    int counter = 8;
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        EnrichmentDomain *ed = ei->giveEnrichmentDomain(1);
        // Enrichment function and its gradient evaluated at the gauss point
        EnrichmentFunction *ef = ei->giveEnrichmentFunction(1);
        double efgp = ef->evaluateFunctionAt(gp, ed);
        FloatArray efgpD;
        ef->evaluateDerivativeAt(efgpD, gp, ed);

        // adds up the number of the dofs from an enrichment item
        // this part is used for the construction of a shifted enrichment
        for ( int j = 1; j <= this->giveNumberOfDofManagers(); j++ ) {

            DofManager *dMan = this->giveDofManager(j);
            if ( ei->isDofManEnriched( dMan ) ) {
                FloatMatrix &BdNode = Bd[j-1];
                BdNode.resize(3, 2);
                BdNode.zero();

                FloatArray *nodecoords = dMan->giveCoordinates();
                // should ask after specific EF in a loop
                double efnode = ef->evaluateFunctionAt(nodecoords, ed);

                // matrix to be added anytime a node is enriched
                // Creates nabla*(ef*N)
                FloatArray help;
                help.resize(2);
                for ( int p = 1; p <= 2; p++ ) {
                    help.at(p) = dNdx.at(j, p) * ( efgp - efnode ) + N.at(j) * efgpD.at(p);
                }

                for ( int k = 1; k <= 2; k++ ) {
                    BdNode.at(k, k) = help.at(k);
                    if ( k == 1 ) {
                        BdNode.at(3, k) = help.at(2);
                    } else if ( k == 2 ) {
                        BdNode.at(3, k) = help.at(1);
                    }
                }
                counter += 2;
            }
        }

        // Create the total B-matrix by appending each contribution to B after one another.
        answer.resize(3, counter);
        answer.zero();
        int column = 1;
        for ( int i = 0; i < 4; i++ ) {
            answer.setSubMatrix(Bc[i],1,column);
            column += 2;
            if ( Bd[i].isNotEmpty() ) {
                answer.setSubMatrix(Bd[i],1,column);
                column += 2;
            }
        }
    }
}


void PlaneStress2dXfem :: computeNmatrixAt(FloatArray &lcoords, FloatMatrix &answer)
{

    FloatArray Nc;
    interpolation.evalN( Nc, lcoords, FEIElementGeometryWrapper(this) );
    // assemble xfem part of strain-displacement matrix
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    FloatArray Nd;
    Nd.resize(4);
    IntArray mask(4);
    int counter = 4;
    FloatArray N, coords;
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        EnrichmentDomain *ed = ei->giveEnrichmentDomain(1);    

        // Enrichment function and its gradient evaluated at the gauss point     
        EnrichmentFunction *ef = ei->giveEnrichmentFunction(1);
        this->computeGlobalCoordinates(coords,lcoords);
        double efgp = ef->evaluateFunctionAt(&coords, ed);

        // adds up the number of the dofs from an enrichment item
        // this part is used for the construction of a shifted enrichment
        for ( int j = 1; j <= this->giveNumberOfDofManagers(); j++ ) {
            DofManager *dMan = this->giveDofManager(j);
            if ( ei->isDofManEnriched( dMan ) ) {
                
                FloatArray *nodecoords = dMan->giveCoordinates();
                double efnode = ef->evaluateFunctionAt(nodecoords, ed);
                Nd.at(j) = ( efgp - efnode ) * Nc.at(j) ;
                
                counter++;
                mask.at(j) = 1;
            } else {
                mask.at(j) = 0;
            }
        }

        // Create the total B-matrix by appending each contribution to B after one another.
        N.resize(counter);
        int column = 1;

        for ( int i = 1; i <= 4; i++ ) {
            N.at(column) = Nc.at(i);
            column ++;
            if ( mask.at(i) ) {
                N.at(column) = Nd.at(i);
                column++;
            }
        }
    }
    answer.beNMatrixOf(N,2);
}



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
    // Returns the total id mask of the dof manager = regular id's + enriched id's
    this->giveDofManager(inode)->giveCompleteMasterDofIDArray(answer);
}



void PlaneStress2dXfem :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    ///@todo: only works for circles
    EDBGCircle *edc = static_cast< EDBGCircle * > ( xMan->giveEnrichmentItem(1)->giveEnrichmentDomain(1) );
    
    FloatArray coords;
    this->computeGlobalCoordinates(coords, *(gp->giveCoordinates()) );
    if ( edc->bg->isInside( coords ) ) {
        Inclusion *ei = static_cast< Inclusion * > ( xMan->giveEnrichmentItem(1) );    
        StructuralMaterial *sm = static_cast< StructuralMaterial * >( ei->giveMaterial() );
        sm->giveCharacteristicMatrix(answer, ReducedForm, rMode, gp, tStep);
    } else {
        PlaneStress2d :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
    }
}

void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
    FloatArray Epsilon;
    this->computeStrainVector(Epsilon, gp, stepN);
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    EDBGCircle *edc = dynamic_cast< EDBGCircle * > ( xMan->giveEnrichmentItem(1)->giveEnrichmentDomain(1) );
    FloatArray coords;
    this->computeGlobalCoordinates(coords, *(gp->giveCoordinates()) );
    if ( edc->bg->isInside( coords ) ) {
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
PlaneStress2dXfem :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) 
{
    this->giveInternalForcesVector_withIRulesAsSubcells(answer, tStep, useUpdatedGpRecord);
}

Element_Geometry_Type 
PlaneStress2dXfem :: giveGeometryType() const 
{ 
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    if ( xMan->isElementEnriched(this) ) {
        //return EGT_Composite;
        return EGT_quad_1; 
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

    XfemManager *xf = this->giveDomain()->giveXfemManager(1);
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

    XfemManager *xf = this->giveDomain()->giveXfemManager(1);
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
