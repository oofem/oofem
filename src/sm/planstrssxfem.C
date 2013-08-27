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

REGISTER_Element( PlaneStress2dXfem );

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
    XfemManager *xMan = this->giveDomain()->giveXfemManager();

    for(int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++)
    {
    	std::vector<FloatArray> intersecPoints;
    	EnrichmentItem *ei = xMan->giveEnrichmentItem(i);

        std::vector< int > intersecEdgeInd;
    	ei->computeIntersectionPoints(intersecPoints, intersecEdgeInd, this);
    	int numIntersecPoints = intersecPoints.size();

        if ( numIntersecPoints > 0 )
        {
            this->XfemElementInterface_updateIntegrationRule();
        } else
        {
            PlaneStress2d ::computeGaussPoints();
        }

    }

}

///@todo: Computation of N and B can be moved to the xfem interface and thus handle all continuum elements in the same way
void PlaneStress2dXfem :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
	XfemElementInterface_createEnrBmatrixAt(answer, *gp, *this);
}


void PlaneStress2dXfem :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
	OOFEM_ERROR("PlaneStress2dXfem :: computeNmatrixAt is not yet implemented.");
/*
    FloatArray Nc;
    interpolation.evalN( Nc, *gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    // assemble xfem part of strain-displacement matrix
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
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
        this->computeGlobalCoordinates(coords, *gp->giveCoordinates());
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
*/
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
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    int nEI = xMan->giveNumberOfEnrichmentItems();
    StructuralMaterial *sm = NULL;

    for(int i = 1; i <= nEI; i++)
    {
    	EnrichmentItem &ei = *(xMan->giveEnrichmentItem(i));
    	if( ei.isMaterialModified(*gp, *this, sm) )
    	{

    		if(sm != NULL)
    		{
    	        sm->giveStiffnessMatrix(answer, rMode, gp, tStep);
    			return;
    		}
    		else
    		{
    			OOFEM_ERROR("PlaneStress2dXfem :: computeConstitutiveMatrixAt: failed to fetch StructuralMaterial\n");
    		}
    	}

    }


    // If no enrichment modifies the material,
    // compute stiffness based on the bulk material.
    PlaneStress2d :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *stepN)
{
//    FloatArray Epsilon;
//    this->computeStrainVector(Epsilon, gp, stepN);


    //////////////////
    // Necessary for postprocessing
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( this->giveCrossSection() );
    cs->giveRealStresses(answer, gp, strain, stepN);
    //////////////////


    XfemManager *xMan = this->giveDomain()->giveXfemManager();

    int nEI = xMan->giveNumberOfEnrichmentItems();

    StructuralMaterial *sm = NULL;
    for(int i = 1; i <= nEI; i++)
    {
    	EnrichmentItem &ei = *(xMan->giveEnrichmentItem(i));
    	if( ei.isMaterialModified(*gp, *this, sm) )
    	{
//    	    printf("In PlaneStress2dXfem :: computeStressVector: Epsilon: "); Epsilon.printYourself();

    		if(sm != NULL)
    		{
    	        sm->giveRealStressVector(answer, gp, strain, stepN);
    			return;
    		}
    		else
    		{
    			OOFEM_ERROR("PlaneStress2dXfem :: computeStressVector: failed to fetch StructuralMaterial\n");
    		}
    	}

    }


    // If no enrichment modifies the material,
    // compute stress based on the bulk material.
//    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( this->giveCrossSection() );
//    cs->giveRealStresses(answer, gp, Epsilon, stepN);
}
/*
void PlaneStress2dXfem :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    this->computeStiffnessMatrix_withIRulesAsSubcells(answer, rMode, tStep);
}

void
PlaneStress2dXfem :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) 
{
    this->giveInternalForcesVector_withIRulesAsSubcells(answer, tStep, useUpdatedGpRecord);
}
*/
Element_Geometry_Type 
PlaneStress2dXfem :: giveGeometryType() const
{ 
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
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

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        PlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( numberOfIntegrationRules > 1 ) {
            PatchIntegrationRule *iRule;
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
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

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        PlaneStress2d :: drawScalar(context);
    } else {
        if ( context.giveIntVarMode() == ISM_local ) {
            int indx;
            double val;
            FloatArray s(3), v;

            indx = context.giveIntVarIndx();

            TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
            PatchIntegrationRule *iRule;
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(0);
                    giveIPValue(v, gp, context.giveIntVarType(), tStep);
                    val += v.at(indx);
                }

                val /= iRule->giveNumberOfIntegrationPoints();
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
