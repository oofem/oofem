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

#include "planstrssxfem.h"
#include "structuralmaterial.h"
#include "xfemelementinterface.h"
#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "enrichmentdomain.h"
#include "structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#include "dynamicinputrecord.h"
#ifdef __OOFEG
 #include "patchintegrationrule.h"
#endif

namespace oofem {
REGISTER_Element(PlaneStress2dXfem)

void PlaneStress2dXfem :: updateYourself(TimeStep *tStep)
{
    PlaneStress2d :: updateYourself(tStep);
    XfemElementInterface :: updateYourselfCZ(tStep);
}

void PlaneStress2dXfem :: postInitialize()
{
    PlaneStress2d :: postInitialize();
    XfemElementInterface :: initializeCZMaterial();
}

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
	if( this->giveDomain()->hasXfemManager() ) {
		XfemManager *xMan = this->giveDomain()->giveXfemManager();

		if ( xMan->isElementEnriched(this) ) {
			if ( !this->XfemElementInterface_updateIntegrationRule() ) {
				PlaneStress2d :: computeGaussPoints();
			}
		} else   {
			PlaneStress2d :: computeGaussPoints();
		}
	}
	else {
		PlaneStress2d :: computeGaussPoints();
	}
}

void PlaneStress2dXfem :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    XfemElementInterface_createEnrBmatrixAt(answer, * gp, * this);
}

void PlaneStress2dXfem :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    XfemElementInterface_createEnrBHmatrixAt(answer, * gp, * this);
}


void PlaneStress2dXfem :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this);
}



int PlaneStress2dXfem :: computeNumberOfDofs()
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
    XfemElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void PlaneStress2dXfem :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    PlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    XfemElementInterface :: computeCohesiveTangent(answer, tStep);
}

void
PlaneStress2dXfem :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    PlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemElementInterface :: computeCohesiveForces(answer, tStep);
}

Element_Geometry_Type
PlaneStress2dXfem :: giveGeometryType() const
{
	if( this->giveDomain()->hasXfemManager() ) {
		XfemManager *xMan = this->giveDomain()->giveXfemManager();
		if ( xMan->isElementEnriched(this) ) {
			//return EGT_Composite;
			return EGT_quad_1;
		} else {
			return EGT_quad_1;
		}
	}
	else {
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
            // TODO: Implement visualization
            /*
             *          PatchIntegrationRule *iRule;
             *          for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
             *              iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );
             *              if ( iRule ) {
             *                  iRule->givePatch()->draw(context);
             *              }
             *          }
             */
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
                // TODO: Implement visualization
                //                iRule->givePatch()->drawWD(context, s);
            }
        } else {
            PlaneStress2d :: drawScalar(context);
        }
    }
}
#endif


IRResultType
PlaneStress2dXfem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    result = PlaneStress2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = XfemElementInterface :: initializeCZFrom(ir);
    return result;
}

MaterialMode PlaneStress2dXfem :: giveMaterialMode()
{
    return XfemElementInterface :: giveMaterialMode();
}

void PlaneStress2dXfem :: giveInputRecord(DynamicInputRecord &input)
{
    PlaneStress2d :: giveInputRecord(input);
    XfemElementInterface :: giveCZInputRecord(input);
}
} // end namespace oofem
