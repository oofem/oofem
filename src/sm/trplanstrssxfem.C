/*
 * trplanstrssxfem.C
 *
 *  Created on: Jun 3, 2013
 *      Author: svennine
 */

#include "trplanstrssxfem.h"

#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "rcm2.h"
#endif


#include "structuralmaterial.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentitem.h"
#include "xfem/enrichmentdomain.h"
#include "structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#ifdef __OOFEG
 #include "xfem/patchintegrationrule.h"
#endif
#include "xfem/delaunay.h"

#include "xfem/XFEMDebugTools.h"
#include <string>
#include <sstream>


namespace oofem {
REGISTER_Element(TrPlaneStress2dXFEM);

void TrPlaneStress2dXFEM :: updateYourself(TimeStep *tStep)
{
    TrPlaneStress2d :: updateYourself(tStep);
    XfemElementInterface :: updateYourselfCZ(tStep);
}

void TrPlaneStress2dXFEM :: postInitialize()
{
    TrPlaneStress2d :: postInitialize();
    XfemElementInterface :: initializeCZMaterial();
}


TrPlaneStress2dXFEM :: ~TrPlaneStress2dXFEM() { }

int TrPlaneStress2dXFEM :: checkConsistency()
{
    TrPlaneStress2d :: checkConsistency();
    return 1;
}

Interface *
TrPlaneStress2dXFEM :: giveInterface(InterfaceType it)
{
    if ( it == XfemElementInterfaceType ) {
        return static_cast< XfemElementInterface * >(this);
    } else if ( it == VTKXMLExportModuleElementInterfaceType ) {
        return static_cast< VTKXMLExportModuleElementInterface * >(this);
    } else {
        return TrPlaneStress2d :: giveInterface(it);
    }
}

int TrPlaneStress2dXFEM :: computeNumberOfDofs()
{
    int ndofs = 0;

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        ndofs += this->giveDofManager(i)->giveNumberOfDofs();
    }

    return ndofs;
}

void TrPlaneStress2dXFEM :: computeGaussPoints()
{
    if ( giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = giveDomain()->giveXfemManager();

        if ( xMan->isElementEnriched(this) ) {
            if ( !this->XfemElementInterface_updateIntegrationRule() ) {
                TrPlaneStress2d :: computeGaussPoints();
            }
        } else {
            TrPlaneStress2d :: computeGaussPoints();
        }
    } else   {
        TrPlaneStress2d :: computeGaussPoints();
    }
}

void TrPlaneStress2dXFEM :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    XfemElementInterface_createEnrBmatrixAt(answer, * gp, * this);
}

void TrPlaneStress2dXFEM :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    XfemElementInterface_createEnrBHmatrixAt(answer, * gp, * this);
}

void TrPlaneStress2dXFEM :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this);
}



void
TrPlaneStress2dXFEM :: giveDofManDofIDMask(int inode, EquationID iEqnId, IntArray &answer) const
{
    // Continuous part
	TrPlaneStress2d::giveDofManDofIDMask(inode, iEqnId, answer);

    // Discontinuous part
	if( this->giveDomain()->hasXfemManager() ) {
		DofManager *dMan = giveDofManager(inode);
		XfemManager *xMan = giveDomain()->giveXfemManager();

        const std::vector<int> &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices( dMan->giveGlobalNumber() );
        for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(nodeEiIndices[i]);
			if ( ei->isDofManEnriched(* dMan) ) {
				IntArray eiDofIdArray;
				ei->computeEnrichedDofManDofIdArray(eiDofIdArray, *dMan);
				answer.followedBy(eiDofIdArray);
			}
		}
	}
}

void
TrPlaneStress2dXFEM :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
TrPlaneStress2dXFEM :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void TrPlaneStress2dXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    TrPlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    XfemElementInterface :: computeCohesiveTangent(answer, tStep);
}

void
TrPlaneStress2dXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    TrPlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemElementInterface :: computeCohesiveForces(answer, tStep);
}


#ifdef __OOFEG
// TODO: FIX OOFEG implementation
void TrPlaneStress2dXFEM :: drawRawGeometry(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        TrPlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( numberOfIntegrationRules > 1 ) {
            int i;
            //            PatchIntegrationRule *iRule;
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                // TODO: Implement visualization.
                /*
                 *              iRule = dynamic_cast< PatchIntegrationRule * >( integrationRulesArray [ i ] );
                 *              if ( iRule ) {
                 *                  iRule->givePatch()->draw(context);
                 *              }
                 */
            }
        } else {
            TrPlaneStress2d :: drawRawGeometry(context);
        }
    }
}

void TrPlaneStress2dXFEM :: drawScalar(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        TrPlaneStress2d :: drawScalar(context);
    } else {
        if ( context.giveIntVarMode() == ISM_local ) {
            int indx;
            double val;
            FloatArray s(3), v;

            indx = context.giveIntVarIndx();

            TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
            PatchIntegrationRule *iRule;
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                iRule = dynamic_cast< PatchIntegrationRule * >(integrationRulesArray [ i ]);

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( GaussPoint *gp: *iRule ) {
                    giveIPValue(v, gp, context.giveIntVarType(), tStep);
                    val += v.at(indx);
                }

                val /= iRule->giveNumberOfIntegrationPoints();
 #endif
                s.at(1) = s.at(2) = s.at(3) = val;
                // TODO: Implement visualization.
                //                iRule->givePatch()->drawWD(context, s);
            }
        } else {
            TrPlaneStress2d :: drawScalar(context);
        }
    }
}
#endif

IRResultType
TrPlaneStress2dXFEM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    result = TrPlaneStress2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = XfemElementInterface :: initializeCZFrom(ir);
    return result;
}

MaterialMode TrPlaneStress2dXFEM :: giveMaterialMode()
{
    return XfemElementInterface :: giveMaterialMode();
}

void TrPlaneStress2dXFEM :: giveInputRecord(DynamicInputRecord &input)
{
    TrPlaneStress2d :: giveInputRecord(input);
    XfemElementInterface :: giveCZInputRecord(input);
}


int
TrPlaneStress2dXFEM :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                        TimeStep *tStep, const FloatArray &coords,
                                                                        FloatArray &answer)
{
    // TODO: Validate implementation.

    FloatArray lcoords, u;
    FloatMatrix n;
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    XfemElementInterface_createEnrNmatrixAt(n, lcoords, * this);

    this->computeVectorOf(EID_MomentumBalance, mode, tStep, u);
    answer.beProductOf(n, u);

    return result;
}

void
TrPlaneStress2dXFEM :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    //    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
    // TODO: For now, take only the continuous part
    int nodeInd = 1;
    TrPlaneStress2d :: giveDofManDofIDMask(nodeInd, EID_MomentumBalance, answer);
}
} /* namespace oofem */
