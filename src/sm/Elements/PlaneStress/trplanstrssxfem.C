/*
 * trplanstrssxfem.C
 *
 *  Created on: Jun 3, 2013
 *      Author: svennine
 */

#include "sm/Elements/PlaneStress/trplanstrssxfem.h"

#include "sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentitem.h"
#include "xfem/delaunay.h"
#include "xfem/XFEMDebugTools.h"
#include "feinterpol.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "parametermanager.h"
#include "paramkey.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "xfem/patchintegrationrule.h"
#endif



#include <string>
#include <sstream>


namespace oofem {
REGISTER_Element(TrPlaneStress2dXFEM);
ParamKey TrPlaneStress2dXFEM::IPK_TrPlaneStress2dXFEM_RegCoeff("reg_coeff");
ParamKey TrPlaneStress2dXFEM::IPK_TrPlaneStress2dXFEM_RegCoeffTol("reg_coeff_tol");

void TrPlaneStress2dXFEM :: updateYourself(TimeStep *tStep)
{
    TrPlaneStress2d :: updateYourself(tStep);
    XfemStructuralElementInterface :: updateYourselfCZ(tStep);
}

void TrPlaneStress2dXFEM :: postInitialize()
{
    TrPlaneStress2d :: postInitialize();
    XfemStructuralElementInterface :: initializeCZMaterial();
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
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this, false);
}



void
TrPlaneStress2dXFEM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Continuous part
    TrPlaneStress2d :: giveDofManDofIDMask(inode, answer);

    // Discontinuous part
    if( this->giveDomain()->hasXfemManager() ) {
        DofManager *dMan = giveDofManager(inode);
        XfemManager *xMan = giveDomain()->giveXfemManager();

        int placeInArray = domain->giveDofManPlaceInArray(dMan->giveGlobalNumber());
        const std::vector<int> &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices( placeInArray );
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
    XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
TrPlaneStress2dXFEM :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void TrPlaneStress2dXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    TrPlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    XfemStructuralElementInterface :: computeCohesiveTangent(answer, tStep);


    const double tol = mRegCoeffTol;
    const double regularizationCoeff = mRegCoeff;
    int numRows = answer.giveNumberOfRows();
    for(int i = 0; i < numRows; i++) {
        if( fabs(answer(i,i)) < tol ) {
            answer(i,i) += regularizationCoeff;
//          printf("Found zero on diagonal.\n");
        }
    }
}

void TrPlaneStress2dXFEM :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface::XfemElementInterface_computeDeformationGradientVector(answer, gp, tStep);
}

void
TrPlaneStress2dXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    TrPlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemStructuralElementInterface :: computeCohesiveForces(answer, tStep);


#if 0
    ////////////////////////////////////////////
    // Contribution from regularization
    FloatMatrix K;
	computeStiffnessMatrix(K, TangentStiffness, tStep);

    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    const double tol = mRegCoeffTol+mRegCoeff;
    const double regularizationCoeff = mRegCoeff;
    int numRows = K.giveNumberOfRows();
    for(int i = 0; i < numRows; i++) {
        if( fabs(K(i,i)) < tol ) {
            answer(i) += regularizationCoeff*u(i);
//          printf("Found zero on diagonal.\n");
        }
    }
#endif
}


Element_Geometry_Type
TrPlaneStress2dXFEM :: giveGeometryType() const
{
    if ( this->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        if ( xMan->isElementEnriched(this) ) {
            return EGT_Composite;
        } else {
            return EGT_Composite;
        }
    } else   {
        return EGT_triangle_1;
    }
}

#ifdef __OOFEG
// TODO: FIX OOFEG implementation
void TrPlaneStress2dXFEM :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        TrPlaneStress2d :: drawRawGeometry(gc, tStep);
    } else {
        if ( integrationRulesArray.size() > 1 ) {
#if 0
            for ( auto &ir: integrationRulesArray ) {
                // TODO: Implement visualization.
                PatchIntegrationRule *iRule = dynamic_cast< PatchIntegrationRule * >( ir );
                if ( iRule ) {
                    iRule->givePatch()->draw(gc);
                }
            }
#endif
        } else {
            TrPlaneStress2d :: drawRawGeometry(gc, tStep);
        }
    }
}

void TrPlaneStress2dXFEM :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        TrPlaneStress2d :: drawScalar(gc, tStep);
    } else {
        if ( gc.giveIntVarMode() == ISM_local ) {
            int indx;
            double val;
            FloatArray s(3), v;

            indx = gc.giveIntVarIndx();

            for ( auto &ir: integrationRulesArray ) {
                PatchIntegrationRule *iRule = dynamic_cast< PatchIntegrationRule * >(ir.get());

 #if 0
                val = iRule->giveMaterial();
 #else
                val = 0.0;
                for ( GaussPoint *gp: *iRule ) {
                    giveIPValue(v, gp, gc.giveIntVarType(), tStep);
                    val += v.at(indx);
                }

                val /= iRule->giveNumberOfIntegrationPoints();
 #endif
                s.at(1) = s.at(2) = s.at(3) = val;
                // TODO: Implement visualization.
                //                iRule->givePatch()->drawWD(gc, s);
            }
        } else {
            TrPlaneStress2d :: drawScalar(gc, tStep);
        }
    }
}
#endif

void
TrPlaneStress2dXFEM :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm = giveDomain()->elementPPM;
    TrPlaneStress2d :: initializeFrom(ir, priority);
    XfemStructuralElementInterface :: initializeCZFrom(ir, priority);

    PM_UPDATE_PARAMETER(mRegCoeff, ppm, ir, this->giveNumber(), IPK_TrPlaneStress2dXFEM_RegCoeff, priority);
    PM_UPDATE_PARAMETER(mRegCoeffTol, ppm, ir, this->giveNumber(), IPK_TrPlaneStress2dXFEM_RegCoeffTol, priority);
}

MaterialMode TrPlaneStress2dXFEM :: giveMaterialMode()
{
    return XfemStructuralElementInterface :: giveMaterialMode();
}

void TrPlaneStress2dXFEM :: giveInputRecord(DynamicInputRecord &input)
{
    TrPlaneStress2d :: giveInputRecord(input);
    XfemStructuralElementInterface :: giveCZInputRecord(input);

    input.setField(mRegCoeff, _IFT_TrPlaneStress2dXFEM_RegCoeff);
    input.setField(mRegCoeffTol, _IFT_TrPlaneStress2dXFEM_RegCoeffTol);

}


void
TrPlaneStress2dXFEM :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    // TODO: Validate implementation.

    FloatArray u;
    FloatMatrix n;

    XfemElementInterface_createEnrNmatrixAt(n, lcoords, * this, false);

    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}


void
TrPlaneStress2dXFEM :: giveElementDofIDMask(IntArray &answer) const
{
    // TODO: For now, take only the continuous part
    TrPlaneStress2d :: giveElementDofIDMask(answer);
}

void
TrPlaneStress2dXFEM :: giveCompositeExportData(std::vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    vtkPieces.resize(1);

    const int numCells = mSubTri.size();
    if(numCells == 0) {
        // Enriched but uncut element
        // Visualize as a triangle
        vtkPieces[0].setNumberOfCells(1);

        int numTotalNodes = 3;
        vtkPieces[0].setNumberOfNodes(numTotalNodes);

        // Node coordinates
        std :: vector< FloatArray >nodeCoords;
        for(int i = 1; i <= 3; i++) {
            const auto &x = giveDofManager(i)->giveCoordinates();
            nodeCoords.push_back(x);

            vtkPieces[0].setNodeCoords(i, x);
        }

        // Connectivity
        IntArray nodes1 = {1, 2, 3};
        vtkPieces[0].setConnectivity(1, nodes1);

        // Offset
        int offset = 3;
        vtkPieces[0].setOffset(1, offset);

        // Cell types
        vtkPieces[0].setCellType(1, 5); // Linear triangle




        // Export nodal variables from primary fields
        vtkPieces[0].setNumberOfPrimaryVarsToExport(primaryVarsToExport,  numTotalNodes);

        for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
            UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);

            for ( int nodeInd = 1; nodeInd <= numTotalNodes; nodeInd++ ) {

                if ( type == DisplacementVector ) { // compute displacement

                        FloatArray u = {0.0, 0.0, 0.0};

                        // Fetch global coordinates (in undeformed configuration)
                        const FloatArray &x = nodeCoords[nodeInd-1];

                        // Compute local coordinates
                        FloatArray locCoord;
                        computeLocalCoordinates(locCoord, x);

                        // Compute displacement in point
                        FloatMatrix NMatrix;
                        computeNmatrixAt(locCoord, NMatrix);
                        FloatArray solVec;
                        computeVectorOf(VM_Total, tStep, solVec);
                        FloatArray uTemp;
                        uTemp.beProductOf(NMatrix, solVec);

                        if(uTemp.giveSize() == 3) {
                            u = uTemp;
                        }
                        else {
                            u = {uTemp[0], uTemp[1], 0.0};
                        }

                        vtkPieces[0].setPrimaryVarInNode(type, nodeInd, u);
                } else {
                    printf("fieldNum: %d\n", fieldNum);
                    // TODO: Implement
//                    ZZNodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
//                    for ( int j = 1; j <= numCellNodes; j++ ) {
//                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values [ j - 1 ]);
//                        nodeNum += 1;
//                    }
                }
            }
        }


        // Export nodal variables from internal fields
        vtkPieces[0].setNumberOfInternalVarsToExport(internalVarsToExport, numTotalNodes);


        // Export cell variables
        vtkPieces[0].setNumberOfCellVarsToExport(cellVarsToExport, 1);
        for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
            FloatArray average;
            std :: unique_ptr< IntegrationRule > &iRule = integrationRulesArray [ 0 ];
            VTKXMLExportModule :: computeIPAverage(average, iRule.get(), this, type, tStep);
            if(average.giveSize() == 0) {
            	average = {0., 0., 0., 0., 0., 0.};
            }

            if(average.giveSize() == 1) {
				vtkPieces[0].setCellVar( type, 1, average );
            }

            if(average.giveSize() == 6) {
				FloatArray averageV9(9);
				averageV9.at(1) = average.at(1);
				averageV9.at(5) = average.at(2);
				averageV9.at(9) = average.at(3);
				averageV9.at(6) = averageV9.at(8) = average.at(4);
				averageV9.at(3) = averageV9.at(7) = average.at(5);
				averageV9.at(2) = averageV9.at(4) = average.at(6);

				vtkPieces[0].setCellVar( type, 1, averageV9 );
            }
        }


        // Export of XFEM related quantities
        if ( domain->hasXfemManager() ) {
            XfemManager *xMan = domain->giveXfemManager();

            int nEnrIt = xMan->giveNumberOfEnrichmentItems();
            vtkPieces[0].setNumberOfInternalXFEMVarsToExport(xMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);

            const int nDofMan = giveNumberOfDofManagers();


            for ( int field = 1; field <= xMan->vtkExportFields.giveSize(); field++ ) {
                XFEMStateType xfemstype = ( XFEMStateType ) xMan->vtkExportFields [ field - 1 ];

                for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
                    EnrichmentItem *ei = xMan->giveEnrichmentItem(enrItIndex);
                    for ( int nodeInd = 1; nodeInd <= numTotalNodes; nodeInd++ ) {

                        const FloatArray &x = nodeCoords[nodeInd-1];
                        FloatArray locCoord;
                        computeLocalCoordinates(locCoord, x);

                        FloatArray N;
                        FEInterpolation *interp = giveInterpolation();
                        interp->evalN( N, locCoord, FEIElementGeometryWrapper(this) );


                        if ( xfemstype == XFEMST_LevelSetPhi ) {
                            double levelSet = 0.0, levelSetInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber(), dMan->giveCoordinates() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPieces[0].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                            double levelSet = 0.0, levelSetInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber(), dMan->giveCoordinates() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPieces[0].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
                            double nodeEnrMarker = 0.0, nodeEnrMarkerInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalNodeEnrMarkerInNode(nodeEnrMarkerInNode, dMan->giveGlobalNumber() );

                                nodeEnrMarker += N.at(elNodeInd)*nodeEnrMarkerInNode;
                            }


                            FloatArray valueArray = {nodeEnrMarker};
                            vtkPieces[0].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);
                        }

                    }
                }
            }
        }

    }
    else {
        // Enriched and cut element

        XfemStructuralElementInterface::giveSubtriangulationCompositeExportData(vtkPieces, primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep);


    }

}

} /* namespace oofem */
