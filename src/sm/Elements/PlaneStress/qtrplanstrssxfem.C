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

#include "qtrplanstrssxfem.h"

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
#include "domain.h"
#include "dynamicinputrecord.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "xfem/patchintegrationrule.h"
#endif



#include <string>
#include <sstream>

namespace oofem {
REGISTER_Element(QTrPlaneStress2dXFEM);

QTrPlaneStress2dXFEM::~QTrPlaneStress2dXFEM() {
	// TODO Auto-generated destructor stub
}

Interface *
QTrPlaneStress2dXFEM :: giveInterface(InterfaceType it)
{
    if ( it == XfemElementInterfaceType ) {
        return static_cast< XfemElementInterface * >(this);
    } else if ( it == VTKXMLExportModuleElementInterfaceType ) {
        return static_cast< VTKXMLExportModuleElementInterface * >(this);
    } else {
        return QTrPlaneStress2d :: giveInterface(it);
    }
}

void QTrPlaneStress2dXFEM :: updateYourself(TimeStep *tStep)
{
    QTrPlaneStress2d :: updateYourself(tStep);
    XfemStructuralElementInterface :: updateYourselfCZ(tStep);
}

void QTrPlaneStress2dXFEM :: postInitialize()
{
    QTrPlaneStress2d :: postInitialize();
    XfemStructuralElementInterface :: initializeCZMaterial();
}

int QTrPlaneStress2dXFEM :: computeNumberOfDofs()
{
    int ndofs = 0;

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        ndofs += this->giveDofManager(i)->giveNumberOfDofs();
    }

    return ndofs;
}

void QTrPlaneStress2dXFEM :: computeGaussPoints()
{
    if ( giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = giveDomain()->giveXfemManager();

        if ( xMan->isElementEnriched(this) ) {
            if ( !this->XfemElementInterface_updateIntegrationRule() ) {
                QTrPlaneStress2d :: computeGaussPoints();
            }
        } else {
            QTrPlaneStress2d :: computeGaussPoints();
        }
    } else   {
        QTrPlaneStress2d :: computeGaussPoints();
    }
}

void QTrPlaneStress2dXFEM :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this, false);
}

void QTrPlaneStress2dXFEM :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    XfemElementInterface_createEnrBmatrixAt(answer, * gp, * this);
}

void QTrPlaneStress2dXFEM :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    XfemElementInterface_createEnrBHmatrixAt(answer, * gp, * this);
}

void
QTrPlaneStress2dXFEM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Continuous part
    QTrPlaneStress2d :: giveDofManDofIDMask(inode, answer);

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
QTrPlaneStress2dXFEM :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
QTrPlaneStress2dXFEM :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void QTrPlaneStress2dXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    QTrPlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
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

void QTrPlaneStress2dXFEM :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface::XfemElementInterface_computeDeformationGradientVector(answer, gp, tStep);
}

void
QTrPlaneStress2dXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    QTrPlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemStructuralElementInterface :: computeCohesiveForces(answer, tStep);
}

Element_Geometry_Type
QTrPlaneStress2dXFEM :: giveGeometryType() const
{
    if ( this->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        if ( xMan->isElementEnriched(this) ) {
            return EGT_Composite;
        } else {
            return EGT_Composite;
        }
    } else   {
        return EGT_triangle_2;
    }
}

void
QTrPlaneStress2dXFEM :: initializeFrom(InputRecord &ir)
{
    QTrPlaneStress2d :: initializeFrom(ir);
    XfemStructuralElementInterface :: initializeCZFrom(ir);

    if ( ir.hasField(_IFT_QTrPlaneStress2dXFEM_RegCoeff) ) {
        ir.giveOptionalField(mRegCoeff, _IFT_QTrPlaneStress2dXFEM_RegCoeff);
    }

    if ( ir.hasField(_IFT_QTrPlaneStress2dXFEM_RegCoeffTol) ) {
        ir.giveOptionalField(mRegCoeffTol, _IFT_QTrPlaneStress2dXFEM_RegCoeffTol);
    }
}

MaterialMode QTrPlaneStress2dXFEM :: giveMaterialMode()
{
    return XfemStructuralElementInterface :: giveMaterialMode();
}

void QTrPlaneStress2dXFEM :: giveInputRecord(DynamicInputRecord &input)
{
    QTrPlaneStress2d :: giveInputRecord(input);
    XfemStructuralElementInterface :: giveCZInputRecord(input);

    input.setField(mRegCoeff, _IFT_QTrPlaneStress2dXFEM_RegCoeff);
    input.setField(mRegCoeffTol, _IFT_QTrPlaneStress2dXFEM_RegCoeffTol);

}

void
QTrPlaneStress2dXFEM :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    // TODO: Validate implementation.

    FloatArray u;
    FloatMatrix n;

    XfemElementInterface_createEnrNmatrixAt(n, lcoords, * this, false);

    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}

void
QTrPlaneStress2dXFEM :: giveCompositeExportData(std::vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    vtkPieces.resize(1);

    const int numCells = mSubTri.size();
    if(numCells == 0) {
        // Enriched but uncut element
        // Visualize as a triangle
        vtkPieces[0].setNumberOfCells(1);

        int numTotalNodes = 6;
        vtkPieces[0].setNumberOfNodes(numTotalNodes);

        // Node coordinates
        std :: vector< FloatArray >nodeCoords;
        for(int i = 1; i <= 6; i++) {
            const auto &x = giveDofManager(i)->giveCoordinates();
            nodeCoords.push_back(x);

            vtkPieces[0].setNodeCoords(i, x);
        }

        // Connectivity
        IntArray nodes1 = {1, 2, 3, 4, 5, 6};
        vtkPieces[0].setConnectivity(1, nodes1);

        // Offset
        int offset = 6;
        vtkPieces[0].setOffset(1, offset);

        // Cell types
        vtkPieces[0].setCellType(1, 22); // Quadratic triangle




        // Export nodal variables from primary fields
        vtkPieces[0].setNumberOfPrimaryVarsToExport(primaryVarsToExport, numTotalNodes);

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


} /* namespace OOFEM */
