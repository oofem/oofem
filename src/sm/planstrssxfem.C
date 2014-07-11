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
#include "xfem/xfemelementinterface.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentitem.h"
#include "xfem/enrichmentdomain.h"
#include "structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(PlaneStress2dXfem)

void PlaneStress2dXfem :: updateYourself(TimeStep *tStep)
{
    PlaneStress2d :: updateYourself(tStep);
    XfemStructuralElementInterface :: updateYourselfCZ(tStep);
}

void PlaneStress2dXfem :: postInitialize()
{
    PlaneStress2d :: postInitialize();
    XfemStructuralElementInterface :: initializeCZMaterial();
}

Interface *
PlaneStress2dXfem :: giveInterface(InterfaceType it)
{
    if ( it == XfemElementInterfaceType ) {
        return static_cast< XfemElementInterface * >(this);
    } else if ( it == VTKXMLExportModuleElementInterfaceType ) {
        return static_cast< VTKXMLExportModuleElementInterface * >(this);
    } else {
        return PlaneStress2d :: giveInterface(it);
    }
}


void
PlaneStress2dXfem :: computeGaussPoints()
{
    if ( this->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();

        if ( xMan->isElementEnriched(this) ) {
            if ( !this->XfemElementInterface_updateIntegrationRule() ) {
                PlaneStress2d :: computeGaussPoints();
            }
        } else {
            PlaneStress2d :: computeGaussPoints();
        }
    } else   {
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
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this, false);
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
PlaneStress2dXfem :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Continuous part
    PlaneStress2d :: giveDofManDofIDMask(inode, answer);

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



void PlaneStress2dXfem :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void PlaneStress2dXfem :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    PlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    XfemStructuralElementInterface :: computeCohesiveTangent(answer, tStep);
}

void
PlaneStress2dXfem :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    PlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemStructuralElementInterface :: computeCohesiveForces(answer, tStep);
}

Element_Geometry_Type
PlaneStress2dXfem :: giveGeometryType() const
{
    if ( this->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        if ( xMan->isElementEnriched(this) ) {
            return EGT_Composite;
            //return EGT_quad_1;
        } else {
            return EGT_Composite;
            //return EGT_quad_1;
        }
    } else   {
        return EGT_quad_1;
    }
}



#ifdef __OOFEG
void PlaneStress2dXfem :: drawRawGeometry(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }
#if 0
    XfemManager *xf = this->giveDomain()->giveXfemManager();
    if ( !xf->isElementEnriched(this) ) {
        PlaneStress2d :: drawRawGeometry(context);
    } else {
        if ( integrationRulesArray.size() > 1 ) {
            // TODO: Implement visualization
            /*
             *          PatchIntegrationRule *iRule;
             *          for ( auto &iRule: integrationRulesArray ) {
             *              iRule = dynamic_cast< PatchIntegrationRule * >( ir );
             *              if ( iRule ) {
             *                  iRule->givePatch()->draw(context);
             *              }
             *          }
             */
        } else {
            PlaneStress2d :: drawRawGeometry(context);
        }
    }
#endif
}

void PlaneStress2dXfem :: drawScalar(oofegGraphicContext &context)
{
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }
#if 0
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
            for ( auto &ir: integrationRulesArray ) {
                PatchIntegrationRule *iRule = dynamic_cast< PatchIntegrationRule * >(ir);

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
                // TODO: Implement visualization
                //                iRule->givePatch()->drawWD(context, s);
            }
        } else {
            PlaneStress2d :: drawScalar(context);
        }
    }
#endif
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

    result = XfemStructuralElementInterface :: initializeCZFrom(ir);
    return result;
}

MaterialMode PlaneStress2dXfem :: giveMaterialMode()
{
    return XfemStructuralElementInterface :: giveMaterialMode();
}

void PlaneStress2dXfem :: giveInputRecord(DynamicInputRecord &input)
{
    PlaneStress2d :: giveInputRecord(input);
    XfemStructuralElementInterface :: giveCZInputRecord(input);
}

void
PlaneStress2dXfem :: giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    const int numCells = mSubTri.size();

    if(numCells == 0) {
        // Enriched but uncut element
        // Visualize as a quad
        vtkPiece.setNumberOfCells(1);

        int numTotalNodes = 4;
        vtkPiece.setNumberOfNodes(numTotalNodes);

        // Node coordinates
        std :: vector< FloatArray >nodeCoords;
        for(int i = 1; i <= 4; i++) {
            FloatArray &x = *(giveDofManager(i)->giveCoordinates());
            nodeCoords.push_back(x);

            vtkPiece.setNodeCoords(i, x);
        }

        // Connectivity
        IntArray nodes1 = {1, 2, 3, 4};
        vtkPiece.setConnectivity(1, nodes1);

        // Offset
        int offset = 4;
        vtkPiece.setOffset(1, offset);

        // Cell types
        vtkPiece.setCellType(1, 9); // Linear quad




        // Export nodal variables from primary fields
        vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

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

                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeInd, u);
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
        vtkPiece.setNumberOfInternalVarsToExport(0, numTotalNodes);


        // Export cell variables
        vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), 1);
        for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
            FloatArray average;
            IntegrationRule *iRule = integrationRulesArray [ 0 ];
            VTKXMLExportModule :: computeIPAverage(average, iRule, this, type, tStep);

            FloatArray averageV9(9);
            averageV9.at(1) = average.at(1);
            averageV9.at(5) = average.at(2);
            averageV9.at(9) = average.at(3);
            averageV9.at(6) = averageV9.at(8) = average.at(4);
            averageV9.at(3) = averageV9.at(7) = average.at(5);
            averageV9.at(2) = averageV9.at(4) = average.at(6);

            vtkPiece.setCellVar( i, 1, averageV9 );
        }


        // Export of XFEM related quantities
        if ( domain->hasXfemManager() ) {
            XfemManager *xMan = domain->giveXfemManager();

            int nEnrIt = xMan->giveNumberOfEnrichmentItems();
            vtkPiece.setNumberOfInternalXFEMVarsToExport(xMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);

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
                                ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                            double levelSet = 0.0, levelSetInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
                            double nodeEnrMarker = 0.0, nodeEnrMarkerInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalNodeEnrMarkerInNode(nodeEnrMarkerInNode, dMan->giveGlobalNumber() );

                                nodeEnrMarker += N.at(elNodeInd)*nodeEnrMarkerInNode;
                            }


                            FloatArray valueArray = {nodeEnrMarker};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);
                        }

                    }
                }
            }
        }

    }
    else {
        // Enriched and cut element

        vtkPiece.setNumberOfCells(numCells);

        int numTotalNodes = numCells*3;
        vtkPiece.setNumberOfNodes(numTotalNodes);

        // Node coordinates
        std :: vector< FloatArray >nodeCoords;
        int nodesPassed = 1;
        for(Triangle tri: mSubTri ) {
            for(int i = 1; i <= 3; i++) {
                FloatArray x = tri.giveVertex(i);
                nodeCoords.push_back(x);
                vtkPiece.setNodeCoords(nodesPassed, x);
                nodesPassed++;
            }
        }

        // Connectivity, offset and cell type
        nodesPassed = 1;
        int offset = 3;
        for(size_t i = 1; i <= mSubTri.size(); i++ ) {
            IntArray nodes = {nodesPassed, nodesPassed+1, nodesPassed+2};
            nodesPassed += 3;
            vtkPiece.setConnectivity(i, nodes);

            vtkPiece.setOffset(i, offset);
            offset += 3;

            vtkPiece.setCellType(i, 5); // Linear triangle
        }



        // Export nodal variables from primary fields
        vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

        for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
            UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);


            nodesPassed = 1;

            for(Triangle tri: mSubTri ) {

                FloatArray triCenter = tri.giveVertex(1);
                triCenter.add( tri.giveVertex(2) );
                triCenter.add( tri.giveVertex(3) );
                triCenter.times(1.0/3.0);

                double meanEdgeLength = 0.0;
                meanEdgeLength += (1.0/3.0)*tri.giveVertex(1).distance( tri.giveVertex(2) );
                meanEdgeLength += (1.0/3.0)*tri.giveVertex(2).distance( tri.giveVertex(3) );
                meanEdgeLength += (1.0/3.0)*tri.giveVertex(3).distance( tri.giveVertex(1) );

                const double relPertLength = 1.0e-5;

                for(int i = 1; i <= 3; i++) {

                    if ( type == DisplacementVector ) { // compute displacement
                        FloatArray u = {0.0, 0.0, 0.0};


                        // Fetch global coordinates (in undeformed configuration)
                        const FloatArray &x = tri.giveVertex(i);

                        FloatArray locCoordNode;
                        computeLocalCoordinates(locCoordNode, x);


                        FloatArray pertVec;
                        FloatArray locCoord;
                        double pertLength = relPertLength;

                        int numTries = 1;
                        for(int j = 0; j < numTries; j++) {
                            // Perturb towards triangle center, to ensure that
                            // we end up on the correct side of the crack
                            pertVec.beDifferenceOf(triCenter, x);
                            pertVec.times(pertLength);
                            FloatArray xPert = x;
                            xPert.add(pertVec);


                            // Compute local coordinates
                            computeLocalCoordinates(locCoord, xPert);

                        }

                        // Compute displacement in point
                        FloatMatrix NMatrix;
                        FloatArray solVec;

                        // Use only regular basis functions for edge nodes to get linear interpolation
                        // (to make the visualization look nice)
                        FloatArray N;
                        FEInterpolation *interp = giveInterpolation();
                        interp->evalN( N, locCoordNode, FEIElementGeometryWrapper(this) );
                        const int nDofMan = giveNumberOfDofManagers();

                        XfemManager *xMan = domain->giveXfemManager();
                        int numEI =  xMan->giveNumberOfEnrichmentItems();

                        bool joinNodes = true;

                        for(int eiIndex = 1; eiIndex <= numEI; eiIndex++) {
                            EnrichmentItem *ei = xMan->giveEnrichmentItem(eiIndex);

                            double levelSetTang = 0.0, levelSetNormal = 0.0, levelSetInNode = 0.0;

                            bool evaluationSucceeded = true;
                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);

                                if(!ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber() )) {
                                    evaluationSucceeded = false;
                                }
                                levelSetTang += N.at(elNodeInd)*levelSetInNode;

                                if(!ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber() )) {
                                    evaluationSucceeded = false;
                                }
                                levelSetNormal += N.at(elNodeInd)*levelSetInNode;

                            }

                            double tangSignDist = 0.0, arcPos = 0.0;
                            ei->giveEnrichmentDomain()->computeTangentialSignDist(tangSignDist, x, arcPos);

                            if(tangSignDist > (1.0e-3)*meanEdgeLength && fabs(levelSetNormal) < (1.0e-2)*meanEdgeLength && evaluationSucceeded) {
                                joinNodes = false;
                            }

                        }

                        if(joinNodes) {
                            // if point on edge
                            XfemElementInterface_createEnrNmatrixAt(NMatrix, locCoord, * this, true);
                            computeVectorOf(VM_Total, tStep, solVec);
                        }
                        else {
                            XfemElementInterface_createEnrNmatrixAt(NMatrix, locCoord, * this, false);
                            computeVectorOf(VM_Total, tStep, solVec);
                        }


                        FloatArray uTemp;
                        uTemp.beProductOf(NMatrix, solVec);

                        if(uTemp.giveSize() == 3) {
                            u = uTemp;
                        }
                        else {
                            u = {uTemp[0], uTemp[1], 0.0};
                        }


                        FloatArray valuearray = u;
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodesPassed, valuearray);
                    } else {
                        // TODO: Implement
                        printf("fieldNum: %d\n", fieldNum);

                    }

                    nodesPassed++;
                }
            }


         }


        // Export nodal variables from internal fields
        vtkPiece.setNumberOfInternalVarsToExport(0, numTotalNodes);


        vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
        for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
            InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);
            FloatArray average;
            IntegrationRule *iRule = integrationRulesArray [ 0 ];
            VTKXMLExportModule :: computeIPAverage(average, iRule, this, type, tStep);

            FloatArray averageV9(9);
            averageV9.at(1) = average.at(1);
            averageV9.at(5) = average.at(2);
            averageV9.at(9) = average.at(3);
            averageV9.at(6) = averageV9.at(8) = average.at(4);
            averageV9.at(3) = averageV9.at(7) = average.at(5);
            averageV9.at(2) = averageV9.at(4) = average.at(6);

            for(size_t triInd = 1; triInd <= mSubTri.size(); triInd++) {
                vtkPiece.setCellVar( i, triInd, averageV9 );
            }
        }



        // Export of XFEM related quantities
        if ( domain->hasXfemManager() ) {
            XfemManager *xMan = domain->giveXfemManager();

            int nEnrIt = xMan->giveNumberOfEnrichmentItems();
            vtkPiece.setNumberOfInternalXFEMVarsToExport(xMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);

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
                                ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                            double levelSet = 0.0, levelSetInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber() );

                                levelSet += N.at(elNodeInd)*levelSetInNode;
                            }


                            FloatArray valueArray = {levelSet};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        } else if ( xfemstype == XFEMST_NodeEnrMarker ) {

                            double nodeEnrMarker = 0.0, nodeEnrMarkerInNode = 0.0;

                            for(int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++) {
                                DofManager *dMan = giveDofManager(elNodeInd);
                                ei->evalNodeEnrMarkerInNode(nodeEnrMarkerInNode, dMan->giveGlobalNumber() );

                                nodeEnrMarker += N.at(elNodeInd)*nodeEnrMarkerInNode;
                            }



                            FloatArray valueArray = {nodeEnrMarker};
                            vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);

                        }

                    }
                }
            }
        }



    }
}

} // end namespace oofem
