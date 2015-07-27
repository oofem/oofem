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

#include "xfemstructuralelementinterface.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "engngm.h"
#include "Elements/nlstructuralelement.h"

#include "xfem/patchintegrationrule.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/XFEMDebugTools.h"
#include "xfem/xfemtolerances.h"

#include "xfem/enrichmentfronts/enrichmentfrontintersection.h"

#include "vtkxmlexportmodule.h"

#include <string>
#include <sstream>

namespace oofem {
XfemStructuralElementInterface :: XfemStructuralElementInterface(Element *e) :
    XfemElementInterface(e),
    mpCZMat(NULL),
    mCZMaterialNum(-1),
    mCSNumGaussPoints(4)
{}

XfemStructuralElementInterface :: ~XfemStructuralElementInterface() {}

bool XfemStructuralElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    const double tol2 = 1.0e-18;

    bool partitionSucceeded = false;


    if ( mpCZMat != NULL ) {
        mpCZIntegrationRules.clear();
        mCZEnrItemIndices.clear();
        mCZTouchingEnrItemIndices.clear();
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        if ( mpCZMat == NULL && mCZMaterialNum > 0 ) {
            initializeCZMaterial();
        }


        MaterialMode matMode = element->giveMaterialMode();

        bool firstIntersection = true;

        std :: vector< std :: vector< FloatArray > >pointPartitions;
        mSubTri.clear();

        std :: vector< int >enrichingEIs;
        int elPlaceInArray = xMan->giveDomain()->giveElementPlaceInArray( element->giveGlobalNumber() );
        xMan->giveElementEnrichmentItemIndices(enrichingEIs, elPlaceInArray);


        for ( size_t p = 0; p < enrichingEIs.size(); p++ ) {
            // Index of current ei
            int eiIndex = enrichingEIs [ p ];

            // Indices of other ei interaction with this ei through intersection enrichment fronts.
            std :: vector< int >touchingEiIndices;
            giveIntersectionsTouchingCrack(touchingEiIndices, enrichingEIs, eiIndex, * xMan);

            if ( firstIntersection ) {
                // Get the points describing each subdivision of the element
                double startXi, endXi;
                bool intersection = false;
                this->XfemElementInterface_prepareNodesForDelaunay(pointPartitions, startXi, endXi, eiIndex, intersection);

                if ( intersection ) {
                    firstIntersection = false;

                    // Use XfemElementInterface_partitionElement to subdivide the element
                    for ( int i = 0; i < int ( pointPartitions.size() ); i++ ) {
                        // Triangulate the subdivisions
                        this->XfemElementInterface_partitionElement(mSubTri, pointPartitions [ i ]);
                    }


                    if ( mpCZMat != NULL ) {
                        Crack *crack = dynamic_cast< Crack * >( xMan->giveEnrichmentItem(eiIndex) );
                        if ( crack == NULL ) {
                            OOFEM_ERROR("Cohesive zones are only available for cracks.")
                        }

                        // We have xi_s and xi_e. Fetch sub polygon.
                        std :: vector< FloatArray >crackPolygon;
                        crack->giveSubPolygon(crackPolygon, startXi, endXi);

                        ///////////////////////////////////
                        // Add cohesive zone Gauss points
                        size_t numSeg = crackPolygon.size() - 1;

                        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                            int czRuleNum = 1;
                            mpCZIntegrationRules.emplace_back( new GaussIntegrationRule(czRuleNum, element) );

                            // Add index of current ei
                            mCZEnrItemIndices.push_back(eiIndex);

                            // Add indices of other ei, that cause interaction through
                            // intersection enrichment fronts
                            mCZTouchingEnrItemIndices.push_back(touchingEiIndices);

                            // Compute crack normal
                            FloatArray crackTang;
                            crackTang.beDifferenceOf(crackPolygon [ segIndex + 1 ], crackPolygon [ segIndex ]);

                            if ( crackTang.computeSquaredNorm() > tol2 ) {
                                crackTang.normalize();
                            }

                            FloatArray crackNormal = {
                                -crackTang.at(2), crackTang.at(1)
                            };

                            mpCZIntegrationRules [ segIndex ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                           crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);

                            for ( GaussPoint *gp: *mpCZIntegrationRules [ segIndex ] ) {
                                double gw = gp->giveWeight();
                                double segLength = crackPolygon [ segIndex ].distance(crackPolygon [ segIndex + 1 ]);
                                gw *= 0.5 * segLength;
                                gp->setWeight(gw);

                                // Fetch material status and set normal
                                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                                if ( ms == NULL ) {
                                    OOFEM_ERROR("Failed to fetch material status.");
                                }

                                ms->letNormalBe(crackNormal);

                                // Give Gauss point reference to the enrichment item
                                // to simplify post processing.
                                crack->AppendCohesiveZoneGaussPoint(gp);
                            }
                        }
                    }



                    partitionSucceeded = true;
                }
            } // if(firstIntersection)
            else {
                // Loop over triangles
                std :: vector< Triangle >allTriCopy;
                for ( size_t triIndex = 0; triIndex < mSubTri.size(); triIndex++ ) {
                    // Call alternative version of XfemElementInterface_prepareNodesForDelaunay
                    std :: vector< std :: vector< FloatArray > >pointPartitionsTri;
                    double startXi, endXi;
                    bool intersection = false;
                    XfemElementInterface_prepareNodesForDelaunay(pointPartitionsTri, startXi, endXi, mSubTri [ triIndex ], eiIndex, intersection);

                    if ( intersection ) {
                        // Use XfemElementInterface_partitionElement to subdivide triangle j
                        for ( int i = 0; i < int ( pointPartitionsTri.size() ); i++ ) {
                            this->XfemElementInterface_partitionElement(allTriCopy, pointPartitionsTri [ i ]);
                        }


                        // Add cohesive zone Gauss points

                        if ( mpCZMat != NULL ) {
                            Crack *crack = dynamic_cast< Crack * >( xMan->giveEnrichmentItem(eiIndex) );
                            if ( crack == NULL ) {
                                OOFEM_ERROR("Cohesive zones are only available for cracks.")
                            }

                            // We have xi_s and xi_e. Fetch sub polygon.
                            std :: vector< FloatArray >crackPolygon;
                            crack->giveSubPolygon(crackPolygon, startXi, endXi);

                            int numSeg = crackPolygon.size() - 1;

                            for ( int segIndex = 0; segIndex < numSeg; segIndex++ ) {
                                int czRuleNum = 1;
                                mpCZIntegrationRules.emplace_back( new GaussIntegrationRule(czRuleNum, element) );
                                size_t newRuleInd = mpCZIntegrationRules.size() - 1;
                                mCZEnrItemIndices.push_back(eiIndex);

                                mCZTouchingEnrItemIndices.push_back(touchingEiIndices);

                                // Compute crack normal
                                FloatArray crackTang;
                                crackTang.beDifferenceOf(crackPolygon [ segIndex + 1 ], crackPolygon [ segIndex ]);

                                if ( crackTang.computeSquaredNorm() > tol2 ) {
                                    crackTang.normalize();
                                }

                                FloatArray crackNormal = {
                                    -crackTang.at(2), crackTang.at(1)
                                };

                                mpCZIntegrationRules [ newRuleInd ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                                 crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);

                                for ( GaussPoint *gp: *mpCZIntegrationRules [ newRuleInd ] ) {
                                    double gw = gp->giveWeight();
                                    double segLength = crackPolygon [ segIndex ].distance(crackPolygon [ segIndex + 1 ]);
                                    gw *= 0.5 * segLength;
                                    gp->setWeight(gw);

                                    // Fetch material status and set normal
                                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                                    if ( ms == NULL ) {
                                        OOFEM_ERROR("Failed to fetch material status.");
                                    }

                                    ms->letNormalBe(crackNormal);

                                    // Give Gauss point reference to the enrichment item
                                    // to simplify post processing.
                                    crack->AppendCohesiveZoneGaussPoint(gp);
                                }
                            }
                        }
                    } else {
                        allTriCopy.push_back(mSubTri [ triIndex ]);
                    }
                }

                mSubTri = allTriCopy;
            }
        }

        // Refine triangles if desired
        int numRefs = xMan->giveNumTriRefs();

        for(int i = 0; i < numRefs; i++) {

            std :: vector< Triangle > triRef;

            for(const Triangle &tri : mSubTri) {
                Triangle::refineTriangle(triRef, tri);
            }

            mSubTri = triRef;
        }

        ////////////////////////////////////////
        // When we reach this point, we have a
        // triangulation that is adapted to all
        // cracks passing through the element.
        // Therefore, we can set up integration
        // points on each triangle.

        if ( xMan->giveVtkDebug() ) {
            std :: stringstream str3;
            int elIndex = this->element->giveGlobalNumber();
            str3 << "TriEl" << elIndex << ".vtk";
            std :: string name3 = str3.str();

            if ( mSubTri.size() > 0 ) {
                XFEMDebugTools :: WriteTrianglesToVTK(name3, mSubTri);
            }
        }


        int ruleNum = 1;

        if ( partitionSucceeded ) {
            std :: vector< std :: unique_ptr< IntegrationRule > >intRule;
            intRule.emplace_back( new PatchIntegrationRule(ruleNum, element, mSubTri) );
            intRule [ 0 ]->SetUpPointsOnTriangle(xMan->giveNumGpPerTri(), matMode);
            element->setIntegrationRules( std :: move(intRule) );
        }


        if ( xMan->giveVtkDebug() ) {
            ////////////////////////////////////////////////////////////////////////
            // Write CZ GP to VTK

            std :: vector< FloatArray >czGPCoord;

            for ( size_t czRulInd = 0; czRulInd < mpCZIntegrationRules.size(); czRulInd++ ) {
                for ( GaussPoint *gp: *mpCZIntegrationRules [ czRulInd ] ) {
                    czGPCoord.push_back( gp->giveGlobalCoordinates() );
                }
            }

            double time = 0.0;

            Domain *dom = element->giveDomain();
            if ( dom != NULL ) {
                EngngModel *em = dom->giveEngngModel();
                if ( em != NULL ) {
                    TimeStep *ts = em->giveCurrentStep();
                    if ( ts != NULL ) {
                        time = ts->giveTargetTime();
                    }
                }
            }

            std :: stringstream str;
            int elIndex = this->element->giveGlobalNumber();
            str << "CZGaussPointsTime" << time << "El" << elIndex << ".vtk";
            std :: string name = str.str();

            XFEMDebugTools :: WritePointsToVTK(name, czGPCoord);
            ////////////////////////////////////////////////////////////////////////
        }
    }

    return partitionSucceeded;
}

void XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = element->giveDomain()->giveXfemManager();
        CrossSection *cs = NULL;

        const std :: vector< int > &materialModifyingEnrItemIndices = xMan->giveMaterialModifyingEnrItemIndices();
        for ( size_t i = 0; i < materialModifyingEnrItemIndices.size(); i++ ) {
            EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(materialModifyingEnrItemIndices [ i ]) );

            if ( ei.isMaterialModified(* gp, * element, cs) ) {
                StructuralCrossSection *structCS = dynamic_cast< StructuralCrossSection * >( cs );

                if ( structCS != NULL ) {
                    if ( mUsePlaneStrain ) {
                        structCS->giveStiffnessMatrix_PlaneStrain(answer, rMode, gp, tStep);
                    } else {
                        structCS->giveStiffnessMatrix_PlaneStress(answer, rMode, gp, tStep);
                    }
                    return;
                } else {
                    OOFEM_ERROR("failed to fetch StructuralMaterial");
                }
            }
        }
    }

    // If no enrichment modifies the material,
    // compute stiffness based on the bulk material.
    StructuralCrossSection *cs = dynamic_cast< StructuralCrossSection * >( element->giveCrossSection() );
    if ( mUsePlaneStrain ) {
        cs->giveStiffnessMatrix_PlaneStrain(answer, rMode, gp, tStep);
    } else {
        cs->giveStiffnessMatrix_PlaneStress(answer, rMode, gp, tStep);
    }
}

void XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = dynamic_cast< StructuralCrossSection * >( element->giveCrossSection() );
    if ( cs == NULL ) {
        OOFEM_ERROR("cs == NULL.");
    }

    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = element->giveDomain()->giveXfemManager();


        CrossSection *csInclusion = NULL;
        const std :: vector< int > &materialModifyingEnrItemIndices = xMan->giveMaterialModifyingEnrItemIndices();
        for ( size_t i = 0; i < materialModifyingEnrItemIndices.size(); i++ ) {
            EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(materialModifyingEnrItemIndices [ i ]) );

            if ( ei.isMaterialModified(* gp, * element, csInclusion) ) {
                StructuralCrossSection *structCSInclusion = dynamic_cast< StructuralCrossSection * >( csInclusion );

                if ( structCSInclusion != NULL ) {
                    if ( mUsePlaneStrain ) {
                        structCSInclusion->giveRealStress_PlaneStrain(answer, gp, strain, tStep);
                    } else {
                        structCSInclusion->giveRealStress_PlaneStress(answer, gp, strain, tStep);
                    }

                    return;
                } else {
                    OOFEM_ERROR("failed to fetch StructuralCrossSection");
                }
            }
        }
    }

    // If no enrichment modifies the material:
    if ( mUsePlaneStrain ) {
        cs->giveRealStress_PlaneStrain(answer, gp, strain, tStep);
    } else {
        cs->giveRealStress_PlaneStress(answer, gp, strain, tStep);
    }
}

void XfemStructuralElementInterface :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        FloatArray solVec;
        element->computeVectorOf(VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();
        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            for ( GaussPoint *gp: *mpCZIntegrationRules [ segIndex ] ) {
                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ], mCZTouchingEnrItemIndices [ segIndex ]);
                ////////////////////////////////////////////////////////


                // Traction
                FloatArray T2D;



                // Fetch material status and get normal
                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                if ( ms == NULL ) {
                    OOFEM_ERROR("Failed to fetch material status.");
                }

                ms->setNewlyInserted(false); //TODO: Do this in a better place. /ES

                FloatArray crackNormal( ms->giveNormal() );

                // Compute jump vector
                FloatArray jump2D;
                computeDisplacementJump(* gp, jump2D, solVec, NMatrix);

                computeGlobalCohesiveTractionVector(T2D, jump2D, crackNormal, NMatrix, * gp, tStep);

                // Add to internal force
                FloatArray NTimesT;

                NTimesT.beTProductOf(NMatrix, T2D);
                CrossSection *cs  = element->giveCrossSection();
                double thickness = cs->give(CS_Thickness, gp);
                double dA = thickness * gp->giveWeight();
                answer.add(dA, NTimesT);
            }
        }
    }
}

void XfemStructuralElementInterface :: computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep)
{
    FloatMatrix F;
    F.resize(3, 3);
    F.beUnitMatrix();         // TODO: Compute properly


    FloatArray jump3D = {
        iJump.at(1), iJump.at(2), 0.0
    };


    FloatArray crackNormal3D = {
        iCrackNormal.at(1), iCrackNormal.at(2), 0.0
    };

    FloatArray ez = {
        0.0, 0.0, 1.0
    };
    FloatArray crackTangent3D;
    crackTangent3D.beVectorProductOf(crackNormal3D, ez);

    FloatMatrix locToGlob(3, 3);
    locToGlob.setColumn(crackTangent3D, 1);
    locToGlob.setColumn(crackNormal3D, 2);
    locToGlob.setColumn(ez, 3);

    FloatArray TLoc, jump3DLoc, TLocRenumbered(3);
    jump3DLoc.beTProductOf(locToGlob, jump3D);

    FloatArray jump3DLocRenumbered = {
        jump3DLoc.at(3), jump3DLoc.at(1), jump3DLoc.at(2)
    };

    mpCZMat->giveFirstPKTraction_3d(TLocRenumbered, & iGP, jump3DLocRenumbered, F, tStep);

    TLoc = {
        TLocRenumbered.at(2), TLocRenumbered.at(3), TLocRenumbered.at(1)
    };


    FloatArray T;
    T.beProductOf(locToGlob, TLoc);

    oT = {
        T.at(1), T.at(2)
    };
}

void XfemStructuralElementInterface :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        FloatArray solVec;
        element->computeVectorOf(VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();

        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            for ( GaussPoint *gp: *mpCZIntegrationRules [ segIndex ] ) {
                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ], mCZTouchingEnrItemIndices [ segIndex ]);

                ////////////////////////////////////////////////////////
                // Compute jump vector
                FloatArray jump2D;
                computeDisplacementJump(* gp, jump2D, solVec, NMatrix);

                FloatArray jump3D = {
                    0.0, jump2D.at(1), jump2D.at(2)
                };

                // Compute traction
                FloatMatrix F;
                F.resize(3, 3);
                F.beUnitMatrix();                     // TODO: Compute properly

                FloatMatrix K3DRenumbered, K3DGlob;


                FloatMatrix K2D;
                K2D.resize(2, 2);
                K2D.zero();

                if ( mpCZMat->hasAnalyticalTangentStiffness() ) {
                    ///////////////////////////////////////////////////
                    // Analytical tangent

                    FloatMatrix K3D;
                    mpCZMat->give3dStiffnessMatrix_dTdj(K3DRenumbered, TangentStiffness, gp, tStep);

                    K3D.resize(3, 3);
                    K3D.zero();
                    K3D.at(1, 1) = K3DRenumbered.at(2, 2);
                    K3D.at(1, 2) = K3DRenumbered.at(2, 3);
                    K3D.at(1, 3) = K3DRenumbered.at(2, 1);

                    K3D.at(2, 1) = K3DRenumbered.at(3, 2);
                    K3D.at(2, 2) = K3DRenumbered.at(3, 3);
                    K3D.at(2, 3) = K3DRenumbered.at(3, 1);

                    K3D.at(3, 1) = K3DRenumbered.at(1, 2);
                    K3D.at(3, 2) = K3DRenumbered.at(1, 3);
                    K3D.at(3, 3) = K3DRenumbered.at(1, 1);


                    // Fetch material status and get normal
                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                    if ( ms == NULL ) {
                        OOFEM_ERROR("Failed to fetch material status.");
                    }

                    FloatArray crackNormal( ms->giveNormal() );

                    FloatArray crackNormal3D = {
                        crackNormal.at(1), crackNormal.at(2), 0.0
                    };

                    FloatArray ez = {
                        0.0, 0.0, 1.0
                    };
                    FloatArray crackTangent3D;
                    crackTangent3D.beVectorProductOf(crackNormal3D, ez);

                    FloatMatrix locToGlob(3, 3);
                    locToGlob.setColumn(crackTangent3D, 1);
                    locToGlob.setColumn(crackNormal3D, 2);
                    locToGlob.setColumn(ez, 3);


                    FloatMatrix tmp3(3, 3);
                    tmp3.beProductTOf(K3D, locToGlob);
                    K3DGlob.beProductOf(locToGlob, tmp3);

                    K2D.at(1, 1) = K3DGlob.at(1, 1);
                    K2D.at(1, 2) = K3DGlob.at(1, 2);
                    K2D.at(2, 1) = K3DGlob.at(2, 1);
                    K2D.at(2, 2) = K3DGlob.at(2, 2);
                } else {
                    ///////////////////////////////////////////////////
                    // Numerical tangent
                    double eps = 1.0e-9;

                    FloatArray T, TPert;

                    // Fetch material status and get normal
                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                    if ( ms == NULL ) {
                        OOFEM_ERROR("Failed to fetch material status.");
                    }

                    FloatArray crackNormal( ms->giveNormal() );

                    computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, * gp, tStep);


                    FloatArray jump2DPert;


                    jump2DPert = jump2D;
                    jump2DPert.at(1) += eps;
                    computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, * gp, tStep);

                    K2D.at(1, 1) = ( TPert.at(1) - T.at(1) ) / eps;
                    K2D.at(2, 1) = ( TPert.at(2) - T.at(2) ) / eps;

                    jump2DPert = jump2D;
                    jump2DPert.at(2) += eps;
                    computeGlobalCohesiveTractionVector(TPert, jump2DPert, crackNormal, NMatrix, * gp, tStep);

                    K2D.at(1, 2) = ( TPert.at(1) - T.at(1) ) / eps;
                    K2D.at(2, 2) = ( TPert.at(2) - T.at(2) ) / eps;

                    computeGlobalCohesiveTractionVector(T, jump2D, crackNormal, NMatrix, * gp, tStep);
                }

                FloatMatrix tmp, tmp2;
                tmp.beProductOf(K2D, NMatrix);
                tmp2.beTProductOf(NMatrix, tmp);

                CrossSection *cs  = element->giveCrossSection();
                double thickness = cs->give(CS_Thickness, gp);
                double dA = thickness * gp->giveWeight();
                answer.add(dA, tmp2);
            }
        }
    }
}

void XfemStructuralElementInterface :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        printf("Entering XfemElementInterface :: computeCohesiveTangentAt().\n");
    }
}

void XfemStructuralElementInterface :: XfemElementInterface_computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    StructuralElement *structEl = dynamic_cast< StructuralElement * >( element );
    if ( structEl == NULL ) {
        OOFEM_ERROR("Not a structural element");
    }

    int ndofs = structEl->computeNumberOfDofs();
    double density, dV;
    FloatMatrix n;
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !structEl->isActivated(tStep) ) {
        return;
    }

    structEl->giveMassMtrxIntegrationgMask(mask);

    mass = 0.;

    for ( GaussPoint *gp: *element->giveIntegrationRule(0) ) {
        structEl->computeNmatrixAt(gp->giveNaturalCoordinates(), n);
        density = structEl->giveStructuralCrossSection()->give('d', gp);

        if ( ipDensity != NULL ) {
            // Override density if desired
            density = * ipDensity;
        }

        dV = structEl->computeVolumeAround(gp);
        mass += density * dV;

        if ( mask.isEmpty() ) {
            answer.plusProductSymmUpper(n, n, density * dV);
        } else {
            for ( int i = 1; i <= ndofs; i++ ) {
                for ( int j = i; j <= ndofs; j++ ) {
                    double summ = 0.;
                    for ( int k = 1; k <= n.giveNumberOfRows(); k++ ) {
                        if ( mask.at(k) == 0 ) {
                            continue;
                        }

                        summ += n.at(k, i) * n.at(k, j);
                    }

                    answer.at(i, j) += summ * density * dV;
                }
            }
        }
    }

    answer.symmetrized();

    const double tol = 1.0e-9;
    const double regularizationCoeff = 1.0e-6;
    int numRows = answer.giveNumberOfRows();
    for ( int i = 0; i < numRows; i++ ) {
        if ( fabs( answer(i, i) ) < tol ) {
            answer(i, i) += regularizationCoeff;
            //          printf("Found zero on diagonal.\n");
        }
    }
}

IRResultType
XfemStructuralElementInterface :: initializeCZFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    int material = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_XfemElementInterface_CohesiveZoneMaterial);
    mCZMaterialNum = material;
    //    printf("In XfemElementInterface :: initializeCZFrom(): mCZMaterialNum: %d\n", mCZMaterialNum );


    // Number of Gauss points used when integrating the cohesive zone
    IR_GIVE_OPTIONAL_FIELD(ir, mCSNumGaussPoints, _IFT_XfemElementInterface_NumIntPointsCZ);
    //    printf("mCSNumGaussPoints: %d\n", mCSNumGaussPoints );

    int planeStrainFlag = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, planeStrainFlag, _IFT_XfemElementInterface_PlaneStrain);
    if ( planeStrainFlag == 1 ) {
        mUsePlaneStrain = true;
    }

    return IRRT_OK;
}

void XfemStructuralElementInterface :: giveCZInputRecord(DynamicInputRecord &input)
{
    if ( mCZMaterialNum > 0 ) {
        input.setField(mCZMaterialNum, _IFT_XfemElementInterface_CohesiveZoneMaterial);
    }

    if ( mUsePlaneStrain ) {
        input.setField(1, _IFT_XfemElementInterface_PlaneStrain);
    }

    input.setField(mCSNumGaussPoints, _IFT_XfemElementInterface_NumIntPointsCZ);
}

void XfemStructuralElementInterface :: initializeCZMaterial()
{
    if ( mCZMaterialNum > 0 ) {
        mpCZMat = dynamic_cast< StructuralInterfaceMaterial * >( this->element->giveDomain()->giveMaterial(mCZMaterialNum) );

        if ( mpCZMat == NULL ) {
            OOFEM_ERROR("Failed to fetch pointer for mpCZMat.");
        }
    }
}

void XfemStructuralElementInterface :: XfemElementInterface_computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    NLStructuralElement *nlStructEl = static_cast<NLStructuralElement*>( element );

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    nlStructEl->computeVectorOf(VM_Total, tStep, u); // solution vector
    if ( nlStructEl->initialDisplacements ) {
        u.subtract(* nlStructEl->initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    nlStructEl->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}

void XfemStructuralElementInterface :: giveIntersectionsTouchingCrack(std :: vector< int > &oTouchingEnrItemIndices, const std :: vector< int > &iCandidateIndices, int iEnrItemIndex, XfemManager &iXMan)
{
    EnrichmentItem *ei = iXMan.giveEnrichmentItem(iEnrItemIndex);


    for ( int candidateIndex : iCandidateIndices ) {
        if ( candidateIndex != iEnrItemIndex ) {
            // Fetch candidate enrichment item
            EnrichmentItem *eiCandidate = iXMan.giveEnrichmentItem(candidateIndex);

            // This treatment is only necessary if the enrichment front
            // is an EnrFrontIntersection. Therefore, start by trying a
            // dynamic cast.

            // Check start tip
            EnrFrontIntersection *efStart = dynamic_cast< EnrFrontIntersection * >( eiCandidate->giveEnrichmentFrontStart() );
            if ( efStart != NULL ) {
                const TipInfo &tipInfo = efStart->giveTipInfo();

                if ( ei->tipIsTouchingEI(tipInfo) ) {
                    //printf("Crack %d is touched by a tip on crack %d.\n", iEnrItemIndex, candidateIndex);
                    oTouchingEnrItemIndices.push_back(candidateIndex);
                }
            }

            // Check end tip
            EnrFrontIntersection *efEnd = dynamic_cast< EnrFrontIntersection * >( eiCandidate->giveEnrichmentFrontEnd() );
            if ( efEnd != NULL ) {
                const TipInfo &tipInfo = efEnd->giveTipInfo();

                if ( ei->tipIsTouchingEI(tipInfo) ) {
                    //printf("Crack %d is touched by a tip on crack %d.\n", iEnrItemIndex, candidateIndex);
                    oTouchingEnrItemIndices.push_back(candidateIndex);
                }
            }
        }
    }
}

void XfemStructuralElementInterface :: giveSubtriangulationCompositeExportData(std :: vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep)
{
    const int numCells = mSubTri.size();
    vtkPieces [ 0 ].setNumberOfCells(numCells);

    int numTotalNodes = numCells * 3;
    vtkPieces [ 0 ].setNumberOfNodes(numTotalNodes);

    // Node coordinates
    std :: vector< FloatArray >nodeCoords;
    int nodesPassed = 1;
    for ( auto &tri: mSubTri ) {
        for ( int i = 1; i <= 3; i++ ) {
            FloatArray x = tri.giveVertex(i);
            nodeCoords.push_back(x);
            vtkPieces [ 0 ].setNodeCoords(nodesPassed, x);
            nodesPassed++;
        }
    }

    // Connectivity, offset and cell type
    nodesPassed = 1;
    int offset = 3;
    for ( size_t i = 1; i <= mSubTri.size(); i++ ) {
        IntArray nodes = {
            nodesPassed, nodesPassed + 1, nodesPassed + 2
        };
        nodesPassed += 3;
        vtkPieces [ 0 ].setConnectivity(i, nodes);

        vtkPieces [ 0 ].setOffset(i, offset);
        offset += 3;

        vtkPieces [ 0 ].setCellType(i, 5); // Linear triangle
    }



    // Export nodal variables from primary fields
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);


        nodesPassed = 1;

        for ( auto &tri: mSubTri ) {
            FloatArray triCenter = tri.giveVertex(1);
            triCenter.add( tri.giveVertex(2) );
            triCenter.add( tri.giveVertex(3) );
            triCenter.times(1.0 / 3.0);

            double meanEdgeLength = 0.0;
            meanEdgeLength += ( 1.0 / 3.0 ) * tri.giveVertex(1).distance( tri.giveVertex(2) );
            meanEdgeLength += ( 1.0 / 3.0 ) * tri.giveVertex(2).distance( tri.giveVertex(3) );
            meanEdgeLength += ( 1.0 / 3.0 ) * tri.giveVertex(3).distance( tri.giveVertex(1) );

            const double relPertLength = XfemTolerances :: giveRelLengthTolTight();

            for ( int i = 1; i <= 3; i++ ) {
                if ( type == DisplacementVector ) { // compute displacement
                    FloatArray u = {
                        0.0, 0.0, 0.0
                    };


                    // Fetch global coordinates (in undeformed configuration)
                    const FloatArray &x = tri.giveVertex(i);

                    FloatArray locCoordNode;
                    element->computeLocalCoordinates(locCoordNode, x);


                    FloatArray pertVec;
                    FloatArray locCoord;
                    double pertLength = relPertLength;

                    int numTries = 1;
                    for ( int j = 0; j < numTries; j++ ) {
                        // Perturb towards triangle center, to ensure that
                        // we end up on the correct side of the crack
                        pertVec.beDifferenceOf(triCenter, x);
                        pertVec.times(pertLength);
                        FloatArray xPert = x;
                        xPert.add(pertVec);


                        // Compute local coordinates
                        element->computeLocalCoordinates(locCoord, xPert);
                    }

                    // Compute displacement in point
                    FloatMatrix NMatrix;
                    FloatArray solVec;

                    // Use only regular basis functions for edge nodes to get linear interpolation
                    // (to make the visualization look nice)
                    FloatArray N;
                    FEInterpolation *interp = element->giveInterpolation();
                    interp->evalN( N, locCoordNode, FEIElementGeometryWrapper(element) );
                    const int nDofMan = element->giveNumberOfDofManagers();

                    XfemManager *xMan = element->giveDomain()->giveXfemManager();
                    int numEI =  xMan->giveNumberOfEnrichmentItems();

                    bool joinNodes = false;

                    for ( int eiIndex = 1; eiIndex <= numEI; eiIndex++ ) {
                        EnrichmentItem *ei = xMan->giveEnrichmentItem(eiIndex);

                        double levelSetTang = 0.0, levelSetNormal = 0.0, levelSetInNode = 0.0;

                        bool evaluationSucceeded = true;
                        for ( int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++ ) {
                            DofManager *dMan = element->giveDofManager(elNodeInd);
                            const FloatArray &nodeCoord = * ( dMan->giveCoordinates() );

                            if ( !ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber(), nodeCoord) ) {
                                evaluationSucceeded = false;
                            }
                            levelSetTang += N.at(elNodeInd) * levelSetInNode;

                            if ( !ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber(), nodeCoord) ) {
                                evaluationSucceeded = false;
                            }
                            levelSetNormal += N.at(elNodeInd) * levelSetInNode;
                        }

                        double tangSignDist = levelSetTang, arcPos = 0.0;

                        GeometryBasedEI *geoEI = dynamic_cast< GeometryBasedEI * >( ei );
                        if ( geoEI != NULL ) {
                            // TODO: Consider removing this special treatment. /ES
                            geoEI->giveGeometry()->computeTangentialSignDist(tangSignDist, x, arcPos);
                        }


//                        if ( ( tangSignDist > ( 1.0e-3 ) * meanEdgeLength && fabs(levelSetNormal) < ( 1.0e-2 ) * meanEdgeLength ) && evaluationSucceeded ) {
//                            joinNodes = false;
//                        }

                        if ( ( tangSignDist < ( 1.0e-3 ) * meanEdgeLength || fabs(levelSetNormal) > ( 1.0e-2 ) * meanEdgeLength ) || !evaluationSucceeded ) {
                            joinNodes = false;
                        }
                    }

                    if ( joinNodes ) {
                        // if point on edge
                        XfemElementInterface_createEnrNmatrixAt(NMatrix, locCoord, * element, true);
                        element->computeVectorOf(VM_Total, tStep, solVec);
                    } else   {
                        XfemElementInterface_createEnrNmatrixAt(NMatrix, locCoord, * element, false);
                        element->computeVectorOf(VM_Total, tStep, solVec);
                    }


                    FloatArray uTemp;
                    uTemp.beProductOf(NMatrix, solVec);

                    if ( uTemp.giveSize() == 3 ) {
                        u = uTemp;
                    } else   {
                        u = {
                            uTemp [ 0 ], uTemp [ 1 ], 0.0
                        };
                    }


                    FloatArray valuearray = u;
                    vtkPieces [ 0 ].setPrimaryVarInNode(fieldNum, nodesPassed, valuearray);
                } else {
                    // TODO: Implement
                    printf("fieldNum: %d\n", fieldNum);
                }

                nodesPassed++;
            }
        }
    }


    // Export nodal variables from internal fields
    vtkPieces [ 0 ].setNumberOfInternalVarsToExport(0, numTotalNodes);


    vtkPieces [ 0 ].setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);

        for ( size_t triInd = 1; triInd <= mSubTri.size(); triInd++ ) {

            FloatArray average;
            IntegrationRule *iRule = element->giveIntegrationRule(0);
            computeIPAverageInTriangle(average, iRule, element, type, tStep, mSubTri[triInd-1]);

            if(average.giveSize() == 0) {
                VTKXMLExportModule :: computeIPAverage(average, iRule, element, type, tStep);
            }


            FloatArray averageVoigt;

            if( average.giveSize() == 6 ) {

                averageVoigt.resize(9);

                averageVoigt.at(1) = average.at(1);
                averageVoigt.at(5) = average.at(2);
                averageVoigt.at(9) = average.at(3);
                averageVoigt.at(6) = averageVoigt.at(8) = average.at(4);
                averageVoigt.at(3) = averageVoigt.at(7) = average.at(5);
                averageVoigt.at(2) = averageVoigt.at(4) = average.at(6);
            }
            else {
                if(average.giveSize() == 1) {
                    averageVoigt.resize(1);
                    averageVoigt.at(1) = average.at(1);
                }
            }

            vtkPieces [ 0 ].setCellVar(i, triInd, averageVoigt);
        }
    }



    // Export of XFEM related quantities
    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = element->giveDomain()->giveXfemManager();

        int nEnrIt = xMan->giveNumberOfEnrichmentItems();
        vtkPieces [ 0 ].setNumberOfInternalXFEMVarsToExport(xMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);

        const int nDofMan = element->giveNumberOfDofManagers();


        for ( int field = 1; field <= xMan->vtkExportFields.giveSize(); field++ ) {
            XFEMStateType xfemstype = ( XFEMStateType ) xMan->vtkExportFields [ field - 1 ];

            for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
                EnrichmentItem *ei = xMan->giveEnrichmentItem(enrItIndex);
                for ( int nodeInd = 1; nodeInd <= numTotalNodes; nodeInd++ ) {
                    const FloatArray &x = nodeCoords [ nodeInd - 1 ];
                    FloatArray locCoord;
                    element->computeLocalCoordinates(locCoord, x);

                    FloatArray N;
                    FEInterpolation *interp = element->giveInterpolation();
                    interp->evalN( N, locCoord, FEIElementGeometryWrapper(element) );


                    if ( xfemstype == XFEMST_LevelSetPhi ) {
                        double levelSet = 0.0, levelSetInNode = 0.0;

                        for ( int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++ ) {
                            DofManager *dMan = element->giveDofManager(elNodeInd);
                            const FloatArray &nodeCoord = * ( dMan->giveCoordinates() );
                            ei->evalLevelSetNormalInNode(levelSetInNode, dMan->giveGlobalNumber(), nodeCoord);

                            levelSet += N.at(elNodeInd) * levelSetInNode;
                        }


                        FloatArray valueArray = {
                            levelSet
                        };
                        vtkPieces [ 0 ].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);
                    } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                        double levelSet = 0.0, levelSetInNode = 0.0;

                        for ( int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++ ) {
                            DofManager *dMan = element->giveDofManager(elNodeInd);
                            const FloatArray &nodeCoord = * ( dMan->giveCoordinates() );
                            ei->evalLevelSetTangInNode(levelSetInNode, dMan->giveGlobalNumber(), nodeCoord);

                            levelSet += N.at(elNodeInd) * levelSetInNode;
                        }


                        FloatArray valueArray = {
                            levelSet
                        };
                        vtkPieces [ 0 ].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);
                    } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
                        double nodeEnrMarker = 0.0, nodeEnrMarkerInNode = 0.0;

                        for ( int elNodeInd = 1; elNodeInd <= nDofMan; elNodeInd++ ) {
                            DofManager *dMan = element->giveDofManager(elNodeInd);
                            ei->evalNodeEnrMarkerInNode( nodeEnrMarkerInNode, dMan->giveGlobalNumber() );

                            nodeEnrMarker += N.at(elNodeInd) * nodeEnrMarkerInNode;
                        }



                        FloatArray valueArray = {
                            nodeEnrMarker
                        };
                        vtkPieces [ 0 ].setInternalXFEMVarInNode(field, enrItIndex, nodeInd, valueArray);
                    }
                }
            }
        }
    }
}

void XfemStructuralElementInterface :: computeIPAverageInTriangle(FloatArray &answer, IntegrationRule *iRule, Element *elem, InternalStateType isType, TimeStep *tStep, const Triangle &iTri)
{
    // Computes the volume average (over an element) for the quantity defined by isType
    double gptot = 0.0;
    answer.clear();
    FloatArray temp;
    if ( iRule ) {
        for ( IntegrationPoint *ip: *iRule ) {

            FloatArray globCoord = ip->giveGlobalCoordinates();
            globCoord.resizeWithValues(2);

            if( iTri.pointIsInTriangle(globCoord) ) {
                elem->giveIPValue(temp, ip, isType, tStep);
                gptot += ip->giveWeight();
                answer.add(ip->giveWeight(), temp);
            }
        }

        answer.times(1. / gptot);
    }

}

} /* namespace oofem */
