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

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"
#include "structuralelement.h"
#include "structuralcrosssection.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "dynamicinputrecord.h"

#include "xfem/patchintegrationrule.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/XFEMDebugTools.h"

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
        for ( size_t i = 0; i < mpCZIntegrationRules.size(); i++ ) {
            if ( mpCZIntegrationRules [ i ] != NULL ) {
                delete mpCZIntegrationRules [ i ];
            }
        }

        mpCZIntegrationRules.clear();
        mCZEnrItemIndices.clear();
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        if ( mpCZMat == NULL && mCZMaterialNum > 0 ) {
            initializeCZMaterial();
        }


        MaterialMode matMode = element->giveMaterialMode();

        bool firstIntersection = true;

        std :: vector< std :: vector< FloatArray > >pointPartitions;
        std :: vector< Triangle >allTri;

        std :: vector< int >enrichingEIs;
        int elPlaceInArray = xMan->giveDomain()->giveElementPlaceInArray( element->giveGlobalNumber() );
        xMan->giveElementEnrichmentItemIndices(enrichingEIs, elPlaceInArray);


        for ( size_t p = 0; p < enrichingEIs.size(); p++ ) {
            int eiIndex = enrichingEIs [ p ];

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
                        this->XfemElementInterface_partitionElement(allTri, pointPartitions [ i ]);
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
                            mpCZIntegrationRules.push_back( new GaussIntegrationRule(czRuleNum, element) );
                            mCZEnrItemIndices.push_back(eiIndex);

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
                                                                        crackPolygon[segIndex],crackPolygon[segIndex + 1]);

                            for ( GaussPoint *gp: *mpCZIntegrationRules [ segIndex ] ) {
                                double gw = gp->giveWeight();
                                double segLength = crackPolygon [ segIndex ].distance(crackPolygon [ segIndex + 1 ]);
                                gw *= 0.5 * segLength;
                                gp->setWeight(gw);

                                FloatArray locCoord;
                                element->computeLocalCoordinates(locCoord, *(gp->giveCoordinates()) );
                                gp->setLocalCoordinates(locCoord);

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
                for ( size_t triIndex = 0; triIndex < allTri.size(); triIndex++ ) {
                    // Call alternative version of XfemElementInterface_prepareNodesForDelaunay
                    std :: vector< std :: vector< FloatArray > >pointPartitionsTri;
                    double startXi, endXi;
                    bool intersection = false;
                    XfemElementInterface_prepareNodesForDelaunay(pointPartitionsTri, startXi, endXi, allTri [ triIndex ], eiIndex, intersection);

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
                                mpCZIntegrationRules.push_back( new GaussIntegrationRule(czRuleNum, element) );
                                size_t newRuleInd = mpCZIntegrationRules.size() - 1;
                                mCZEnrItemIndices.push_back(eiIndex);

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
                                                                            crackPolygon[segIndex],crackPolygon[segIndex + 1]);

                                for ( GaussPoint *gp: * mpCZIntegrationRules [ newRuleInd ] ) {
                                    double gw = gp->giveWeight();
                                    double segLength = crackPolygon [ segIndex ].distance(crackPolygon [ segIndex + 1 ]);
                                    gw *= 0.5 * segLength;
                                    gp->setWeight(gw);

                                    FloatArray locCoord;
                                    element->computeLocalCoordinates(locCoord, *(gp->giveCoordinates()) );
                                    gp->setLocalCoordinates(locCoord);

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
                        allTriCopy.push_back(allTri [ triIndex ]);
                    }
                }

                allTri = allTriCopy;
            }
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

            XFEMDebugTools :: WriteTrianglesToVTK(name3, allTri);
        }


        int ruleNum = 1;

        if ( partitionSucceeded ) {
            IntegrationRule *intRule = new PatchIntegrationRule(ruleNum, element, allTri);
            intRule->SetUpPointsOnTriangle(xMan->giveNumGpPerTri(), matMode);
            element->setIntegrationRules({intRule});
        }


        if ( xMan->giveVtkDebug() ) {
            ////////////////////////////////////////////////////////////////////////
            // Write CZ GP to VTK

            std :: vector< FloatArray >czGPCoord;

            for ( size_t czRulInd = 0; czRulInd < mpCZIntegrationRules.size(); czRulInd++ ) {
                for ( GaussPoint *gp: * mpCZIntegrationRules [ czRulInd ] ) {
                    czGPCoord.push_back( * ( gp->giveCoordinates() ) );
                }
            }

            double time = 0.0;

            Element *el = element;

            Domain *dom = el->giveDomain();
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
                    structCS->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
                    return;
                } else {
                    OOFEM_ERROR("failed to fetch StructuralMaterial");
                }
            }
        }
    }

    // If no enrichment modifies the material,
    // compute stiffness based on the bulk material.
    StructuralElement &structEl = dynamic_cast< StructuralElement & >( ( * element ) );
    structEl.StructuralElement :: computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
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
                    } else   {
                        structCSInclusion->giveRealStress_PlaneStress(answer, gp, strain, tStep);
                    }

                    return;
                } else {
                    OOFEM_ERROR("failed to fetch StructuralCrossSection");
                }
            }
        }
    }


    if ( mUsePlaneStrain ) {
        cs->giveRealStress_PlaneStrain(answer, gp, strain, tStep);
    } else   {
        cs->giveRealStress_PlaneStress(answer, gp, strain, tStep);
    }
}

void XfemStructuralElementInterface :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep)
{
    if ( hasCohesiveZone() ) {
        FloatArray solVec;
        element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();
        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {

            for ( GaussPoint *gp: * mpCZIntegrationRules [ segIndex ] ) {
                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ]);
                ////////////////////////////////////////////////////////


                // Traction
                FloatArray T2D;



                // Fetch material status and get normal
                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                if ( ms == NULL ) {
                    OOFEM_ERROR("Failed to fetch material status.");
                }

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

    FloatArray TLoc(3), jump3DLoc, TLocRenumbered(3);
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
        element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();

        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            for ( GaussPoint *gp: * mpCZIntegrationRules [ segIndex ] ) {

                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ]);

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
    IntegrationRule *iRule = element->giveIntegrationRule(0);
    IntArray mask;

    answer.resize(ndofs, ndofs);
    answer.zero();
    if ( !structEl->isActivated(tStep) ) {
        return;
    }

    structEl->giveMassMtrxIntegrationgMask(mask);

    mass = 0.;

    for ( GaussPoint *gp: *iRule ) {
        structEl->computeNmatrixAt(* ( gp->giveLocalCoordinates() ), n);
        density = structEl->giveMaterial()->give('d', gp);

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
} /* namespace oofem */
