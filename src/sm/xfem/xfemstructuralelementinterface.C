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
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerial.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "engngm.h"
#include "sm/Elements/nlstructuralelement.h"
#include "mathfem.h"

#include "sm/Materials/structuralfe2material.h"
#include "prescribedgradienthomogenization.h"

#include "xfem/patchintegrationrule.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/XFEMDebugTools.h"
#include "xfem/xfemtolerances.h"

#include "xfem/enrichmentfronts/enrichmentfrontintersection.h"
#include "xfem/xfemstructuremanager.h"

#include "vtkxmlexportmodule.h"

#include "prescribedgradientbcweak.h"
#include "mathfem.h"
#include <cmath>

#include <string>
#include <sstream>


////////////////////////////////////
//
// Alt. I  : include_bulk_jump = 1 and include_bulk_corr = 1
//
// Alt. II : include_bulk_jump = 1 and include_bulk_corr = 0
//
// Alt. III: include_bulk_jump = 0 and include_bulk_corr = 0


namespace oofem {
XfemStructuralElementInterface :: XfemStructuralElementInterface(Element *e) :
    XfemElementInterface(e)
{}

bool XfemStructuralElementInterface :: XfemElementInterface_updateIntegrationRule()
{
    const double tol2 = 1.0e-18;

    bool partitionSucceeded = false;


    if ( mpCZMat != nullptr ) {
        mpCZIntegrationRules_tmp.clear();
        mpCZExtraIntegrationRules_tmp.clear();
        mIntRule_tmp.clear();
        mCZEnrItemIndices.clear();
        mCZTouchingEnrItemIndices.clear();
    }

    XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
    if ( xMan->isElementEnriched(element) ) {
        if ( mpCZMat == nullptr && mCZMaterialNum > 0 ) {
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
                    for ( auto &partition : pointPartitions ) {
                        // Triangulate the subdivisions
                        this->XfemElementInterface_partitionElement(mSubTri, partition);
                    }


                    if ( mpCZMat != nullptr ) {
                        Crack *crack = dynamic_cast< Crack * >( xMan->giveEnrichmentItem(eiIndex) );
                        if ( crack == nullptr ) {
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
                            mpCZIntegrationRules_tmp.push_back( std::make_unique<GaussIntegrationRule>(czRuleNum, element) );

                            if ( mIncludeBulkCorr ) {
                                mpCZExtraIntegrationRules_tmp.push_back( std::make_unique<GaussIntegrationRule>(czRuleNum, element) );
                            }

                            size_t cz_rule_ind = mpCZIntegrationRules_tmp.size() - 1;

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
                            } else {
                                // Oops, we got a segment of length zero.
                                // These Gauss weights will be zero, so we can
                                // set the tangent to anything reasonable
                                crackTang = {0.0, 1.0};
                            }

                            FloatArray crackNormal = {
                                -crackTang.at(2), crackTang.at(1), 0.
                            };

                            mpCZIntegrationRules_tmp [ cz_rule_ind ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                           crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);

                            if ( mIncludeBulkCorr ) {
                                mpCZExtraIntegrationRules_tmp [ cz_rule_ind ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                           crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                            }

                            for ( auto &gp: *mpCZIntegrationRules_tmp [ cz_rule_ind ] ) {
                                double gw = gp->giveWeight();
                                double segLength = distance(crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                                gw *= 0.5 * segLength;
                                gp->setWeight(gw);

                                // Fetch material status and set normal
                                StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                                if ( ms ) {
                                    ms->letNormalBe(crackNormal);
                                } else {
                                    StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*>( mpCZMat->giveStatus(gp) );

                                    if ( fe2ms ) {
                                        fe2ms->letNormalBe(crackNormal);

                                        PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>( fe2ms->giveBC() );

                                        if ( bc ) {
                                            FloatArray periodicityNormal = crackNormal;

                                            periodicityNormal.normalize();

                                            if ( periodicityNormal(0) < 0.0 && periodicityNormal(1) < 0.0 ) {
                                                // Rotate 90 degrees (works equally well for periodicity)
                                                periodicityNormal = {periodicityNormal(1), -periodicityNormal(0)};
                                            }

                                            bc->setPeriodicityNormal(periodicityNormal);
                                            bc->recomputeTractionMesh();
                                        }
                                    } else {
                                        OOFEM_ERROR("Failed to fetch material status.");
                                    }
                                }

                                // Give Gauss point reference to the enrichment item
                                // to simplify post processing.
                                crack->AppendCohesiveZoneGaussPoint(gp);
                            }


                            if ( mIncludeBulkCorr ) {
                                for ( auto &gp: *mpCZExtraIntegrationRules_tmp [ cz_rule_ind ] ) {
                                    double gw = gp->giveWeight();
                                    double segLength = distance(crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                                    gw *= 0.5 * segLength;
                                    gp->setWeight(gw);

                                    // Fetch material status and set normal
                                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( gp->giveMaterialStatus() );
                                    if ( ms ) {
                                        ms->letNormalBe(crackNormal);
                                    } else {
                                        StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*>( gp->giveMaterialStatus() );

                                        if ( fe2ms ) {
                                            fe2ms->letNormalBe(crackNormal);

                                            PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>( fe2ms->giveBC() );

                                            if ( bc) {
                                                FloatArray periodicityNormal = crackNormal;
                                                periodicityNormal.normalize();

                                                if ( periodicityNormal(0) < 0.0 && periodicityNormal(1) < 0.0 ) {
                                                    // Rotate 90 degrees (works equally well for periodicity)
                                                    periodicityNormal = {periodicityNormal(1), -periodicityNormal(0)};
                                                }

                                                bc->setPeriodicityNormal(periodicityNormal);
                                                bc->recomputeTractionMesh();

                                            }
                                        } else {
                                            // Macroscale material model: nothing needs to be done.
                                        }
                                    }
                                }
                            }
                        }
                    }

                    partitionSucceeded = true;
                }
            } else { // if ( firstIntersection )
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

                        if ( mpCZMat != nullptr ) {
                            Crack *crack = dynamic_cast< Crack * >( xMan->giveEnrichmentItem(eiIndex) );
                            if ( crack == nullptr ) {
                                OOFEM_ERROR("Cohesive zones are only available for cracks.")
                            }

                            // We have xi_s and xi_e. Fetch sub polygon.
                            std :: vector< FloatArray >crackPolygon;
                            crack->giveSubPolygon(crackPolygon, startXi, endXi);

                            int numSeg = crackPolygon.size() - 1;

                            for ( int segIndex = 0; segIndex < numSeg; segIndex++ ) {
                                int czRuleNum = 1;
                                mpCZIntegrationRules_tmp.emplace_back( std::make_unique<GaussIntegrationRule>(czRuleNum, element) );

                                if ( mIncludeBulkCorr ) {
                                    mpCZExtraIntegrationRules_tmp.emplace_back( std::make_unique<GaussIntegrationRule>(czRuleNum, element) );
                                }

                                size_t newRuleInd = mpCZIntegrationRules_tmp.size() - 1;
                                mCZEnrItemIndices.push_back(eiIndex);

                                mCZTouchingEnrItemIndices.push_back(touchingEiIndices);

                                // Compute crack normal
                                FloatArray crackTang;
                                crackTang.beDifferenceOf(crackPolygon [ segIndex + 1 ], crackPolygon [ segIndex ]);

                                if ( crackTang.computeSquaredNorm() > tol2 ) {
                                    crackTang.normalize();
                                } else {
                                    // Oops, we got a segment of length zero.
                                    // These Gauss weights will be zero, so we can
                                    // set the tangent to anything reasonable
                                    crackTang = {0.0, 1.0};
                                }

                                FloatArrayF<3> crackNormal = { -crackTang.at(2), crackTang.at(1), 0.};

                                mpCZIntegrationRules_tmp [ newRuleInd ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                                 crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);

                                if ( mIncludeBulkCorr ) {
                                    mpCZExtraIntegrationRules_tmp [ newRuleInd ]->SetUpPointsOn2DEmbeddedLine(mCSNumGaussPoints, matMode,
                                                                                                 crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                                }

                                for ( auto &gp: *mpCZIntegrationRules_tmp [ newRuleInd ] ) {
                                    double gw = gp->giveWeight();
                                    double segLength = distance(crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                                    gw *= 0.5 * segLength;
                                    gp->setWeight(gw);

                                    // Fetch material status and set normal
                                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                                    if ( ms ) {
                                        ms->letNormalBe(crackNormal);
                                    } else {
                                        StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*>( mpCZMat->giveStatus(gp) );

                                        if ( fe2ms ) {
                                            fe2ms->letNormalBe(crackNormal);

                                            PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>( fe2ms->giveBC() );

                                            if ( bc ) {
                                                FloatArray periodicityNormal = crackNormal;

                                                periodicityNormal.normalize();

                                                if ( periodicityNormal(0) < 0.0 && periodicityNormal(1) < 0.0 ) {
                                                    // Rotate 90 degrees (works equally well for periodicity)
                                                    periodicityNormal = {periodicityNormal(1), -periodicityNormal(0)};
                                                }

                                                bc->setPeriodicityNormal(periodicityNormal);
                                                bc->recomputeTractionMesh();
                                            }
                                        } else {
                                            OOFEM_ERROR("Failed to fetch material status.");
                                        }
                                    }

                                    // Give Gauss point reference to the enrichment item
                                    // to simplify post processing.
                                    crack->AppendCohesiveZoneGaussPoint(gp);
                                }

                                if ( mIncludeBulkCorr ) {

                                    for ( auto &gp: *mpCZExtraIntegrationRules_tmp [ newRuleInd ] ) {
                                        double gw = gp->giveWeight();
                                        double segLength = distance(crackPolygon [ segIndex ], crackPolygon [ segIndex + 1 ]);
                                        gw *= 0.5 * segLength;
                                        gp->setWeight(gw);


                                        // Fetch material status and set normal
                                        StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                                        if ( ms ) {
                                            ms->letNormalBe(crackNormal);
                                        } else {
                                            StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*>( mpCZMat->giveStatus(gp) );

                                            if ( fe2ms ) {
                                                fe2ms->letNormalBe(crackNormal);

                                                PrescribedGradientBCWeak *bc = dynamic_cast<PrescribedGradientBCWeak*>( fe2ms->giveBC() );

                                                if ( bc ) {
                                                    FloatArray periodicityNormal = crackNormal;

                                                    periodicityNormal.normalize();

                                                    if ( periodicityNormal(0) < 0.0 && periodicityNormal(1) < 0.0 ) {
                                                        // Rotate 90 degrees (works equally well for periodicity)
                                                        periodicityNormal = {periodicityNormal(1), -periodicityNormal(0)};
                                                    }

                                                    bc->setPeriodicityNormal(periodicityNormal);
                                                    bc->recomputeTractionMesh();
                                                }
                                            } else {
                                                OOFEM_ERROR("Failed to fetch material status.");
                                            }
                                        }

                                        // Give Gauss point reference to the enrichment item
                                        // to simplify post processing.
                                        crack->AppendCohesiveZoneGaussPoint(gp);
                                    }
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

        for ( int i = 0; i < numRefs; i++ ) {
            std :: vector< Triangle > triRef;
            for ( const auto &tri : mSubTri ) {
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

//        printf("totalCrackLengthInEl: %e\n", totalCrackLengthInEl);

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
            mIntRule_tmp.push_back( std::make_unique<PatchIntegrationRule>(ruleNum, element, mSubTri) );
            mIntRule_tmp [ 0 ]->SetUpPointsOnTriangle(xMan->giveNumGpPerTri(), matMode);
        }


        if ( xMan->giveVtkDebug() ) {
            ////////////////////////////////////////////////////////////////////////
            // Write CZ GP to VTK

            std :: vector< FloatArray >czGPCoord;

            for ( size_t czRulInd = 0; czRulInd < mpCZIntegrationRules_tmp.size(); czRulInd++ ) {
                for ( auto &gp: *mpCZIntegrationRules_tmp [ czRulInd ] ) {
                    czGPCoord.push_back( gp->giveGlobalCoordinates() );
                }
            }

            double time = 0.0;

            Domain *dom = element->giveDomain();
            if ( dom != nullptr ) {
                EngngModel *em = dom->giveEngngModel();
                if ( em != nullptr ) {
                    TimeStep *ts = em->giveCurrentStep();
                    if ( ts != nullptr ) {
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


        bool map_state_variables = true;
        if ( map_state_variables ) {

            ////////////////////////////////////////////////////////////////////////
            // Copy bulk GPs
            int num_ir_new = mIntRule_tmp.size();
            for ( int i = 0; i < num_ir_new; i++ ) {

                for ( auto &gp_new : *mIntRule_tmp[i] ) {


                    // Fetch new material status. Create it if it does not exist.
                    MaterialStatus *ms_new = dynamic_cast<MaterialStatus*>( gp_new->giveMaterialStatus() );
                    if ( !ms_new ) {
                        StructuralElement *s_el = dynamic_cast<StructuralElement*>(element);
                        StructuralCrossSection *cs = s_el->giveStructuralCrossSection();
                        cs->createMaterialStatus(*gp_new);
                        ms_new = dynamic_cast<MaterialStatus*>( gp_new->giveMaterialStatus() );
                    }

                    // Find closest old GP.
                    double closest_dist = 0.0;
                    MaterialStatus *ms_old = giveClosestGP_MatStat( closest_dist, element->giveIntegrationRulesArray() , gp_new->giveGlobalCoordinates());

                    if ( ms_old ) {
                        // Copy state variables.
                        MaterialStatusMapperInterface *mapper_interface = dynamic_cast<MaterialStatusMapperInterface*>( ms_new );

                        if ( mapper_interface ) {
                            //printf("Successfully casted to MaterialStatusMapperInterface.\n");
                            mapper_interface->copyStateVariables(*ms_old);
                        } else {
                            OOFEM_ERROR("Failed casting to MaterialStatusMapperInterface.")
                        }
                    }

                }

            }

            ////////////////////////////////////////////////////////////////////////
            // Copy CZ GPs
            if ( mCZMaterialNum > 0 ) {

                num_ir_new = mpCZIntegrationRules_tmp.size();
                for ( int i = 0; i < num_ir_new; i++ ) {

                    for ( auto &gp_new : *(mpCZIntegrationRules_tmp[i]) ) {

                        // Fetch new material status. Create it if it does not exist.
                        MaterialStatus *ms_new = mpCZMat->giveStatus(gp_new);

                        // Find closest old GP.
                        double closest_dist_cz = 0.0;
                        MaterialStatus *ms_old = giveClosestGP_MatStat(closest_dist_cz, mpCZIntegrationRules , gp_new->giveGlobalCoordinates());

                        // If we are using the nonstandard cz formulation, we
                        // also consider mapping from bulk GPs.
                        XfemStructureManager *xsMan = dynamic_cast<XfemStructureManager*>( element->giveDomain()->giveXfemManager() );
                        bool non_std_cz = false;
                        if ( xsMan ) {
                            non_std_cz = xsMan->giveUseNonStdCz();
                        }

                        if ( non_std_cz ) {
                            double closest_dist_bulk = 0.0;
                            MaterialStatus *ms_old_bulk = giveClosestGP_MatStat(closest_dist_bulk, element->giveIntegrationRulesArray() , gp_new->giveGlobalCoordinates());
                            if ( closest_dist_bulk < closest_dist_cz ) {
                                printf("Bulk is closest. Dist: %e\n", closest_dist_bulk);
                                ms_old = ms_old_bulk;
                            }
                        }

                        if ( ms_old ) {

                            // Copy state variables.
                            MaterialStatusMapperInterface *mapper_interface = dynamic_cast<MaterialStatusMapperInterface*>( ms_new );

                            if ( mapper_interface ) {
                                mapper_interface->copyStateVariables(*ms_old);
                            } else {
                                OOFEM_ERROR("Failed casting to MaterialStatusMapperInterface.")
                            }
                        }

                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////
            // Copy "extra" CZ GPs (Used for non-standard CZ model)
            if ( mCZMaterialNum > 0 && mIncludeBulkCorr ) {

                num_ir_new = mpCZExtraIntegrationRules_tmp.size();
                for ( int i = 0; i < num_ir_new; i++ ) {

                    for ( auto &gp_new : *(mpCZExtraIntegrationRules_tmp[i]) ) {

                        // Fetch new material status. Create it if it does not exist.
                        MaterialStatus *ms_new = mpCZMat->giveStatus(gp_new);

                        // Find closest old GP.
                        double closest_dist_cz = 0.0;
                        MaterialStatus *ms_old = giveClosestGP_MatStat(closest_dist_cz, mpCZExtraIntegrationRules , gp_new->giveGlobalCoordinates());

                        // If we are using the nonstandard cz formulation, we
                        // also consider mapping from bulk GPs.
                        XfemStructureManager *xsMan = dynamic_cast<XfemStructureManager*>( element->giveDomain()->giveXfemManager() );
                        bool non_std_cz = false;
                        if ( xsMan ) {
                            non_std_cz = xsMan->giveUseNonStdCz();
                        }

                        if ( non_std_cz ) {
                            double closest_dist_bulk = 0.0;
                            MaterialStatus *ms_old_bulk = giveClosestGP_MatStat(closest_dist_bulk, element->giveIntegrationRulesArray() , gp_new->giveGlobalCoordinates());
                            if ( closest_dist_bulk < closest_dist_cz ) {
                                printf("Bulk is closest. Dist: %e\n", closest_dist_bulk);
                                ms_old = ms_old_bulk;
                            }
                        }

                        if ( ms_old ) {

                            // Copy state variables.
                            MaterialStatusMapperInterface *mapper_interface = dynamic_cast<MaterialStatusMapperInterface*>( ms_new );

                            if ( mapper_interface ) {
                                mapper_interface->copyStateVariables(*ms_old);
                            } else {
                                OOFEM_ERROR("Failed casting to MaterialStatusMapperInterface.")
                            }
                        }

                    }
                }
            } else {
                //printf("No extra CZ GPs to map.\n");
            }
        }
    }

    if ( partitionSucceeded ) {
        mpCZIntegrationRules = std::move(mpCZIntegrationRules_tmp);

        if ( mIncludeBulkCorr ) {
            mpCZExtraIntegrationRules = std::move(mpCZExtraIntegrationRules_tmp);
        }

        element->setIntegrationRules( std :: move(mIntRule_tmp) );
    }

    return partitionSucceeded;
}

MaterialStatus* XfemStructuralElementInterface :: giveClosestGP_MatStat(double &oClosestDist, std :: vector< std :: unique_ptr< IntegrationRule > > &iRules, const FloatArray &iCoord)
{
    double min_dist2 = std::numeric_limits<double>::max();
    MaterialStatus *closest_ms = nullptr;

    for ( size_t i = 0; i < iRules.size(); i++ ) {
        for ( auto &gp : *(iRules[i]) ) {

            const FloatArray &x = gp->giveGlobalCoordinates();
            double d2 = distance_square(x, iCoord);
            if ( d2 < min_dist2 ) {
                min_dist2 = d2;
                closest_ms = dynamic_cast<MaterialStatus*>( gp->giveMaterialStatus() );
            }

        }
    }

    oClosestDist = sqrt(min_dist2);
    return closest_ms;
}

double XfemStructuralElementInterface :: computeEffectiveSveSize(StructuralFE2MaterialStatus *iFe2Ms)
{

#if 0
    //return 1.0*sqrt( iFe2Ms->giveBC()->domainSize() );
    //return 2.0*sqrt( iFe2Ms->giveBC()->domainSize() );

#else
    // TODO: Cover also angle < 0 and angle > 90.

    const FloatArray &n = iFe2Ms->giveNormal();
    double l_box = sqrt( iFe2Ms->giveBC()->domainSize() );

    const FloatArray t = {n(1), -n(0)};
    double angle = atan2( t(1), t(0) );

    if ( angle < 0.25*M_PI ) {
        // angle < 45 degrees

        double l_s = l_box*(cos(angle));
        return l_s;
    } else {
        // angle >= 45 degrees

        double l_s = l_box*(sin(angle));
        return l_s;
    }
#endif
}

void XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = element->giveDomain()->giveXfemManager();
        CrossSection *cs = nullptr;

        const std :: vector< int > &materialModifyingEnrItemIndices = xMan->giveMaterialModifyingEnrItemIndices();
        for ( size_t i = 0; i < materialModifyingEnrItemIndices.size(); i++ ) {
            EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(materialModifyingEnrItemIndices [ i ]) );

            if ( ei.isMaterialModified(* gp, * element, cs) ) {
                StructuralCrossSection *structCS = dynamic_cast< StructuralCrossSection * >( cs );

                if ( structCS != nullptr ) {
                    if ( mUsePlaneStrain ) {
                        answer = structCS->giveStiffnessMatrix_PlaneStrain(rMode, gp, tStep);
                    } else {
                        answer = structCS->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
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
    auto cs = dynamic_cast< StructuralCrossSection * >( element->giveCrossSection() );
    if ( mUsePlaneStrain ) {
        answer = cs->giveStiffnessMatrix_PlaneStrain(rMode, gp, tStep);
    } else {
        answer = cs->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
    }
}

void XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    auto cs = dynamic_cast< StructuralCrossSection * >( element->giveCrossSection() );
    if ( cs == nullptr ) {
        OOFEM_ERROR("cs == nullptr.");
    }

    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = element->giveDomain()->giveXfemManager();


        CrossSection *csInclusion = nullptr;
        const std :: vector< int > &materialModifyingEnrItemIndices = xMan->giveMaterialModifyingEnrItemIndices();
        for ( size_t i = 0; i < materialModifyingEnrItemIndices.size(); i++ ) {
            EnrichmentItem &ei = * ( xMan->giveEnrichmentItem(materialModifyingEnrItemIndices [ i ]) );

            if ( ei.isMaterialModified(* gp, * element, csInclusion) ) {
                auto structCSInclusion = dynamic_cast< StructuralCrossSection * >( csInclusion );

                if ( structCSInclusion != nullptr ) {
                    if ( mUsePlaneStrain ) {
                        answer = structCSInclusion->giveRealStress_PlaneStrain(strain, gp, tStep);
                    } else {
                        answer = structCSInclusion->giveRealStress_PlaneStress(strain, gp, tStep);
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
        answer = cs->giveRealStress_PlaneStrain(strain, gp, tStep);
    } else {
        answer = cs->giveRealStress_PlaneStress(strain, gp, tStep);
    }
}

void XfemStructuralElementInterface :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep)
{
    if ( !useNonStdCz() ) {

        if ( hasCohesiveZone() ) {
            FloatArray solVec;
            element->computeVectorOf(VM_Total, tStep, solVec);

            size_t numSeg = mpCZIntegrationRules.size();
            for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                for ( auto &gp: *mpCZIntegrationRules [ segIndex ] ) {
                    ////////////////////////////////////////////////////////
                    // Compute a (slightly modified) N-matrix

                    FloatMatrix NMatrix;
                    computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ], mCZTouchingEnrItemIndices [ segIndex ]);
                    ////////////////////////////////////////////////////////


                    // Traction
                    FloatArray T2D;


                    // Fetch material status and get normal
                    StructuralInterfaceMaterialStatus *ms = dynamic_cast< StructuralInterfaceMaterialStatus * >( mpCZMat->giveStatus(gp) );
                    if ( ms == nullptr ) {
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
    } else {
        // Non-standard cz formulation.
        //printf("Using non-standard cz formulation.\n");

        if ( hasCohesiveZone() ) {
            FloatArray solVec;
            element->computeVectorOf(VM_Total, tStep, solVec);

            StructuralFE2Material *fe2Mat = dynamic_cast<StructuralFE2Material*>(mpCZMat);
            if ( !fe2Mat ) {
                OOFEM_ERROR("Failed to cast StructuralFE2Material*.")
            }

            size_t numSeg = mpCZIntegrationRules.size();
            for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                for ( int gpInd = 0; gpInd < mpCZIntegrationRules[ segIndex ]->giveNumberOfIntegrationPoints(); gpInd++ ) {

                    GaussPoint *gp = mpCZIntegrationRules[ segIndex ]->getIntegrationPoint(gpInd);

                    GaussPoint *bulk_gp = nullptr;
                    if ( mIncludeBulkCorr ) {
                        bulk_gp = mpCZExtraIntegrationRules[ segIndex ]->getIntegrationPoint(gpInd);
                    }

                    StructuralMaterial *bulkMat = dynamic_cast<StructuralMaterial*>( element->giveCrossSection()->giveMaterial(bulk_gp) );
                    if ( !bulkMat ) {
                        OOFEM_ERROR("Failed to fetch bulk material.")
                    }

                    StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*> ( gp->giveMaterialStatus() );

                    if ( fe2ms == nullptr ) {
                        OOFEM_ERROR("The material status is not of an allowed type.")
                    }

                    ////////////////////////////////////////////////////////
                    // Compute jump

                    // Compute a (slightly modified) N-matrix
                    FloatMatrix NMatrix;
                    computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ], mCZTouchingEnrItemIndices [ segIndex ]);

                    FloatArray jump2D;
                    computeDisplacementJump(* gp, jump2D, solVec, NMatrix);

                    ////////////////////////////////////////////////////////
                    // Fetch normal
                    FloatArray crackNormal( fe2ms->giveNormal() );


                    ////////////////////////////////////////////////////////
                    // Fetch L_s
                    double l_s = computeEffectiveSveSize( fe2ms );


                    ////////////////////////////////////////////////////////
                    // Construct strain (only consider the smeared jump for now)
                    FloatArray smearedJumpStrain = {jump2D(0)*crackNormal(0)/l_s, jump2D(1)*crackNormal(1)/l_s, 0.0, 0.0, 0.0, (1.0/l_s)*( jump2D(0)*crackNormal(1) + jump2D(1)*crackNormal(0) )};

                    FloatArray smearedBulkStrain(6);
                    smearedBulkStrain.zero();
                    FloatMatrix BAvg;

                    if ( mIncludeBulkJump ) {
                        ////////////////////////////////////////////////////////
                        // Bulk contribution to SVE strain

                        // Crack gp coordinates
                        const FloatArray &xC = gp->giveGlobalCoordinates();

                        // For now, we will just perturb the coordinates of the GP to compute B^- and B^+ numerically.
                        double eps = 1.0e-6;
                        FloatArray xPert = xC;

                        xPert.add(eps, crackNormal);
                        FloatArray locCoordPert;
                        element->computeLocalCoordinates(locCoordPert, xPert);

                        FloatMatrix BPlus;
                        this->ComputeBOrBHMatrix(BPlus, *gp, *element, false, locCoordPert);

                        xPert = xC;
                        xPert.add(-eps, crackNormal);
                        element->computeLocalCoordinates(locCoordPert, xPert);

                        FloatMatrix BMinus;
                        this->ComputeBOrBHMatrix(BMinus, *gp, *element, false, locCoordPert);

                        BAvg = BPlus;
                        BAvg.add(1.0, BMinus);
                        BAvg.times(0.5);

                        smearedBulkStrain.beProductOf(BAvg, solVec);

                        if ( smearedBulkStrain.giveSize() == 4 ) {
                            smearedBulkStrain = {smearedBulkStrain(0), smearedBulkStrain(1), smearedBulkStrain(2), 0.0, 0.0, smearedBulkStrain(3)};
                        }

                        FloatArray smearedJumpStrainTemp = smearedJumpStrain;

                        smearedJumpStrain.add(smearedBulkStrain);
                    }


                    ////////////////////////////////////////////////////////
                    // Compute homogenized stress
                    StructuralElement *se = dynamic_cast<StructuralElement*>(this->element);
                    if ( !se ) {
                        OOFEM_ERROR("Failed to cast StructuralElement.")
                    }

                    auto stressVec = fe2Mat->giveRealStressVector_3d(smearedJumpStrain, gp, tStep);
                    //printf("stressVec: "); stressVec.printYourself();

                    FloatArray trac = {stressVec[0]*crackNormal[0]+stressVec[5]*crackNormal[1], stressVec[5]*crackNormal[0]+stressVec[1]*crackNormal[1]};

                    ////////////////////////////////////////////////////////
                    // Standard part

                    // Add to internal force
                    FloatArray NTimesT;

                    NTimesT.beTProductOf(NMatrix, trac);
                    CrossSection *cs  = element->giveCrossSection();
                    double thickness = cs->give(CS_Thickness, gp);
                    double dA = thickness * gp->giveWeight();
                    answer.add(dA, NTimesT);

                    if ( mIncludeBulkCorr ) {
                        auto stressVecBulk = bulkMat->giveRealStressVector_3d(smearedBulkStrain, bulk_gp, tStep);

                        ////////////////////////////////////////////////////////
                        // Non-standard jump part
                        FloatArray stressV4 = {stressVec[0]-stressVecBulk[0], stressVec[1]-stressVecBulk[1], stressVec[2]-stressVecBulk[2], stressVec[5]-stressVecBulk[5]};
                        FloatArray BTimesT;
                        BTimesT.beTProductOf(BAvg, stressV4);
                        answer.add(1.0*dA*l_s, BTimesT);
                    }
                }
            }
        }
    }
}

void XfemStructuralElementInterface :: computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep)
{
    auto F = eye<3>(); // TODO: Compute properly

    FloatArrayF<3> jump3D = {iJump.at(1), iJump.at(2), 0.0};

    FloatArrayF<3> crackNormal3D = {iCrackNormal.at(1), iCrackNormal.at(2), 0.0};

    FloatArrayF<3> ez = {0.0, 0.0, 1.0};
    auto crackTangent3D = cross(crackNormal3D, ez);

    FloatMatrixF<3,3> locToGlob;
    locToGlob.setColumn(crackTangent3D, 0);
    locToGlob.setColumn(crackNormal3D, 1);
    locToGlob.setColumn(ez, 2);

    FloatArrayF<3> jump3DLoc = Tdot(locToGlob, jump3D);

    FloatArrayF<3> jump3DLocRenumbered = {jump3DLoc.at(3), jump3DLoc.at(1), jump3DLoc.at(2)};

    StructuralInterfaceMaterial *intMat = dynamic_cast<StructuralInterfaceMaterial*>(mpCZMat);
    if ( intMat == nullptr ) {
        OOFEM_ERROR("Failed to cast StructuralInterfaceMaterial*.")
    }
    auto TLocRenumbered = intMat->giveFirstPKTraction_3d(jump3DLocRenumbered, F, & iGP, tStep);

    FloatArrayF<3> TLoc = {TLocRenumbered.at(2), TLocRenumbered.at(3), TLocRenumbered.at(1)};

    auto T = dot(locToGlob, TLoc);

    oT = {T.at(1), T.at(2)};
}

void XfemStructuralElementInterface :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
{
    if ( !useNonStdCz() ) {

        if ( hasCohesiveZone() ) {
            FloatArray solVec;
            element->computeVectorOf(VM_Total, tStep, solVec);

            size_t numSeg = mpCZIntegrationRules.size();

            StructuralInterfaceMaterial *intMat = dynamic_cast<StructuralInterfaceMaterial*>(mpCZMat);
            if ( !intMat ) {
                OOFEM_ERROR("Failed to cast StructuralInterfaceMaterial*.")
            }

            for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
                for ( auto &gp: *mpCZIntegrationRules [ segIndex ] ) {
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

                    FloatMatrix K3DGlob;

                    FloatMatrix K2D;
                    K2D.resize(2, 2);
                    K2D.zero();

                    if ( intMat->hasAnalyticalTangentStiffness() ) {
                        ///////////////////////////////////////////////////
                        // Analytical tangent

                        auto K3DRenumbered = intMat->give3dStiffnessMatrix_dTdj(TangentStiffness, gp, tStep);

                        FloatMatrix K3D(3,3);
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
                        if ( ms == nullptr ) {
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
                        if ( ms == nullptr ) {
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
    } else {
        // Non-standard cz formulation.


        FloatArray solVec;
        element->computeVectorOf(VM_Total, tStep, solVec);

        size_t numSeg = mpCZIntegrationRules.size();

        StructuralFE2Material *fe2Mat = dynamic_cast<StructuralFE2Material*>(mpCZMat);
        if ( !fe2Mat ) {
            OOFEM_ERROR("Failed to cast StructuralFE2Material*.")
        }

        for ( size_t segIndex = 0; segIndex < numSeg; segIndex++ ) {
            //for ( auto &gp: *mpCZIntegrationRules [ segIndex ] ) {
            for ( int gpInd = 0; gpInd < mpCZIntegrationRules[ segIndex ]->giveNumberOfIntegrationPoints(); gpInd++ ) {

                GaussPoint *gp = mpCZIntegrationRules[ segIndex ]->getIntegrationPoint(gpInd);

                GaussPoint *bulk_gp = nullptr;
                StructuralMaterial *bulkMat = nullptr;
                if ( mIncludeBulkCorr ) {
                    bulk_gp = mpCZExtraIntegrationRules[ segIndex ]->getIntegrationPoint(gpInd);
                    bulkMat = dynamic_cast<StructuralMaterial*>( element->giveCrossSection()->giveMaterial(bulk_gp) );

                    if ( !bulkMat ) {
                        OOFEM_ERROR("Failed to fetch bulk material.")
                    }

                }


                StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*> ( gp->giveMaterialStatus() );

                if ( fe2ms == nullptr ) {
                    OOFEM_ERROR("The material status is not of an allowed type.")
                }

                ////////////////////////////////////////////////////////
                // Compute a (slightly modified) N-matrix

                FloatMatrix NMatrix;
                computeNCohesive(NMatrix, * gp, mCZEnrItemIndices [ segIndex ], mCZTouchingEnrItemIndices [ segIndex ]);


                ////////////////////////////////////////////////////////
                // Fetch normal
                FloatArray n( fe2ms->giveNormal() );

                // Traction part of tangent
                auto C = fe2Mat->give3dMaterialStiffnessMatrix(TangentStiffness, gp, tStep);

                ////////////////////////////////////////////////////////
                // Fetch L_s
                double l_s = computeEffectiveSveSize( fe2ms );

                FloatMatrix Ka(2,2);
                double a1 = 1.0;
                double a2 = 1.0;
                double a3 = 1.0;
                Ka(0,0) = (0.5/l_s)*(    C(0,0)*n(0)*n(0) + a2*C(0,5)*n(0)*n(1) + a1*C(5,0)*n(1)*n(0) + a3*C(5,5)*n(1)*n(1) ) +
                        (0.5/l_s)*(    C(0,0)*n(0)*n(0) + a2*C(0,5)*n(0)*n(1) + a1*C(5,0)*n(1)*n(0) + a3*C(5,5)*n(1)*n(1) );

                Ka(0,1) = (0.5/l_s)*( a2*C(0,5)*n(0)*n(0) +    C(0,1)*n(0)*n(1) + a3*C(5,5)*n(1)*n(0) + a1*C(5,1)*n(1)*n(1) ) +
                        (0.5/l_s)*( a2*C(0,5)*n(0)*n(0) +    C(0,1)*n(0)*n(1) + a3*C(5,5)*n(1)*n(0) + a1*C(5,1)*n(1)*n(1) );


                Ka(1,0) = (0.5/l_s)*( a1*C(5,0)*n(0)*n(0) + a3*C(5,5)*n(0)*n(1) +    C(1,0)*n(1)*n(0) + a2*C(1,5)*n(1)*n(1) ) +
                        (0.5/l_s)*( a1*C(5,0)*n(0)*n(0) + a3*C(5,5)*n(0)*n(1) +    C(1,0)*n(1)*n(0) + a2*C(1,5)*n(1)*n(1) );


                Ka(1,1) = (0.5/l_s)*( a3*C(5,5)*n(0)*n(0) + a1*C(5,1)*n(0)*n(1) + a2*C(1,5)*n(1)*n(0) +    C(1,1)*n(1)*n(1) ) +
                        (0.5/l_s)*( a3*C(5,5)*n(0)*n(0) + a1*C(5,1)*n(0)*n(1) + a2*C(1,5)*n(1)*n(0) +    C(1,1)*n(1)*n(1) );


                FloatMatrix tmp, tmp2;
                tmp.beProductOf(Ka, NMatrix);
                tmp2.beTProductOf(NMatrix, tmp);

                CrossSection *cs  = element->giveCrossSection();
                double thickness = cs->give(CS_Thickness, gp);
                double dA = thickness * gp->giveWeight();
                answer.add(dA, tmp2);



                FloatMatrix BAvg;

                if ( mIncludeBulkJump ) {
                    ////////////////////////////////////////////////////////
                    // Bulk contribution to SVE strain

                    // Crack gp coordinates
                    const FloatArray &xC = gp->giveGlobalCoordinates();

                    // For now, we will just perturb the coordinates of the GP to compute B^- and B^+ numerically.
                    double eps = 1.0e-6;
                    FloatArray xPert = xC;

                    xPert.add(eps, n);
                    FloatArray locCoordPert;
                    element->computeLocalCoordinates(locCoordPert, xPert);

                    FloatMatrix BPlus;
                    this->ComputeBOrBHMatrix(BPlus, *gp, *element, false, locCoordPert);

                    xPert = xC;
                    xPert.add(-eps, n);
                    element->computeLocalCoordinates(locCoordPert, xPert);

                    FloatMatrix BMinus;
                    this->ComputeBOrBHMatrix(BMinus, *gp, *element, false, locCoordPert);

                    BAvg = BPlus;
                    BAvg.add(1.0, BMinus);
                    BAvg.times(0.5);


                    FloatMatrix Kb(2,4); // Implicitly assumes plane strain. Fix later.

                    Kb(0,0) = C(0,0)*n(0) + C(5,0)*n(1);
                    Kb(0,1) = C(0,1)*n(0) + C(5,1)*n(1);
                    Kb(0,3) = C(0,5)*n(0) + C(5,5)*n(1);

                    Kb(1,0) = C(5,0)*n(0) + C(1,0)*n(1);
                    Kb(1,1) = C(5,1)*n(0) + C(1,1)*n(1);
                    Kb(1,3) = C(5,5)*n(0) + C(1,5)*n(1);

                    tmp.beProductOf(Kb, BAvg);
                    tmp2.beTProductOf(NMatrix, tmp);
                    answer.add(dA, tmp2);
                }

                ////////////////////////////////////////////////////////
                // Non-standard bulk contribution

                if ( mIncludeBulkCorr ) {

                    auto CBulk = bulkMat->give3dMaterialStiffnessMatrix(TangentStiffness, bulk_gp, tStep);

                    tmp.beTranspositionOf(tmp2);
                    answer.add(1.0*dA, tmp);


                    FloatMatrix C4(4,4);
                    C4(0,0) = C(0,0);
                    C4(0,1) = C(0,1);
                    C4(0,2) = C(0,2);
                    C4(0,3) = C(0,5);

                    C4(1,0) = C(1,0);
                    C4(1,1) = C(1,1);
                    C4(1,2) = C(1,2);
                    C4(1,3) = C(1,5);

                    C4(2,0) = C(2,0);
                    C4(2,1) = C(2,1);
                    C4(2,2) = C(2,2);
                    C4(2,3) = C(2,5);

                    C4(3,0) = C(5,0);
                    C4(3,1) = C(5,1);
                    C4(3,2) = C(5,2);
                    C4(3,3) = C(5,5);

                    tmp.beProductOf(C4, BAvg);
                    tmp2.beTProductOf(BAvg, tmp);
                    answer.add(1.0*dA*l_s, tmp2);

                    FloatMatrix C4Bulk(4,4);
                    C4Bulk(0,0) = CBulk(0,0);
                    C4Bulk(0,1) = CBulk(0,1);
                    C4Bulk(0,2) = CBulk(0,2);
                    C4Bulk(0,3) = CBulk(0,5);

                    C4Bulk(1,0) = CBulk(1,0);
                    C4Bulk(1,1) = CBulk(1,1);
                    C4Bulk(1,2) = CBulk(1,2);
                    C4Bulk(1,3) = CBulk(1,5);

                    C4Bulk(2,0) = CBulk(2,0);
                    C4Bulk(2,1) = CBulk(2,1);
                    C4Bulk(2,2) = CBulk(2,2);
                    C4Bulk(2,3) = CBulk(2,5);

                    C4Bulk(3,0) = CBulk(5,0);
                    C4Bulk(3,1) = CBulk(5,1);
                    C4Bulk(3,2) = CBulk(5,2);
                    C4Bulk(3,3) = CBulk(5,5);

                    tmp.beProductOf(C4Bulk, BAvg);
                    tmp2.beTProductOf(BAvg, tmp);
                    answer.add(-1.0*dA*l_s, tmp2);
                }

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
    if ( structEl == nullptr ) {
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

    for ( auto &gp: *element->giveIntegrationRule(0) ) {
        structEl->computeNmatrixAt(gp->giveNaturalCoordinates(), n);
        density = structEl->giveStructuralCrossSection()->give('d', gp);

        if ( ipDensity != nullptr ) {
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

void
XfemStructuralElementInterface :: initializeCZFrom(InputRecord &ir)
{
    int material = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_XfemElementInterface_CohesiveZoneMaterial);
    mCZMaterialNum = material;


    // Number of Gauss points used when integrating the cohesive zone
    IR_GIVE_OPTIONAL_FIELD(ir, mCSNumGaussPoints, _IFT_XfemElementInterface_NumIntPointsCZ);
    //    printf("mCSNumGaussPoints: %d\n", mCSNumGaussPoints );

    int planeStrainFlag = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, planeStrainFlag, _IFT_XfemElementInterface_PlaneStrain);
    if ( planeStrainFlag == 1 ) {
        mUsePlaneStrain = true;
    }
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

        mpCZMat = this->element->giveDomain()->giveMaterial(mCZMaterialNum);

        if ( mpCZMat == nullptr ) {
            OOFEM_ERROR("Failed to fetch pointer for mpCZMat.");
        }
    }
}

bool XfemStructuralElementInterface :: useNonStdCz()
{
    if ( element->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->element->giveDomain()->giveXfemManager();
        XfemStructureManager *xsMan = dynamic_cast<XfemStructureManager*>( xMan );

        if ( xsMan ) {
            if ( xsMan->giveUseNonStdCz() ) {
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        return false;
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
            if ( efStart != nullptr ) {
                const TipInfo &tipInfo = efStart->giveTipInfo();

                if ( ei->tipIsTouchingEI(tipInfo) ) {
                    //printf("Crack %d is touched by a tip on crack %d.\n", iEnrItemIndex, candidateIndex);
                    oTouchingEnrItemIndices.push_back(candidateIndex);
                }
            }

            // Check end tip
            EnrFrontIntersection *efEnd = dynamic_cast< EnrFrontIntersection * >( eiCandidate->giveEnrichmentFrontEnd() );
            if ( efEnd != nullptr ) {
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
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(primaryVarsToExport, numTotalNodes);

    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);


        nodesPassed = 1;

        for ( auto &tri: mSubTri ) {
            FloatArray triCenter = tri.giveVertex(1);
            triCenter.add( tri.giveVertex(2) );
            triCenter.add( tri.giveVertex(3) );
            triCenter.times(1.0 / 3.0);

            double meanEdgeLength = 0.0;
            meanEdgeLength += ( 1.0 / 3.0 ) * distance(tri.giveVertex(1), tri.giveVertex(2) );
            meanEdgeLength += ( 1.0 / 3.0 ) * distance(tri.giveVertex(2), tri.giveVertex(3) );
            meanEdgeLength += ( 1.0 / 3.0 ) * distance(tri.giveVertex(3), tri.giveVertex(1) );

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
//                    interp->evalN( N, locCoord, FEIElementGeometryWrapper(element) );
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
                            const auto &nodeCoord = dMan->giveCoordinates();

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
                        if ( geoEI != nullptr ) {
                            // TODO: Consider removing this special treatment. /ES
                            geoEI->giveGeometry()->computeTangentialSignDist(tangSignDist, x, arcPos);
                        }


                        if ( !evaluationSucceeded ) {
//                            printf("!evaluationSucceeded.\n");
                        }

//                        if ( ( tangSignDist > ( 1.0e-3 ) * meanEdgeLength && fabs(levelSetNormal) < ( 1.0e-2 ) * meanEdgeLength ) && evaluationSucceeded ) {
//                            joinNodes = false;
//                        }

//                        if ( ( tangSignDist < ( 1.0e-3 ) * meanEdgeLength || fabs(levelSetNormal) > ( 1.0e-2 ) * meanEdgeLength ) || !evaluationSucceeded ) {
                        if ( ( tangSignDist < ( 1.0e-3 ) * meanEdgeLength || fabs(levelSetNormal) > ( 1.0e-2 ) * meanEdgeLength ) && false ) {
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
                    vtkPieces [ 0 ].setPrimaryVarInNode(type, nodesPassed, valuearray);
                } else {
                    // TODO: Implement
                    printf("fieldNum: %d\n", fieldNum);
                }

                nodesPassed++;
            }
        }
    }


    // Export nodal variables from internal fields
    vtkPieces [ 0 ].setNumberOfInternalVarsToExport(internalVarsToExport, numTotalNodes);


    vtkPieces [ 0 ].setNumberOfCellVarsToExport(cellVarsToExport, numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);

        for ( size_t triInd = 1; triInd <= mSubTri.size(); triInd++ ) {

            FloatArray average;
            IntegrationRule *iRule = element->giveIntegrationRule(0);
            computeIPAverageInTriangle(average, iRule, element, type, tStep, mSubTri[triInd-1]);

            if ( average.giveSize() == 0 ) {
                VTKXMLExportModule :: computeIPAverage(average, iRule, element, type, tStep);
            }


            FloatArray averageVoigt;

            if ( average.giveSize() == 6 ) {

                averageVoigt.resize(9);

                averageVoigt.at(1) = average.at(1);
                averageVoigt.at(5) = average.at(2);
                averageVoigt.at(9) = average.at(3);
                averageVoigt.at(6) = averageVoigt.at(8) = average.at(4);
                averageVoigt.at(3) = averageVoigt.at(7) = average.at(5);
                averageVoigt.at(2) = averageVoigt.at(4) = average.at(6);
            } else {
                if ( average.giveSize() == 1 ) {
                    averageVoigt.resize(1);
                    averageVoigt.at(1) = average.at(1);
                }
            }

            vtkPieces [ 0 ].setCellVar(type, triInd, averageVoigt);
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
                            const auto &nodeCoord = dMan->giveCoordinates();
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
                            const auto &nodeCoord = dMan->giveCoordinates();
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
//            globCoord.resizeWithValues(2);

            if ( iTri.pointIsInTriangle(globCoord) ) {
                elem->giveIPValue(temp, ip, isType, tStep);
                gptot += ip->giveWeight();
                answer.add(ip->giveWeight(), temp);
            }
        }

        answer.times(1. / gptot);
    }

}

} /* namespace oofem */
