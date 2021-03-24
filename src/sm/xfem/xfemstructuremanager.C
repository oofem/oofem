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

#include "xfemstructuremanager.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "datareader.h"
#include "mathfem.h"

#include "sm/xfem/enrichmentitems/crack.h"
#include "xfem/enrichmentfronts/enrichmentfrontintersection.h"
#include "geometry.h"

#include "matforceevaluator.h"

namespace oofem {
REGISTER_XfemManager(XfemStructureManager)

XfemStructureManager :: XfemStructureManager(Domain *domain) :
    XfemManager(domain),
    mSplitCracks(false),
    mNonstandardCz(false),
    mMinCrackLength(0.0),
    mCrackMergeTol(0.0),
    mpMatForceEvaluator( new MaterialForceEvaluator() )
{}

XfemStructureManager :: ~XfemStructureManager()
{}

void XfemStructureManager :: initializeFrom(InputRecord &ir)
{
    XfemManager :: initializeFrom(ir);

    int splitCracks = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, splitCracks, _IFT_XfemStructureManager_splitCracks);
    if ( splitCracks == 1 ) {
        mSplitCracks = true;
    }

    int nonStdCz = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nonStdCz, _IFT_XfemStructureManager_nonstandardCZ);
    if ( nonStdCz == 1 ) {
        mNonstandardCz = true;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, mMinCrackLength, _IFT_XfemStructureManager_minCrackLength);

#if 0
    if ( mMinCrackLength < 1.0e-12 ) {
        printf("mMinCrackLength: %e\n", mMinCrackLength);
    }
#endif

    IR_GIVE_OPTIONAL_FIELD(ir, mCrackMergeTol, _IFT_XfemStructureManager_crackMergeTol);

    if ( mCrackMergeTol > 1.0e-12 ) {
        printf("mCrackMergeTol: %e\n", mCrackMergeTol);
    }
}

void XfemStructureManager :: postInitialize() {
    XfemManager :: postInitialize();

    if ( mSplitCracks ) {
            splitCracks();
        }
        mergeCloseCracks();
        updateNodeEnrichmentItemMap();

}

void XfemStructureManager :: giveInputRecord(DynamicInputRecord &input)
{
    XfemManager :: giveInputRecord(input);

    if ( mSplitCracks ) {
        input.setField(1, _IFT_XfemStructureManager_splitCracks);
    }

    input.setField(mMinCrackLength, _IFT_XfemStructureManager_minCrackLength);

    input.setField(mCrackMergeTol, _IFT_XfemStructureManager_crackMergeTol);

    if ( mNonstandardCz ) {
        input.setField(1, _IFT_XfemStructureManager_nonstandardCZ);
    }
}

int XfemStructureManager :: instanciateYourself(DataReader &dr)
{
    int result = XfemManager :: instanciateYourself(dr);

    
    return result;
}

void XfemStructureManager :: propagateFronts(bool &oAnyFronHasPropagated)
{
    oAnyFronHasPropagated = false;

    for ( auto &ei: enrichmentItemList ) {

        bool eiHasPropagated = false;
        ei->propagateFronts(eiHasPropagated);

        if ( eiHasPropagated ) {
            oAnyFronHasPropagated = true;
        }
    }

    updateNodeEnrichmentItemMap();
}


void XfemStructureManager :: updateYourself(TimeStep *tStep)
{
    XfemManager::updateYourself(tStep);
}

void XfemStructureManager :: splitCracks()
{
    // Loop over cracks
    for ( int i = 1; i <= giveNumberOfEnrichmentItems(); i++ ) {
        Crack *crack_i = dynamic_cast< Crack * >( this->giveEnrichmentItem(i) );
        if ( crack_i ) {
            // Check if crack i intersects with any of the cracks [1,i-1]:
            for ( int j = 1; j < i; j++ ) {
                // TODO: To improve performance, we may wish to use
                //       a tree structure here.
                bool splittedCrack = false;

                Crack *crack_j = dynamic_cast< Crack * >( this->giveEnrichmentItem(j) );
                if ( crack_j ) {
                    // If so, find the arc length positions of the intersections on crack i ...

                    std :: vector< FloatArray >intersectionPoints;
                    std :: vector< double >arcPositions_i, arcPositions_j;
                    crack_i->computeCrackIntersectionPoints(* crack_j, intersectionPoints, arcPositions_i);
                    crack_j->computeArcPoints(intersectionPoints, arcPositions_j);

                    const double arcLengthTol = 1.0e-6;

                    for ( int k = 0; k < int( arcPositions_i.size() ); k++ ) {
                        if ( arcPositions_i [ k ] < arcLengthTol || arcPositions_i [ k ] > ( 1.0 - arcLengthTol ) || arcPositions_j [ k ] < arcLengthTol || arcPositions_j [ k ] > ( 1.0 - arcLengthTol ) ) {
                            arcPositions_i.erase(arcPositions_i.begin() + k);
                            arcPositions_j.erase(arcPositions_j.begin() + k);
                            k--;
                        }
                    }

                    if ( arcPositions_i.size() > 0 ) {
                        arcPositions_i.insert(arcPositions_i.begin(), 0.0);
                        arcPositions_i.push_back(1.0);
                        arcPositions_j.insert(arcPositions_j.begin(), 0.0);
                        arcPositions_j.push_back(1.0);

                        for ( int k = 1; k < int( arcPositions_i.size() ); k++ ) {
                            // Only include segments of finite length
                            if ( fabs(arcPositions_i [ k ] - arcPositions_i [ k - 1 ]) > arcLengthTol ) {
                                //printf("arcPositions.size(): %lu\n", arcPositions.size() );

                                DynamicDataReader dataReader("xfem-crack");
                                crack_i->appendInputRecords(dataReader);
                                // ... split crack i at the intersection and add intersection enrichment
                                // fronts at the newly detected intersections.

                                int n1 = this->giveNumberOfEnrichmentItems() + 1;
                                //                        EnrichmentItem *newEI_1 = new Crack(n1, this, this->giveDomain() );
                                auto newCrack = std::make_unique<Crack>( n1, this, this->giveDomain() );

                                auto &ir = dataReader.giveInputRecord(DataReader :: IR_enrichItemRec, i);
                                newCrack->initializeFrom(ir);
                                newCrack->instanciateYourself(dataReader);

                                PolygonLine *new_pl = dynamic_cast< PolygonLine * >( newCrack->giveGeometry() );
                                //                                EDCrack *ed = dynamic_cast<EDCrack*>( newEI_1->giveEnrichmentDomain() );

                                if ( !new_pl ) {
                                    OOFEM_ERROR("Failed to cast PolygonLine *new_pl.")
                                } else {
                                    //printf("arcPositions_i[k-1]: %e arcPositions_i[k]: %e\n", arcPositions_i[k-1], arcPositions_i[k] );
                                    new_pl->cropPolygon(arcPositions_i [ k - 1 ], arcPositions_i [ k ]);


                                    PolygonLine *polygonLine_j = dynamic_cast< PolygonLine * >( crack_j->giveGeometry() );

                                    if ( !polygonLine_j ) {
                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_j.")
                                    }

                                    PolygonLine *polygonLine_i = dynamic_cast< PolygonLine * >( crack_i->giveGeometry() );

                                    if ( !polygonLine_i ) {
                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_i.")
                                    }

                                    // Take enrichment front tangent direction
                                    // as the normal direction of crack_j
                                    //                                    EnrichmentDomain_BG *ed_crack_j = dynamic_cast<EnrichmentDomain_BG*>( crack_j->giveEnrichmentDomain() );
                                    //                                    if(ed_crack_j == NULL) {
                                    //                                        OOFEM_ERROR("Failed to cast EnrichmentDomain_BG *ed_crack_j.")
                                    //                                    }
                                    //
                                    //                                    PolygonLine *polygonLine_j = dynamic_cast<PolygonLine*>( ed_crack_j->bg );
                                    //                                    if(polygonLine_j == NULL) {
                                    //                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_j.")
                                    //                                    }
                                    //
                                    //                                    EnrichmentDomain_BG *ed_crack_i = dynamic_cast<EnrichmentDomain_BG*>( crack_i->giveEnrichmentDomain() );
                                    //                                    if(ed_crack_i == NULL) {
                                    //                                        OOFEM_ERROR("Failed to cast EnrichmentDomain_BG *ed_crack_i.")
                                    //                                    }
                                    //
                                    //                                    PolygonLine *polygonLine_i = dynamic_cast<PolygonLine*>( ed_crack_i->bg );
                                    //                                    if(polygonLine_i == NULL) {
                                    //                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_i.")
                                    //                                    }


                                    // Change enrichment fronts
                                    if ( k - 1 > 0 ) {
                                        FloatArray frontTangent1;
                                        polygonLine_j->giveNormal(frontTangent1, arcPositions_j [ k - 1 ]);

                                        FloatArray crackTangent1;
                                        polygonLine_i->giveTangent(crackTangent1, arcPositions_i [ k - 1 ]);
                                        crackTangent1.times(-1.0);

                                        if ( frontTangent1.dotProduct(crackTangent1) < 0.0 ) {
                                            frontTangent1.times(-1.0);
                                        }
                                        auto ef = std::make_unique<EnrFrontIntersection>();
                                        ef->setTangent(frontTangent1);

                                        newCrack->setEnrichmentFrontStart(std::move(ef));
                                    }

                                    if ( k < int( arcPositions_i.size() ) - 1 ) {
                                        FloatArray frontTangent1;
                                        polygonLine_j->giveNormal(frontTangent1, arcPositions_j [ k ]);

                                        FloatArray crackTangent1;
                                        polygonLine_i->giveTangent(crackTangent1, arcPositions_i [ k ]);


                                        if ( frontTangent1.dotProduct(crackTangent1) < 0.0 ) {
                                            frontTangent1.times(-1.0);
                                        }
                                        auto ef = std::make_unique<EnrFrontIntersection>();
                                        ef->setTangent(frontTangent1);

                                        newCrack->setEnrichmentFrontEnd(std::move(ef));
                                    }
                                }
                                //this->enrichmentItemList[i-1] = std :: move(ei);

                                this->enrichmentItemList.push_back(nullptr);
                                newCrack->updateGeometry();
                                this->enrichmentItemList [ enrichmentItemList.size() - 1 ] = std :: move(newCrack);


                                splittedCrack = true;
                            }
                        }
                    }
                }

                if ( splittedCrack ) {
                    enrichmentItemList.erase(enrichmentItemList.begin() + i - 1);
                    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();
                    i--;
                    break;
                }
            }
        }
    }

    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();

    //printf("After splitting: Number of ei: %d\n", giveNumberOfEnrichmentItems() );

    removeShortCracks();

    for ( size_t i = 0; i < enrichmentItemList.size(); i++ ) {
        enrichmentItemList [ i ]->setNumber(i + 1);
    }

    for ( size_t i = 0; i < enrichmentItemList.size(); i++ ) {
        enrichmentItemList [ i ]->updateGeometry();
    }
}

void XfemStructureManager :: removeShortCracks()
{
    //printf("Entering XfemStructureManager :: removeShortCracks()\n");

    //printf("Number of ei before removal: %d\n", giveNumberOfEnrichmentItems());

    double l_tol = mMinCrackLength;

    for ( int i = 1; i <= giveNumberOfEnrichmentItems(); i++ ) {
        Crack *crack = dynamic_cast< Crack * >( this->giveEnrichmentItem(i) );
        if ( crack ) {
            double l = crack->computeLength();

            if ( l < l_tol ) {
                printf("Removing short crack with l: %e\n", l);

                // Explicitly erasing things is a mess...
                //crack->removeEnrichedDofs();
                //enrichmentItemList.erase(enrichmentItemList.begin() + i - 1);
                //i--;

                // ...therefore, just remove the geometry.
                PolygonLine *polygonLine = dynamic_cast< PolygonLine * >( crack->giveGeometry() );
                polygonLine->clear();
                //FloatArray tmp = {0.0, 0.0, 0.0};
                FloatArray tmp = {-1.0e3, -1.0e3, -1.0e3};
                polygonLine->insertVertexBack(tmp);

            }
        }
    }

    //printf("Number of ei after removal: %d\n", giveNumberOfEnrichmentItems());

}

bool XfemStructureManager :: tipsHaveOppositeDirection(EnrichmentFront *iEf1, EnrichmentFront *iEf2)
{
    const TipInfo &t1 = iEf1->giveTipInfo();
    const TipInfo &t2 = iEf2->giveTipInfo();

    return t1.mTangDir.dotProduct( t2.mTangDir ) <= 0.0;
}

void XfemStructureManager :: mergeCloseCracks()
{
    //printf("Entering XfemStructureManager :: mergeCloseCracks().\n");

    const double &dist_tol = mCrackMergeTol;

    // Loop over cracks and check if two crack tips are closer to
    // each other than a predefined distance. If so, merge the cracks.
    // Loop over cracks
    for ( int i = 1; i <= giveNumberOfEnrichmentItems(); i++ ) {
        Crack *crack_i = dynamic_cast< Crack * >( this->giveEnrichmentItem(i) );
        if ( crack_i ) {

            BasicGeometry *bg_i = crack_i->giveGeometry();
            TipInfo startTip_i, endTip_i;
            bg_i->giveTips(startTip_i, endTip_i) ;

            const FloatArray &ps_i =startTip_i.mGlobalCoord;
            const FloatArray &pe_i =endTip_i.mGlobalCoord;

            PolygonLine *polygonLine_i = dynamic_cast< PolygonLine * >( crack_i->giveGeometry() );

            if ( !polygonLine_i ) {
                OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_i.")
            }

            for ( int j = i+1; j <= giveNumberOfEnrichmentItems(); j++ ) {
                // TODO: To improve performance, we may wish to use
                //       a tree structure here.
//                bool mergedCrack = false;

                Crack *crack_j = dynamic_cast< Crack * >( this->giveEnrichmentItem(j) );

                if ( crack_i->computeLength() < 1.0e-18 || crack_j->computeLength() < 1.0e-18 ) {
                    continue;
                }

                if ( crack_j ) {
                    BasicGeometry *bg_j = crack_j->giveGeometry();
                    TipInfo startTip_j, endTip_j;
                    bg_j->giveTips(startTip_j, endTip_j) ;

                    const FloatArray &ps_j =startTip_j.mGlobalCoord;
                    const FloatArray &pe_j =endTip_j.mGlobalCoord;

                    PolygonLine *polygonLine_j = dynamic_cast< PolygonLine * >( crack_j->giveGeometry() );

                    if ( !polygonLine_j ) {
                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_j.")
                    }

                    ////////////////////////////////////////////////////////////////
                    if ( distance(ps_i, ps_j) < dist_tol ) {
                        printf("distance(ps_i, ps_j) < dist_tol\n");

                        // Merging is not reasonable if the tips are close to parallel
                        if ( !tipsHaveOppositeDirection(crack_i->giveEnrichmentFrontStart(), crack_j->giveEnrichmentFrontStart()) ) {
                            printf("Preventing merge due to parallel tips.\n");
                        } else {
                            // Append points to the start of polygonLine_i
                            int n = polygonLine_j->giveNrVertices();
                            for (int k = 1; k <= n; k++) {
                                polygonLine_i->insertVertexFront( polygonLine_j->giveVertex(k) );
                            }

                            polygonLine_i->removeDuplicatePoints(1.0e-18);


                            polygonLine_j->clear();
                            //FloatArray tmp = { 0.0, 0.0, 0.0};
                            FloatArray tmp = {-1.0e3, -1.0e3, -1.0e3};
                            polygonLine_j->insertVertexBack(tmp);

                            // Fix tips
                            EnrichmentFront *ef_tmp = crack_i->giveEnrichmentFrontStart();
                            crack_i->setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront>(crack_j->giveEnrichmentFrontEnd()), false );
                            crack_j->setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront>(ef_tmp), false);


                            //mergedCrack = true;
                            break;
                        }

                    }

                    ////////////////////////////////////////////////////////////////
                    if ( distance(ps_i, pe_j) < dist_tol ) {
                        printf("distance(ps_i, pe_j) < dist_tol\n");

#if 1
                        if ( !tipsHaveOppositeDirection(crack_i->giveEnrichmentFrontStart(), crack_j->giveEnrichmentFrontEnd()) ) {
                            printf("Preventing merge due to parallel tips.\n");
                        } else {

                            // Append points to the start of polygonLine_i
                            int n = polygonLine_j->giveNrVertices();
                            for(int k = n; k > 0; k--) {
                                polygonLine_i->insertVertexFront( polygonLine_j->giveVertex(k) );
                            }

                            polygonLine_i->removeDuplicatePoints(1.0e-18);


                            polygonLine_j->clear();
                            //FloatArray tmp = {0.0, 0.0, 0.0};
                            FloatArray tmp = {-1.0e3, -1.0e3, -1.0e3};
                            polygonLine_j->insertVertexBack(tmp);

                            // Fix tips
                            ///@todo Quite ugly; come up with bettermethods for swapping without releasing from the ptrs.
                            EnrichmentFront *ef_tmp = crack_i->giveEnrichmentFrontStart();
                            crack_i->setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront>(crack_j->giveEnrichmentFrontStart()), false );
                            crack_j->setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront>(ef_tmp), false);

                            //mergedCrack = true;
                            break;
                        }
#endif
                    }


                    ////////////////////////////////////////////////////////////////
                    if ( distance(pe_i, ps_j) < dist_tol ) {
                        printf("distance(pe_i, ps_j) < dist_tol\n");

                        if ( !tipsHaveOppositeDirection(crack_i->giveEnrichmentFrontEnd(), crack_j->giveEnrichmentFrontStart()) ) {
                            printf("Preventing merge due to parallel tips.\n");
                        } else {

                            // Append points to the end of polygonLine_i
                            int n = polygonLine_j->giveNrVertices();
                            for(int k = 1; k <= n; k++) {
                                polygonLine_i->insertVertexBack( polygonLine_j->giveVertex(k) );
                            }

                            polygonLine_i->removeDuplicatePoints(1.0e-18);


                            polygonLine_j->clear();
                            //FloatArray tmp = {0.0, 0.0, 0.0};
                            FloatArray tmp = {-1.0e3, -1.0e3, -1.0e3};
                            polygonLine_j->insertVertexBack(tmp);


                            // Fix tips
                            EnrichmentFront *ef_tmp = crack_i->giveEnrichmentFrontEnd();
                            crack_i->setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront>(crack_j->giveEnrichmentFrontEnd()), false );
                            crack_j->setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront>(ef_tmp), false);

                            //mergedCrack = true;
                            break;
                        }
                    }

                    ////////////////////////////////////////////////////////////////
                    if ( distance(pe_i, pe_j) < dist_tol ) {
                        printf("distance(pe_i, pe_j) < dist_tol\n");

                        if ( !tipsHaveOppositeDirection(crack_i->giveEnrichmentFrontEnd(), crack_j->giveEnrichmentFrontEnd()) ) {
                            printf("Preventing merge due to parallel tips.\n");
                        } else {

                            // Append points to the end of polygonLine_i
                            int n = polygonLine_j->giveNrVertices();
                            for(int k = n; k > 0; k--) {
                                polygonLine_i->insertVertexBack( polygonLine_j->giveVertex(k) );
                            }

                            polygonLine_i->removeDuplicatePoints(1.0e-18);

                            polygonLine_j->clear();
                            //FloatArray tmp = {0.0, 0.0, 0.0};
                            FloatArray tmp = {-1.0e3, -1.0e3, -1.0e3};
                            polygonLine_j->insertVertexBack(tmp);


                            // Fix tips
                            EnrichmentFront *ef_tmp = crack_i->giveEnrichmentFrontEnd();
                            crack_i->setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront>(crack_j->giveEnrichmentFrontStart()), false );
                            crack_j->setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront>(ef_tmp), false);

                            //mergedCrack = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    removeShortCracks();

    for ( size_t i = 0; i < enrichmentItemList.size(); i++ ) {
        enrichmentItemList [ i ]->setNumber(i + 1);
    }


    for ( size_t i = 0; i < enrichmentItemList.size(); i++ ) {
        enrichmentItemList [ i ]->updateGeometry();
    }

    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();

}


double XfemStructureManager :: computeTotalCrackLength()
{
    double l_tot = 0.0;

    for ( int i = 1; i <= giveNumberOfEnrichmentItems(); i++ ) {
        auto crack = dynamic_cast< Crack * >( this->giveEnrichmentItem(i) );
        if ( crack ) {
            l_tot += crack->computeLength();
        }
    }

    return l_tot;
}


} /* namespace oofem */
