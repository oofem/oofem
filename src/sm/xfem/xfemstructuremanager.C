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

#include "xfem/enrichmentdomain.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/enrichmentfronts/enrichmentfrontintersection.h"

namespace oofem {
REGISTER_XfemManager(XfemStructureManager)

XfemStructureManager::XfemStructureManager(Domain *domain):
XfemManager(domain),
mSplitCracks(false)
{

}

XfemStructureManager::~XfemStructureManager() {

}

IRResultType XfemStructureManager :: initializeFrom(InputRecord *ir)
{
    XfemManager :: initializeFrom(ir);

    IRResultType result; // Required by IR_GIVE_FIELD macro
    int splitCracks = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, splitCracks, _IFT_XfemStructureManager_splitCracks);
    if ( splitCracks == 1 ) {
        mSplitCracks = true;
    }

    return IRRT_OK;
}

void XfemStructureManager :: giveInputRecord(DynamicInputRecord &input)
{
    XfemManager :: giveInputRecord(input);

    if ( mSplitCracks ) {
        input.setField(1, _IFT_XfemStructureManager_splitCracks);
    }
}

int XfemStructureManager :: instanciateYourself(DataReader *dr)
{
    int result = XfemManager::instanciateYourself(dr);

    if(mSplitCracks) {
        splitCracks();
    }

    updateNodeEnrichmentItemMap();

    return result;
}

void XfemStructureManager :: splitCracks()
{
    std::vector<size_t> eiToRemove;

    // Loop over cracks
    for(int i = 1; i <= giveNumberOfEnrichmentItems(); i++) {

        Crack *crack_i = dynamic_cast<Crack*>( this->giveEnrichmentItem(i) );
        if( crack_i != NULL ) {

            // Check if crack i intersects with any of the cracks [1,i-1]:
            for(int j = 1; j < i; j++) {
                // TODO: To improve performance, we may wish to use
                //       a tree structure here.
                bool splittedCrack = false;

                Crack *crack_j = dynamic_cast<Crack*>( this->giveEnrichmentItem(j) );
                if( crack_j != NULL ) {

                    // If so, find the arc length positions of the intersections on crack i ...

                    std::vector<FloatArray> intersectionPoints;
                    std::vector<double> arcPositions_i, arcPositions_j;
                    crack_i->computeIntersectionPoints(*crack_j, intersectionPoints, arcPositions_i);
                    crack_j->computeArcPoints(intersectionPoints, arcPositions_j);

                    const double arcLengthTol = 1.0e-9;

                    for(int k = 0; k < int(arcPositions_i.size()); k++) {

                        if( arcPositions_i[k] < arcLengthTol || arcPositions_i[k] > (1.0-arcLengthTol) || arcPositions_j[k] < arcLengthTol || arcPositions_j[k] > (1.0-arcLengthTol) ) {
                            arcPositions_i.erase( arcPositions_i.begin()+k );
                            arcPositions_j.erase( arcPositions_j.begin()+k );
                            k--;
                        }
                    }

                    if(arcPositions_i.size() > 0 ) {

                        arcPositions_i.insert(arcPositions_i.begin(), 0.0);
                        arcPositions_i.push_back(1.0);
                        arcPositions_j.insert(arcPositions_j.begin(), 0.0);
                        arcPositions_j.push_back(1.0);

                        for(int k = 1; k < int(arcPositions_i.size()); k++) {

                            // Only include segments of finite length
                            if( fabs(arcPositions_i[k] - arcPositions_i[k-1]) > arcLengthTol ) {

                                //printf("arcPositions.size(): %lu\n", arcPositions.size() );

                                DynamicDataReader dataReader;
                                crack_i->appendInputRecords(dataReader);
                                // ... split crack i at the intersection and add intersection enrichment
                                // fronts at the newly detected intersections.

                                int n1 = this->giveNumberOfEnrichmentItems()+1;
        //                        EnrichmentItem *newEI_1 = new Crack(n1, this, this->giveDomain() );
                                std :: unique_ptr< EnrichmentItem > newEI_1( new Crack(n1, this, this->giveDomain() ) );

                                InputRecord *ir = dataReader.giveInputRecord(DataReader :: IR_enrichItemRec, i);
                                newEI_1->initializeFrom( ir );
                                newEI_1->instanciateYourself(&dataReader);

                                EDCrack *ed = dynamic_cast<EDCrack*>( newEI_1->giveEnrichmentDomain() );
                                if(ed != NULL) {
                                    //printf("arcPositions[k-1]: %e arcPositions[k]: %e\n", arcPositions[k-1], arcPositions[k] );
                                    ed->cropPolygon(arcPositions_i[k-1], arcPositions_i[k]);


                                    // Take enrichment front tangent direction
                                    // as the normal direction of crack_j
                                    EnrichmentDomain_BG *ed_crack_j = dynamic_cast<EnrichmentDomain_BG*>( crack_j->giveEnrichmentDomain() );
                                    if(ed_crack_j == NULL) {
                                        OOFEM_ERROR("Failed to cast EnrichmentDomain_BG *ed_crack_j.")
                                    }

                                    PolygonLine *polygonLine_j = dynamic_cast<PolygonLine*>( ed_crack_j->bg );
                                    if(polygonLine_j == NULL) {
                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_j.")
                                    }

                                    EnrichmentDomain_BG *ed_crack_i = dynamic_cast<EnrichmentDomain_BG*>( crack_i->giveEnrichmentDomain() );
                                    if(ed_crack_i == NULL) {
                                        OOFEM_ERROR("Failed to cast EnrichmentDomain_BG *ed_crack_i.")
                                    }

                                    PolygonLine *polygonLine_i = dynamic_cast<PolygonLine*>( ed_crack_i->bg );
                                    if(polygonLine_i == NULL) {
                                        OOFEM_ERROR("Failed to cast PolygonLine *polygonLine_i.")
                                    }


                                    // Change enrichment fronts
                                    if( k-1 > 0 ) {

                                        FloatArray frontTangent1;
                                        polygonLine_j->giveNormal(frontTangent1, arcPositions_j[k-1]);

                                        FloatArray crackTangent1;
                                        polygonLine_i->giveTangent(crackTangent1, arcPositions_i[k-1]);
                                        crackTangent1.times(-1.0);

                                        EnrFrontIntersection *ef = new EnrFrontIntersection();

                                        if( frontTangent1.dotProduct(crackTangent1) < 0.0 ) {
                                            frontTangent1.times(-1.0);
                                        }
                                        ef->setTangent(frontTangent1);

                                        newEI_1->setEnrichmentFrontStart(ef);
                                    }

                                    if( k < int(arcPositions_i.size())-1) {

                                        FloatArray frontTangent1;
                                        polygonLine_j->giveNormal(frontTangent1, arcPositions_j[k]);

                                        FloatArray crackTangent1;
                                        polygonLine_i->giveTangent(crackTangent1, arcPositions_i[k]);

                                        EnrFrontIntersection *ef = new EnrFrontIntersection();

                                        if( frontTangent1.dotProduct(crackTangent1) < 0.0 ) {
                                            frontTangent1.times(-1.0);
                                        }
                                        ef->setTangent(frontTangent1);

                                        newEI_1->setEnrichmentFrontEnd(ef);
                                    }



                                }
                                //this->enrichmentItemList[i-1] = std :: move(ei);

                                this->enrichmentItemList.push_back(NULL);
                                newEI_1->updateGeometry();
                                this->enrichmentItemList[ enrichmentItemList.size()-1 ] = std :: move(newEI_1);


                                splittedCrack = true;
                            }
                        }
                    }
                }

                if(splittedCrack) {
                    enrichmentItemList.erase( enrichmentItemList.begin() + i-1 );
                    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();
                    i--;
                    break;
                }
            }

        }
    }

    numberOfEnrichmentItems = giveNumberOfEnrichmentItems();

    for( size_t i = 0; i < enrichmentItemList.size(); i++ ) {
        enrichmentItemList[i]->setNumber(i+1);
        enrichmentItemList[i]->writeVtkDebug();
        enrichmentItemList[i]->updateGeometry();
    }

}

} /* namespace oofem */
