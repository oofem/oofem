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

#include "ncprincipalstrain.h"

#include "error.h"
#include "xfem/enrichmentitem.h"
#include "domain.h"
#include "element.h"
#include "gausspoint.h"

#include "Materials/structuralms.h"
#include "Materials/structuralmaterial.h"

#include "xfem/enrichmentitems/crack.h"
#include "xfem/xfemmanager.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentfronts/enrichmentfrontlinbranchfunconeel.h"
#include "xfem/enrichmentfronts/enrichmentfrontcohesivebranchfunconeel.h"
#include "xfem/propagationlaws/plhoopstresscirc.h"
#include "xfem/propagationlaws/plprincipalstrain.h"
#include "dynamicdatareader.h"
#include "dynamicinputrecord.h"
#include "geometry.h"
#include "classfactory.h"
#include "spatiallocalizer.h"
#include "crosssection.h"

#include <memory>

namespace oofem {
REGISTER_NucleationCriterion(NCPrincipalStrain)

NCPrincipalStrain::NCPrincipalStrain(Domain *ipDomain):
NucleationCriterion(ipDomain),
mStrainThreshold(0.0),
mInitialCrackLength(0.0),
mCutOneEl(false),
mCrossSectionInd(1),
mIncrementLength(1.0)
{

}

NCPrincipalStrain::~NCPrincipalStrain() {

}

std::vector<std::unique_ptr<EnrichmentItem>> NCPrincipalStrain::nucleateEnrichmentItems() {


	SpatialLocalizer *octree = this->mpDomain->giveSpatialLocalizer();
	XfemManager *xMan = mpDomain->giveXfemManager();

	std::vector<std::unique_ptr<EnrichmentItem>> eiList;

	// Center coordinates of newly inserted cracks
	std::vector<FloatArray> center_coord_inserted_cracks;

	// Loop over all elements and all bulk GP.
	for(auto &el : mpDomain->giveElements() ) {

		int numIR = el->giveNumberOfIntegrationRules();

		int csNum = el->giveCrossSection()->giveNumber();

		if(csNum == mCrossSectionInd) {

			for(int irInd = 0; irInd < numIR; irInd++) {
				IntegrationRule *ir = el->giveIntegrationRule(irInd);


				int numGP = ir->giveNumberOfIntegrationPoints();

				for(int gpInd = 0; gpInd < numGP; gpInd++) {
					GaussPoint *gp = ir->getIntegrationPoint(gpInd);


						StructuralMaterialStatus *ms = dynamic_cast<StructuralMaterialStatus*>(gp->giveMaterialStatus());

						if(ms != NULL) {

							const FloatArray &strain = ms->giveTempStrainVector();

							FloatArray principalVals;
							FloatMatrix principalDirs;
							StructuralMaterial::computePrincipalValDir(principalVals, principalDirs, strain, principal_strain);

							if(principalVals[0] > mStrainThreshold) {

								FloatArray crackNormal;
								crackNormal.beColumnOf(principalDirs, 1);
		//						printf("crackNormal: "); crackNormal.printYourself();

								FloatArray crackTangent = {-crackNormal(1), crackNormal(0)};
								crackTangent.normalize();
		//						printf("crackTangent: "); crackTangent.printYourself();


								// Create geometry
								FloatArray pc = {gp->giveGlobalCoordinates()(0), gp->giveGlobalCoordinates()(1)};
		//						printf("Global coord: "); pc.printYourself();


								FloatArray ps = pc;
								ps.add(-0.5*mInitialCrackLength, crackTangent);

								FloatArray pe = pc;
								pe.add(0.5*mInitialCrackLength, crackTangent);

								if(mCutOneEl) {
									// If desired, ensure that the crack cuts exactly one element.
									Line line(ps, pe);
									std::vector<FloatArray> intersecPoints;
		//							line.computeIntersectionPoints(el.get(), intersecPoints);

									for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
										int n1 = i;
										int n2 = 0;
										if ( i < el->giveNumberOfDofManagers() ) {
											n2 = i + 1;
										} else {
											n2 = 1;
										}
									}

									if(intersecPoints.size() == 2) {
										ps = std::move(intersecPoints[0]);
										pe = std::move(intersecPoints[1]);
									}
									else {
										OOFEM_ERROR("intersecPoints.size() != 2")
									}
								}

								FloatArray points = {ps(0), ps(1), pc(0), pc(1), pe(0), pe(1)};


								// Check if nucleation is allowed, by checking for already existing cracks close to the GP.
								// Idea: Nucleation is not allowed if we are within an enriched element. In this way, branching is not
								// completely prohibited, but we avoid initiating multiple similar cracks.
								bool insertionAllowed = true;

								Element *el_s = octree->giveElementContainingPoint(ps);
								if(el_s) {
									if( xMan->isElementEnriched(el_s) ) {
										insertionAllowed = false;
									}
								}

								Element *el_c = octree->giveElementContainingPoint(pc);
								if(el_c) {
									if( xMan->isElementEnriched(el_c) ) {
										insertionAllowed = false;
									}
								}

								Element *el_e = octree->giveElementContainingPoint(pe);
								if(el_e) {
									if( xMan->isElementEnriched(el_e) ) {
										insertionAllowed = false;
									}
								}

								for(const auto &x: center_coord_inserted_cracks) {
									if( x.distance(pc) <  2.0*mInitialCrackLength) {
										insertionAllowed = false;
										break;
										printf("Preventing insertion.\n");
									}
								}

								if(insertionAllowed) {
									int n = xMan->giveNumberOfEnrichmentItems() + 1;
									std::unique_ptr<Crack> crack(new Crack(n, xMan, mpDomain));


									// Geometry
									std::unique_ptr<BasicGeometry> geom = std::unique_ptr<BasicGeometry>(new PolygonLine());
									geom->insertVertexBack(ps);
									geom->insertVertexBack(pc);
									geom->insertVertexBack(pe);
									crack->setGeometry(std::move(geom));

									// Enrichment function
									EnrichmentFunction *ef = new HeavisideFunction(1, mpDomain);
									crack->setEnrichmentFunction(ef);

									// Enrichment fronts
//									EnrichmentFront *efStart = new EnrFrontLinearBranchFuncOneEl();
									EnrichmentFront *efStart = new EnrFrontCohesiveBranchFuncOneEl();
									crack->setEnrichmentFrontStart(efStart);

//									EnrichmentFront *efEnd = new EnrFrontLinearBranchFuncOneEl();
									EnrichmentFront *efEnd = new EnrFrontCohesiveBranchFuncOneEl();
									crack->setEnrichmentFrontEnd(efEnd);




									///////////////////////////////////////
									// Propagation law

									// Options
									PLPrincipalStrain *pl = new PLPrincipalStrain();
									pl->setRadius(1.0);
									pl->setIncrementLength(mIncrementLength);
									pl->setStrainThreshold(mStrainThreshold);

									crack->setPropagationLaw(pl);

									crack->updateDofIdPool();

									center_coord_inserted_cracks.push_back(pc);
									eiList.push_back( std::unique_ptr<EnrichmentItem>(std::move(crack)) );

									printf("NCPrincipalStrain: Nucleating a crack. principalVals[0]: %e\n", principalVals[0] );

									// We only introduce one crack per element in a single time step.
									break;
								}
							}
						}

				}
			}
		} // If correct csNum
	}


	return std::move( eiList );
}


IRResultType NCPrincipalStrain::initializeFrom(InputRecord *ir) {

    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mStrainThreshold, _IFT_NCPrincipalStrain_StrainThreshold);
    printf("mStrainThreshold: %e\n", mStrainThreshold);

    IR_GIVE_FIELD(ir, mInitialCrackLength, _IFT_NCPrincipalStrain_InitialCrackLength);
    printf("mInitialCrackLength: %e\n", mInitialCrackLength);

    IR_GIVE_FIELD(ir, mIncrementLength, _IFT_NCPrincipalStrain_IncrementLength);
    printf("mIncrementLength: %e\n", mIncrementLength);

    return NucleationCriterion::initializeFrom(ir);
}

void NCPrincipalStrain :: appendInputRecords(DynamicDataReader &oDR)
{
    DynamicInputRecord *ir = new DynamicInputRecord();

    ir->setRecordKeywordField( this->giveInputRecordName(), 1 );

    ir->setField(mStrainThreshold, _IFT_NCPrincipalStrain_StrainThreshold);
    ir->setField(mInitialCrackLength, _IFT_NCPrincipalStrain_InitialCrackLength);
    ir->setField(mIncrementLength, _IFT_NCPrincipalStrain_IncrementLength);

    oDR.insertInputRecord(DataReader :: IR_crackNucleationRec, ir);

    // Enrichment function
    DynamicInputRecord *efRec = new DynamicInputRecord();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, efRec);
}

} /* namespace oofem */
