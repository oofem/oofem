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

/*
 * gnuplotexportmodule.C
 *
 *  Created on: Jan 29, 2014
 *      Author: svennine
 */

#include "gnuplotexportmodule.h"
#include "classfactory.h"
#include "valuemodetype.h"
#include "structengngmodel.h"
#include "outputmanager.h"
#include "dofmanager.h"
#include "boundarycondition.h"
#include "xfem/enrichmentitem.h"
#include "xfem/xfemmanager.h"
#include "structuralinterfacematerialstatus.h"
#include "xfem/enrichmentdomain.h"
#include "xfem/XFEMDebugTools.h"
#include "prescribedgradient.h"
#include "prescribedgradientbcneumann.h"
#include "gausspoint.h"
#include "timestep.h"
#include "xfem/enrichmentitems/crack.h"

#include <sstream>

namespace oofem {
REGISTER_ExportModule(GnuplotExportModule)

GnuplotExportModule::GnuplotExportModule(int n, EngngModel *e):
ExportModule(n, e),
mExportReactionForces(false),
mExportBoundaryConditions(false)
{

}

GnuplotExportModule::~GnuplotExportModule() {

}

IRResultType GnuplotExportModule::initializeFrom(InputRecord *ir)
{
    mExportReactionForces = ir->hasField(_IFT_GnuplotExportModule_ReactionForces);
    mExportBoundaryConditions = ir->hasField(_IFT_GnuplotExportModule_BoundaryConditions);
    return ExportModule::initializeFrom(ir);
}

void GnuplotExportModule::doOutput(TimeStep *tStep, bool forcedOutput)
{
    if (!(testTimeStepOutput(tStep) || forcedOutput)) {
        return;
    }

	// Export the sum of reaction forces for each Dirichlet BC
	if(mExportReactionForces) {
		outputReactionForces(tStep);
	}

	Domain *domain = emodel->giveDomain(1);

	// Export output from boundary conditions
	if(mExportBoundaryConditions) {
		int numBC = domain->giveNumberOfBoundaryConditions();

		for(int i = 1; i <= numBC; i++) {

			PrescribedGradient *presGradBC = dynamic_cast<PrescribedGradient*>( domain->giveBc(i) );
			if(presGradBC != NULL) {
				outputBoundaryCondition(*presGradBC, tStep);
			}


			PrescribedGradientBCNeumann *presGradBCNeumann = dynamic_cast<PrescribedGradientBCNeumann*>( domain->giveBc(i) );
			if(presGradBCNeumann != NULL) {
				outputBoundaryCondition(*presGradBCNeumann, tStep);
			}

		}
	}

	if(domain->hasXfemManager()) {
		XfemManager *xMan = domain->giveXfemManager();

		int numEI = xMan->giveNumberOfEnrichmentItems();

		for(int i = 1; i <= numEI; i++) {
			EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
			ei->callGnuplotExportModule(*this);
		}
	}

}

void GnuplotExportModule::initialize()
{

}

void GnuplotExportModule::terminate()
{

}

/////////////////////////////////////////////////
// Help functions
void GnuplotExportModule::outputReactionForces(TimeStep *tStep)
{
	// Add sum of reaction forces to arrays
    // Compute sum of reaction forces for each BC number
	Domain *domain = emodel->giveDomain(1);
	StructuralEngngModel *seMod = dynamic_cast<StructuralEngngModel* >(emodel);
	if(seMod == NULL) {
		OOFEM_ERROR("failed to cast to StructuralEngngModel.\n");
	}

    IntArray ielemDofMask;
    FloatArray reactions;
    IntArray dofManMap, dofMap, eqnMap;

    // test if solution step output is active
    if ( !domain->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    // map contains corresponding dofmanager and dofs numbers corresponding to prescribed equations
    // sorted according to dofmanger number and as a minor crit. according to dof number
    // this is necessary for extractor, since the sorted output is expected
    seMod->buildReactionTable(dofManMap, dofMap, eqnMap, tStep, 1);

    // compute reaction forces
    seMod->computeReaction(reactions, tStep, 1);

    // Find highest index of prescribed dofs
    int maxIndPresDof = 0;
    for ( int i = 1; i <= dofManMap.giveSize(); i++ ) {
        maxIndPresDof = std::max(maxIndPresDof, dofMap.at(i));
    }

    int numBC = domain->giveNumberOfBoundaryConditions();

    while ( mReactionForceHistory.size() < size_t(numBC) ) {
    	std::vector<FloatArray> emptyArray;
    	mReactionForceHistory.push_back( emptyArray );
    }

    while ( mDispHist.size() < size_t(numBC) ) {
    	std::vector<double> emptyArray;
    	mDispHist.push_back( emptyArray );
    }

    for(int bcInd = 0; bcInd < numBC; bcInd++) {
    	FloatArray fR(maxIndPresDof), disp(numBC);
    	fR.zero();


        for ( int i = 1; i <= dofManMap.giveSize(); i++ ) {
        	DofManager *dMan = domain->giveDofManager( dofManMap.at(i) );
        	Dof *dof = dMan->giveDof( dofMap.at(i) );

        	if( dof->giveBcId() == bcInd+1) {
        		fR.at( dofMap.at(i) ) += reactions.at( eqnMap.at(i) );

        		// Slightly dirty
        		BoundaryCondition *bc = dynamic_cast<BoundaryCondition*> (domain->giveBc(bcInd+1));
        		if(bc != NULL) {
        			disp.at(bcInd+1) = bc->give(dof, VM_Total, tStep);
        		}
        	}
        }

        mDispHist[bcInd].push_back(disp.at(bcInd+1));
        mReactionForceHistory[bcInd].push_back(fR);



        // X
        FILE * pFileX;
        char fileNameX[100];
        sprintf(fileNameX, "ReactionForceGnuplotBC%dX.dat", bcInd+1);
        pFileX = fopen ( fileNameX , "wb" );

        fprintf(pFileX, "# %s\n", fileNameX);
        fprintf(pFileX, "#u Fx\n");
        for ( size_t j = 0; j < mDispHist[bcInd].size(); j++ ) {
    		fprintf(pFileX, "%e %e\n", mDispHist[bcInd][j], mReactionForceHistory[bcInd][j].at(1) );
        }

        fclose(pFileX);

        // Y
        FILE * pFileY;
        char fileNameY[100];
        sprintf(fileNameY, "ReactionForceGnuplotBC%dY.dat", bcInd+1);
        pFileY = fopen ( fileNameY , "wb" );

        fprintf(pFileY, "# %s\n", fileNameY);
        fprintf(pFileY, "#u Fx\n");
        for ( size_t j = 0; j < mDispHist[bcInd].size(); j++ ) {
        	if( mReactionForceHistory[bcInd][j].giveSize() >= 2 ) {
        		fprintf(pFileY, "%e %e\n", mDispHist[bcInd][j], mReactionForceHistory[bcInd][j].at(2) );
        	}
        }

        fclose(pFileY);

    }
}

void GnuplotExportModule::outputXFEM(EnrichmentItem &iEI)
{

}

void GnuplotExportModule::outputXFEM(Crack &iCrack)
{
	const std::vector<GaussPoint*> &czGaussPoints = iCrack.giveCohesiveZoneGaussPoints();
	size_t numPoints = czGaussPoints.size();

	std::vector<double> arcLengthPositions, normalJumps, tangJumps;

	const EnrichmentDomain *ed = iCrack.giveEnrichmentDomain();

	for(size_t i = 0; i < numPoints; i++) {
		GaussPoint *gp = czGaussPoints[i];

		StructuralInterfaceMaterialStatus *matStat = dynamic_cast<StructuralInterfaceMaterialStatus*> ( gp->giveMaterialStatus() );
		if(matStat != NULL) {

			// Compute arc length position of the Gauss point
			const FloatArray &coord = *(gp->giveCoordinates());
			double tangDist = 0.0, arcPos = 0.0;
			ed->computeTangentialSignDist(tangDist, coord, arcPos);
			arcLengthPositions.push_back(arcPos);

			// Compute displacement jump in normal and tangential direction
			// Local numbering: (tang_z, tang, normal)
			const FloatArray &jumpLoc 		= matStat->giveJump();

			double normalJump = jumpLoc.at(3);
			normalJumps.push_back(normalJump);


			tangJumps.push_back( jumpLoc.at(2) );
		}
	}



	Domain *domain = emodel->giveDomain(1);
    XfemManager *xMan = domain->giveXfemManager();
    if ( xMan != NULL ) {
    	double time = 0.0;

    	TimeStep *ts = emodel->giveCurrentStep();
    	if ( ts != NULL ) {
    		time = ts->giveTargetTime();
    	}

    	int eiIndex = iCrack.giveNumber();

    	std :: stringstream strNormalJump;
    	strNormalJump << "NormalJumpGnuplotEI" << eiIndex << "Time" << time << ".dat";
    	std :: string nameNormalJump = strNormalJump.str();
    	XFEMDebugTools::WriteArrayToGnuplot(nameNormalJump, arcLengthPositions, normalJumps);

    	std :: stringstream strTangJump;
    	strTangJump << "TangJumpGnuplotEI" << eiIndex << "Time" << time << ".dat";
    	std :: string nameTangJump = strTangJump.str();
    	XFEMDebugTools::WriteArrayToGnuplot(nameTangJump, arcLengthPositions, tangJumps);
    }
}

void GnuplotExportModule::outputBoundaryCondition(PrescribedGradient &iBC, TimeStep *tStep)
{
	FloatArray stress;
	iBC.computeField(stress, EID_MomentumBalance, tStep);
	printf("Mean stress computed in Gnuplot export module: "); stress.printYourself();

	double time = 0.0;

	TimeStep *ts = emodel->giveCurrentStep();
	if ( ts != NULL ) {
		time = ts->giveTargetTime();
	}

	int bcIndex = iBC.giveNumber();

	std :: stringstream strMeanStress;
	strMeanStress << "PrescribedGradientGnuplotMeanStress" << bcIndex << "Time" << time << ".dat";
	std :: string nameMeanStress = strMeanStress.str();
	std::vector<double> componentArray, stressArray;

	for(int i = 1; i <= stress.giveSize(); i++) {
		componentArray.push_back(i);
		stressArray.push_back(stress.at(i));
	}

	XFEMDebugTools::WriteArrayToGnuplot(nameMeanStress, componentArray, stressArray);

}

void GnuplotExportModule::outputBoundaryCondition(PrescribedGradientBCNeumann &iBC, TimeStep *tStep)
{
	int numDofs = iBC.giveInternalDofManager(1)->giveNumberOfDofs();
	FloatArray stress(numDofs);
	for(int dofInd = 1; dofInd <= numDofs; dofInd++) {
		stress.at(dofInd) = iBC.giveInternalDofManager(1)->giveDof(dofInd)->giveUnknown(VM_Total, tStep);
	}
	
		printf("Mean stress computed in Gnuplot export module: "); stress.printYourself();

	double time = 0.0;

	TimeStep *ts = emodel->giveCurrentStep();
	if ( ts != NULL ) {
		time = ts->giveTargetTime();
	}

	int bcIndex = iBC.giveNumber();

	std :: stringstream strMeanStress;
	strMeanStress << "PrescribedGradientGnuplotMeanStress" << bcIndex << "Time" << time << ".dat";
	std :: string nameMeanStress = strMeanStress.str();
	std::vector<double> componentArray, stressArray;

	for(int i = 1; i <= stress.giveSize(); i++) {
		componentArray.push_back(i);
		stressArray.push_back(stress.at(i));
	}

	XFEMDebugTools::WriteArrayToGnuplot(nameMeanStress, componentArray, stressArray);
	
}


} // end namespace oofem
