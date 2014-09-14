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
#include "../sm/EngineeringModels/structengngmodel.h"
#include "outputmanager.h"
#include "dofmanager.h"
#include "boundarycondition.h"
#include "xfem/enrichmentitem.h"
#include "xfem/xfemmanager.h"
#include "../sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "xfem/enrichmentdomain.h"
#include "xfem/XFEMDebugTools.h"
#include "prescribedgradient.h"
#include "prescribedgradientbcneumann.h"
#include "prescribedgradientbcweak.h"
#include "gausspoint.h"
#include "timestep.h"
#include "feinterpol.h"
#include "xfem/enrichmentitems/crack.h"

#include <sstream>

namespace oofem {
REGISTER_ExportModule(GnuplotExportModule)

GnuplotExportModule::GnuplotExportModule(int n, EngngModel *e):
ExportModule(n, e),
mExportReactionForces(false),
mExportBoundaryConditions(false),
mExportMesh(false),
mExportXFEM(false)
{

}

GnuplotExportModule::~GnuplotExportModule() {

}

IRResultType GnuplotExportModule::initializeFrom(InputRecord *ir)
{
    mExportReactionForces = ir->hasField(_IFT_GnuplotExportModule_ReactionForces);
    mExportBoundaryConditions = ir->hasField(_IFT_GnuplotExportModule_BoundaryConditions);
    mExportMesh = ir->hasField(_IFT_GnuplotExportModule_mesh);
    mExportXFEM = ir->hasField(_IFT_GnuplotExportModule_xfem);
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

            PrescribedGradientBCWeak *presGradBCWeak = dynamic_cast<PrescribedGradientBCWeak*>( domain->giveBc(i) );
            if(presGradBCWeak != NULL) {
                outputBoundaryCondition(*presGradBCWeak, tStep);
            }

		}
	}

	if(mExportXFEM) {
        if(domain->hasXfemManager()) {
            XfemManager *xMan = domain->giveXfemManager();

            int numEI = xMan->giveNumberOfEnrichmentItems();

            std::vector< std::vector<FloatArray> > points;

            for(int i = 1; i <= numEI; i++) {
                EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
                ei->callGnuplotExportModule(*this);

                std::vector<FloatArray> eiPoints;
                ei->giveSubPolygon(eiPoints, 0.0, 1.0);
                points.push_back(eiPoints);
            }

            outputXFEMGeometry(points);
        }
	}

	if(mExportMesh) {
	    outputMesh(*domain);
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
		OOFEM_ERROR("failed to cast to StructuralEngngModel.");
	}

    IntArray ielemDofMask;
    FloatArray reactions;
    IntArray dofManMap, dofidMap, eqnMap;

    // test if solution step output is active
    if ( !domain->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    // map contains corresponding dofmanager and dofs numbers corresponding to prescribed equations
    // sorted according to dofmanger number and as a minor crit. according to dof number
    // this is necessary for extractor, since the sorted output is expected
    seMod->buildReactionTable(dofManMap, dofidMap, eqnMap, tStep, 1);

    // compute reaction forces
    seMod->computeReaction(reactions, tStep, 1);

    // Find highest index of prescribed dofs
    int maxIndPresDof = 0;
    for ( int i = 1; i <= dofManMap.giveSize(); i++ ) {
        maxIndPresDof = std::max(maxIndPresDof, dofidMap.at(i));
    }

    int numBC = domain->giveNumberOfBoundaryConditions();

    while ( mReactionForceHistory.size() < size_t(numBC) ) {
    	std::vector<FloatArray> emptyArray;
    	mReactionForceHistory.push_back( emptyArray );
    }

    maxIndPresDof = domain->giveNumberOfSpatialDimensions();

    while ( mDispHist.size() < size_t(numBC) ) {
    	std::vector<double> emptyArray;
    	mDispHist.push_back( emptyArray );
    }

    for(int bcInd = 0; bcInd < numBC; bcInd++) {
    	FloatArray fR(maxIndPresDof), disp(numBC);
    	fR.zero();


        for ( int i = 1; i <= dofManMap.giveSize(); i++ ) {
        	DofManager *dMan = domain->giveDofManager( dofManMap.at(i) );
        	Dof *dof = dMan->giveDofWithID( dofidMap.at(i) );

        	if( dof->giveBcId() == bcInd+1) {
        		fR.at( dofidMap.at(i) ) += reactions.at( eqnMap.at(i) );

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

	std::vector<double> arcLengthPositions, normalJumps, tangJumps;

	const EnrichmentDomain *ed = iCrack.giveEnrichmentDomain();

	for( GaussPoint *gp: czGaussPoints ) {

		StructuralInterfaceMaterialStatus *matStat = dynamic_cast<StructuralInterfaceMaterialStatus*> ( gp->giveMaterialStatus() );
		if(matStat != NULL) {

			// Compute arc length position of the Gauss point
			const FloatArray &coord = *(gp->giveNaturalCoordinates());
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

void GnuplotExportModule::outputXFEMGeometry(const std::vector< std::vector<FloatArray> > &iEnrItemPoints)
{
    double time = 0.0;

    TimeStep *ts = emodel->giveCurrentStep();
    if ( ts != NULL ) {
        time = ts->giveTargetTime();
    }

    std :: stringstream strCracks;
    strCracks << "CracksTime" << time << ".dat";
    std :: string nameCracks = strCracks.str();
    WritePointsToGnuplot(nameCracks, iEnrItemPoints);
}

void GnuplotExportModule::outputBoundaryCondition(PrescribedGradient &iBC, TimeStep *tStep)
{
	FloatArray stress;
	iBC.computeField(stress, tStep);
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
    FloatArray stress;
    iBC.computeField(stress, tStep);

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

void GnuplotExportModule::outputBoundaryCondition(PrescribedGradientBCWeak &iBC, TimeStep *tStep)
{
    FloatArray stress;
    iBC.computeField(stress, tStep);

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



    // Traction node coordinates
    std::vector< std::vector<FloatArray> > nodePointArray;
    size_t numTracEl = iBC.giveNumberOfTractionElements();
    for(size_t i = 0; i < numTracEl; i++) {

        std::vector<FloatArray> points;
        FloatArray xS, xE;
        iBC.giveTractionElCoord(i, xS, xE);
        points.push_back(xS);
        points.push_back(xE);

        nodePointArray.push_back(points);
    }

    std :: stringstream strTractionNodes;
    strTractionNodes << "TractionNodesGnuplotTime" << time << ".dat";
    std :: string nameTractionNodes = strTractionNodes.str();

    WritePointsToGnuplot(nameTractionNodes, nodePointArray);



    // Traction element normal direction
    std::vector< std::vector<FloatArray> > nodeNormalArray;
    for(size_t i = 0; i < numTracEl; i++) {

        std::vector<FloatArray> points;
        FloatArray n,t;
        iBC.giveTractionElNormal(i, n,t);
        points.push_back(n);
        points.push_back(n);

        nodeNormalArray.push_back(points);
    }

    std :: stringstream strTractionNodeNormals;
    strTractionNodeNormals << "TractionNodeNormalsGnuplotTime" << time << ".dat";
    std :: string nameTractionNodeNormals = strTractionNodeNormals.str();

    WritePointsToGnuplot(nameTractionNodeNormals, nodeNormalArray);



    // Traction (x,y)
    std::vector< std::vector<FloatArray> > nodeTractionArray;
    for(size_t i = 0; i < numTracEl; i++) {

        std::vector<FloatArray> tractions;
        FloatArray tS, tE;

        iBC.giveTraction(i, tS, tE, VM_Total, tStep);

        tractions.push_back(tS);
        tractions.push_back(tE);
        nodeTractionArray.push_back(tractions);
    }

    std :: stringstream strTractions;
    strTractions << "TractionsGnuplotTime" << time << ".dat";
    std :: string nameTractions = strTractions.str();

    WritePointsToGnuplot(nameTractions, nodeTractionArray);



    // Arc position along the boundary
    std::vector< std::vector<FloatArray> > arcPosArray;
    for(size_t i = 0; i < numTracEl; i++) {
        std::vector<FloatArray> arcPos;
        double xiS = 0.0, xiE = 0.0;
        iBC.giveTractionElArcPos(i, xiS, xiE);
        arcPos.push_back( FloatArray{xiS} );
        arcPos.push_back( FloatArray{xiE} );

        arcPosArray.push_back(arcPos);
    }

    std :: stringstream strArcPos;
    strArcPos << "ArcPosGnuplotTime" << time << ".dat";
    std :: string nameArcPos = strArcPos.str();

    WritePointsToGnuplot(nameArcPos, arcPosArray);


    // Traction (normal, tangent)
    std::vector< std::vector<FloatArray> > nodeTractionNTArray;
    for(size_t i = 0; i < numTracEl; i++) {

        std::vector<FloatArray> tractions;
        FloatArray tS, tE;

        iBC.giveTraction(i, tS, tE, VM_Total, tStep);
        FloatArray n,t;
        iBC.giveTractionElNormal(i, n, t);


        double tSn = tS.dotProduct(n,2);
        double tSt = tS.dotProduct(t,2);
        tractions.push_back( {tSn ,tSt} );

        double tEn = tE.dotProduct(n,2);
        double tEt = tE.dotProduct(t,2);
        tractions.push_back( {tEn, tEt} );
        nodeTractionNTArray.push_back(tractions);
    }

    std :: stringstream strTractionsNT;
    strTractionsNT << "TractionsNormalTangentGnuplotTime" << time << ".dat";
    std :: string nameTractionsNT = strTractionsNT.str();

    WritePointsToGnuplot(nameTractionsNT, nodeTractionNTArray);



    // Boundary points and displacements
    IntArray boundaries, bNodes;
    iBC.giveBoundaries(boundaries);

    std::vector< std::vector<FloatArray> > bndNodes;

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {

        Element *e = iBC.giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        std::vector<FloatArray> bndSegNodes;

        // Add the start and end nodes of the segment
        DofManager *startNode   = e->giveDofManager( bNodes[0] );
        FloatArray xS    = *(startNode->giveCoordinates());

        Dof *dSu = startNode->giveDofWithID(D_u);
        double dU = dSu->giveUnknown(VM_Total, tStep);
        xS.push_back(dU);

        Dof *dSv = startNode->giveDofWithID(D_v);
        double dV = dSv->giveUnknown(VM_Total, tStep);
        xS.push_back(dV);

        bndSegNodes.push_back(xS);

        DofManager *endNode     = e->giveDofManager( bNodes[1] );
        FloatArray xE    = *(endNode->giveCoordinates());

        Dof *dEu = endNode->giveDofWithID(D_u);
        dU = dEu->giveUnknown(VM_Total, tStep);
        xE.push_back(dU);

        Dof *dEv = endNode->giveDofWithID(D_v);
        dV = dEv->giveUnknown(VM_Total, tStep);
        xE.push_back(dV);

        bndSegNodes.push_back(xE);

        bndNodes.push_back(bndSegNodes);
    }

    std :: stringstream strBndNodes;
    strBndNodes << "BndNodesGnuplotTime" << time << ".dat";
    std :: string nameBndNodes = strBndNodes.str();

    WritePointsToGnuplot(nameBndNodes, bndNodes);

}

void GnuplotExportModule::outputMesh(Domain &iDomain)
{
    std::vector< std::vector<FloatArray> > pointArray;

    if(iDomain.giveNumberOfSpatialDimensions() == 2) {
        // Write all element edges to gnuplot
        int numEl = iDomain.giveNumberOfElements();
        for(int elInd = 1; elInd <= numEl; elInd++) {
            Element *el = iDomain.giveElement(elInd);

            int numEdges = el->giveNumberOfNodes();


            for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
                std::vector<FloatArray> points;

                IntArray bNodes;
                el->giveInterpolation()->boundaryGiveNodes(bNodes, edgeIndex);

                int niLoc = bNodes.at(1);
                const FloatArray &xS = *(el->giveNode(niLoc)->giveCoordinates() );
                points.push_back(xS);

                int njLoc = bNodes.at( bNodes.giveSize() );
                const FloatArray &xE = *(el->giveNode(njLoc)->giveCoordinates() );
                points.push_back(xE);

                pointArray.push_back(points);
            }

        }


        double time = 0.0;

        TimeStep *ts = emodel->giveCurrentStep();
        if ( ts != NULL ) {
            time = ts->giveTargetTime();
        }

        std :: stringstream strMesh;
        strMesh << "MeshGnuplotTime" << time << ".dat";
        std :: string nameMesh = strMesh.str();

        WritePointsToGnuplot(nameMesh, pointArray);
    }
}

void GnuplotExportModule :: WritePointsToGnuplot(const std :: string &iName, const std :: vector< std::vector<FloatArray> > &iPoints)
{
    std :: ofstream file;
    file.open( iName.data() );

    // Set some output options
    file << std :: scientific;

    file << "# x y\n";

    for(auto posVec: iPoints) {
        for(auto pos: posVec) {

            for(int i = 0; i < pos.giveSize(); i++) {
                file << pos[i] << " ";
            }
            file << "\n";

//            file << pos[0] << " " << pos[1] << "\n";
        }
        file << "\n";
    }

    file.close();

}


} // end namespace oofem
