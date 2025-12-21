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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#include "sm/EngineeringModels/structengngmodel.h"
#include "outputmanager.h"
#include "dofmanager.h"
#include "boundarycondition.h"
#include "xfem/enrichmentitem.h"
#include "xfem/xfemmanager.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "xfem/XFEMDebugTools.h"
#include "prescribedgradient.h"
#include "prescribedgradientbcneumann.h"
#include "prescribedgradientbcweak.h"
#include "gausspoint.h"
#include "timestep.h"
#include "feinterpol.h"
#include "xfem/enrichmentitems/crack.h"
#include "dofmanager.h"
#include "xfem/matforceevaluator.h"
#include "function.h"
#include "sm/Elements/Interfaces/structuralinterfaceelement.h"
#include "sm/Materials/structuralfe2material.h"

#include <sstream>

namespace oofem {
REGISTER_ExportModule(GnuplotExportModule)

GnuplotExportModule::GnuplotExportModule(int n, EngngModel *e):
    ExportModule(n, e),
    mExportReactionForces(false),
    mExportBoundaryConditions(false),
    mExportBoundaryConditionsExtra(false),
    mExportMesh(false),
    mExportXFEM(false),
    mExportCrackLength(false),
	mExportInterfaceEl(false),
    mMonitorNodeIndex(-1),
    mpMatForceEvaluator( new  MaterialForceEvaluator() )
{}

GnuplotExportModule::~GnuplotExportModule()
{}

void GnuplotExportModule::initializeFrom(InputRecord &ir)
{
    ExportModule::initializeFrom(ir);

    mExportReactionForces = ir.hasField(_IFT_GnuplotExportModule_ReactionForces);
    mExportBoundaryConditions = ir.hasField(_IFT_GnuplotExportModule_BoundaryConditions);
    mExportBoundaryConditionsExtra = ir.hasField(_IFT_GnuplotExportModule_BoundaryConditionsExtra);
    mExportMesh = ir.hasField(_IFT_GnuplotExportModule_mesh);
    mExportXFEM = ir.hasField(_IFT_GnuplotExportModule_xfem);
    mExportCrackLength = ir.hasField(_IFT_GnuplotExportModule_cracklength);
    mExportInterfaceEl = ir.hasField(_IFT_GnuplotExportModule_interface_el);

    ir.giveOptionalField(mMonitorNodeIndex, _IFT_GnuplotExportModule_monitornode);

    ir.giveOptionalField(mMatForceRadii, _IFT_GnuplotExportModule_materialforceradii);
//    printf("mMatForceRadii: "); mMatForceRadii.printYourself();
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

    mTimeHist.push_back( tStep->giveTargetTime() );

    if(mExportXFEM) {
        if(domain->hasXfemManager()) {
            XfemManager *xMan = domain->giveXfemManager();

            int numEI = xMan->giveNumberOfEnrichmentItems();

            std::vector< std::vector<FloatArray> > points;

            for(int i = 1; i <= numEI; i++) {
                EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
                ei->callGnuplotExportModule(*this, tStep);

                GeometryBasedEI *geoEI = dynamic_cast<GeometryBasedEI*>(ei);
                if(geoEI != NULL) {
                    std::vector<FloatArray> eiPoints;
                    geoEI->giveSubPolygon(eiPoints, 0.0, 1.0);
                    points.push_back(eiPoints);
                }
            }

            outputXFEMGeometry(points);
        }
    }

    if(mExportMesh) {
        outputMesh(*domain);
    }

    if(mMonitorNodeIndex != -1) {
        DofManager *dMan = domain->giveDofManager(mMonitorNodeIndex);
        outputNodeDisp(*dMan, tStep);
    }

    if(mExportInterfaceEl) {
    	outputInterfaceEl(*domain, tStep);
    }
}

void GnuplotExportModule::initialize()
{
    ExportModule :: initialize();
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

    FloatArray reactions;
    IntArray dofManMap, dofidMap, eqnMap;

    // test if solution step output is active
    if ( !testTimeStepOutput(tStep) ) {
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

            if ( dof->giveBcId() == bcInd+1 ) {
                fR.at( dofidMap.at(i) ) += reactions.at( eqnMap.at(i) );

                // Slightly dirty
                BoundaryCondition *bc = dynamic_cast<BoundaryCondition*> (domain->giveBc(bcInd+1));
                if ( bc != NULL ) {
                    disp.at(bcInd+1) = std::max( disp.at(bcInd+1), bc->give(dof, VM_Total, tStep->giveTargetTime()) );
                }
                ///@todo This function should be using the primaryfield instead of asking BCs directly. / Mikael
            }
        }

        mDispHist[bcInd].push_back(disp.at(bcInd+1));
        mReactionForceHistory[bcInd].push_back(fR);



        // X
        FILE * pFileX;
        char fileNameX[100];
        sprintf(fileNameX, "ReactionForceGnuplotBC%dX.dat", bcInd+1);
        pFileX = fopen ( fileNameX , "wb" );

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

        fprintf(pFileY, "#u Fx\n");
        for ( size_t j = 0; j < mDispHist[bcInd].size(); j++ ) {
            if( mReactionForceHistory[bcInd][j].giveSize() >= 2 ) {
                fprintf(pFileY, "%e %e\n", mDispHist[bcInd][j], mReactionForceHistory[bcInd][j].at(2) );
            }
        }

        fclose(pFileY);

    }
}

void GnuplotExportModule::outputXFEM(EnrichmentItem &iEI, TimeStep *tStep)
{

}

void GnuplotExportModule::outputXFEM(Crack &iCrack, TimeStep *tStep)
{
    const std::vector<GaussPoint*> &czGaussPoints = iCrack.giveCohesiveZoneGaussPoints();

    std::vector<double> arcLengthPositions, normalJumps, tangJumps, normalTractions, tangTractions;

    std::vector<double> arcLengthPositionsB = iCrack.giveCohesiveZoneArcPositions();

    const BasicGeometry *bg = iCrack.giveGeometry();

//    for( auto *gp: czGaussPoints ) {
//
//        StructuralInterfaceMaterialStatus *matStat = dynamic_cast<StructuralInterfaceMaterialStatus*> ( gp->giveMaterialStatus() );
//        if(matStat != NULL) {
//
//            // Compute arc length position of the Gauss point
//            const FloatArray &coord = (gp->giveGlobalCoordinates());
//            printf("coord: "); coord.printYourself();
//            double tangDist = 0.0, arcPos = 0.0;
//            bg->computeTangentialSignDist(tangDist, coord, arcPos);
//            arcLengthPositions.push_back(arcPos);
//
////            // Compute displacement jump in normal and tangential direction
////            // Local numbering: (tang_z, tang, normal)
////            const FloatArray &jumpLoc = matStat->giveJump();
////
////            double normalJump = jumpLoc.at(3);
////            normalJumps.push_back(normalJump);
////
////
////            tangJumps.push_back( jumpLoc.at(2) );
////
////
////            const FloatArray &trac = matStat->giveFirstPKTraction();
////            normalTractions.push_back(trac.at(3));
//        }
//    }

#if 1
    size_t num_cz_gp = czGaussPoints.size();

//    for( auto *gp: czGaussPoints ) {
    for(size_t gp_ind = 0; gp_ind < num_cz_gp; gp_ind++) {
//    	printf("gp_ind: %lu\n", gp_ind );
    	GaussPoint *gp = czGaussPoints[gp_ind];

    	if( gp != NULL ) {

			StructuralInterfaceMaterialStatus *matStat = dynamic_cast<StructuralInterfaceMaterialStatus*> ( gp->giveMaterialStatus() );
			if(matStat != NULL) {

				// Compute arc length position of the Gauss point
				const FloatArray &coord = (gp->giveGlobalCoordinates());
				double tangDist = 0.0, arcPos = 0.0;
				bg->computeTangentialSignDist(tangDist, coord, arcPos);
//				printf("arcPos: %e\n", arcPos );
				arcLengthPositions.push_back(arcPos);

				// Compute displacement jump in normal and tangential direction
				// Local numbering: (tang_z, tang, normal)
				const FloatArray &jumpLoc = matStat->giveJump();

				double normalJump = jumpLoc.at(3);
				normalJumps.push_back(normalJump);


				tangJumps.push_back( jumpLoc.at(2) );


				const FloatArray &trac = matStat->giveFirstPKTraction();
				normalTractions.push_back(trac.at(3));

				tangTractions.push_back(trac.at(2));
			}
			else {
			    StructuralFE2MaterialStatus *fe2ms = dynamic_cast<StructuralFE2MaterialStatus*>(gp->giveMaterialStatus());

			    if(fe2ms != NULL) {
//			    	printf("Casted to StructuralFE2MaterialStatus.\n");

					const FloatArray &coord = (gp->giveGlobalCoordinates());
					double tangDist = 0.0, arcPos = 0.0;
					bg->computeTangentialSignDist(tangDist, coord, arcPos);
	//				printf("arcPos: %e\n", arcPos );
					arcLengthPositions.push_back(arcPos);

			    	const FloatArray &n = fe2ms->giveNormal();
//			    	printf("n: "); n.printYourself();

			    	const FloatArray &sig_v = fe2ms->giveStressVector();
//			    	printf("sig_v: "); sig_v.printYourself();

			    	FloatMatrix sig_m(2,2);
			    	sig_m(0,0) = sig_v(0);
			    	sig_m(1,1) = sig_v(1);
			    	sig_m(0,1) = sig_m(1,0) = sig_v( sig_v.giveSize()-1 );

			    	FloatArray trac(2);
			    	trac(0) = sig_m(0,0)*n(0) + sig_m(0,1)*n(1);
			    	trac(1) = sig_m(1,0)*n(0) + sig_m(1,1)*n(1);
			    	double trac_n = trac(0)*n(0) + trac(1)*n(1);
					normalTractions.push_back(trac_n);

					const FloatArray t = Vec2(-n(1), n(0));
			    	double trac_t = trac(0)*t(0) + trac(1)*t(1);
					tangTractions.push_back(trac_t);


			    }

			}
    	}
    }
#endif



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

        std :: stringstream strNormalTrac;
        strNormalTrac << "NormalTracGnuplotEI" << eiIndex << "Time" << time << ".dat";
        std :: string nameNormalTrac = strNormalTrac.str();
        XFEMDebugTools::WriteArrayToGnuplot(nameNormalTrac, arcLengthPositions, normalTractions);

        std :: stringstream strTangTrac;
        strTangTrac << "TangTracGnuplotEI" << eiIndex << "Time" << time << ".dat";
        std :: string nameTangTrac = strTangTrac.str();
        XFEMDebugTools::WriteArrayToGnuplot(nameTangTrac, arcLengthPositions, tangTractions);


        std::vector<FloatArray> matForcesStart, matForcesEnd;
        std::vector<double> radii;

        // Material forces
        for(double matForceRadius : mMatForceRadii) {

            radii.push_back(matForceRadius);

            EnrichmentFront *efStart = iCrack.giveEnrichmentFrontStart();
            const TipInfo &tipInfoStart = efStart->giveTipInfo();

            FloatArray matForceStart;
            mpMatForceEvaluator->computeMaterialForce(matForceStart, *domain, tipInfoStart, tStep, matForceRadius);

            if(matForceStart.giveSize() > 0) {
                matForcesStart.push_back(matForceStart);
            }
            else {
                matForcesStart.push_back(Vec2(0.0,0.0));
            }


            EnrichmentFront *efEnd = iCrack.giveEnrichmentFrontEnd();
            const TipInfo &tipInfoEnd = efEnd->giveTipInfo();

            FloatArray matForceEnd;
            mpMatForceEvaluator->computeMaterialForce(matForceEnd, *domain, tipInfoEnd, tStep, matForceRadius);

            if(matForceEnd.giveSize() > 0) {
                matForcesEnd.push_back(matForceEnd);
            }
            else {
                matForcesEnd.push_back(Vec2(0.0,0.0));
            }

        }

        std::vector< std::vector<FloatArray> > matForcesStartArray, matForcesEndArray;
        matForcesStartArray.push_back(matForcesStart);
        matForcesEndArray.push_back(matForcesEnd);


        std :: stringstream strRadii;
        strRadii << "MatForceRadiiGnuplotTime" << time << "Crack" << iCrack.giveNumber() << ".dat";
        XFEMDebugTools::WriteArrayToGnuplot(strRadii.str(), radii, radii);


        std :: stringstream strMatForcesStart;
        strMatForcesStart << "MatForcesStartGnuplotTime" << time << "Crack" << iCrack.giveNumber() << ".dat";
        WritePointsToGnuplot(strMatForcesStart.str(), matForcesStartArray);


        std :: stringstream strMatForcesEnd;
        strMatForcesEnd << "MatForcesEndGnuplotTime" << time << "Crack" << iCrack.giveNumber() << ".dat";
        WritePointsToGnuplot(strMatForcesEnd.str(), matForcesEndArray);

        double crackLength = iCrack.computeLength();
    //    printf("crackLength: %e\n", crackLength );
        mCrackLengthHist[eiIndex].push_back(crackLength);
        std :: stringstream strCrackLength;
        strCrackLength << "CrackLengthGnuplotEI" << eiIndex << "Time" << time << ".dat";
        XFEMDebugTools::WriteArrayToGnuplot(strCrackLength.str(), mTimeHist, mCrackLengthHist[eiIndex]);
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

void GnuplotExportModule :: outputBoundaryCondition(PrescribedGradient &iBC, TimeStep *tStep)
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

    // Homogenized stress
    std :: stringstream strMeanStress;
    strMeanStress << "PrescribedGradientGnuplotMeanStress" << bcIndex << "Time" << time << ".dat";
    std :: string nameMeanStress = strMeanStress.str();
    std::vector<double> componentArray, stressArray;

    for(int i = 1; i <= stress.giveSize(); i++) {
        componentArray.push_back(i);
        stressArray.push_back(stress.at(i));
    }

    XFEMDebugTools::WriteArrayToGnuplot(nameMeanStress, componentArray, stressArray);

    FloatArray grad;
    iBC.giveGradientVoigt(grad);
    outputGradient(iBC.giveNumber(), *iBC.giveDomain(), grad, tStep);
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

    // Homogenized strain
    
    FloatArray grad;
    iBC.giveGradientVoigt(grad);
    outputGradient(iBC.giveNumber(), *iBC.giveDomain(), grad, tStep);
}

void GnuplotExportModule::outputBoundaryCondition(PrescribedGradientBCWeak &iBC, TimeStep *tStep)
{
    FloatArray stress;
    iBC.computeField(stress, tStep);

    printf("Mean stress computed in Gnuplot export module: "); stress.printYourself();
    printf("sigXX: %.12e\n", stress(0) );

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


    // Homogenized strain
    FloatArray grad;
    iBC.giveGradientVoigt(grad);
    outputGradient(iBC.giveNumber(), *iBC.giveDomain(), grad, tStep);

#if 0
    FloatArray grad;
    iBC.giveGradientVoigt(grad);
    double timeFactor = iBC.giveTimeFunction()->evaluate(ts, VM_Total);
    printf("timeFactor: %e\n", timeFactor );
    grad.times(timeFactor);
    printf("Mean grad computed in Gnuplot export module: "); grad.printYourself();

    std :: stringstream strMeanGrad;
    strMeanGrad << "PrescribedGradientGnuplotMeanGrad" << bcIndex << "Time" << time << ".dat";
    std :: string nameMeanGrad = strMeanGrad.str();
    std::vector<double> componentArrayGrad, gradArray;

    for(int i = 1; i <= grad.giveSize(); i++) {
        componentArrayGrad.push_back(i);
        gradArray.push_back(grad.at(i));
    }

    XFEMDebugTools::WriteArrayToGnuplot(nameMeanGrad, componentArrayGrad, gradArray);
#endif

    if(mExportBoundaryConditionsExtra) {

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
            arcPos.push_back( Vec1(xiS) );
            arcPos.push_back( Vec1(xiE) );

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
            tractions.push_back(Vec2(tSn ,tSt));

            double tEn = tE.dotProduct(n,2);
            double tEt = tE.dotProduct(t,2);
            tractions.push_back(Vec2(tEn, tEt));
            nodeTractionNTArray.push_back(tractions);
        }

        std :: stringstream strTractionsNT;
        strTractionsNT << "TractionsNormalTangentGnuplotTime" << time << ".dat";
        std :: string nameTractionsNT = strTractionsNT.str();

        WritePointsToGnuplot(nameTractionsNT, nodeTractionNTArray);



        // Boundary points and displacements
        IntArray boundaries;
        iBC.giveBoundaries(boundaries);

        std::vector< std::vector<FloatArray> > bndNodes;

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {

            Element *e = iBC.giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary, e->giveGeometryType());

            std::vector<FloatArray> bndSegNodes;

            // Add the start and end nodes of the segment
            DofManager *startNode = e->giveDofManager( bNodes[0] );

            Dof *dSu = startNode->giveDofWithID(D_u);
            double dU = dSu->giveUnknown(VM_Total, tStep);

            Dof *dSv = startNode->giveDofWithID(D_v);
            double dV = dSv->giveUnknown(VM_Total, tStep);

            bndSegNodes.push_back(FloatArray::fromConcatenated({startNode->giveCoordinates(),Vec2(dU,dV)}));

            DofManager *endNode = e->giveDofManager( bNodes[1] );

            Dof *dEu = endNode->giveDofWithID(D_u);
            dU = dEu->giveUnknown(VM_Total, tStep);

            Dof *dEv = endNode->giveDofWithID(D_v);
            dV = dEv->giveUnknown(VM_Total, tStep);

            bndSegNodes.push_back(FloatArray::fromConcatenated({endNode->giveCoordinates(),Vec2(dU,dV)}));

            bndNodes.push_back(bndSegNodes);
        }

        std :: stringstream strBndNodes;
        strBndNodes << "BndNodesGnuplotTime" << time << ".dat";
        std :: string nameBndNodes = strBndNodes.str();

        WritePointsToGnuplot(nameBndNodes, bndNodes);

    }
}

void GnuplotExportModule::outputGradient(int bc, Domain &d, FloatArray &grad, TimeStep *tStep)
{
    // Homogenized strain
    double timeFactor = d.giveBc(bc)->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
    printf("timeFactor: %e\n", timeFactor );
    grad.times(timeFactor);
    printf("Mean grad computed in Gnuplot export module: "); grad.printYourself();

    double time = tStep->giveTargetTime();


    std :: stringstream strMeanGrad;
    strMeanGrad << "PrescribedGradientGnuplotMeanGrad" << bc << "Time" << time << ".dat";
    std :: string nameMeanGrad = strMeanGrad.str();
    std::vector<double> componentArrayGrad, gradArray;

    for(int i = 1; i <= grad.giveSize(); i++) {
        componentArrayGrad.push_back(i);
        gradArray.push_back(grad.at(i));
    }

    XFEMDebugTools::WriteArrayToGnuplot(nameMeanGrad, componentArrayGrad, gradArray);
}

void GnuplotExportModule::outputMesh(Domain &iDomain)
{
    std::vector< std::vector<FloatArray> > pointArray;

    if(iDomain.giveNumberOfSpatialDimensions() == 2) {
        // Write all element edges to gnuplot
        for ( auto &el : iDomain.giveElements() ) {
            int numEdges = el->giveNumberOfNodes();


            for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
                std::vector<FloatArray> points;

                const auto &bNodes = el->giveInterpolation()->boundaryGiveNodes(edgeIndex, el->giveGeometryType());

                int niLoc = bNodes.at(1);
                const auto &xS = el->giveNode(niLoc)->giveCoordinates();
                points.push_back(xS);

                int njLoc = bNodes.at( bNodes.giveSize() );
                const auto &xE = el->giveNode(njLoc)->giveCoordinates();
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

void GnuplotExportModule :: outputNodeDisp(DofManager &iDMan, TimeStep *tStep)
{
    // Append node solution to history
    FloatArray nodeSol;
    iDMan.giveCompleteUnknownVector(nodeSol, VM_Total, tStep);

    mMonitorNodeDispHist.push_back(nodeSol);


    // Write to file
    std::vector< std::vector<FloatArray> > nodeDispHist;
    nodeDispHist.push_back(mMonitorNodeDispHist);

    std :: string name = "MonitorNodeSolGnuplot.dat";
    WritePointsToGnuplot(name, nodeDispHist);
}

void GnuplotExportModule :: outputInterfaceEl(Domain &d, TimeStep *tStep) {
//	printf("Exporting interface el.\n");

	std::vector<FloatArray> points, tractions, tractionsProj;

	int numEl = d.giveNumberOfElements();
	for(int i = 1; i <= numEl; i++) {
		Element *el = d.giveElement(i);

		StructuralInterfaceElement *intEl = dynamic_cast<StructuralInterfaceElement*>(el);

		if(intEl) {
			printf("Found StructuralInterfaceElement.\n");

			IntegrationRule *ir = intEl->giveDefaultIntegrationRulePtr();

			int numGP = ir->giveNumberOfIntegrationPoints();
//			printf("numGP: %d\n", numGP );

			for(int gpInd = 0; gpInd < numGP; gpInd++ ) {
				GaussPoint *gp = ir->getIntegrationPoint(gpInd);
//				gp->giveGlobalCoordinates().printYourself();

				points.push_back( (gp->giveGlobalCoordinates()) );

			    MaterialStatus *ms = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
				StructuralInterfaceMaterialStatus *ims = static_cast<StructuralInterfaceMaterialStatus*>(ms);

				FloatArray traction = ims->giveTraction();
//				printf("traction: "); traction.printYourself();
				tractions.push_back(traction);


				FloatArray tractionProj = ims->giveProjectedTraction();
				tractionsProj.push_back(tractionProj);

			}

		}
	}

    // Export x vs normal traction

    double time = 0.0;

    TimeStep *ts = emodel->giveCurrentStep();
    if ( ts != NULL ) {
        time = ts->giveTargetTime();
    }

    FILE * pFileX;
    std :: stringstream strFileNameX;
    strFileNameX << "NormalTractionVsXTime" << time << ".dat";
    std :: string nameStringX = strFileNameX.str();

    pFileX = fopen ( nameStringX.c_str() , "wb" );

    fprintf(pFileX, "#x tn\n");
    for ( size_t j = 0; j < points.size(); j++ ) {
    	fprintf(pFileX, "%e %e\n", points[j][0], tractions[j].at(3) );
    }

    fclose(pFileX);


//    FILE * pFileXProj;
//    std :: stringstream strFileNameXProj;
//    strFileNameXProj << "NormalTractionProjVsXTime" << time << ".dat";
//    std :: string nameStringXProj = strFileNameXProj.str();
//
//    pFileXProj = fopen ( nameStringXProj.c_str() , "wb" );
//
//    fprintf(pFileXProj, "#x tn\n");
//    for ( size_t j = 0; j < points.size(); j++ ) {
////    	printf("tractionsProj[j]: "); tractionsProj[j].printYourself();
//    	fprintf(pFileXProj, "%e %e\n", points[j][0], tractionsProj[j].at(3) );
//    }
//
//    fclose(pFileXProj);


    // Export x vs shear traction

    FILE * pFileXshear;
    std :: stringstream strFileNameXshear;
    strFileNameXshear << "ShearTractionVsXTime" << time << ".dat";
    std :: string nameStringXshear = strFileNameXshear.str();

    pFileXshear = fopen ( nameStringXshear.c_str() , "wb" );

    fprintf(pFileXshear, "#x tn\n");
    for ( size_t j = 0; j < points.size(); j++ ) {
    	fprintf(pFileXshear, "%e %e\n", points[j][0], tractions[j].at(1) );
    }

    fclose(pFileXshear);




    // Export y vs normal traction
    FILE * pFileY;
    std :: stringstream strFileNameY;
    strFileNameY << "NormalTractionVsYTime" << time << ".dat";
    std :: string nameStringY = strFileNameY.str();

    pFileY = fopen ( nameStringY.c_str() , "wb" );

    fprintf(pFileY, "#y tn\n");
    for ( size_t j = 0; j < points.size(); j++ ) {
    	fprintf(pFileY, "%e %e\n", points[j][1], tractions[j].at(3) );
    }

    fclose(pFileY);




    // Export y vs shear traction
    FILE * pFileYshear;
    std :: stringstream strFileNameYshear;
    strFileNameYshear << "ShearTractionVsYTime" << time << ".dat";
    std :: string nameStringYshear = strFileNameYshear.str();

    pFileYshear = fopen ( nameStringYshear.c_str() , "wb" );

    fprintf(pFileYshear, "#y tn\n");
    for ( size_t j = 0; j < points.size(); j++ ) {
    	fprintf(pFileYshear, "%e %e\n", points[j][1], tractions[j].at(1) );
    }

    fclose(pFileYshear);

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
