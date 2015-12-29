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

#include "fe2structuralmaterial.h"

#include "structuralmaterial.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "oofemtxtdatareader.h"
#include "domain.h"
#include "contextioerr.h"
#include "gausspoint.h"

#include "engngm.h"
#include "util.h"

#include <sstream>

namespace oofem {
REGISTER_Material(FE2StructuralMaterial);
int FE2StructuralMaterial :: n = 1;

FE2StructuralMaterialStatus :: FE2StructuralMaterialStatus(int n, Domain *d, GaussPoint *gp, const std :: string &inputfile) :
    StructuralMaterialStatus(n, d, gp)
{

    if ( !this->createRVE(n, gp, inputfile) ) {
        OOFEM_ERROR("Couldn't create RVE");
    }
}

FE2StructuralMaterialStatus :: ~FE2StructuralMaterialStatus()
{

}

bool FE2StructuralMaterialStatus :: createRVE(int n, GaussPoint *gp, const std :: string &inputfile)
{
    OOFEMTXTDataReader dr( inputfile.c_str() );
    EngngModel *em = InstanciateProblem(& dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->setProblemScale(microScale);
    em->checkProblemConsistency();
    em->initMetaStepAttributes( em->giveMetaStep(1) );
    em->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    em->init();

    this->rve.reset( em );

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }

    this->rve->letOutputBaseFileNameBe( name.str() );

    this->bc = dynamic_cast< PrescribedGradientBCWeak * >( this->rve->giveDomain(1)->giveBc(1) );
    if ( !this->bc ) {
        OOFEM_ERROR("RVE doesn't have necessary boundary condition; should have MixedGradientPressure as first b.c. (in first domain)");
    }

    return true;
}

void FE2StructuralMaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
}

void FE2StructuralMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
}

void FE2StructuralMaterialStatus :: initTempStatus()
{
	StructuralMaterialStatus :: initTempStatus();
}

void FE2StructuralMaterialStatus :: updateYourself(TimeStep *tStep)
{
	StructuralMaterialStatus :: updateYourself(tStep);

    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);
}

contextIOResultType FE2StructuralMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return this->rve->saveContext(&stream, mode, obj);
}

contextIOResultType FE2StructuralMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
//    this->markOldTangents();

    return this->rve->restoreContext(&stream, mode, obj);
}

FE2StructuralMaterial::FE2StructuralMaterial(int n, Domain * d): StructuralMaterial(n,d) {


}

FE2StructuralMaterial::~FE2StructuralMaterial() {

}

IRResultType FE2StructuralMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    IR_GIVE_FIELD(ir, this->inputfile, _IFT_FE2StructuralMaterial_fileName);
    return StructuralMaterial :: initializeFrom(ir);
}

void FE2StructuralMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
    input.setField(this->inputfile, _IFT_FE2StructuralMaterial_fileName);
}

MaterialStatus *FE2StructuralMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new FE2StructuralMaterialStatus(n++, this->giveDomain(), gp, this->inputfile);
}

void FE2StructuralMaterial ::giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
{
//	printf("Entering FE2StructuralMaterial ::giveRealStressVector_3d.\n");

    FE2StructuralMaterialStatus *ms = static_cast< FE2StructuralMaterialStatus * >( this->giveStatus(gp) );

    ms->setTimeStep(tStep);

    PrescribedGradientBCWeak *bc = ms->giveBC();

    // Set strain
//    printf("reducedE: "); reducedE.printYourself();

    bc->setPrescribedGradientVoigt(reducedE);

    // Solve RVE problem
    ms->giveRVE()->solveYourselfAt(ms->giveRVE()->giveCurrentStep());

    // Compute homogenized stress
    bc->computeField(answer, tStep);
//    printf("Homogenized stress: "); answer.printYourself();

    ms->letTempStressVectorBe(answer);

//	OOFEM_ERROR("Not supported.")
}

void FE2StructuralMaterial ::give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
//	printf("Entering FE2StructuralMaterial ::give3dMaterialStiffnessMatrix.\n");


	// Numerical tangent
	FE2StructuralMaterialStatus *status = static_cast<FE2StructuralMaterialStatus*>( this->giveStatus( gp ) );
    double h = 1.0e-9;

    const FloatArray &epsRed = status->giveTempStrainVector();
    FloatArray eps;
    StructuralMaterial::giveFullSymVectorForm(eps, epsRed, gp->giveMaterialMode() );


    int dim = eps.giveSize();
    answer.resize(dim, dim);
    answer.zero();

    FloatArray sig, sigPert, epsPert;

    for(int i = 1; i <= dim; i++) {
    	// Add a small perturbation to the strain
    	epsPert = eps;
    	epsPert.at(i) += h;

    	giveRealStressVector_3d(sigPert, gp, epsPert, tStep);
    	answer.setColumn(sigPert, i);
    }

	giveRealStressVector_3d(sig, gp, eps, tStep);

    for(int i = 1; i <= dim; i++) {
        for(int j = 1; j <= dim; j++) {
        	answer.at(j,i) -= sig.at(j);
        	answer.at(j,i) /= h;
        }
    }

}


} /* namespace oofem */
