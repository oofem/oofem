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

#include "bondceb.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_Material(BondCEBMaterial);

BondCEBMaterial :: BondCEBMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    //
    // constructor
    //
{
    tauf = 0.;
    alpha = 0.4;
}


BondCEBMaterial :: ~BondCEBMaterial()
//
// destructor
//
{ }

int
BondCEBMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _2dInterface || mode == _3dInterface;
}


void
BondCEBMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                            MatResponseMode mode,
                                                            GaussPoint *gp,
                                                            TimeStep *tStep)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
    OOFEM_ERROR("not implemented");
}


void
BondCEBMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &totalStrain,
                                                   TimeStep *tStep)
//
// returns the stress (traction) at the end of the step
// corresponding to the given strain at the end of the step
//
{
    MaterialMode mMode = gp->giveMaterialMode();    
    int ntc = 2; // number of traction vector components (2 or 3)
    if ( mMode == _3dInterface ) {
      ntc = 3;
    }
    int i;
    BondCEBMaterialStatus *status = static_cast< BondCEBMaterialStatus * >( this->giveStatus(gp) );

    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    
    // normal traction evaluated elastically
    answer.resize(ntc);
    answer.at(1) = kn * totalStrain.at(1);

    // trial values of shear tractions evaluated elastically
    double s = 0., dKappa = 0.;
    for (i=2; i<=ntc; i++) {
      double depsi = totalStrain.at(i) - status->giveStrainVector().at(i);
      answer.at(i) = status->giveStressVector().at(i) + ks * depsi;
      s += answer.at(i)*answer.at(i);
      dKappa += depsi*depsi;
     }

    // norm of trial shear traction
    s = sqrt(s);

    // cumulative slip increment
    dKappa = sqrt(dKappa);

    // cumulative slip at the end of the step
    double tempKappa = status->giveKappa() + dKappa;

    // maximum allowed norm of shear traction
    double smax = evaluateBondStress(tempKappa);

    // reduce shear tractions, if needed
    if (s>smax){
      for (i=2; i<=ntc; i++) {
	answer.at(i) *= (smax/s);
      }
    }
 
    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempKappa(tempKappa);
}

double
BondCEBMaterial :: evaluateBondStress(const double kappa)
{
  if (kappa<=0.)
    return 0.;
  if (kappa<=s1)
    return taumax*pow(kappa/s1,alpha);
  if (kappa<=s2)
    return taumax;
  if (kappa<=s3)
    return taumax - (taumax-tauf) * (kappa-s2) / (s3-s2);
  return tauf;
}

void
BondCEBMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                                  MatResponseMode rMode,
                                                  GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _2dInterface:
        answer.resize(2, 2);
	answer.zero();
        answer.at(1, 1) = kn;
	if ( rMode == ElasticStiffness ) {
	  answer.at(2, 2) = ks;
	} else {
	  answer.at(2, 2) = ks;
	}
        break;
    case _3dInterface:
        answer.resize(3, 3);
	answer.zero();
	answer.at(1, 1) = kn;
	if ( rMode == ElasticStiffness ) {
	  answer.at(2, 2) = answer.at(3, 3) = ks;
	} else {
	  answer.at(2, 2) = answer.at(3, 3) = ks;
	}
        break;
    default:
        StructuralMaterial :: giveStiffnessMatrix(answer, rMode, gp, tStep);
    }
}

int
BondCEBMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    BondCEBMaterialStatus *status = static_cast< BondCEBMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

IRResultType
BondCEBMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // mandatory parameters
    IR_GIVE_FIELD(ir, kn, _IFT_BondCEBMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_BondCEBMaterial_ks);
    IR_GIVE_FIELD(ir, s1, _IFT_BondCEBMaterial_s1);
    IR_GIVE_FIELD(ir, s2, _IFT_BondCEBMaterial_s2);
    IR_GIVE_FIELD(ir, s3, _IFT_BondCEBMaterial_s3);
    IR_GIVE_FIELD(ir, taumax, _IFT_BondCEBMaterial_taumax);

    // optional parameters
    IR_GIVE_OPTIONAL_FIELD(ir, tauf, _IFT_BondCEBMaterial_tauf);
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_BondCEBMaterial_al);

    // dependent parameter
    s0 = pow(pow(s1,-alpha)*taumax/ks,1./(1.-alpha));
    if (s0>s1) {
      s0 = s1;
      ks = taumax/s1;
      OOFEM_WARNING("Parameter ks adjusted");
    }

    return StructuralMaterial :: initializeFrom(ir);
}


void
BondCEBMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_BondCEBMaterial_kn);
    input.setField(this->ks, _IFT_BondCEBMaterial_ks);
    input.setField(this->s1, _IFT_BondCEBMaterial_s1);
    input.setField(this->s2, _IFT_BondCEBMaterial_s2);
    input.setField(this->s3, _IFT_BondCEBMaterial_s3);
    input.setField(this->taumax, _IFT_BondCEBMaterial_taumax);
}



BondCEBMaterialStatus :: BondCEBMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
}


BondCEBMaterialStatus :: ~BondCEBMaterialStatus()
{ }


void
BondCEBMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->kappa > 0.0 ) {
        fprintf(file, "kappa %g ", this->kappa);
    }

    fprintf(file, "}\n");
}


void
BondCEBMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
}

void
BondCEBMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
}


contextIOResultType
BondCEBMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
BondCEBMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}
} // end namespace oofem
