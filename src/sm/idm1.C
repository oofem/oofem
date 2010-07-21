/* $Header: /home/cvs/bp/oofem/sm/src/idm1.C,v 1.8.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

// file: idm1.C


#include "idm1.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mmaclosestiptransfer.h"
#include "datastream.h"
#include "contextioerr.h"

#ifndef __MAKEDEPEND
#include <math.h>
#endif

namespace oofem {

#ifdef IDM_USE_MMAClosestIPTransfer
MMAClosestIPTransfer IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMAContainingElementProjection
MMAContainingElementProjection IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMAShapeFunctProjection
MMAShapeFunctProjection IsotropicDamageMaterial1 :: mapper;
#endif

#ifdef IDM_USE_MMALeastSquareProjection
MMALeastSquareProjection IsotropicDamageMaterial1 :: mapper;
#endif

IsotropicDamageMaterial1 :: IsotropicDamageMaterial1(int n, Domain *d) : IsotropicDamageMaterial(n, d),
									 RandomMaterialExtensionInterface()
    //
    // constructor
    //
{
    // deleted by paren, where linearElasticMaterial instance declared
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    equivStrainType = EST_Unknown;
    softType = ST_Unknown;
    k = 0.;
    md = 1.;
}


IsotropicDamageMaterial1 :: ~IsotropicDamageMaterial1()
//
// destructor
//
{ }

IRResultType
IsotropicDamageMaterial1 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int equivStrainType, soft;
    IsotropicDamageMaterial :: initializeFrom(ir);
    RandomMaterialExtensionInterface :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    // specify the type of formula for equivalent strain
    IR_GIVE_OPTIONAL_FIELD(ir, equivStrainType, IFT_IsotropicDamageMaterial1_equivstraintype, "equivstraintype"); // Macro
    if ( equivStrainType == 1 ) {
        this->equivStrainType = EST_Rankine;
    } else if ( equivStrainType == 2 ) {
        this->equivStrainType = EST_ElasticEnergy;
    } else if (equivStrainType == 3) { 
      this->equivStrainType = EST_Mises;
      IR_GIVE_FIELD (ir, k, IFT_IsotropicDamageMaterial1_k, "k"); 
    } else {
        this->equivStrainType = EST_Mazars; // default
    }

    // specify the type of formula for damage evolution law
    soft = 0;
    IR_GIVE_OPTIONAL_FIELD (ir, soft, IFT_IsotropicDamageMaterial1_softeningtype, "damlaw"); 
    IR_GIVE_FIELD (ir, e0, IFT_IsotropicDamageMaterial1_e0, "e0"); 
    switch(soft){
    case 1: // linear law
      if (ir->hasField(IFT_IsotropicDamageMaterial1_wf, "wf")) {
	this->softType = ST_Linear_Cohesive_Crack;
	IR_GIVE_FIELD (ir, ef, IFT_IsotropicDamageMaterial1_wf, "wf");
      } else {
	this->softType = ST_Linear;
	IR_GIVE_FIELD (ir, ef, IFT_IsotropicDamageMaterial1_ef, "ef");
      } 
      break;
      //case 2: // reserved for bilinear softening
      //case 3: // reserved for Hordijk's law
    case 4:
      this->softType = ST_Mazars;
      IR_GIVE_FIELD (ir, At, IFT_IsotropicDamageMaterial1_At, "at"); 
      IR_GIVE_FIELD (ir, Bt, IFT_IsotropicDamageMaterial1_Bt, "bt"); 
      break;
    case 5:
      this->softType = ST_Smooth;
      md = 1.;
      IR_GIVE_OPTIONAL_FIELD (ir, md, IFT_IsotropicDamageMaterial1_md, "md"); 
      break;
    default:
      if (ir->hasField(IFT_IsotropicDamageMaterial1_wf, "wf")) {
	this->softType = ST_Exponential_Cohesive_Crack;
	IR_GIVE_FIELD (ir, ef, IFT_IsotropicDamageMaterial1_wf, "wf");
      } else { 
	this->softType = ST_Exponential;
	IR_GIVE_FIELD (ir, ef, IFT_IsotropicDamageMaterial1_ef, "ef"); 
      }
    }
    
    this->mapper.initializeFrom(ir);
    
    return IRRT_OK;
}


int
IsotropicDamageMaterial1 :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    IsotropicDamageMaterial :: giveInputRecordString(str, keyword);
    linearElasticMaterial->giveInputRecordString(str, false);
    sprintf(buff, " e0 %e ef %e equivstraintype %d", this->e0, this->ef, ( int ) this->equivStrainType);
    str += buff;

    return 1;
}




void
IsotropicDamageMaterial1 :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime)
{
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }

    if ( this->equivStrainType == EST_Mazars ) {
        double posNorm = 0.0;
        FloatArray principalStrains, fullstrain;

        crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);

        // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            double nu = lmat->give(NYxz,gp);
            fullstrain.at(3) = -nu * ( fullstrain.at(1) + fullstrain.at(2) ) / ( 1. - nu );
        } else if ( gp->giveMaterialMode() == _1dMat ) {
            double nu = lmat->give(NYxz,gp);
            fullstrain.at(2) = -nu *fullstrain.at(1);
            fullstrain.at(3) = -nu *fullstrain.at(1);
        }

        this->computePrincipalValues(principalStrains, fullstrain, principal_strain);

        for ( int i = 1; i <= 3; i++ ) {
            if ( principalStrains.at(i) > 0.0 ) {
                posNorm += principalStrains.at(i) * principalStrains.at(i);
            }
        }

        kappa = sqrt(posNorm);
    } else if ( this->equivStrainType == EST_Rankine ) {
        // EST_Rankine equiv strain measure
        int i;
        FloatMatrix de;
        FloatArray stress, fullStress, principalStress;
        double sum = 0.;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        crossSection->giveFullCharacteristicVector(fullStress, gp, stress);
        this->computePrincipalValues(principalStress, fullStress, principal_stress);
        for ( i = 1; i <= 3; i++ ) {
            if ( principalStress.at(i) > 0.0 ) {
                sum += principalStress.at(i) * principalStress.at(i);
            }
        }

        kappa = sqrt(sum) / lmat->give('E',gp);
    } else if ( this->equivStrainType == EST_ElasticEnergy ) {
        FloatMatrix de;
        FloatArray stress;
        double sum;

        lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
        stress.beProductOf(de, strain);
        sum = dotProduct( strain, stress, strain.giveSize() );

        kappa = sqrt( sum / lmat->give('E',gp) );
    } else if (this->equivStrainType == EST_Mises){
      double nu=lmat->give (NYxz,NULL);
      FloatArray principalStrains, fullstrain;  
      crossSection->giveFullCharacteristicVector(fullstrain, gp, strain);
      if (gp->giveMaterialMode() == _PlaneStress) {
	fullstrain.at(3) = -nu * (fullstrain.at(1)+fullstrain.at(2))/(1.-nu);
      } else if (gp->giveMaterialMode() == _1dMat) {
	fullstrain.at(2) = -nu * fullstrain.at(1);
	fullstrain.at(3) = -nu * fullstrain.at(1);
      }
      this->computePrincipalValues (principalStrains, fullstrain, principal_strain);
      double I1e,J2e;
      this->computeStrainInvariants(&principalStrains,&I1e,&J2e);
      double a,b,c; 
      a=(k-1)*I1e/(2*k*(1-2*nu));
      b=(k-1)*(k-1)*I1e*I1e/((1-2*nu)*(1-2*nu));
      c=12*k*J2e/((1+nu)*(1+nu));
      kappa = a + 1/(2*k)*sqrt(b+c);
    } else {
      _error("computeEquivalentStrain: unknown EquivStrainType");
    }
}

void
IsotropicDamageMaterial1:: computeStrainInvariants(FloatArray* strainVector, double* I1e, double* J2e) 
{
  *I1e = strainVector->at(1)+strainVector->at(2)+strainVector->at(3);
  double s1=strainVector->at(1)*strainVector->at(1);
  double s2 = strainVector->at(2)*strainVector->at(2);
  double s3 =strainVector->at(3)*strainVector->at(3);
  *J2e=1./2.*(s1+s2+s3)-1./6.*((*I1e)*(*I1e));
}

  /* OLD VERSION, ABANDONED ON 20 JULY 2010
void
IsotropicDamageMaterial1 :: computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp)
{
    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);

    if ( kappa > e0 ) {
        // omega = 1.0-(e0/kappa)*exp(-(kappa-e0)/(ef-e0));
        IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
        int nite = 0;
        double R, Lhs, E, Ft, help;

        //
        // iterations with Newton-Raphson method to achieve mesh objectivity
        // we are finding a state, where the elastic stress is equal to
        // the stress from crack-opening relation (ef = wf characterizes the crack opening diagram)

        omega = 0.0;
        E = this->giveLinearElasticMaterial()->give('E',gp);
        Ft = E * e0;
        //printf ("\nle=%f, kappa=%f", status->giveLe(), kappa);
        do {
            nite++;
            help = status->giveLe() * omega * kappa / ef;
            R = ( 1. - omega ) * E * kappa - Ft *exp(-help);
            Lhs = E * kappa - Ft *exp(-help) * status->giveLe() * kappa / ef;
            omega += R / Lhs;
            // printf ("\n%d: R=%f, omega=%f",nite, R, omega);
            if ( nite > 40 ) {
                _error("computeDamageParam: algorithm for the crack-opening objectivity not converging");
            }
        } while ( fabs(R) >= IDM1_ITERATION_LIMIT );

        if ( ( omega > 1.0 ) || ( omega < 0.0 ) ) {
            _error("Damage parameter out of range <0.0;1.0>");
        }
    } else {
        omega = 0.0;
    }
}
  */

void 
IsotropicDamageMaterial1 :: computeDamageParam (double& omega, double kappa, const FloatArray& strain, GaussPoint* gp) 
{
  if (this->softType==ST_Exponential_Cohesive_Crack || this->softType==ST_Linear_Cohesive_Crack) {
    // adjustment of softening law according to the element size
    computeDamageParamForCohesiveCrack (omega, kappa, gp);
  } else {
    // no adjustment according to element size
    omega = damageFunction(kappa,gp);
  }
}

void 
IsotropicDamageMaterial1 :: computeDamageParamForCohesiveCrack (double& omega, double kappa, GaussPoint* gp) 
{
  omega = 0.0; 
  const double e0 = this->give(e0_ID, gp); // e0 can be a random property
  if (kappa > e0) {
    double ef = this->give(ef_ID, gp); // ef can be a random property
    IsotropicDamageMaterial1Status *status = (IsotropicDamageMaterial1Status*) this -> giveStatus (gp);
    double Le = status->giveLe();
    if (this->softType==ST_Linear_Cohesive_Crack){
      ef /= Le; // crack opening divided by element size gives strain
      if (kappa < ef) {
	omega = (ef/kappa) * (kappa-e0) / (ef-e0);
      } else { 
	omega = 1.0;
      }
      return;
    }
    // exponential cohesive crack - iteration needed
    double R, Lhs, help;
    int nite = 0;
    // iteration to achieve objectivity
    // we are looking for a state in which the elastic stress is equal to
    // the stress from crack-opening relation 
    // (ef = wf has the meaning of crack opening, not of strain)
    do {
      nite++;
      help = Le*omega*kappa/ef;
      R = (1.-omega)*kappa - e0*exp(-help);
      Lhs = kappa-e0*exp(-help)*Le*kappa/ef;
      omega += R/Lhs;
      if (nite>40) _error ("computeDamageParamForCohesiveCrack: algorithm not converging");
    } while (fabs(R) >= e0*IDM1_ITERATION_LIMIT);
    if ((omega > 1.0) || (omega < 0.0)) _error ("computeDamageParamForCohesiveCrack: internal error");
  }
}

double
IsotropicDamageMaterial1 :: damageFunction(double kappa, GaussPoint* gp)
{
  const double e0 = this->give(e0_ID, gp); // e0 can be a random property
  double ef = 0.;
  if (softType==ST_Linear || softType==ST_Exponential) {
    ef = this->give(ef_ID, gp); // ef can be a random property
  }

  switch (softType){

  case ST_Linear:  
    if (kappa <= e0)
      return 0.0;
    else if (kappa < ef)
      return (ef/kappa) * (kappa-e0) / (ef-e0);
    else 
      return 1.0;

  case ST_Exponential:
    if (kappa > e0)
      return 1.0 - (e0/kappa) * exp(-(kappa-e0)/(ef-e0));
    else 
      return 0.0;

  case ST_Mazars:
    return 1.0 - (1.0-At)*e0/kappa - At*exp(-Bt*(kappa-e0));

  case ST_Smooth:	
    return 1.0-exp(-pow(kappa/e0,md));

    /*
  case ST_SmoothExtended: 
  // this was used in fitting of a specific load-displacement diagram
  // constructed by the particle model (collaboration with Peter Grassl)
    if (kappa <= e1)
      return 1.0 - exp(-pow(kappa/e0,md));
    else
      return 1.0 - s1*exp(-(kappa-e1)/(ef*(1.+pow((kappa-e1)/e2,nd))))/kappa;
    */

  default:
    printf("IsotropicDamageMaterial1::damageFunction ... undefined softening type %d\n",softType);
  }
}

double
IsotropicDamageMaterial1 :: complianceFunction(double kappa, GaussPoint* gp)
{
  double om = damageFunction(kappa,gp);
  return om/(1.-om);
}

void
IsotropicDamageMaterial1 :: initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
    FloatMatrix principalDir(3, 3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);
    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    const double e0 = this->give(e0_ID, gp);
    const double ef = this->give(ef_ID, gp);


    crossSection->giveFullCharacteristicVector(fullstrain, gp, strainVector);

    if ( ( kappa > e0 ) && ( status->giveDamage() == 0. ) ) {
        this->computePrincipalValDir(principalStrains, principalDir, fullstrain, principal_strain);
        // find index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        le = gp->giveElement()->giveCharacteristicLenght(gp, crackPlaneNormal);
        // remember le in corresponding status
        status->setLe(le);

        if ( e0 >= ef / le ) {
            _warning4("Fracturing strain ef/le=%f is lower than the elastic strain e0=%f, possible snap-back. Increase crack opening ef to %f", ef / le, e0, e0 * le);
        }
    }
}

double
IsotropicDamageMaterial1 :: give(int aProperty, GaussPoint* gp)
{
  double answer;
  if (RandomMaterialExtensionInterface::give(aProperty, gp, answer)) {
    return answer;
  } else if (aProperty == e0_ID) {
    return this->e0;
  } else if (aProperty == ef_ID) {
    return this->ef;
  } else {
    return IsotropicDamageMaterial::give(aProperty, gp);
  }
}


Interface *
IsotropicDamageMaterial1 :: giveInterface(InterfaceType type)
{
    if ( type == MaterialModelMapperInterfaceType ) {
        return ( MaterialModelMapperInterface * ) this;
    } else {
        return NULL;
    }
}


MaterialStatus*
IsotropicDamageMaterial1::CreateStatus(GaussPoint *gp) const
{
  IsotropicDamageMaterial1Status* answer = new IsotropicDamageMaterial1Status(1, IsotropicDamageMaterial1 :: domain, gp);
  return answer;
}

MaterialStatus *
IsotropicDamageMaterial1 :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status;

    status = gp->giveMaterialStatus();
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus(status);
	    this->_generateStatusVariables (gp);
        }
    }

    return status;
}


int
IsotropicDamageMaterial1 :: MMI_map(GaussPoint *gp, Domain *oldd, TimeStep *tStep)
{
    int result;
    FloatArray intVal, strainIncr(3);
    IntArray toMap(3);
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);


    toMap.at(1) = ( int ) IST_MaxEquivalentStrainLevel;
    toMap.at(2) = ( int ) IST_DamageTensor;
    toMap.at(3) = ( int ) IST_StrainTensor;
    this->mapper.init(oldd, toMap, gp, tStep);

    result = mapper.mapVariable(intVal, gp, IST_MaxEquivalentStrainLevel, tStep);
    if ( result ) {
        status->setTempKappa( intVal.at(1) );
    }

    result = mapper.mapVariable(intVal, gp, IST_DamageTensor, tStep);
    if ( result ) {
        status->setTempDamage( intVal.at(1) );
    }

#ifdef IDM_USE_MAPPEDSTRAIN
    result = mapper.mapVariable(intVal, gp, IST_StrainTensor, tStep);
    if ( result ) {
        status->letTempStrainVectorBe(intVal);
    }

#endif
    status->updateYourself(tStep);

    return result;
}




int
IsotropicDamageMaterial1 :: MMI_update(GaussPoint *gp,  TimeStep *tStep, FloatArray *estrain)
{
    int result = 1;
    FloatArray intVal, strain;
    IsotropicDamageMaterial1Status *status = ( IsotropicDamageMaterial1Status * ) this->giveStatus(gp);

    // now update all internal vars accordingly
    strain = status->giveStrainVector();
#ifdef IDM_USE_MAPPEDSTRAIN
    this->giveRealStressVector(intVal, ReducedForm, gp, strain, tStep);
#else
    this->giveRealStressVector(intVal, ReducedForm, gp, * estrain, tStep);
#endif
    this->updateYourself(gp, tStep);
    return result;
}


int
IsotropicDamageMaterial1 :: MMI_finish(TimeStep *tStep)
{
    this->mapper.finish(tStep);
    return 1;
}


IsotropicDamageMaterial1Status :: IsotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g) :
  IsotropicDamageMaterialStatus(n, d, g), RandomMaterialStatusExtensionInterface()
{
    le = 0.0;
}

void
IsotropicDamageMaterial1Status :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterialStatus :: initTempStatus();
}



void
IsotropicDamageMaterial1Status :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterialStatus :: updateYourself(atTime);
}

Interface *
IsotropicDamageMaterial1Status :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return ( RandomMaterialStatusExtensionInterface * ) this;
    } else {
        return NULL;
    }
}


contextIOResultType
IsotropicDamageMaterial1Status :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
IsotropicDamageMaterial1Status :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

} // end namespace oofem
