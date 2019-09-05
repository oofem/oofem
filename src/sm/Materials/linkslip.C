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

#include "linkslip.h"
#include "linearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"
// #ifdef __TM_MODULE
//  #include "latticetransportelement.h"
//  #include "latticetransmat.h"
//  #include "latticelammat.h"
//  #include "pore.h"
//  #include "discretetransportproblem.h"
// #endif

namespace oofem {
  REGISTER_Material(LinkSlip);

  /// constructor which creates a dummy material without a status and without random extension interface
  LinkSlip :: LinkSlip(int n, Domain *d) : LinearElasticMaterial(n, d)
  {
  }



  LinkSlip :: ~LinkSlip()
  //
  // destructor
  //
  {}

  int
  LinkSlip :: hasMaterialModeCapability(MaterialMode mode)
  {
    if ( mode == _3dMat ) {
      return 1;
    }

    return 0;
  }


  IRResultType
  LinkSlip :: initializeFrom(InputRecord *ir)
  {
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    LinearElasticMaterial :: initializeFrom(ir);

    //axial stiffness
    IR_GIVE_FIELD(ir, kNormal, _IFT_LinkSlip_kn); // Macro

    //Ratio of lateral to axial stiffness
    alphaOne = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LinkSlip_a1); // Macro

    //Two models available
    //0 default. Linear elastic-perfect plastic
    //1 bond slip according first part of CEB model and then constant
    //2 bond slip according to CEB model
    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_LinkSlip_type); // Macro
    
    IR_GIVE_FIELD(ir, tauMax, _IFT_LinkSlip_t0); // Macro

    if(type == 1 || type == 2){
      IR_GIVE_FIELD(ir, s1, _IFT_LinkSlip_s1);
      IR_GIVE_FIELD(ir, alpha, _IFT_LinkSlip_alpha);
    }
    
    if(type == 2){
      IR_GIVE_FIELD(ir, s2, _IFT_LinkSlip_s2);
      IR_GIVE_FIELD(ir, s3, _IFT_LinkSlip_s3);
      IR_GIVE_FIELD(ir, tauFinal, _IFT_LinkSlip_tf);
    }

    if(type == 1 || type == 2){
      if ((type == 1 || type == 2) && this->kNormal < this->tauMax/this->s1) {
	this->kNormal = this->tauMax/this->s1;
	OOFEM_WARNING("Parameter kN adjusted");
      }
    }
    
    return IRRT_OK;
  }


  MaterialStatus *
  LinkSlip :: CreateStatus(GaussPoint *gp) const
  {
    LinkSlipStatus *answer = new LinkSlipStatus(gp);

    return answer;
  }
  
  double
  LinkSlip :: evaluateBondStress(const double kappa) const
  {
    if(this->type == 0){//elastic-perfectly plastic
      return this->tauMax;
    }
    else if(this->type == 1){//modified bond model      
      if ( kappa <= 0. ){
        return 0.;
      }
      if ( kappa <= s1 ){
        return tauMax * pow(kappa/s1, alpha);
      }     
      return tauMax;
    }
    else if(this->type == 2){
      if ( kappa <= 0. ){
        return 0.;
      }
      if ( kappa <= s1 ){
        return tauMax * pow(kappa/s1, alpha);
      }
      if ( kappa <= s2 ){
        return tauMax;
      }
      if ( kappa <= s3 ){
        return tauMax - (tauMax-tauFinal) * (kappa-s2) / (s3-s2);
      }
      return tauFinal;
    }
    else{//unknown type
      OOFEM_ERROR("Unknown bond model type. Type should be 0, 1 or 2.");
    }
    
    return 0.; //Should not be here. 
  }



  void
  LinkSlip :: giveRealStressVector_3d(FloatArray &answer,
				   GaussPoint *gp,
				   const FloatArray &totalStrain,
				   TimeStep *atTime)
  {
    LinkSlipStatus *status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );

    //strain has the meanig of slip. Stress is force
    
    FloatArray strainVector;

    double tempPlasticStrain = status->givePlasticStrain();
    double tempKappa = status->giveKappa();

    FloatMatrix stiffnessMatrix;
    this->giveStiffnessMatrix(stiffnessMatrix, ElasticStiffness, gp, atTime);

    answer.resize(3);
    answer.zero();

    /*First component is the slip one for which the stress should be limited using plasiticity (frictional slip between fibre and matrix). The other components are kept elastic. */
    answer.at(1) = (totalStrain.at(1)-tempPlasticStrain)*stiffnessMatrix.at(1, 1);
    
    
    double f = fabs(answer.at(1)) - evaluateBondStress(tempKappa);

    if(f>0){//plastic response.
      //Reduced stress by increasing plastic strain.
      tempPlasticStrain = tempPlasticStrain + sgn(answer.at(1))*f/stiffnessMatrix.at(1, 1);
      tempKappa = tempKappa + f/stiffnessMatrix.at(1,1);
      answer.at(1) = (totalStrain.at(1)-tempPlasticStrain)*stiffnessMatrix.at(1, 1);
    }

    //Compute the final stress components
    for ( int i = 2; i <= 3; i++ ) { // only diagonal terms matter
      answer.at(i) =  stiffnessMatrix.at(i, i) * totalStrain.at(i);
    }

    //Set temp values in status needed for dissipation
    status->letTempKappaBe(tempKappa);
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

    return;
  }

  Interface *
  LinkSlip :: giveInterface(InterfaceType type)
  {
    return NULL;
  }



  void
  LinkSlip :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
  {
 
    /* Returns elastic moduli in reduced stress-strain space*/
    answer.resize(3, 3);
    answer.zero();

    answer.at(1, 1) = 1.;
    answer.at(2, 2) = this->alphaOne; // shear
    answer.at(3, 3) = this->alphaOne; // shear    
    answer.times(this->kNormal);    
  }

  LinkSlipStatus :: LinkSlipStatus(GaussPoint *g) :  StructuralMaterialStatus(g)
  {

  }
  

  void
  LinkSlipStatus :: initTempStatus()
  //
  // initializes temp variables according to variables form previous equlibrium state.
  // builds new crackMap
  //
  {
    StructuralMaterialStatus :: initTempStatus();
    this->tempPlasticStrain = this->plasticStrain;
  }
  

  void
  LinkSlip :: giveThermalDilatationVector(FloatArray &answer,
					  GaussPoint *gp,  TimeStep *tStep)
  //
  // returns a FloatArray(3) of initial strain vector
  // caused by unit temperature in direction of
  // gp (element) local axes
  //
  {
    double alpha = this->give(tAlpha, gp);
  
    answer.resize(3);
    answer.zero();

    answer.at(1) = alpha;
  }

  void
  LinkSlipStatus :: printOutputAt(FILE *file, TimeStep *tStep)
  {
    MaterialStatus :: printOutputAt(file, tStep);

        fprintf(file, "  jump ");
    for ( auto &val : this->strainVector ) {
        fprintf(file, " %.4e", val );
    }

    fprintf(file, "\n              traction ");
    for ( auto &val : this->stressVector ) {
        fprintf(file, " %.4e", val );
    }
    fprintf(file, "\n");
    
    fprintf(file, "plasticStrain %.8e, kappa %.8e\n", this->plasticStrain, this->kappa);
    return;
  }
  
  void
  LinkSlipStatus :: saveContext(DataStream &stream, ContextMode mode)
  //
  // saves full information stored in this Status
  // no temp variables stored
  //
  {
    // save parent class status
  
    StructuralMaterialStatus :: saveContext(stream, mode);
  
    // write a raw data
    if ( !stream.write(plasticStrain) ) {
      THROW_CIOERR(CIO_IOERR);
    }

    // write a raw data
    if ( !stream.write(kappa) ) {
      THROW_CIOERR(CIO_IOERR);
    }
    
  }


  void
  LinkSlipStatus :: restoreContext(DataStream &stream, ContextMode mode)
  //
  // restores full information stored in stream to this Status
  //
  {
    StructuralMaterialStatus :: restoreContext(stream, mode);
    
    // read raw data
    if ( !stream.read(plasticStrain) ) {
      THROW_CIOERR(CIO_IOERR);
    }
  
    // read raw data
    if ( !stream.read(kappa) ) {
      THROW_CIOERR(CIO_IOERR);
    }
   
  }
 
  void
  LinkSlipStatus :: updateYourself(TimeStep *atTime)
  //
  // updates variables (nonTemp variables describing situation at previous equilibrium state)
  // after a new equilibrium state has been reached
  // temporary variables are having values corresponding to newly reached equilibrium.
  //
  {
    StructuralMaterialStatus :: updateYourself(atTime);
    this->plasticStrain = this->tempPlasticStrain;
    this->kappa = this->tempKappa;
  }

 
  int
  LinkSlip :: giveIPValue(FloatArray &answer,
			  GaussPoint *gp,
			  InternalStateType type,
			  TimeStep *atTime)
  {
    // return LinearElasticMaterial :: giveIPValue(answer, gp, type, atTime);
    return Material :: giveIPValue(answer, gp, type, atTime);
  }
}
