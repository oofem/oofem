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

#include "fcm.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "stressvector.h"

#include <cstring>

namespace oofem {
FCMMaterial :: FCMMaterial(int n, Domain *d) : StructuralMaterial(n, d),
    linearElasticMaterial(n, d),
    ecsMethod(ECSM_Unknown)
{}


bool
FCMMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain;
}


void
FCMMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                    const FloatArray &totalStrain,
                                    TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
    // g = global direction
    // l = local (crack) direction
    // normal = normal components only (in direction of crack planes, "l")

    // stress and strain transformation matrices
    FloatMatrix epsG2L, sigmaL2G, sigmaG2L;
    // stiffness matrices (total, elastic, cracking)
    FloatMatrix D, De, Dcr, DeRed;
    FloatMatrix principalDirs;

    // STRESSES
    FloatArray stressIncrement_l, stressVector_g, trialStress_g, trialStress_l, sigmaResid, sigmaElast_l, sigmaCrack_l, sigmaResid_old;
    FloatArray stressVector_l;
   
    // STRAINS
    FloatArray reducedStrain_g, strainIncrement_g, strainIncrement_l, crackStrain, crackStrainIncrement, elasticStrain_l, elasticStrainIncrement_l;     //    FloatArray reducedStrain_l

    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    MaterialMode mMode = gp->giveMaterialMode();

    int nCr = status->giveNumberOfCracks();
    int nMaxCr = status->giveMaxNumberOfCracks(gp);

    double maxErr = 0.;
    double maxTau = 0.;

    double G;

    double DII;

    double gamma_cr, d_gamma_cr;
    double tau_el, tau_cr;
    double d_tau;
    double d_tau_old = 0.;
    bool illinoisFlag = false;

    int iterLimitGlobal = 20;
    int iterLimitGradient = 20;
    int iterLimitNormal = 100;
    int iterLimitShear = 100;
    int indexLimit = 10;

    FloatArray tempNormalCrackStrain, normalCrackStrain;
    //    FloatArray normalStrainVector_l;
    FloatArray normalStrainIncrement_l;
    FloatArray normalStressVector_l;
    FloatArray normalCrackStrainIncrement;
    FloatArray normalElasticStrainIncrement_l;
    FloatArray normalElasticStrain_l;


    // for the bisection method:
    int index, indexOld, indexCount;

    bool changeIndexFlag = false;

    bool plus, minus;
    FloatArray crackStrainPlus, crackStrainMinus, residPlus, residMinus; //local

    FloatArray helpRF;

    int iterN;
    
    plus = false;
    minus = false;

    sigmaResid_old.resize(nMaxCr);
    
    index = indexOld = indexCount = 0;

    
    this->initTempStatus(gp);

    if ( !this->isActivated(tStep) ) {
        FloatArray zeros;
        zeros.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        zeros.zero();
        status->letTempStrainVectorBe(zeros);
        status->letTempStressVectorBe(zeros);
        answer = zeros;
        return;
    }

    // get elastic stiffness matrix
    this->giveStiffnessMatrix(De, ElasticStiffness, gp, tStep);

    // equilibrated global stress
    stressVector_g = status->giveStressVector();


    this->giveStressDependentPartOfStrainVector(reducedStrain_g, gp, totalStrain, tStep, VM_Incremental);
    strainIncrement_g.beDifferenceOf( reducedStrain_g, status->giveStrainVector() );

    trialStress_g.beProductOf(De, strainIncrement_g);
    trialStress_g.add(stressVector_g);
  
    // ELASTIC CASE - NO CRACKING SO FAR (in prevs steps)
    if ( nCr == 0 ) {

      this->computePrincipalValDir(trialStress_l, principalDirs, trialStress_g, principal_stress);

        // REMAINS ELASTIC or element cracking is prevented from above
        if ( !this->isStrengthExceeded( principalDirs, gp, tStep, 1, trialStress_l.at(1) ) ) {
            answer = trialStress_g;
            status->letTempStrainVectorBe(totalStrain);
            status->letTempStressVectorBe(trialStress_g);
            return;

            // STARTS CRACKING - 1st crack
        } else {
	  this->initializeCrack(gp, tStep, principalDirs, 1);

            for ( int iCrack = 2; iCrack <= min(nMaxCr, this->nAllowedCracks); iCrack++ ) {
                if ( this->isStrengthExceeded( principalDirs, gp, tStep, iCrack, trialStress_l.at(iCrack) ) ) { // for "nonlocal model"
		  this->initializeCrack(gp, tStep, principalDirs, iCrack);
                }
            }
        }

        // CRACKING OCCURED IN PREVIOUS STEPS
    } else {
        // check other principal directions - if strength is not exceeded, take transformation matrices from previous analysi
        if ( nCr < nMaxCr ) {

            principalDirs = status->giveCrackDirs();

            for ( int iCrack = nCr + 1; iCrack <= min(nMaxCr, this->nAllowedCracks); iCrack++ ) {
                // condition to guarantee that third  crack won't be initiated if second crack does not exist
                if ( ( status->giveNumberOfTempCracks() + 1 ) < iCrack ) {
                    break;
                }

                if ( !this->checkStrengthCriterion(principalDirs, trialStress_g, gp, tStep, iCrack) ) { // = strength of composite
                    // second and third crack is now initialized
		  this->initializeCrack(gp, tStep, principalDirs, iCrack);
                }
            } // end for
        } // end nMaxCr
    } // cracking in previous steps


    epsG2L = status->giveG2LStrainVectorTransformationMtrx();
    sigmaG2L = status->giveG2LStressVectorTransformationMtrx();
    sigmaL2G = status->giveL2GStressVectorTransformationMtrx();


    // PREPARE STRAIN AND STRESS INCREMENT IN CRACK DIRECTIONS

    // strain and strain increment in crack direction
    strainIncrement_l.beProductOf(epsG2L, strainIncrement_g); // from global to local
    
    // stress increment in crack direction
    stressIncrement_l.beProductOf(De, strainIncrement_l);

    // (equilibrated) stress in crack direction

    stressVector_l.beProductOf(sigmaG2L, status->giveStressVector() );
    

    // SIMILARLY TO RCM2: THE STRESS INCREMENT IN ELASTIC
    // AND CRACKING UNIT MUST BE EQUAL - ITERATE

    // average shear modulus
    G = computeOverallElasticShearModulus(gp, tStep);
    
    DeRed = De;
    DeRed.resizeWithData(nMaxCr,nMaxCr);
    
    crackStrain = status->giveCrackStrainVector();
    
    tempNormalCrackStrain = status->giveCrackStrainVector();
    tempNormalCrackStrain.resize(nMaxCr);
    
    normalCrackStrain = tempNormalCrackStrain;
    
    
    normalStressVector_l = stressVector_l;
    normalStressVector_l.resize(nMaxCr);
    
    normalStrainIncrement_l = strainIncrement_l;
    normalStrainIncrement_l.resizeWithValues(nMaxCr);
    
    for ( int iterG = 1; iterG <= iterLimitGlobal; iterG++ ) {
      
      if ( iterG == 1 ) {
	// residuum
	sigmaResid = stressIncrement_l;
      }
      
      for  ( int indexCount = 1; indexCount <= indexLimit; indexCount++ ) {

	if (indexCount > 1 ) {
	  if ( !changeIndexFlag ) {// same index as before - do not even try to fix it
	    OOFEM_WARNING("FCM: maximum error on the same component, max. stress error %f, index %d, indexCount %d, element %d\n", maxErr, index, indexCount, gp->giveElement()->giveNumber() );
	    
	    break;
	  }
	}

	minus = plus = false;
	illinoisFlag = false;
	changeIndexFlag = false;
	iterN = 0;
	sigmaResid_old.zero();
	

	do { // find convergence in the crack directions
	  
	  iterN++;
	
	  if ( (iterN == 1) && (indexCount == 1) ) {
	    // shrink the total residuum
	    sigmaResid.resizeWithValues(nMaxCr);  
	  }
	
	  this->giveNormalLocalCrackedStiffnessMatrix(Dcr, TangentStiffness, gp, tStep);
	  
	  D = DeRed;
	  D.add(Dcr);
	  
	  D.solveForRhs(sigmaResid, normalCrackStrainIncrement);
	  
	  // update and store cracking strain
	  tempNormalCrackStrain.add(normalCrackStrainIncrement);
	  status->setTempNormalCrackStrainVector(tempNormalCrackStrain);
	  
	  // update statuses
	  this->updateCrackStatus(gp);
	  
	  // compute and compare stress in elastic and cracking unit:
	  // CR
	  this->giveNormalLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);
	  
	  sigmaCrack_l.resize( nMaxCr );
	  for ( int i = 1; i <= nMaxCr; i++ ) {
	  sigmaCrack_l.at(i) = Dcr.at(i, i) * tempNormalCrackStrain.at(i);
	  }
	  
	  // compute and compare stress in elastic and cracking unit:     
	  // increment of elastic strain 
	  // d_eps_el = d_eps - d_eps_cr
	  // d_eps_cr = eps_cr_n - eps_cr_n-1
	  normalElasticStrainIncrement_l = normalStrainIncrement_l;
	  normalElasticStrainIncrement_l.subtract(tempNormalCrackStrain);
	  normalElasticStrainIncrement_l.add(normalCrackStrain);
	  
	  // d_sigma = Del * d_eps_el
	  sigmaElast_l.beProductOf(DeRed, normalElasticStrainIncrement_l);
	  
	  // sigma_n = sigma_n-1 + d_sigma
	  sigmaElast_l.add(normalStressVector_l);
	  
	  // residuum
	  sigmaResid = sigmaElast_l;
	  sigmaResid.subtract(sigmaCrack_l);
	  
	  maxErr = 0.;
	  indexOld = index;
	  
	  for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
	    if ( fabs( sigmaResid.at(i) ) > maxErr ) {
	      maxErr = fabs( sigmaResid.at(i) );

	      //	      if (indexCount == 1) {
		index = i;
		//	      }
	    }
	  }

	  

	  if ( maxErr <= fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) {
	    break; // convergence criteria satisfied
	  }
	  
	  // try only gradient method and do not slow down with the preparation for the bisection method
	  if ( indexCount == 1 ) {
	    
	    if (iterN < iterLimitGradient/2.) {
	      continue;
	      
	    }
	  }
	  // identify the biggest error only in the first run of the gradient method
	  //} else {
	  
	  if (index != indexOld) {
	    minus = plus = false;
		//	crackStrainPlus.resize(nMaxCr);
		/*crackStrainPlus.zero();
		//crackStrainMinus.resize(nMaxCr);
		crackStrainMinus.zero();
		//	residPlus.resize(nMaxCr);
		residPlus.zero();
		//residMinus.resize(nMaxCr);
		residMinus.zero();*/
	  }
	  
	
	
	  if ( sigmaResid.at(index) > 0 ) {
	    if ( !plus) { // first positive error
	      plus = true;
	      residPlus = sigmaResid;
	      crackStrainPlus = tempNormalCrackStrain;
	    } else if ( fabs(sigmaResid.at(index)) < fabs(residPlus.at(index)) ) { // smaller error
	      residPlus = sigmaResid;
	      crackStrainPlus = tempNormalCrackStrain;
	    }
	    
	  } else {
	    if ( !minus) { // first negative error
	      minus = true;
	      residMinus = sigmaResid;
	      crackStrainMinus = tempNormalCrackStrain;
	    }
	    
	    else if ( fabs(sigmaResid.at(index)) < fabs(residMinus.at(index)) ) { // smaller error
	      residMinus = sigmaResid;
	      crackStrainMinus = tempNormalCrackStrain;
	    }
	  }
	  
	
	  //	  if ( (! plus) || (! minus) ) {

	  if ( ( plus) && ( minus) ) {
	    if ( ( indexCount == 1 ) && ( iterN < iterLimitGradient ) ) {
	      continue;
	    }

	    

	    //( iterN < iterLimitGradient ) || (! plus) || (! minus) ) {

	  } else {
	    continue;
	  }
	  
	  do { // perform bisection method
	    
	    // check if the cracking strain at index is negative and try to fix it
	    tempNormalCrackStrain.zero();
	    
	    if(  (crackStrainMinus.at(index) < 0.) || (crackStrainPlus.at(index) < 0.) ) {
	      tempNormalCrackStrain.add(crackStrainMinus);
	      tempNormalCrackStrain.add(crackStrainPlus);
	      tempNormalCrackStrain.times(0.5);
	      
	    } else {
	      // OOFEM_WARNING("Trying regula-falsi iterations\n");
		  
	      // try regula-falsi
	      // eps_cr_new = (eps_cr_min  * resid_plus - eps_cr_plus * resid_minus) / (resid_plus - resid_minus)
	      // https://en.wikipedia.org/wiki/False_position_method#The_Regula_Falsi_.28False_Position.29_Method

	      if (illinoisFlag) { // same sign of the error (at least) twice in a row

		if ( sigmaResid.at(index) == residPlus.at(index) ) {  // error is positive

		  tempNormalCrackStrain.add(crackStrainMinus);
		  tempNormalCrackStrain.times(residPlus.at(index) );
		  
		  helpRF = crackStrainPlus;
		  helpRF.times( 0.5 * residMinus.at(index) );
		  
		  tempNormalCrackStrain.subtract(helpRF);
		  tempNormalCrackStrain.times( 1. / ( residPlus.at(index) - 0.5 * residMinus.at(index) ) );

		} else { // error is negative

		  tempNormalCrackStrain.add(crackStrainMinus);
		  tempNormalCrackStrain.times( 0.5 * residPlus.at(index) );
		  
		  helpRF = crackStrainPlus;
		  helpRF.times( residMinus.at(index) );
		  
		  tempNormalCrackStrain.subtract(helpRF);
		  tempNormalCrackStrain.times( 1. / ( 0.5 * residPlus.at(index) - residMinus.at(index) ) );
	  
		}

	      } else { // different sign of the error (consecutive errors)
		// perform standard regula-falsi
		
		tempNormalCrackStrain.add(crackStrainMinus);
		tempNormalCrackStrain.times(residPlus.at(index) );
		
		helpRF = crackStrainPlus;
		helpRF.times(residMinus.at(index) );
		
		tempNormalCrackStrain.subtract(helpRF);
		tempNormalCrackStrain.times( 1. / ( residPlus.at(index) - residMinus.at(index) ) );
	      }
	    }
	    
	    
	    // too small strain and the crack has not existed before
	    if ( ( fabs (crackStrainMinus.at(index) ) <= fcm_THRESHOLD_CRACK_STRAIN ) &&
		 ( fabs (crackStrainPlus.at(index) ) <= fcm_THRESHOLD_CRACK_STRAIN ) &&
		 ( status->giveCrackStatus(index) == 0 ) ) {
	      tempNormalCrackStrain.at(index) = 0.;
	    }
	    
	    // replace negative cracking strain by a zero strain
	    if (tempNormalCrackStrain.at(index) < 0.) {
	      tempNormalCrackStrain.at(index) = 0.;
	    }
	    
	    status->setTempNormalCrackStrainVector(tempNormalCrackStrain);
	    // update statuses - to have correct stiffnesses (unlo-relo vs. softening)
	    this->updateCrackStatus(gp);
	    
	    // compute and compare stress in elastic and cracking unit:
	    // CR
	    this->giveNormalLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);
	    sigmaCrack_l.resize( nMaxCr );
	    for ( int i = 1; i <= nMaxCr; i++ ) {
	      sigmaCrack_l.at(i) = Dcr.at(i, i) * tempNormalCrackStrain.at(i);
	    }
	    
	    // increment of elastic strain 
	    // d_eps_el = d_eps - d_eps_cr
	    // d_eps_cr = eps_cr_n - eps_cr_n-1
	    normalElasticStrainIncrement_l = normalStrainIncrement_l;
	    normalElasticStrainIncrement_l.subtract(tempNormalCrackStrain);
	    normalElasticStrainIncrement_l.add(normalCrackStrain);
	    
	    // d_sigma = Del * d_eps_el
	    sigmaElast_l.beProductOf(DeRed, normalElasticStrainIncrement_l);
	    
	    //	sigma_n = sigma_n-1 + d_sigma
	    sigmaElast_l.add(normalStressVector_l);
	    
	    // residuum
	    sigmaResid = sigmaElast_l;
	    sigmaResid.subtract(sigmaCrack_l);


	    /*  if ( (status->giveTempCrackStatus(index) == 4) && (tempNormalCrackStrain.at(index) < fcm_SMALL_STRAIN) ) {
	      sigmaResid.at(index) = 0.;
	      }*/
	    
	    
	    
	    // adjust residuum if the crack does not exist
	    /* if ( status->giveTempCrackStatus(index) == 0 ) {
	      sigmaResid.at(index) = 0.;
	      }*/
	    
	    // condition for the given stress component
	    if  ( fabs( sigmaResid.at(index) ) < fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) {
	      changeIndexFlag = true;
	      break;
	    }
	    
	    
	    if ( (sigmaResid.at(index) >= 0.) && ( sigmaResid.at(index) < residPlus.at(index) ) ) {
	      residPlus = sigmaResid;
	      crackStrainPlus = tempNormalCrackStrain;
	      
	    } else if ( (sigmaResid.at(index) <= 0.) && ( fabs(sigmaResid.at(index)) < fabs(residMinus.at(index) ) ) ) {

	      residMinus = sigmaResid;
	      crackStrainMinus = tempNormalCrackStrain;
	      
	    } else { // error is increasing -> change index
	      changeIndexFlag = true;
	    }

	    
	    // cracking strain at index component is close to zero
	    if ( ( ( fabs( crackStrainMinus.at(index) - crackStrainPlus.at(index) ) < 1.e-16 ) &&
		 ( max( fabs( crackStrainMinus.at(index) ), fabs( crackStrainPlus.at(index) ) ) < 1.e-14 ) ) ||
		 ( status->giveTempCrackStatus(index) == pscm_NONE ) ) {
	      tempNormalCrackStrain.zero();
	      tempNormalCrackStrain.add(crackStrainMinus);
	      tempNormalCrackStrain.add(crackStrainPlus);
	      tempNormalCrackStrain.times(0.5);
	      
	      tempNormalCrackStrain.at(index) = 0.;

#if DEBUG		  		    
	      //	      OOFEM_WARNING("Fixed crack model: cracking strain component %d set to zero", index);
#endif
		
	      
	      status->setTempNormalCrackStrainVector(tempNormalCrackStrain);
	      
	      // update statuses
	      this->updateCrackStatus(gp);
	      
	      // compute and compare stress in elastic and cracking unit:
	      for ( int i_count = 1; i_count <= 2; i_count++ ) {
		
		// EL
		// increment of elastic strain 
		// d_eps_el = d_eps - d_eps_cr
		// d_eps_cr = eps_cr_n - eps_cr_n-1
		normalElasticStrainIncrement_l = normalStrainIncrement_l;
		normalElasticStrainIncrement_l.subtract(tempNormalCrackStrain);
		normalElasticStrainIncrement_l.add(normalCrackStrain);
		// d_sigma = Del * d_eps_el
		sigmaElast_l.beProductOf(DeRed, normalElasticStrainIncrement_l);
		
		// eliminating stress error at index component
		tempNormalCrackStrain.at(index) = sigmaElast_l.at(index) / this->giveCrackingModulus(SecantStiffness, gp, tStep, index);
		status->setTempNormalCrackStrainVector(tempNormalCrackStrain);
		    
		// update statuses
		this->updateCrackStatus(gp);
	      }
	      
	      // stress in cracking unit
	      this->giveNormalLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);			  
	      sigmaCrack_l.resize( nMaxCr );
	      for ( int i = 1; i <= nMaxCr; i++ ) {
		sigmaCrack_l.at(i) = Dcr.at(i, i) * tempNormalCrackStrain.at(i);
	      }
	      
	      // reevaluate the residuum 
	      sigmaResid = sigmaElast_l;
	      sigmaResid.subtract(sigmaCrack_l);
	      
	    }  // end - crack strain is close to zero

	    // evaluate maximum error 
	    indexOld = index;
	    
	    maxErr = 0.;
	    for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
	      if ( fabs( sigmaResid.at(i) ) > maxErr ) {
		maxErr = fabs( sigmaResid.at(i) );
		index = i;
	      }
	    }

	    // biggest error on different index
	    if (indexOld != index) {
	      changeIndexFlag = true;
	      illinoisFlag = false;

	    } else { // evaluate

	      if ( sgn(sigmaResid.at(index) ) == sgn(sigmaResid_old.at(index) ) ) {
		illinoisFlag = true;
	      } else {
		illinoisFlag = false;
	      }

	      sigmaResid_old = sigmaResid;
	      
	    }

	
	    // end loop bisection method
	  } while ( (!changeIndexFlag) && (iterN <= iterLimitNormal) && ( fabs( sigmaResid.at(index) ) < fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) );
	  

	  // end normal loop
	} while ( ( !changeIndexFlag ) && ( maxErr > fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) && (iterN <= iterLimitNormal) );

	if ( maxErr <= fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) {
	  break; // break index loop
	}

	if ( indexCount >= indexLimit/2.) {
	  OOFEM_WARNING("FCM: slowly converging equilibrium in crack direction!, max. stress error %f, index %d, indexCount %d, element %d\n", maxErr, index, indexCount, gp->giveElement()->giveNumber() );
	}

	if ( indexCount >= indexLimit ) {
	  OOFEM_WARNING("FCM: slowly converging equilibrium in crack direction!, max. stress error %f, index %d, indexCount %d, element %d\n", maxErr, index, indexCount, gp->giveElement()->giveNumber() );
	  break;
	}

	
      }  // end index loop in normal direction
      
      // FIND EQUILIBRIUM FOR SHEAR
      tau_cr = 0.;
      DII = 0.;
      

      double d_gamma_el;
      
      for ( int i = nMaxCr+1; i <= crackStrain.giveSize(); i++ ) {

	plus = false;
	minus = false;
	    
	crackStrainPlus.resize(1);
	crackStrainPlus.zero();
	crackStrainMinus.resize(1);
	crackStrainMinus.zero();
	residPlus.resize(1);
	residPlus.zero();
	residMinus.resize(1);
	residMinus.zero();

	
	d_tau = G * strainIncrement_l.at(i);
	gamma_cr = status->giveCrackStrain(i);
	
	for ( int iterS = 1; iterS <= iterLimitShear; iterS++ ) {
	
	  switch ( mMode ) {
	  case  _PlaneStress:
	    
	    DII = computeTotalD2Modulus(gp, tStep, 6);
	    break;
	    
	  case _PlaneStrain: // check
	    DII = computeTotalD2Modulus(gp, tStep,  6);
	    break;
	    
	  case _3dMat:
	    DII = computeTotalD2Modulus(gp, tStep, i);
	    break;
	    
	  default:
	    OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
	  }

	  
	  if ( ( !plus ) || ( !minus ) ) {
	    
	    if (iterS == 1) {
	      d_gamma_cr = ( stressVector_l.at(i) - DII * status->giveCrackStrain(i) + d_tau ) / ( G + DII );
	    } else {
	      d_gamma_cr = ( d_tau ) / ( G + DII );
	    }
	    
	    gamma_cr += d_gamma_cr;
	  } else {
	    // same sign of the cracking strain
	    if (sgn( crackStrainMinus.at(1) ) == sgn( crackStrainPlus.at(1) ) ) {

	      if (illinoisFlag) { // same sign of the error (at least) twice in a row

		if (d_tau == residPlus.at(1) ) { // error is positive
		  		  gamma_cr = ( crackStrainMinus.at(1) * residPlus.at(1) - 0.5 * crackStrainPlus.at(1) * residMinus.at(1) ) / ( residPlus.at(1) - 0.5 * residMinus.at(1) );

		} else {
		  		  gamma_cr = ( 0.5 * crackStrainMinus.at(1) * residPlus.at(1) - crackStrainPlus.at(1) * residMinus.at(1) ) / ( 0.5 * residPlus.at(1) - residMinus.at(1) );
		}
		
		// different sign of the error (consecutive errors)				
	      } else {
		gamma_cr = ( crackStrainMinus.at(1) * residPlus.at(1) -  crackStrainPlus.at(1) * residMinus.at(1) ) / ( residPlus.at(1) - residMinus.at(1) );
	      }
	    } else { // different sign of the cracking strain
	      gamma_cr = ( crackStrainMinus.at(1) + crackStrainPlus.at(1) ) / 2.;
	    }
	  }
	  
	  
	  status->setTempCrackStrain(i, gamma_cr);
	  this->updateCrackStatus(gp);
	  
	  switch ( mMode ) {
	  case  _PlaneStress:
	    
	    maxTau = this->maxShearStress(gp, tStep, 6);
	    tau_cr = computeTotalD2Modulus(gp, tStep, 6) * gamma_cr;
	    break;
	    
	  case _PlaneStrain: // check
	    
	    maxTau = this->maxShearStress(gp, tStep, 6);
	    tau_cr = computeTotalD2Modulus(gp, tStep, 6) * gamma_cr;
	    break;
	    
	  case _3dMat:
	    
	    maxTau = this->maxShearStress(gp, tStep, i);
	    tau_cr = computeTotalD2Modulus(gp, tStep, i) * gamma_cr;
	    break;
	    
	  default:
	    OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
	  }
	  
	  if( fabs(tau_cr) > maxTau) {

	    if ( ( !plus ) || ( !minus ) ) {
	    
	    // d_eps_el = (tau_max - tau_n-1)/G
	    d_gamma_el = ( sgn(tau_cr)* maxTau - stressVector_l.at(i) ) / G;

	    // eps_cr_n = eps_cr_n-1 + d_eps - d_eps_el
	    gamma_cr = status->giveCrackStrain(i) + strainIncrement_l.at(i) - d_gamma_el;

	    // maxTau can be dependent on the actual cracking strain
	    // (such as in the case of FRC)

	    } else {
	      d_gamma_el =  strainIncrement_l.at(i) - (gamma_cr - status->giveCrackStrain(i) );
	    }

	    status->setTempCrackStrain(i, gamma_cr);
	    this->updateCrackStatus(gp);

	    // reevaluate TauMax
	    switch ( mMode ) {
	    case  _PlaneStress:
	      
	      maxTau = this->maxShearStress(gp, tStep, 6);
	      break;
	      
	    case _PlaneStrain: // check
	    
	      maxTau = this->maxShearStress(gp, tStep, 6);
	      break;
	      
	    case _3dMat:
	      
	      maxTau = this->maxShearStress(gp, tStep, i);
	      break;
	      
	    default:
	      OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
	    }

   	    tau_el =  stressVector_l.at(i) + d_gamma_el * G;

	    tau_cr = min( maxTau, DII * fabs(gamma_cr) );
	    tau_cr *= sgn(gamma_cr);
	    
	  } else { // tau is smaller than tau_max

	    // old shear stress
	    tau_el = stressVector_l.at(i);
	    // shear stress increment
	    tau_el += G * ( strainIncrement_l.at(i) - gamma_cr + status->giveCrackStrain(i) );
	  }
	  
	  // compute residuum = sig_el - sig_cr
	  d_tau = tau_el - tau_cr;

	  if  ( fabs( d_tau ) < fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) {
	    break;
	    
	  } else {
	    
	    if (d_tau < 0.) {
	      minus = true;

	      residMinus.at(1) = d_tau;
	      crackStrainMinus.at(1) = gamma_cr;
	      
	    } else {
	      plus = true;

	      residPlus.at(1) = d_tau;
	      crackStrainPlus.at(1) = gamma_cr;
	    }

	    // compare signs of the error
	    if (sgn(d_tau_old) == sgn(d_tau) ) { // same sign of the error as before
	      illinoisFlag = true;
	      
	    } else {
	      illinoisFlag = false;
	    }
	    
	    d_tau_old = d_tau;
	   
	  }

#if DEBUG		  		    
	  if ( iterS > iterLimitShear*0.3) {
	    OOFEM_WARNING("Fixed crack model: Local equilibrium in shear not converging!, max. stress error %f, index %d, iter %d, element %d\n", fabs(d_tau), i, iterS, gp->giveElement()->giveNumber() );
	  }
#endif

	  
	  if ( iterS == iterLimitShear ) {
	    OOFEM_WARNING("Fixed crack model: Local equilibrium in shear direction(s) not reached!, max. stress error %f, element %d\n", fabs(d_tau), gp->giveElement()->giveNumber() );
	  }
	  
	} // i-th shear component equlibrium loop
	
      } // shear components loop

      // compute final stress residuum
      crackStrain = status->giveTempCrackStrainVector();
      
      // el
      elasticStrainIncrement_l = strainIncrement_l;
      elasticStrainIncrement_l.subtract(crackStrain);
      elasticStrainIncrement_l.add(status->giveCrackStrainVector() );

      sigmaElast_l.beProductOf(De, elasticStrainIncrement_l);
      sigmaElast_l.add(stressVector_l);
      
      this->giveTotalLocalCrackedStiffnessMatrix(Dcr, SecantStiffness, gp, tStep);

      sigmaCrack_l.resize( crackStrain.giveSize() );
      
      for ( int i = 1; i <= crackStrain.giveSize(); i++ ) {
	sigmaCrack_l.at(i) = Dcr.at(i, i) * crackStrain.at(i);

	if (this->isThisShearComponent(gp, i) ) {

	  switch ( mMode ) {
	  case  _PlaneStress:
	    maxTau = this->maxShearStress(gp, tStep, 6);
	    break;
	    
	  case _PlaneStrain:
	    maxTau = this->maxShearStress(gp, tStep, 6);
	    break;
	    
	  case _3dMat:
	    maxTau = this->maxShearStress(gp, tStep, i);
	    break;
	    
	  default:
	    OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
	  }
	  

	  if ( fabs( sigmaCrack_l.at(i) ) > maxTau ) {
	    sigmaCrack_l.at(i) = sgn(sigmaCrack_l.at(i)) * maxTau;
	  }
	}
      }

      // residuum
      sigmaResid = sigmaElast_l;
      sigmaResid.subtract(sigmaCrack_l);
      
      maxErr = 0.;
      for ( int i = 1; i <= sigmaResid.giveSize(); i++ ) {
	if ( fabs( sigmaResid.at(i) ) > maxErr ) {
	  maxErr = fabs( sigmaResid.at(i) );
	}
      }

      if  ( maxErr < fcm_TOLERANCE * this->giveTensileStrength(gp, tStep) ) {
	break;
      }

      if ( iterG == iterLimitGlobal ) {
	OOFEM_WARNING("Fixed crack model: Local equilibrium not reached!, max. stress error %f, element %d\n", maxErr, gp->giveElement()->giveNumber() );
      }

      if ( iterG > iterLimitGlobal/3 ) {
	OOFEM_WARNING("Slowly converning (normal vs. shear interaction)!,  iter %d/ iter limit %d,  max. stress error %f, element %d\n", iterG, iterLimitGlobal, maxErr, gp->giveElement()->giveNumber() );
      }      

    } // global loop


    // correct non-zeros, statuses, etc...
    for ( int i = 1; i <= nMaxCr; i++ ) {
#if DEBUG
        // check for NaNs
        if ( crackStrain.at(i) != crackStrain.at(i) ) {
            OOFEM_ERROR( "crack strain is NaN: %f", crackStrain.at(i) );
        }
#endif

        if ( ( status->giveTempCrackStatus(i) == pscm_NONE ) && ( crackStrain.at(i) <= 0. ) ) {
            status->setTempCrackStrain(i, 0.);
        } else if ( ( status->giveTempCrackStatus(i) == pscm_JUST_INIT ) && ( crackStrain.at(i) <= fcm_SMALL_STRAIN ) ) {
            if ( i + 1 > status->giveMaxNumberOfCracks(gp) ) {
                status->setTempCrackStatus(i, pscm_NONE);
                status->setTempCrackStrain(i, 0.);
                if ( i == 1 ) {
                    principalDirs.zero();
                    status->setCrackDirs(principalDirs);
                }

                // can't change crack if the subsequent crack exists
            } else if ( status->giveTempCrackStatus(i + 1) != pscm_NONE ) {} else {
                status->setTempCrackStatus(i, pscm_NONE);
                status->setTempCrackStrain(i, 0.);
                if ( i == 1 ) {
                    principalDirs.zero();
                    status->setCrackDirs(principalDirs);
                }
            }
        } else if ( crackStrain.at(i) < 0. ) {
            status->setTempCrackStrain(i, 0.);
        }
    }

    // at first crack initiation there should not be any other cracking strains
    if ( ( status->giveTempCrackStatus(1) == pscm_JUST_INIT ) && ( status->giveTempCrackStatus(2) == pscm_NONE ) ) {
        status->setTempCrackStrain(2, 0.);
        status->setTempCrackStrain(3, 0.);
    }

    // convert from local to global coordinates  
    answer.beProductOf(sigmaL2G, sigmaElast_l);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

    return;
}
  

void
FCMMaterial :: initializeCrack(GaussPoint *gp, TimeStep *tStep, FloatMatrix &base, int nCrack)
{
    MaterialMode mMode = gp->giveMaterialMode();
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    FloatArray crackVector(3);
    crackVector.zero();

    FloatMatrix epsG2L, sigmaG2L, epsL2G, sigmaL2G;

    int nMaxCracks;
    double Le;

    for ( int i = 1; i <= base.giveNumberOfRows(); i++ ) {
        crackVector.at(i) = base.at(i, nCrack);
    }

    Le = this->giveCharacteristicElementLength(gp, crackVector);
    status->setCharLength(nCrack, Le);

    this->checkSnapBack(gp, tStep, nCrack);

    status->setTempCrackStatus(nCrack, pscm_JUST_INIT);
    //status->setTempCrackStatus(nCrack, pscm_SOFTENING);

    // change the base crack vector and transformation matrices
    // when the initiated crack is not the last possible one

    // return for the second crack in plane stress - base is set and transformation matrices are given by the first crack
    if ( ( mMode == _PlaneStress ) && ( nCrack == 2 ) ) {
        return;
    }

    nMaxCracks = status->giveMaxNumberOfCracks(gp);

    if ( nCrack <= nMaxCracks ) {
        // make sure that the matrix of base vectors is right-handed
        if ( nMaxCracks == 3 ) {
            base.at(1, 3) = base.at(2, 1) * base.at(3, 2) - base.at(3, 1) * base.at(2, 2);
            base.at(2, 3) = base.at(3, 1) * base.at(1, 2) - base.at(1, 1) * base.at(3, 2);
            base.at(3, 3) = base.at(1, 1) * base.at(2, 2) - base.at(2, 1) * base.at(1, 2);
        }

        //    status->setTempCrackDirs(base);
        status->setCrackDirs(base);

        if ( mMode == _PlaneStress ) {
            sigmaG2L = this->givePlaneStressVectorTranformationMtrx(base, 0);
            epsG2L = this->give2DStrainVectorTranformationMtrx(base, 0);

            sigmaL2G = this->givePlaneStressVectorTranformationMtrx(base, 1);
            epsL2G = this->give2DStrainVectorTranformationMtrx(base, 1);
        } else {
            sigmaG2L = this->giveStressVectorTranformationMtrx(base, 0);
            epsG2L = this->giveStrainVectorTranformationMtrx(base, 0);

            sigmaL2G = this->giveStressVectorTranformationMtrx(base, 1);
            epsL2G = this->giveStrainVectorTranformationMtrx(base, 1);
        }

        status->setG2LStressVectorTransformationMtrx(sigmaG2L);
        status->setG2LStrainVectorTransformationMtrx(epsG2L);

        status->setL2GStressVectorTransformationMtrx(sigmaL2G);
        status->setL2GStrainVectorTransformationMtrx(epsL2G);
    }
}


bool
FCMMaterial :: isIntact(GaussPoint *gp, int icrack) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    if ( icrack >= 4 ) {
        OOFEM_ERROR("Unexpected crack number");
    }

    if (  status->giveCrackStatus(icrack) != pscm_NONE ) {
      return false;
    } else if ( ( status->giveTempCrackStatus(icrack) != pscm_NONE ) && ( status->giveTempCrackStatus(icrack) != pscm_CLOSED ) ) {
      return false;
    } else {
      return true;
    }
}


bool
FCMMaterial :: isIntactForShear(GaussPoint *gp, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    int normal_1, normal_2;

    if ( i == 4 ) { // y-z
        normal_1 = 2;
        normal_2 = 3;
    } else if ( i == 5 ) { // x-z
        normal_1 = 1;
        normal_2 = 3;
    } else if ( i == 6 ) { // x-y
        normal_1 = 1;
        normal_2 = 2;
    } else {
        OOFEM_ERROR("Unexpected number for shear stress (must be either 4, 5 or 6).");
        normal_1 = normal_2 = 0;
    }

    // any crack had existed in previous steps - newly inserted condition
    if  ( (status->giveCrackStatus(normal_1) != pscm_NONE ) || (status->giveCrackStatus(normal_2) != pscm_NONE ) )  { 
        return false;

	//    } else if ( ( status->giveTempCrackStatus(normal_1) != pscm_NONE ) && ( status->giveTempCrackStatus(normal_1) != pscm_CLOSED ) ) {
	// return false;
	//    } else if ( ( status->giveTempCrackStatus(normal_2) != pscm_NONE ) && ( status->giveTempCrackStatus(normal_2) != pscm_CLOSED ) ) {
	// return false;
    } else {
        return true;
    }
}



bool
FCMMaterial :: isThisShearComponent(GaussPoint *gp, int component) {

    MaterialMode mMode = gp->giveMaterialMode();

    switch ( mMode ) {

    case  _PlaneStress:
      if (component == 3) {
	return true;
      }
      break;

    case  _PlaneStrain:
      if (component == 4) {
	return true;
      }
      break;

    case  _3dMat:

      if (component >=4 ) {
	return true;
      }
      break;

    default:
      OOFEM_ERROR("Unsupported material mode");
    }

    return false;
}


  

double
FCMMaterial :: computeNormalCrackOpening(GaussPoint *gp, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double crackOpening, N;

    if ( i > status->giveNumberOfTempCracks() ) {
        crackOpening = 0;
    } else {
        crackOpening = max(status->giveCharLength(i) * status->giveTempCrackStrain(i), 0.);
        N = this->giveNumberOfCracksInDirection(gp, i);
        crackOpening /= N;
    }

    return crackOpening;
}



double
FCMMaterial :: computeMaxNormalCrackOpening(GaussPoint *gp, TimeStep *tStep, int i) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double crackOpening, N;

    if ( i > status->giveNumberOfTempCracks() ) {
        crackOpening = 0;
    } else {
        crackOpening = max(status->giveCharLength(i) * status->giveMaxCrackStrain(i), 0.);
        N = this->giveNumberOfCracksInDirection(gp, i);
        crackOpening /= N;
    }

    return crackOpening;
}



double
FCMMaterial :: computeShearSlipOnCrack(GaussPoint *gp, TimeStep *tStep, int icrack) {
    MaterialMode mMode = gp->giveMaterialMode();

    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    // total shear slip
    double slip = 0;

    // number of cracks in the direction of iCrack
    int nCracks = this->giveNumberOfCracksInDirection(gp, icrack);

    // shear directions on icrack plane
    int dir_j, dir_k;

    // cracking shear strains
    double gamma_cr_ij, gamma_cr_ik;

    // factor for redistribution of shear crack strain
    double factor_ij, factor_ik;

    double u_ij, u_ik;

    if ( status->giveTempCrackStatus(icrack) == pscm_NONE ) { // crack not initiated
        return slip;
    }

    if ( icrack > 3 ) {
        OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
    }


    if ( ( mMode == _PlaneStress ) || ( mMode == _PlaneStrain ) ) {
        // get shear strain
        if ( mMode == _PlaneStress ) {
            gamma_cr_ij = status->giveTempCrackStrain(3);
        } else {
            gamma_cr_ij = status->giveTempCrackStrain(6);
        }

        factor_ij = 1.;

        // determine if it necessary to split the cracking strain
        // both cracks exist in 2D

        if ( ( icrack == 2 ) || ( status->giveTempCrackStatus(2) != pscm_NONE ) ) {
            if ( icrack == 1 ) {
                dir_j = 2;
            } else { // icrack == 2
                dir_j = 1;
            }


            // well this is quite unfortunate. The problem is that the shear slip should be redistributed according to the D2 moduli for the individual crack. In reality this is a big problem for the FRC-FCM problem beacuse it would have led to recursive calling (D2 modulus evaluates damage and damage is computed from shear slip etc.

            factor_ij = this->computeShearStiffnessRedistributionFactor(gp, tStep, icrack, dir_j);
        }

        slip = factor_ij * fabs(gamma_cr_ij) * status->giveCharLength(icrack) / nCracks;
    } else if ( mMode == _3dMat ) {
        if ( icrack == 1 ) {
            dir_j = 2;
            dir_k = 3;

            gamma_cr_ij = status->giveTempCrackStrain(6);
            gamma_cr_ik = status->giveTempCrackStrain(5);
        } else if ( icrack == 2 ) {
            dir_j = 1;
            dir_k = 3;

            gamma_cr_ij = status->giveTempCrackStrain(6);
            gamma_cr_ik = status->giveTempCrackStrain(4);
        } else { // icrack == 3
            dir_j = 1;
            dir_k = 2;

            gamma_cr_ij = status->giveTempCrackStrain(5);
            gamma_cr_ik = status->giveTempCrackStrain(4);
        }


        factor_ij = 1.;
        factor_ik = 1.;


        if ( ( status->giveTempCrackStatus(dir_j) != pscm_NONE ) || ( status->giveTempCrackStatus(dir_k) != pscm_NONE ) ) {
            if ( status->giveTempCrackStatus(dir_j) != pscm_NONE ) {
	      factor_ij = this->computeShearStiffnessRedistributionFactor(gp, tStep, icrack, dir_j);
            }

            if ( status->giveTempCrackStatus(dir_k) != pscm_NONE ) {
                factor_ik = this->computeShearStiffnessRedistributionFactor(gp, tStep, icrack, dir_k);
            }
        }

        u_ij = factor_ij * fabs(gamma_cr_ij) * status->giveCharLength(icrack) / nCracks;
        u_ik = factor_ik * fabs(gamma_cr_ik) * status->giveCharLength(icrack) / nCracks;

        slip = sqrt( pow(u_ij, 2) + pow(u_ik, 2) );
    } else {
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
    }

    return slip;
}


bool
FCMMaterial :: isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress) {
  if ( trialStress > this->giveTensileStrength(gp, tStep) ) {
        return true;
    } else {
        return false;
    }
}




double
FCMMaterial :: computeShearStiffnessRedistributionFactor(GaussPoint *gp, TimeStep *tStep, int ithCrackPlane, int jthCrackDirection) {
    double factor_ij;
    double D2_i, D2_j;

    D2_i = this->computeD2ModulusForCrack(gp, tStep, ithCrackPlane);
    D2_j = this->computeD2ModulusForCrack(gp, tStep, jthCrackDirection);

    factor_ij = D2_j / ( D2_i + D2_j );

    return factor_ij;
}


bool
FCMMaterial :: checkStrengthCriterion(FloatMatrix &newBase, const FloatArray &globalStress, GaussPoint *gp, TimeStep *tStep, int nCrack) {
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    double sigX, sigY, tau, sig2;
    FloatArray trialStress, planeStress, princStress, crackingStrain, shearStrains, newShearStrains;
    FloatMatrix sigmaG2L, princCrackDir, oldBase, princCrackDirExt;

    sigmaG2L = status->giveG2LStressVectorTransformationMtrx();

    // rotate stress to local coordinates
    //    trialStress = globalStress;
    //    trialStress.rotatedWith(sigmaG2L, 'n');
    trialStress.beProductOf(sigmaG2L, globalStress);

    if ( status->giveMaxNumberOfCracks(gp) < 3 ) { // plane stress
        // test material strength but keep crack coordinates
        if ( this->isStrengthExceeded( newBase, gp, tStep, nCrack, trialStress.at(nCrack) ) ) {
            return false;
        } else {
            return true;
        }
    } else if ( nCrack == 3 ) { // 3D but is the last crack
        if ( this->isStrengthExceeded( newBase, gp, tStep, nCrack, trialStress.at(nCrack) ) ) {
            return false;
        } else {
            return true;
        }
    } else { // second crack in 3D case
        // if the stress-strength criterion is violated
        // the crack coordinates have to be changed
        // in-plane shear strain remains unchanged

        // compute second principal stress cheap - if passed, compute eigenvectors and eigenvalues
        sigX = trialStress.at(2); //1st normal stress in crack plane
        sigY = trialStress.at(3); //2nd normal stress in crack plane
        tau = trialStress.at(4); // shear stress in crack plane

        sig2 = ( sigX + sigY ) / 2. + sqrt( ( sigX - sigY ) * ( sigX - sigY ) / 4. + tau * tau );

        if ( this->isStrengthExceeded(newBase, gp, tStep, 2, sig2) ) {
            oldBase = newBase;
            planeStress.resize(3);
            planeStress.zero();

            planeStress.at(1) = trialStress.at(2); //1st normal stress in crack plane
            planeStress.at(2) = trialStress.at(3); //2nd normal stress in crack plane
            planeStress.at(3) = trialStress.at(4); //shear in crack plane


            // compute principal stresses on the already existing crack plane
            this->computePrincipalValDir(princStress, princCrackDir, planeStress, principal_stress);

            // right-hand orientation
            princCrackDir.at(1, 2) = -1. * princCrackDir.at(2, 1);
            princCrackDir.at(2, 2) = princCrackDir.at(1, 1);

            // establish vector of crack directions on the crack plane
            // the first vector is the normal direction, [1, 0, 0]

            princCrackDirExt.resize(3, 3);
            princCrackDirExt.zero();
            princCrackDirExt.at(1, 1) = 1.;

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    princCrackDirExt.at(i + 1, j + 1) = princCrackDir.at(i, j);
                }
            }

            // new crack base vector in global coordinates
            newBase.beProductOf(oldBase, princCrackDirExt);

            /*
             * // algorithm which gives the same results - rotation of axes in space
             * double theta;
             *
             * if ( princCrackDir.at(1,1) > 0. ) {
             * theta = asin( princCrackDir.at(2,1) );
             * } else {
             * theta = M_PI - asin ( princCrackDir.at(2,1) );
             * }
             *
             * if (theta < 0.) {
             * theta += 2 * M_PI;
             * }
             *
             *
             * FloatArray baseVector2;
             * FloatArray baseVector3;
             *
             * FloatArray newBaseVector2;
             * FloatArray newBaseVector3;
             *
             * baseVector2.resize(3);
             * baseVector3.resize(3);
             *
             * newBaseVector2.resize(3);
             * newBaseVector3.resize(3);
             *
             * for (int i = 1; i <= 3; i++) {
             * baseVector2.at(i) = oldBase.at(i,2);
             * baseVector3.at(i) = oldBase.at(i,3);
             * }
             *
             * FloatMatrix K;
             * FloatMatrix R;
             * FloatMatrix help;
             *
             * K.resize(3,3);
             * K.zero();
             *
             * K.at(1,2) = -1. * oldBase.at(3,1);
             * K.at(1,3) = oldBase.at(2,1);
             * K.at(2,1) = oldBase.at(3,1);
             * K.at(2,3) = -1. * oldBase.at(1,1);
             * K.at(3,1) = -1. * oldBase.at(2,1);
             * K.at(3,2) = oldBase.at(1,1);
             *
             * R.resize(3,3);
             * R.beUnitMatrix();
             *
             * help = K;
             * help.times( sin(theta) );
             *
             * R.add(help);
             *
             * help.beProductOf (K, K);
             * help.times( 1.-cos(theta) );
             *
             * R.add(help);
             *
             * newBaseVector2.beProductOf(R, baseVector2);
             * newBaseVector3.beProductOf(R, baseVector3);
             */


            /*
             * //	the same way - R evaluated symbolically in Matlab
             * //	c = cos, s = sin
             * //	k = axis of rotation
             * R.at(1,1) = 1. + (c-1.) * k2*k2 + (c-1.) * k3*k3;
             * R.at(1,2) = -k3*s - k1*k2 * (c-1.);
             * R.at(1,3) =  k2*s - k1*k3 * (c-1.);
             *
             * R.at(2,1) = k3*s - k1*k2 * (c-1.);
             * R.at(2,2) = 1. + (c-1.) * k1*k1 + (c-1.) * k3*k3;
             * R.at(2,3) = -k1*s - k2*k3 * (c-1.);
             *
             * R.at(3,1) = -k2*s - k1*k3 * (c-1.);
             * R.at(3,2) = k1*s - k2*k3 * (c-1.);
             * R.at(3,3) = 1. + (c-1.) * k1*k1 + (c-1.) * k2*k2;
             */


            // modify shear strains

            crackingStrain = status->giveCrackStrainVector();

            shearStrains.resize(2);
            shearStrains.at(1) = crackingStrain.at(6); //aka x-y
            shearStrains.at(2) = crackingStrain.at(5); //aka x-z

            // gamma_xy_l   [ c   s ]    gamma_xy_g
            //            = [       ] *
            // gamma_xz_l   [ -s  c ]    gamma_xz_g

            newShearStrains.beTProductOf(princCrackDir, shearStrains);

            crackingStrain.at(5) =  newShearStrains.at(2); // aka x-z
            crackingStrain.at(6) =  newShearStrains.at(1); // aka x-y

            status->setCrackStrainVector(crackingStrain);
            status->setTempCrackStrainVector(crackingStrain);

            // modify maximum shear strains
            // get max shear strains in the original configuration
            shearStrains.at(1) = status->giveMaxCrackStrain(6); //aka x-y
            shearStrains.at(2) = status->giveMaxCrackStrain(5); //aka x-z

            // gamma_xy_l   [ c   s ]    gamma_xy_g
            //            = [       ] *
            // gamma_xz_l   [ -s  c ]    gamma_xz_g

            newShearStrains.beTProductOf(princCrackDir, shearStrains);

            status->setMaxCrackStrain( 5, max( newShearStrains.at(2), crackingStrain.at(5) ) ); // aka x-z
            status->setMaxCrackStrain( 6, max( newShearStrains.at(1), crackingStrain.at(6) ) ); // aka x-y

            status->setTempMaxCrackStrain( 5, max( newShearStrains.at(2), crackingStrain.at(5) ) ); // aka x-z
            status->setTempMaxCrackStrain( 6, max( newShearStrains.at(1), crackingStrain.at(6) ) ); // aka x-y

            return false;
        } else { // strength not reached
            return true;
        }
    } // number of crack in given stress-state
}


double
FCMMaterial :: giveCharacteristicElementLength(GaussPoint *gp, const FloatArray &crackPlaneNormal)
{
    return gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
}


void
FCMMaterial :: updateCrackStatus(GaussPoint *gp)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );

    int nMaxCracks = status->giveMaxNumberOfCracks(gp);

    FloatArray crackStrain, maxCrackStrain;
    crackStrain = status->giveTempCrackStrainVector();
    maxCrackStrain = status->giveMaxCrackStrainVector();


    // loop over NORMAL components of the crack strain vector
    // statuses are only for normal components
    for ( int i = 1; i <= nMaxCracks; i++ ) {
        // changes in status only in previously cracked (or at least initiated in this step) material
        if ( status->giveTempCrackStatus(i) != pscm_NONE ) {
            // temp crack bigger or equal than max strain so far
            if ( crackStrain.at(i) >= maxCrackStrain.at(i) ) {
                status->setTempMaxCrackStrain( i, crackStrain.at(i) );

                if ( status->giveCrackStatus(i) != pscm_NONE ) { //already existing crack in prev step
                    status->setTempCrackStatus(i, pscm_SOFTENING);
                } else if ( ( fabs( crackStrain.at(i) ) <= 1.e-18 ) && ( status->giveCrackStatus(i) == pscm_NONE ) ) { // just initiated crack with zero strain
                    if ( i == nMaxCracks ) {
                        status->setTempCrackStatus(i, pscm_NONE);
                    } else if ( nMaxCracks > i ) { // curent crack is not the last possible one
                        if ( status->giveTempCrackStatus(i + 1) == pscm_NONE ) { // the next crack is not initiated
                            status->setTempCrackStatus(i, pscm_NONE); // reset the current one to NONE
                        } else {
                            status->setTempCrackStatus(i, pscm_JUST_INIT); // else keep as just initiated
                        }
                    } else {
                        status->setTempCrackStatus(i, pscm_JUST_INIT);
                    }
                } else {
                    status->setTempCrackStatus(i, pscm_JUST_INIT);
                }

                // closed previously
            } else if ( crackStrain.at(i) <= 0. ) {

	      if (status->giveCrackStatus(i) != pscm_NONE ) {
                status->setTempCrackStatus(i, pscm_CLOSED);
	      } else {
                status->setTempCrackStatus(i, pscm_NONE);
		status->setTempMaxCrackStrain(i, 0.);
	      }


		//                if ( status->giveMaxCrackStrain(i) == 0. ) { //not existing crack in previous step that wants to close?


		
              //      status->setTempMaxCrackStrain(i, 0.);

		    //		    status->setTempCrackStatus(i, pscm_NONE);
		    
                    //	  status->setTempCrackStatus(i, pscm_NONE);
                    // CLOSED status should have similar behavior and it is safer. The status can be changed within one iteration loop on the local scale
                    // the crack can be changed from CLOSED to NONE during overall updating
                    //	  status->setTempCrackStatus(i, pscm_CLOSED);
	      //}

                // unloading--reloading
            } else if ( crackStrain.at(i) < maxCrackStrain.at(i) ) {
                status->setTempCrackStatus(i, pscm_UNLO_RELO);
            } else {
	      OOFEM_ERROR("Unexpected value of %d-th cracking strain %f ", i, crackStrain.at(i) );
            }
        }
    }

    // SHEAR components
    for ( int i = nMaxCracks + 1; i <=  crackStrain.giveSize(); i++ ) {
        if ( fabs( crackStrain.at(i) ) > maxCrackStrain.at(i) ) {
            status->setTempMaxCrackStrain( i, fabs( crackStrain.at(i) ) );
        }
    }
}





void
FCMMaterial :: giveTotalLocalCrackedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    int dim, j;
    FloatMatrix Dcr;
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );

    dim = indx.giveSize();

    Dcr.resize(dim, dim);

    for ( int i = 1; i <= dim; i++ ) {
        j = indx.at(i);

        if ( j <= 3 ) {
	  Dcr.at(i, i) = this->giveCrackingModulus(rMode, gp, tStep, j);
        } else {
	  Dcr.at(i, i) = this->computeTotalD2Modulus(gp, tStep, j);
        }
    }

    answer = Dcr;
}


  
void
FCMMaterial :: giveNormalLocalCrackedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{

    FloatMatrix Dcr;
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    int nMaxCr = status->giveMaxNumberOfCracks(gp);

    Dcr.resize(nMaxCr, nMaxCr);

    for ( int i = 1; i <= nMaxCr; i++ ) {
      Dcr.at(i, i) = this->giveCrackingModulus(rMode, gp, tStep, i);
    }

    answer = Dcr;
}  



void
FCMMaterial :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    //    StructuralMaterial *lMat = static_cast< StructuralMaterial * >( this->giveLinearElasticMaterial() );

    MaterialMode mMode = gp->giveMaterialMode();

    int numberOfActiveCracks, nMaxCracks;
    double overallElasticStiffness;
    double G = this->computeOverallElasticShearModulus(gp, tStep);

    FloatMatrix D, De, C, stiffnessL2G;

    FloatMatrix DeHelp, Dcr, DcrHelp, DcrHelp2;

    numberOfActiveCracks = status->giveNumberOfTempCracks();

    // ELASTIC MATRIX
    linearElasticMaterial.giveStiffnessMatrix(De, rMode, gp, tStep);
    //    lMat->giveStiffnessMatrix(De, rMode, gp, tStep);
    
    overallElasticStiffness = this->computeOverallElasticStiffness(gp, tStep);
    if ( overallElasticStiffness != ( this->give('E', gp) ) ) {
      De.times( overallElasticStiffness / ( this->give('E', gp) ) );
    }

    if ( ( rMode == ElasticStiffness ) || ( numberOfActiveCracks == 0 ) ) {	
        answer = De;
        return;
    }

    // SECANT OR TANGENT MATRIX - first in local direction
    nMaxCracks = status->giveMaxNumberOfCracks(gp);
    De.resizeWithData(nMaxCracks,nMaxCracks);


    double E_crack;
    
    if ( rMode == SecantStiffness ) {

      C.beInverseOf(De);
      
      for ( int i = 1; i <= numberOfActiveCracks; i++ ) {

	E_crack =  this->giveCrackingModulus(rMode, gp, tStep, i);

	C.at(i, i) += 1. / E_crack;

      }
	

      D.beInverseOf(C);
      
    } else if ( rMode == TangentStiffness ) {

      //    if ( ( rMode == SecantStiffness ) || ( rMode == TangentStiffness ) ) {
      
      // shrink to square matrix number of cracks x number of cracks
      DeHelp = De;
      DeHelp.resizeWithData(numberOfActiveCracks, numberOfActiveCracks);
      
      Dcr.resize(numberOfActiveCracks, numberOfActiveCracks);
      Dcr.zero();
      
      for ( int i = 1; i <= numberOfActiveCracks; i++ ) {
	if ( status->giveTempCrackStrain(i) > fcm_THRESHOLD_CRACK_STRAIN ) {
	  Dcr.at(i, i) = this->giveCrackingModulus(rMode, gp, tStep, i);
	} else {
	  Dcr.at(i, i) = fcm_BIGNUMBER * overallElasticStiffness;
	}
      }
      
      Dcr.add(DeHelp);
      C.beInverseOf(Dcr);
      C.resizeWithData(nMaxCracks, nMaxCracks);
      
      DcrHelp.beProductOf(De, C); // De (De + Dcr)^-1
      DcrHelp2.beProductOf(DcrHelp, De); // De (De + Dcr)^-1 De
      
      D.zero();
      D.resize(nMaxCracks, nMaxCracks);
      D.add(De);
      D.subtract(DcrHelp2); // De - De (De + Dcr)^-1 De
    }


    if ( this->normalCoeffNumer > 0. ) {
      // apply minimum stiffness
      double nu = this->givePoissonsRatio();
      
      for ( int i = 1; i <= nMaxCracks; i++ ) {
	for ( int j = 1; j <= nMaxCracks; j++ ) {
	  if (i == j) { // diagonal terms
	    D.at(i, j) = max( D.at(i,j), this->normalCoeffNumer * overallElasticStiffness );
	  } else {
	    D.at(i, j) = max( D.at(i,j), this->normalCoeffNumer * overallElasticStiffness * nu/(1-nu) );
	  }
	}
      }
    }




    
    // resize add shear moduli on diagonal
    switch ( mMode ) {
    case _PlaneStress:
        D.resizeWithData(3, 3);

	if ( max( status->giveCrackStrain(1), status->giveCrackStrain(2) ) > fcm_THRESHOLD_CRACK_STRAIN ) { 
	  D.at(3, 3) = this->computeEffectiveShearModulus(gp, tStep, 6);
	  D.at(3, 3) = max( D.at(3,3), this->shearCoeffNumer * G );
	} else {
	  D.at(3, 3) = G;
	}

	// debug:
	/*double help;
	help = this->shearCoeffNumer * G;
	if (help > D.at(3, 3) ) {
	  D.at(3, 3) = help;
	  }*/

        break;

    case _3dMat:
        D.resizeWithData(6, 6);
        for ( int i = 4; i <= 6; i++ ) {
	  if (i == 4) {
	    if ( max( status->giveCrackStrain(2), status->giveCrackStrain(3) ) > fcm_THRESHOLD_CRACK_STRAIN ) {
	      D.at(i, i) = this->computeEffectiveShearModulus(gp, tStep, i);
	      D.at(i, i) = max( D.at(i,i), this->shearCoeffNumer * G );
	    } else {
	      D.at(i, i) = G;
	    }
	    
	  } else if ( i == 5 ) {
	    if ( max( status->giveCrackStrain(1), status->giveCrackStrain(3) ) > fcm_THRESHOLD_CRACK_STRAIN ) {
	      D.at(i, i) = this->computeEffectiveShearModulus(gp, tStep, i);
	      D.at(i, i) = max( D.at(i,i), this->shearCoeffNumer * G );
	    } else {
	      D.at(i, i) = G;
	    }
	    
	  } else { // i == 6
	    if ( max( status->giveCrackStrain(1), status->giveCrackStrain(2) ) > fcm_THRESHOLD_CRACK_STRAIN ) {
	      D.at(i, i) = this->computeEffectiveShearModulus(gp, tStep, i);
	      D.at(i, i) = max( D.at(i,i), this->shearCoeffNumer * G );
	      } else {
	      D.at(i, i) = G;
	    }
	  }
	} // end for

        break;

    case _PlaneStrain:
        D.resizeWithData(6, 6);

	if ( max( status->giveCrackStrain(1), status->giveCrackStrain(2) ) > fcm_THRESHOLD_CRACK_STRAIN ) { 
	  D.at(6, 6) = this->computeEffectiveShearModulus(gp, tStep, 6);
	  D.at(6, 6) = max( D.at(6,6), this->shearCoeffNumer * G );
	} else {
	  D.at(6, 6) = G;
	}

        break;

    default:
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
    }

    // transform stiffnes to global c.s
    // to transform stiffness from LOCAL to GLOBAL coordinate system
    // the transformation matrix is the same as for strain transformation
    // from GLOBAL to LOCAL coordinate system
    stiffnessL2G = status->giveG2LStrainVectorTransformationMtrx();
    D.rotatedWith(stiffnessL2G, 'n');


    switch ( mMode ) {
    case _PlaneStress:
        answer = D;
        return;

        break;
    case _3dMat:
        answer = D;
        return;

        break;
    case _PlaneStrain:
        StructuralMaterial :: giveReducedSymMatrixForm( answer, D, gp->giveMaterialMode() );
        return;

        break;
    default:
        OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
        return;
    }
}



double
FCMMaterial :: computeTotalD2Modulus(GaussPoint *gp, TimeStep *tStep, int shearDirection)
{
    double D2 = 0.;
    int crackA, crackB;
    double D2_1, D2_2;

    if ( this->isIntactForShear(gp, shearDirection) ) {
      D2 = this->computeOverallElasticStiffness(gp, tStep) * fcm_BIGNUMBER;
    } else {
        if ( shearDirection == 4 ) {
            crackA = 2;
            crackB = 3;
        } else if ( shearDirection == 5 ) {
            crackA = 1;
            crackB = 3;
        }  else if ( shearDirection == 6 ) {
            crackA = 1;
            crackB = 2;
        } else {
            OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
            crackA = crackB = 0;
        }

        if ( ( this->isIntact(gp, crackA) ) || ( this->isIntact(gp, crackB) ) ) {
            if ( this->isIntact(gp, crackA) ) {
	      D2 = this->computeD2ModulusForCrack(gp, tStep, crackB);
            } else {
	      D2 = this->computeD2ModulusForCrack(gp, tStep, crackA);
            }
        } else {
	  D2_1 = this->computeD2ModulusForCrack(gp, tStep, crackA);
	  D2_2 = this->computeD2ModulusForCrack(gp, tStep, crackB);

            if ( multipleCrackShear ) {
                // serial coupling of stiffnesses 1/D = 1/D1 + 1/D2
                D2 = D2_1 * D2_2 / ( D2_1 + D2_2 );
            } else {
                D2 = min(D2_1, D2_2);
            }
        }

	// adjust stiffness according to the maximum allowable stress
	/*FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
	MaterialMode mMode = gp->giveMaterialMode();
	double D2_help, gamma_max;
	  switch ( mMode ) {
	  case  _PlaneStress:
	    
	    gamma_max = status->giveMaxCrackStrain(3);
	    break;
	    
	  case _PlaneStrain: // check
	    gamma_max = status->giveMaxCrackStrain(3);
	    break;
	    
	  case _3dMat:
	    gamma_max = status->giveMaxCrackStrain(shearDirection);
	    break;
	    
	  default:
	    OOFEM_ERROR( "Material mode %s not supported", __MaterialModeToString(mMode) );
	  }
	D2_help  = this->maxShearStress(gp, tStep, shearDirection) / gamma_max;
	D2 = min(D2, D2_help);
	*/
    }

    return D2;
}


double
FCMMaterial :: computeNumerD2Modulus(GaussPoint *gp, TimeStep *tStep, int shearDirection)
{
    double D2 = 0.;
    int crackA, crackB;
    double D2_1, D2_2;

    if ( this->isIntactForShear(gp, shearDirection) ) {
      D2 = this->computeOverallElasticStiffness(gp, tStep) * fcm_BIGNUMBER;
    } else {
        if ( shearDirection == 4 ) {
            crackA = 2;
            crackB = 3;
        } else if ( shearDirection == 5 ) {
            crackA = 1;
            crackB = 3;
        }  else if ( shearDirection == 6 ) {
            crackA = 1;
            crackB = 2;
        } else {
            OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
            crackA = crackB = 0;
        }

        if ( ( this->isIntact(gp, crackA) ) || ( this->isIntact(gp, crackB) ) ) {
            if ( this->isIntact(gp, crackA) ) {
	      D2 = this->computeNumerD2ModulusForCrack(gp, tStep, crackB);
            } else {
	      D2 = this->computeNumerD2ModulusForCrack(gp, tStep, crackA);
            }
        } else {
	  D2_1 = this->computeNumerD2ModulusForCrack(gp, tStep, crackA);
	  D2_2 = this->computeNumerD2ModulusForCrack(gp, tStep, crackB);

            if ( multipleCrackShear ) {
                // serial coupling of stiffnesses 1/D = 1/D1 + 1/D2
                D2 = D2_1 * D2_2 / ( D2_1 + D2_2 );
            } else {
                D2 = min(D2_1, D2_2);
            }
        }
    }

    return D2;
}

  

void
FCMMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial.initializeFrom(ir);

    this->nAllowedCracks = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, nAllowedCracks, _IFT_FCM_nAllowedCracks);


    this->crackSpacing = -1.;
    if ( ir.hasField(_IFT_FCM_crackSpacing) ) {
        IR_GIVE_FIELD(ir, crackSpacing, _IFT_FCM_crackSpacing);
    }

    this->multipleCrackShear = false;
    if ( ir.hasField(_IFT_FCM_multipleCrackShear) ) {
        this->multipleCrackShear = true;
    }

    int ecsm = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, ecsm, _IFT_FCM_ecsm);
    switch ( ecsm ) {
    case 1: ecsMethod = ECSM_SquareRootOfArea;
        break;
    case 2: ecsMethod = ECSM_ProjectionCentered;
        break;
    case 3: ecsMethod = ECSM_Oliver1;
        break;
    case 4: ecsMethod = ECSM_Oliver1modified;
        break;
    default: ecsMethod = ECSM_Projection;
    }

    this->shearCoeffNumer = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, shearCoeffNumer, _IFT_FCM_shearCoeffNumer);
    this->normalCoeffNumer = -1. * fcm_BIGNUMBER;
    IR_GIVE_OPTIONAL_FIELD(ir, normalCoeffNumer, _IFT_FCM_normalCoeffNumer);
}


double
FCMMaterial :: give(int aProperty, GaussPoint *gp) const
{
    return linearElasticMaterial.give(aProperty, gp);
}



double
FCMMaterial :: giveCrackSpacing(void)
{
    return crackSpacing;
}


double
FCMMaterial :: giveNumberOfCracksInDirection(GaussPoint *gp, int iCrack)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    double spacing, L, N;

    L = status->giveCharLength(iCrack);
    spacing = this->giveCrackSpacing();

    if ( ( spacing > L ) || ( spacing < 0. ) ) {
        N = 1;
    } else {
        N = ( L / spacing );
    }

    return N;
}


double
FCMMaterial :: giveNumberOfCracksForShearDirection(GaussPoint *gp, int i)
{
    double N;
    int dir_1, dir_2;

    if ( i == 4 ) {
        dir_1 = 2;
        dir_2 = 3;
    } else if ( i == 5 ) {
        dir_1 = 1;
        dir_2 = 3;
    }  else if ( i == 6 ) {
        dir_1 = 1;
        dir_2 = 2;
    } else {
        OOFEM_ERROR("Unexpected value of index i (4, 5, 6 permitted only)");
        dir_1 = dir_2 = 0;
    }

    N = max( this->giveNumberOfCracksInDirection(gp, dir_1), this->giveNumberOfCracksInDirection(gp, dir_2) );

    return N;
}


int
FCMMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FCMMaterialStatus *status = static_cast< FCMMaterialStatus * >( this->giveStatus(gp) );
    double width;
    int index;


    // first/dominant crack vector
    if ( type == IST_CrackVector ) {
        answer.resize(3);
        answer.zero();

        // MAX CRACK
        width = 0.;
        index = 1;
        for ( int i = 1; i <= status->giveNumberOfCracks(); i++ ) {
            if ( status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i)  > width ) {
                width = status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i);
                index = i;
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;
    } else if ( type == IST_2ndCrackVector ) {
        answer.resize(3);
        answer.zero();
        index = 2;

        if ( status->giveNumberOfCracks() <= 1 ) {
            return 1;
        } else if ( status->giveNumberOfCracks() == 2 ) {
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                index = 1;
            }
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 3;
                } else {
                    index = 1;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 3;
                } else {
                    index = 2;
                }
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;
    } else if ( type == IST_3rdCrackVector ) {
        answer.resize(3);
        answer.zero();
        index = 3;

        if ( status->giveNumberOfCracks() <= 2 ) {
            return 1;
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 1;
                } else {
                    index = 3;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 2;
                } else {
                    index = 3;
                }
            }
        }

        for ( int i = 1; i <= status->giveCrackDirs().giveNumberOfRows(); i++ ) {
            answer.at(i) = status->giveCrackDirs().at(i, index);
        }

        return 1;


        // width of a first/dominant crack
    } else if ( type == IST_CrackWidth ) {
        answer.resize(1);
        answer.zero();

        // MAX WIDTH
        width = 0.;
        for ( int i = 1; i <= status->giveNumberOfCracks(); i++ ) {
            width = max( width, status->giveCharLength(i) * status->giveCrackStrain(i) / this->giveNumberOfCracksInDirection(gp, i) );
        }
        answer.at(1) = width;
        return 1;
    } else if ( type == IST_2ndCrackWidth ) {
        answer.resize(1);
        answer.zero();
        index = 2;
        if ( status->giveNumberOfCracks() <= 1 ) {
            return 1;
        } else if ( status->giveNumberOfCracks() == 2 ) {
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2)  > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                index = 1;
            }
        } else { // 3 cracks
            if ( status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) { // #2 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(1) * status->giveCrackStrain(1) / this->giveNumberOfCracksInDirection(gp, 1) ) {
                    index = 3;
                } else {
                    index = 1;
                }
            } else { // #1 is biggest
                if ( status->giveCharLength(3) * status->giveCrackStrain(3) / this->giveNumberOfCracksInDirection(gp, 3) > status->giveCharLength(2) * status->giveCrackStrain(2) / this->giveNumberOfCracksInDirection(gp, 2) ) {
                    index = 3;
                } else {
                    index = 2;
                }
            }
        }

        width = status->giveCharLength(index) * status->giveCrackStrain(index) / this->giveNumberOfCracksInDirection(gp, index);

        answer.at(1) = width;
        return 1;
    } else if ( type == IST_CrackDirs ) {
        const FloatMatrix &help = status->giveCrackDirs();
        answer.resize(9);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) = help.at(1, i);
            answer.at(i + 3) = help.at(2, i);
            answer.at(i + 6) = help.at(3, i);
        }

        return 1;

    } else if ( type == IST_CrackStatuses ) {
        answer.resize(3);
	answer.zero();
        for ( int i = 1; i <= status->giveNumberOfCracks(); i++ ) {
            answer.at(i) = status->giveCrackStatus(i);
        }
        return 1;

    } else if ( type == IST_CrackStatusesTemp ) {
        answer.resize(3);
	answer.zero();
        for ( int i = 1; i <= status->giveNumberOfTempCracks(); i++ ) {
            answer.at(i) = status->giveTempCrackStatus(i);
        }
        return 1;
	
    } else if ( type == IST_CrackStrainTensor ) {
      //        FloatArray crackStrain = status->giveCrackStrainVector();
        FloatArray crackStrain;
        FloatMatrix epsL2G = status->giveL2GStrainVectorTransformationMtrx();
        // from local to global
	//        crackStrain.rotatedWith(epsL2G, 'n');
	crackStrain.beProductOf(epsL2G, status->giveCrackStrainVector() );

        StructuralMaterial :: giveFullSymVectorForm( answer, crackStrain, gp->giveMaterialMode() );

        return 1;

    } else if ( type == IST_CrackSlip ) {
      answer.resize(1);
      answer.zero();

      /*      if ( ( mMode == _PlaneStress ) || ( mMode == _PlaneStrain ) ) {
	answer = this->computeShearSlipOnCrack(gp, tStep, 1);
	} else {*/

	int nMaxCracks = status->giveNumberOfTempCracks();
	for ( int i = 1; i <= nMaxCracks; i++ ) {
	  answer.at(1) = max(answer.at(1), this->computeShearSlipOnCrack(gp, tStep, i) );
	}
	//      }
      return 1;

    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

FloatMatrixF<6,6>
FCMMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *tStep) const
{
    FloatMatrix answer;
    const_cast<FCMMaterial*>(this)->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    return answer;
}


FloatMatrixF<3,3>
FCMMaterial :: givePlaneStressStiffMtrx(MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const
{
    FloatMatrix answer;
    const_cast<FCMMaterial*>(this)->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    return answer;
}


FloatMatrixF<4,4>
FCMMaterial :: givePlaneStrainStiffMtrx(MatResponseMode mode,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const

{
    FloatMatrix answer;
    const_cast<FCMMaterial*>(this)->giveMaterialStiffnessMatrix(answer, mode, gp, tStep);
    return answer;
}


FCMMaterialStatus :: FCMMaterialStatus(GaussPoint *gp) :
    StructuralMaterialStatus(gp),
    crackStatuses(), tempCrackStatuses(),
    maxCrackStrains(), tempMaxCrackStrains(),
    crackStrainVector(), tempCrackStrainVector(),
    crackDirs(),
    charLengths(),
    transMatrix_G2Lstress(), transMatrix_G2Lstrain(),
    transMatrix_L2Gstress(), transMatrix_L2Gstrain()
{
    // resize in constructor according to stress-state
    this->nMaxCracks = giveMaxNumberOfCracks(gp);

    crackStatuses.resize(this->nMaxCracks);
    crackStatuses.zero();
    tempCrackStatuses = crackStatuses;

    charLengths.resize(this->nMaxCracks);
    charLengths.zero();

    crackDirs.resize(this->nMaxCracks, this->nMaxCracks);
    crackDirs.zero();
    for ( int i = 1; i <= this->nMaxCracks; i++ ) {
        crackDirs.at(i, i) = 1.0;
    }


    if ( this->nMaxCracks == 2 ) { //plane stress
        maxCrackStrains.resize(3);
        maxCrackStrains.zero();
        tempMaxCrackStrains = maxCrackStrains;

        crackStrainVector = maxCrackStrains;
        tempCrackStrainVector = maxCrackStrains;


        transMatrix_G2Lstress.resize(3, 3);
        transMatrix_G2Lstress.zero();
        transMatrix_L2Gstrain = transMatrix_L2Gstress = transMatrix_G2Lstrain = transMatrix_G2Lstress;
    } else {
        maxCrackStrains.resize(6);
        maxCrackStrains.zero();
        tempMaxCrackStrains = maxCrackStrains;

        crackStrainVector = maxCrackStrains;
        tempCrackStrainVector = maxCrackStrains;


        transMatrix_G2Lstress.resize(6, 6);
        transMatrix_G2Lstress.zero();
        transMatrix_L2Gstrain = transMatrix_L2Gstress = transMatrix_G2Lstrain = transMatrix_G2Lstress;
    }
}


void
FCMMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    char s [ 11 ];

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->giveNumberOfCracks() > 0 ) {
        for ( int i = 1; i <= crackDirs.giveNumberOfColumns(); i++ ) {
            switch ( crackStatuses.at(i) ) {
            case pscm_NONE:
                strcpy(s, "NONE");
                break;
            case pscm_JUST_INIT:
                strcpy(s, "JUST_INIT");
                break;
            case pscm_SOFTENING:
                strcpy(s, "SOFTENING");
                break;
            case pscm_UNLO_RELO:
                strcpy(s, "UNLO_RELO");
                break;
            case pscm_CLOSED:
                strcpy(s, "CLOSED");
                break;
            default:
                strcpy(s, "UNKNOWN");
                break;
            }

            fprintf(file, "crack %d {status %s, crackplane is normal to { ", i, s);

            for ( int j = 1; j <= crackDirs.giveNumberOfRows(); j++ ) {
                fprintf( file, "%f ", crackDirs.at(j, i) );
            }

            fprintf(file, "} ");
            
	    fprintf(file, "crack_width %.4e", this->giveCharLength(i) * this->giveCrackStrain(i) );

            fprintf(file, " } ");	    
        }
    }

    fprintf(file, " }\n");
}

int
FCMMaterialStatus :: giveNumberOfCracks() const
//
// return number of existing cracks
//
{
    int answer = 0;

    for ( int i = 1; i <= crackStatuses.giveSize(); i++ ) {
        if ( crackStatuses.at(i) != pscm_NONE  ) {
            answer++;
        }
    }

    return answer;
}


int
FCMMaterialStatus :: giveNumberOfTempCracks() const
//
// return number of existing temp cracks
//
{
    int answer = 0;

    for ( int i = 1; i <= tempCrackStatuses.giveSize(); i++ ) {
        if ( tempCrackStatuses.at(i) != pscm_NONE  ) {
            answer++;
        }
    }

    return answer;
}

int
FCMMaterialStatus :: giveMaxNumberOfCracks(GaussPoint *gp)
//
// return number of maximum allowable cracks
//
{
    int nCr = 0;
    IntArray indx;

    if ( this->nMaxCracks == 0 ) {
        StructuralMaterial :: giveVoigtSymVectorMask( indx, gp->giveMaterialMode() );

        for ( int i = 1; i <= 3; i++ ) {
            if ( indx.contains(i) ) {
                nCr++;
            }
        }

        this->nMaxCracks = nCr;
    }

    return this->nMaxCracks;
}


void
FCMMaterialStatus :: setTempNormalCrackStrainVector(FloatArray tempNormalCrackStrain)
//
// return number of maximum allowable cracks
//
{

  for ( int i = 1; i <= tempNormalCrackStrain.giveSize(); i++ ) {
    tempCrackStrainVector.at(i) = tempNormalCrackStrain.at(i);
  }

}


void
FCMMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    StructuralMaterialStatus :: initTempStatus();

    tempCrackStatuses = crackStatuses;
    tempMaxCrackStrains = maxCrackStrains;
    tempCrackStrainVector = crackStrainVector;
}



void
FCMMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(tStep);


    maxCrackStrains = tempMaxCrackStrains;
    crackStrainVector = tempCrackStrainVector;

    //    for ( int i = 1; i <= crackStrainVector.giveSize(); i++ ) {
    for ( int i = 1; i <= this->nMaxCracks; i++ ) { // loop only in normal directions!
        if ( crackStrainVector.at(i) < 0. ) {
            crackStrainVector.at(i) = 0.;
        }
    }

    // updating of statuses has to be done more carefully.
    // consider a crack which does not exist in the previous step and ends as closed in the end of this step
    // this crack is naturally treated as "NONE" in the following steps

    for ( int i = 1; i <= crackStatuses.giveSize(); i++ ) {
        if ( ( tempCrackStatuses.at(i) == pscm_CLOSED ) && ( crackStatuses.at(i) == pscm_NONE ) ) {
            // no other crack so this one can be set as non-existing
            if ( i + 1 > nMaxCracks ) {
                crackStatuses.at(i) = pscm_NONE;


                // be sure that in the second and third crack does not exist, if it does we have to copy CLOSED status
            } else if ( tempCrackStatuses.at(i + 1) != pscm_NONE ) {
                crackStatuses.at(i) = tempCrackStatuses.at(i);
            } else {
                crackStatuses.at(i) = pscm_NONE;
            }
        } else {
            crackStatuses.at(i) = tempCrackStatuses.at(i);
        }
    }
}



void
FCMMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = crackStatuses.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = charLengths.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}

void
FCMMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralMaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = crackStatuses.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = maxCrackStrains.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackDirs.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = charLengths.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = crackStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_G2Lstress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = transMatrix_L2Gstress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}
} // end namespace oofem
