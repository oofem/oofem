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

#include "concrete2.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(Concrete2);

Concrete2 :: Concrete2(int n, Domain *d) : DeformationTheoryMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}


Concrete2 :: ~Concrete2()
{
    delete linearElasticMaterial;
}

IRResultType
Concrete2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, E, _IFT_Concrete2_e);
    IR_GIVE_FIELD(ir, n, _IFT_Concrete2_n);
    IR_GIVE_FIELD(ir, SCCC, _IFT_Concrete2_sccc);
    IR_GIVE_FIELD(ir, SCCT, _IFT_Concrete2_scct);
    IR_GIVE_FIELD(ir, EPP, _IFT_Concrete2_epp);
    IR_GIVE_FIELD(ir, EPU, _IFT_Concrete2_epu);
    IR_GIVE_FIELD(ir, EOPP, _IFT_Concrete2_eopp);
    IR_GIVE_FIELD(ir, EOPU, _IFT_Concrete2_eopu);
    IR_GIVE_FIELD(ir, SHEARTOL, _IFT_Concrete2_sheartol);
    IR_GIVE_FIELD(ir, IS_PLASTIC_FLOW, _IFT_Concrete2_is_plastic_flow);
    IR_GIVE_FIELD(ir, IFAD, _IFT_Concrete2_ifad);

    // stirrups constants
    IR_GIVE_FIELD(ir, stirrE, _IFT_Concrete2_stirr_e);
    IR_GIVE_FIELD(ir, stirrFt, _IFT_Concrete2_stirr_ft);
    IR_GIVE_FIELD(ir, stirrA, _IFT_Concrete2_stirr_a);
    IR_GIVE_FIELD(ir, stirrTOL, _IFT_Concrete2_stirr_tol);
    IR_GIVE_FIELD(ir, stirrEREF, _IFT_Concrete2_stirr_eref);
    IR_GIVE_FIELD(ir, stirrLAMBDA, _IFT_Concrete2_stirr_lambda);

    result = this->linearElasticMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    return Material :: initializeFrom(ir);
}


double
Concrete2 :: give(int aProperty, GaussPoint *gp)
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
{
    double value;

    switch ( aProperty ) {
    case c2_SCCC:
        return this->SCCC;

    case c2_SCCT:
        return this->SCCT;

    case c2_EPP:
        return this->EPP;

    case c2_EPU:
        return this->EPU;

    case c2_EOPP:
        return this->EOPP;

    case c2_EOPU:
        return this->EOPU;

    case c2_SHEARTOL:
        return this->SHEARTOL;

    case c2_E:
        return this->E;

    case c2_n:
        return this->n;

    case stirr_E:
        return this->stirrE;

    case stirr_Ft:
        return this->stirrFt;

    case stirr_A:
        return this->stirrA;

    case stirr_TOL:
        return this->stirrTOL;

    case stirr_EREF:
        return this->stirrEREF;

    case stirr_LAMBDA:
        return this->stirrLAMBDA;

    case c2_IS_PLASTIC_FLOW:
        return this->IS_PLASTIC_FLOW;

    case c2_IFAD:
        return this->IFAD;

    default:
        if ( propertyDictionary.includes(aProperty) ) {
            value = propertyDictionary.at(aProperty);
            return value;
        } else {
            return this->linearElasticMaterial->give(aProperty, gp);
            // error ("give: property not defined");
        }
    }
}


void
Concrete2 :: giveRealStressVector_PlateLayer(FloatArray &answer,
                                             GaussPoint *gp,
                                             const FloatArray &totalStrain,
                                             TimeStep *tStep)
//
// returns total stress vector of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
// Plane stress or uniaxial stress + transverse shear in
// concrete layers with transverse stirrups.
//
//
// gp strains - EXX,EYY,2*EXY,2*EYZ,2*EXZ
//
// material constant of concrete
// stored in this c2_property:
//
// E    - Young modulus
// ny   - Poisson ratio
// Ro   - specific mass
// SCCC<=0 - pressure strength
// SCCT>=0 - tension strength
// EPP>=0  - treshold eff. plastic strain for softening in compress.
// EPU>=epp- ultimate eff. pl. strain
// EOPP>=0 - threshold vlumetric plastic strain for soft. in tension
// EOPU>=EOPP ultimate vol. pl. strain
// SHEARTOL  threshold value of the relative shear deformation
//           (psi**2/eef) at which shear is considered in layers. for
//           lower r.s.d. the transverse shear remains elastic decoupled
//           from bending. default value SHEARTOL = 0.01
// IS_PLASTIC_FLOW   indicates that plastic flow (not deformation theory)
//                   is used in pressure.
// IFAD   IFAD<=0 STATE VARIABLES WILL NOT BE UPDATED
//        >0 UPDATE S.V.
//
//
// state varibles of the layer are stored in corresponding Material Status
//
// SCCM      current pressure strenght
// EPM       max. eff. plastic strain
// SCTM      current tension strenght
// E0PM      max. vol. plastic strain
// SRF       current stress in stirrups
// SEZ       current strain in transverse (z) direction.
// if flow teory of compresion then also available are:
// EXX_PLAST, EYY_PLAST, EZZ_PLAST, EXY_PLAST=2.EXY_PLAST, EYZ_PLAST (2*),EXZ_PLAST
//    - components of plastic strain associated with epp
//
// Note: formulated in Full stress strain space
{
    double us, SCC, SCT, psi2, comp, tol = 0.0, crit = 0.0, nez = 0.0, effstold = 0.0, dezold = 0.0, effst = 0.0, dez = 0.0;
    double epsult = 0.0, ep1 = -1.0, ep2 = -1.0, ep3 = -1.0, relus, G, ez;
    int ifsh, i, ifplas, ifupd;
    FloatArray plasticStrainVector, plasticStrainIncrementVector;

    Concrete2MaterialStatus *status = static_cast< Concrete2MaterialStatus * >( this->giveStatus(gp) );
    FloatArray currentStress, currentStrain, currentEffStrain(6), pVal, ep, pStress, strainIncr, plasticStrain, help, reducedStrain, helpR;
    FloatMatrix pDir;

    pDir.resize(3, 3); // principal directions
    pStress.resize(3); // principal stresses
    strainIncr.resize(1);

    this->initTempStatus(gp);
    
    // subtract stress independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigen strain value
    this->giveStressDependentPartOfStrainVector(reducedStrain, gp,
                                                totalStrain, tStep, VM_Total);

    StructuralMaterial :: giveFullSymVectorForm( currentStrain, reducedStrain, gp->giveMaterialMode() );
    StructuralMaterial :: giveFullSymVectorForm( currentStress, status->giveStressVector(), gp->giveMaterialMode() );

    if ( status->giveTempCurrentTensionStrength() < 0. ) {
        //
        // called for the firs time - init flags
        //
        status->giveTempCurrentPressureStrength() = this->give(c2_SCCC, gp);
        status->giveTempCurrentTensionStrength() = this->give(c2_SCCT, gp);

        status->giveCurrentPressureStrength() = this->give(c2_SCCC, gp);
        status->giveCurrentTensionStrength() = this->give(c2_SCCT, gp);
    }

    plasticStrainVector = status->givePlasticStrainVector();
    StructuralMaterial :: giveFullSymVectorForm( plasticStrain, plasticStrainVector, gp->giveMaterialMode() );

    ep.resize(3);

    //for (i =1; i <= 6; i++) currentStrain.at(i) += strainIncrement.at(i);

    ez  = status->giveTempCurrentStrainInZDir();
    us  = 0.;
    // eep = 0.;
    //
    // Concrete strenghts, SCC and SCT are actual strenghts
    // which may change in the iteration if there is any
    //
    SCC = this->give(c2_SCCC, gp);
    if ( this->give(c2_EPP, gp) != 0. ) {
        SCC = status->giveTempCurrentPressureStrength();
    }

    SCT = this->give(c2_SCCT, gp);
    if ( this->give(c2_EOPP, gp) != 0. ) {
        SCT = status->giveTempCurrentTensionStrength();
    }

    //
    // create full 3d strain tensor
    //

    if ( this->give(c2_IS_PLASTIC_FLOW, gp) > 0. ) {
        currentEffStrain.at(1) = currentStrain.at(1) - plasticStrain.at(1);
        currentEffStrain.at(2) = currentStrain.at(2) - plasticStrain.at(2);
        currentEffStrain.at(3) = ez - plasticStrain.at(3);
        currentEffStrain.at(4) = currentStrain.at(4) - plasticStrain.at(4);
        currentEffStrain.at(5) = currentStrain.at(5) - plasticStrain.at(5);
        currentEffStrain.at(6) = currentStrain.at(6) - plasticStrain.at(6);
    } else {
        currentEffStrain.at(1) = currentStrain.at(1);
        currentEffStrain.at(2) = currentStrain.at(2);
        currentEffStrain.at(3) = ez;
        currentEffStrain.at(4) = currentStrain.at(4);
        currentEffStrain.at(5) = currentStrain.at(5);
        currentEffStrain.at(6) = currentStrain.at(6);
    }

    //
    // check for shear to be considered
    //
    ifsh = 0;
    psi2 = currentEffStrain.at(4) * currentEffStrain.at(4) +
    currentEffStrain.at(5) * currentEffStrain.at(5);
    comp = sqrt(currentEffStrain.at(1) * currentEffStrain.at(1) +
                currentEffStrain.at(2) * currentEffStrain.at(2) +
                currentEffStrain.at(3) * currentEffStrain.at(3) +
                currentEffStrain.at(6) * currentEffStrain.at(6) +
                psi2);

    if ( comp == 0. ) { // strange
        for ( i = 1; i <= 6; i++ ) { // directly modify gp-record
            currentStress.at(i) = 0.;
        }

        // note in status are reduced space variables
        status->letTempStrainVectorBe(totalStrain);

        StructuralMaterial :: giveReducedSymVectorForm( helpR, currentStress, gp->giveMaterialMode() );

        status->letTempStressVectorBe(helpR);

        StructuralMaterial :: giveReducedSymVectorForm( help, plasticStrain, gp->giveMaterialMode() );

        plasticStrainIncrementVector = status->givePlasticStrainIncrementVector();
        plasticStrainIncrementVector.subtract(plasticStrainVector);
        plasticStrainIncrementVector.add(help);
        status->letPlasticStrainIncrementVectorBe(plasticStrainIncrementVector);
        //   plasticStrain->negated()->add (status->givePlasticStrainVector());
        //      status->givePlasticStrainIncrementVector()-> add(plasticStrain);

        StructuralMaterial :: giveReducedSymVectorForm( answer, currentStress, gp->giveMaterialMode() );
        return;
    }

    //
    // Treshold relative shear strain to trigger shear computation
    // is stored in SHEARTOL material constant. Once triggered (EZ.ne 0)
    // then shear computation cannot be abandoned even if the
    // shear strain drops bellow the limit.
    //
    if ( ( ez != 0. ) && ( ( sqrt(psi2) / comp ) < this->give(c2_SHEARTOL, gp) ) ) {
        //
        // no transverse shear - 2d case
        //

        this->computePrincipalValDir(pVal, pDir, currentEffStrain, principal_strain);
        //
        // dtp2 procedure expects first two princ values to be  in xy plane
        // so sort results acordingly
        //
        int j = 1;
        if ( pDir.at(3, 2) > pDir.at(3, j) ) {
            j = 2;                         // find the eigen direction
        }

        if ( pDir.at(3, 3) > pDir.at(3, j) ) {
            j = 3;                         // normal to xy plane
        }

        if ( j != 3 ) {
            // swap  values
            int k;
            double swap;
            for ( k = 1; k <= 3; k++ ) { // swap eigenvectors
                swap = pDir.at(k, 3);
                pDir.at(k, 3) = pDir.at(k, j);
                pDir.at(k, j) = swap;
            }

            swap = pVal.at(3);
            pVal.at(3) = pVal.at(j);
            pVal.at(j) = swap;
        }

        this->dtp2(gp, pVal, pStress, ep, SCC, SCT, & ifplas);
    } else {
        ifsh = 1;
        //
        // optional iteration to achieve equilibrium in the z - direction.
        //   // iteration toleration limit tol (=0. -- no iteration)
        // iteration only possible when stirrups present
        // Number of loops nez
        tol = 0.;
        if ( this->give(stirr_E, gp) > 0 ) {
            tol = this->give(stirr_TOL, gp);
        }

        crit = fabs( this->give(c2_SCCC, gp) ) + this->give(c2_SCCT, gp);
        nez = 0;
        effstold = 0.;
        dezold = 0.;
label18:
        //
        // find principal strains and corresponding principal direction
        //
        this->computePrincipalValDir(pVal, pDir, currentEffStrain, principal_strain);
        //
        // plasticity condition
        //
        this->dtp3(gp, pVal, pStress, ep, SCC, SCT, & ifplas);
        //
        //  transverse equilibrium
        //  transv strain  ez (flags -> at(SEZ))  its incr.  dez  ,
        //  stirrup stress flags->at(SRF)
        //
        //  unbalanced transverse stress us
        //  concrete component
        //
        us = ( pDir.at(3, 1) * pStress.at(1) * pDir.at(3, 1) +
              pDir.at(3, 2) * pStress.at(2) * pDir.at(3, 2) +
              pDir.at(3, 3) * pStress.at(3) * pDir.at(3, 3) );
        //
        // Stirrups component
        //
        if ( this->give(stirr_E, gp) > 0. ) {
            us += status->giveTempCurrentStressInStirrups() * this->give(stirr_A, gp);
        }

        //
        // increment of transverse strain ez
        // effective modulus of concrete in z direction
        // f(1)*red**4  where red is sinus of the angle
        // between z axis and max. principal strain dir.
        // when us becomes positive (may happen if the
        // effective modulus was too small and, hence, the
        // last change dez of the transverse strain ez
        // too large in the previous increment), it loses sense.
        // the Jacobi iteration delivers just approximate princ.
        // vectors which causes spurious pulses in us and fluctuations
        // in sxz,syz,qx,qz. they decrease
        // when toli->0. they are successfully managed by z equil.
        // iteration and do not harm the solution.
        //
        if ( ( ifplas ) && ( us < 0. ) ) {
            //   transverse strain increment if concrete is plastic in tension
            //   and further Z-extension is necessary (US<0).
            //   effst is an iteration factor to get transverse equilibrium
            if ( pVal.at(1) > pVal.at(2) ) {
                if ( pVal.at(1) > pVal.at(3) ) {
                    effst = pDir.at(3, 1);
                } else {
                    effst = pDir.at(3, 3);
                }
            } else {
                if ( pVal.at(2) > pVal.at(3) ) {
                    effst = pDir.at(3, 2);
                } else {
                    effst = pDir.at(3, 3);
                }
            }

            effst = ( 1. - effst * effst );
            effst *= effst;
            if ( effst < 1.e-6 ) {
                effst = 1.e-6;
            }

            effst *= this->give(c2_E, gp);
            if ( this->give(stirr_E, gp) > 0. ) {
                effst += this->give(stirr_E, gp) * this->give(stirr_A, gp);
            }

            effst *= 2.0;
        } else {
            // concrete is elastic
            effst = this->give(c2_E, gp) * 1.1;
        }

        dez = -us / effst;
        // if the increment of ez has changed sign in this z-equil.
        // iteration loop, raise the iteration factor by adding up
        // the two last values
        if ( nez > 1 ) {
            if ( dez * dezold < 0. ) {
                effst += effstold;
                dez = -us / effst;
            }

            dezold = dez;
            effstold = effst;
        }

        //
        ez  +=  dez;
        if ( this->give(stirr_E, gp) > 0. ) {
            strainIncr.at(1) = dez;
            this->updateStirrups(gp, strainIncr, tStep);
            // updates flags->at(SEZ)
            // stirr(dez,sv(l3),fst);
        }
    }

    //
    // strain softening
    //
    ifupd = 0;
    if ( this->give(c2_IFAD, gp) ) {
        if ( this->give(c2_EOPP, gp) >= 0. ) {
            if ( this->give(c2_EOPU, gp) > 0. ) {
                epsult = this->give(c2_EOPU, gp);
            } else {
                // position dependent strain softening
                // NOT IMPLEMENTED
                epsult = this->give(c2_EOPU, gp);
            }
        }

        this->strsoft(gp, epsult, ep, ep1, ep2, ep3, SCC, SCT, ifupd);
    }

    //
    // End of the Ez iteration loop
    //
    if ( ifsh ) {
        relus = fabs(us) / crit;
        if ( ( tol > 0. ) && ( relus > tol ) ) {
            //
            //     add the increment of EZ to the actual strain tenzor
            //
            currentEffStrain.at(1) += pDir.at(3, 1) * dez * pDir.at(3, 1);
            currentEffStrain.at(2) += pDir.at(3, 2) * dez * pDir.at(3, 2);
            currentEffStrain.at(3) += pDir.at(3, 3) * dez * pDir.at(3, 3);
            currentEffStrain.at(4) += pDir.at(3, 2) * dez * pDir.at(3, 3);
            currentEffStrain.at(5) += pDir.at(3, 1) * dez * pDir.at(3, 3);
            currentEffStrain.at(6) += pDir.at(3, 1) * dez * pDir.at(3, 2);
            //     back to EZ iteration start
            //
            if ( this->give(c2_IS_PLASTIC_FLOW, gp) && this->give(c2_IFAD, gp) ) {
                plasticStrain.at(1) += ( pDir.at(1, 1) * ep1 * pDir.at(1, 1) +
                                        pDir.at(1, 2) * ep2 * pDir.at(1, 2) +
                                        pDir.at(1, 3) * ep3 * pDir.at(1, 3) );
                plasticStrain.at(2) += ( pDir.at(2, 1) * ep1 * pDir.at(2, 1) +
                                        pDir.at(2, 2) * ep2 * pDir.at(2, 2) +
                                        pDir.at(2, 3) * ep3 * pDir.at(2, 3) );
                plasticStrain.at(3) += ( pDir.at(3, 1) * ep1 * pDir.at(3, 1) +
                                        pDir.at(3, 2) * ep2 * pDir.at(3, 2) +
                                        pDir.at(3, 3) * ep3 * pDir.at(3, 3) );
                plasticStrain.at(4) += 2. * ( pDir.at(2, 1) * ep1 * pDir.at(3, 1) +
                                             pDir.at(2, 2) * ep2 * pDir.at(3, 2) +
                                             pDir.at(2, 3) * ep3 * pDir.at(3, 3) );
                plasticStrain.at(5) += 2. * ( pDir.at(1, 1) * ep1 * pDir.at(3, 1) +
                                             pDir.at(1, 2) * ep2 * pDir.at(3, 2) +
                                             pDir.at(1, 3) * ep3 * pDir.at(3, 3) );
                plasticStrain.at(6) += 2. * ( pDir.at(1, 1) * ep1 * pDir.at(2, 1) +
                                             pDir.at(1, 2) * ep2 * pDir.at(2, 2) +
                                             pDir.at(1, 3) * ep3 * pDir.at(2, 3) );
            }

            pVal.zero();
            goto label18;
        }
    }

    //
    // update compressive plastic strains.
    //
    if ( ifupd ) {
        if ( this->give(c2_IFAD, gp) ) {
            if ( this->give(c2_IS_PLASTIC_FLOW, gp) ) {
                // -
                if ( ifsh == 0 ) {
                    //
                    //  no transverse shear, 1-2 to x-y
                    //
                    plasticStrain.at(1) += ( pDir.at(1, 1) * ep1 * pDir.at(1, 1) +
                                            pDir.at(1, 2) * ep2 * pDir.at(1, 2) );
                    plasticStrain.at(2) += ( pDir.at(2, 1) * ep1 * pDir.at(2, 1) +
                                            pDir.at(2, 2) * ep2 * pDir.at(2, 2) );
                    plasticStrain.at(6) += 2. * ( pDir.at(1, 1) * ep1 * pDir.at(2, 1) +
                                                 pDir.at(1, 2) * ep2 * pDir.at(2, 2) );
                } else {
                    //
                    //  transverse shear present , transform back to x,y,z
                    //
                    plasticStrain.at(1) += ( pDir.at(1, 1) * ep1 * pDir.at(1, 1) +
                                            pDir.at(1, 2) * ep2 * pDir.at(1, 2) +
                                            pDir.at(1, 3) * ep3 * pDir.at(1, 3) );
                    plasticStrain.at(2) += ( pDir.at(2, 1) * ep1 * pDir.at(2, 1) +
                                            pDir.at(2, 2) * ep2 * pDir.at(2, 2) +
                                            pDir.at(2, 3) * ep3 * pDir.at(2, 3) );
                    plasticStrain.at(3) += ( pDir.at(3, 1) * ep1 * pDir.at(3, 1) +
                                            pDir.at(3, 2) * ep2 * pDir.at(3, 2) +
                                            pDir.at(3, 3) * ep3 * pDir.at(3, 3) );
                    plasticStrain.at(4) += 2. * ( pDir.at(2, 1) * ep1 * pDir.at(3, 1) +
                                                 pDir.at(2, 2) * ep2 * pDir.at(3, 2) +
                                                 pDir.at(2, 3) * ep3 * pDir.at(3, 3) );
                    plasticStrain.at(5) += 2. * ( pDir.at(1, 1) * ep1 * pDir.at(3, 1) +
                                                 pDir.at(1, 2) * ep2 * pDir.at(3, 2) +
                                                 pDir.at(1, 3) * ep3 * pDir.at(3, 3) );
                    plasticStrain.at(6) += 2. * ( pDir.at(1, 1) * ep1 * pDir.at(2, 1) +
                                                 pDir.at(1, 2) * ep2 * pDir.at(2, 2) +
                                                 pDir.at(1, 3) * ep3 * pDir.at(2, 3) );
                }
            }
        }
    }

    //
    // stresses back , 1-2 to x-z , sxz= tau=-sigxz
    //
    // 36
    if ( ifsh == 0 ) {
        G = 0.5 * this->give(c2_E, gp) / ( 1.0 + this->give(c2_n, gp) );

        currentStress.at(1) = ( pDir.at(1, 1) * pStress.at(1) * pDir.at(1, 1) +
                               pDir.at(1, 2) * pStress.at(2) * pDir.at(1, 2) );
        currentStress.at(2) = ( pDir.at(2, 1) * pStress.at(1) * pDir.at(2, 1) +
                               pDir.at(2, 2) * pStress.at(2) * pDir.at(2, 2) );
        currentStress.at(3) = 0.;
        currentStress.at(4) = currentStrain.at(4) * G;
        currentStress.at(5) = currentStrain.at(5) * G;
        currentStress.at(6) = ( pDir.at(1, 1) * pStress.at(1) * pDir.at(2, 1) +
                               pDir.at(1, 2) * pStress.at(2) * pDir.at(2, 2) );
    } else {
        currentStress.at(1) = ( pDir.at(1, 1) * pStress.at(1) * pDir.at(1, 1) +
                               pDir.at(1, 2) * pStress.at(2) * pDir.at(1, 2) +
                               pDir.at(1, 3) * pStress.at(3) * pDir.at(1, 3) );
        currentStress.at(2) = ( pDir.at(2, 1) * pStress.at(1) * pDir.at(2, 1) +
                               pDir.at(2, 2) * pStress.at(2) * pDir.at(2, 2) +
                               pDir.at(2, 3) * pStress.at(3) * pDir.at(2, 3) );
        currentStress.at(3) =  us;
        currentStress.at(4) = ( pDir.at(2, 1) * pStress.at(1) * pDir.at(3, 1) +
                               pDir.at(2, 2) * pStress.at(2) * pDir.at(3, 2) +
                               pDir.at(2, 3) * pStress.at(3) * pDir.at(3, 3) );
        currentStress.at(5) = ( pDir.at(1, 1) * pStress.at(1) * pDir.at(3, 1) +
                               pDir.at(1, 2) * pStress.at(2) * pDir.at(3, 2) +
                               pDir.at(1, 3) * pStress.at(3) * pDir.at(3, 3) );
        currentStress.at(6) = ( pDir.at(1, 1) * pStress.at(1) * pDir.at(2, 1) +
                               pDir.at(1, 2) * pStress.at(2) * pDir.at(2, 2) +
                               pDir.at(1, 3) * pStress.at(3) * pDir.at(2, 3) );
    }

    if ( this->give(c2_IFAD, gp) ) {
        status->giveTempCurrentStrainInZDir() = ez;
    }

    //totalStrainVector = status->giveStrainVector();
    //totalStrainVector.add (totalStrainIncrement);
    status->letTempStrainVectorBe(totalStrain);

    StructuralMaterial :: giveReducedSymVectorForm( helpR, currentStress, gp->giveMaterialMode() );
    status->letTempStressVectorBe(helpR);

    StructuralMaterial :: giveFullSymVectorForm( help, plasticStrain, gp->giveMaterialMode() );
    plasticStrainIncrementVector = status->givePlasticStrainIncrementVector();
    plasticStrainIncrementVector.subtract(plasticStrainVector);
    plasticStrainIncrementVector.add(help);
    status->letPlasticStrainIncrementVectorBe(plasticStrainIncrementVector);

    //plasticStrain->negated()->add (status->givePlasticStrainVector());
    //status->givePlasticStrainIncrementVector()-> add(plasticStrain);

    StructuralMaterial :: giveFullSymVectorForm( answer, currentStress, gp->giveMaterialMode() );
}


void
Concrete2 :: dtp3(GaussPoint *gp, FloatArray &e, FloatArray &s, FloatArray &ep,
                  double SCC, double SCT, int *ifplas)
//
// DEFORMATION THEORY OF PLASTICITY, PRNCIPAL STRESSES
// CONDITION OF PLASTICITY.
//
// e  - principal strains
// s  - principal stresses (output)
// SCT,SCC - tensile and compressive strengths
// ifplas - =0 if material is elastic
//          =1 if plastic in tension
//
//
// state varibles of the layer are stored
// in corresponding gp  mat Status
//
// SCCM      current pressure strenght
// EPM       max. eff. plastic strain
// SCTM      current tension strenght
// E0PM      max. vol. plastic strain
// SRF       current stress in stirrups
// SEZ       current strain in transverse (z) direction.
// if flow teory of compresion then also available are:
// EXX_PLAST, EYY_PLAST, EZZ_PLAST, EXY_PLAST=2.EXY_PLAST, EYZ_PLAST (2*),EXZ_PLAST
//    - components of plastic strain associated with epp in gp->plasticStrainVector
//
{
    int i, ii = 0, j, k, jj = 0, kk = 0;
    double yy, ey, e0, yy1, yy2, s0, sci = 0.0, scj = 0.0, sck = 0.0;

    if ( e.giveSize() != 3 ) {
        OOFEM_ERROR("principal strains size mismatch");
    }

    if ( s.giveSize() != 3 ) {
        OOFEM_ERROR("principal stress size mismatch");
    }

    if ( ep.giveSize() != 3 ) {
        OOFEM_ERROR("plastic strains size mismatch");
    }

    * ifplas = 0;

    yy  = this->give(c2_n, gp);
    ey  = this->give(c2_E, gp);
    e0  = this->give(c2_E, gp) / ( 1. + yy ) / ( 1. - yy - yy );
    yy1 = 1. - yy;
    yy2 = 1. + yy;

    if ( this->give(c2_E, gp) < 1.e-6 ) {
        for ( i = 1; i <= 3; i++ ) {
            s.at(i) = 0.;
            ep.at(i) = 0.;
        }

        return;
    }

    if ( this->give(c2_n, gp) >= 0.01 ) {
        s.at(1) =  e0 * ( yy1 * e.at(1) + yy * ( e.at(2) + e.at(3) ) );
        s.at(2) =  e0 * ( yy1 * e.at(2) + yy * ( e.at(1) + e.at(3) ) );
        s.at(3) =  e0 * ( yy1 * e.at(3) + yy * ( e.at(2) + e.at(1) ) );
        //
        // find intersection with side  i , sig(i)=sci
        //
        s0 = 0.;
        for ( i = 1; i <= 3; i++ ) {
            if ( ( s.at(i) - SCT - s0 ) > 0. ) {
                ii = i;
                sci = SCT;
                s0 = s.at(i) - SCT;
                * ifplas = 1;
            } else {
                if ( ( SCC - s.at(i) - s0 ) > 0. ) {
                    ii = i;
                    sci = SCC;
                    s0 = SCC - s.at(i);
                }
            }
        }

        if ( s0 <=  0. ) {
            //
            //  linear elastic
            //
            ep.at(1) = 0.;
            ep.at(2) = 0.;
            ep.at(3) = 0.;
            return;
        }

        //
        // Strain vector outside elastic cube
        // Plastic state
        //
        j = iperm(ii, 3);
        k = iperm(j, 3);
        i = ii;
        s.at(i) = sci;
        s.at(j) = ( yy * sci + ey / yy2 * ( e.at(j) + yy * e.at(k) ) ) / yy1;
        s.at(k) = ( yy * sci + ey / yy2 * ( e.at(k) + yy * e.at(j) ) ) / yy1;
        s0 = 0.;
        //
        // find intersection with edges j and k
        //
        if ( ( s.at(j) - SCT - s0 ) > 0. ) {
            jj = j;
            kk = k;
            scj = SCT;
            s0  = s.at(j) - SCT;
        }

        if ( ( SCC - s.at(j) - s0 ) > 0. ) {
            jj = j;
            kk = k;
            scj = SCC;
            s0  = SCC - s.at(j);
        }

        if ( ( s.at(k) - SCT - s0 ) > 0. ) {
            jj = k;
            kk = j;
            scj = SCT;
            s0 = s.at(k) - SCT;
        }

        if ( ( SCC - s.at(k) - s0 ) > 0. ) {
            jj = k;
            kk = j;
            scj = SCC;
            s0 = SCC - s.at(k);
        } else {
            //
            // staying on the side  i
            //
            if ( s0 <= 0. ) {
                ep.at(i) = e.at(i) - sci / yy1 / e0 + yy / yy1 * ( e.at(j) + e.at(k) );
                ep.at(j) = 0.;
                ep.at(k) = 0.;
                return;
            }
        }

        //
        // edge  j  on the side  i
        //
        j = jj;
        k = kk;
        s.at(j) = scj;
        s.at(k) = yy * ( sci + scj ) + ey *e.at(k);
        //
        // corner ?
        //
        if ( ( s.at(k) - SCT ) > 0. ) {
            sck = SCT;
        } else {
            if ( ( SCC - s.at(k) ) > 0. ) {
                sck = SCC;
            } else {
                //
                // staying on the edge  j - i
                //
                ep.at(i) = e.at(i) + yy *e.at(k) + yy2 / ey * ( yy * scj - yy1 * sci );
                ep.at(j) = e.at(j) + yy *e.at(k) + yy2 / ey * ( yy * sci - yy1 * scj );
                ep.at(k) = 0.;
                return;
            }
        }

        //
        // corner sck on the edge  j  of the side  i
        //
        s.at(k) = sck;
        ep.at(i) = e.at(i) - ( sci - yy * ( scj + sck ) ) / ey;
        ep.at(j) = e.at(j) - ( scj - yy * ( sck + sci ) ) / ey;
        ep.at(k) = e.at(k) - ( sck - yy * ( sci + scj ) ) / ey;
        return;
    } else {
        //
        // Poisson ratio = 0.
        //
        for ( i = 1; i <= 3; i++ ) {
            sci = 0.;
            ep.at(i) = 0.;
            s.at(i)  = e.at(i) * ey;
            if ( s.at(i) > SCT ) {
                sci = SCT;
                * ifplas = 1;
                s.at(i) = sci;
                ep.at(i) -= sci / ey;
            } else {
                if ( s.at(i) < SCC ) {
                    sci = SCC;
                    s.at(i) = sci;
                    ep.at(i) -= sci / ey;
                }
            }
        }
    }
}

void
Concrete2 :: dtp2(GaussPoint *gp, FloatArray &e, FloatArray &s, FloatArray &ep,
                  double SCC, double SCT, int *ifplas)
//
// DEFORMATION THEORY OF PLASTICITY, PRNCIPAL STRESSES
// CONDITION OF PLASTICITY. - PLANE STRESS
//
// e  - principal strains
// s  - principal stresses (output)
// SCT,SCC - tensile and compressive strengths
// ifplas - =0 if material is elastic
//          =1 if plastic in tension
//
//
// state varibles of the layer are stored
// in gp->matInfo->flagDictionary = flags
//
// SCCM      current pressure strenght
// EPM       max. eff. plastic strain
// SCTM      current tension strenght
// E0PM      max. vol. plastic strain
// SRF       current stress in stirrups
// SEZ       current strain in transverse (z) direction.
// if flow teory of compresion then also available are:
// EXX_PLAST, EYY_PLAST, EZZ_PLAST, EXY_PLAST=2.EXY_PLAST, EYZ_PLAST (2*),EXZ_PLAST
//    - components of plastic strain associated with epp
//
{
    int i, j, k, jj = 0;
    double yy, ey, e0, so, scj = 0.0, sck = 0.0, sci = 0.0;
    // Dictionary *flags = gp->matInfo->flagDictionary;

    if ( e.giveSize() != 3 ) {
        OOFEM_ERROR("principal strains size mismatch");
    }

    if ( s.giveSize() != 3 ) {
        OOFEM_ERROR("principal stress size mismatch");
    }

    if ( ep.giveSize() != 3 ) {
        OOFEM_ERROR("plastic strains size mismatch");
    }

    * ifplas = 0;

    yy  = this->give(c2_n, gp);
    ey  = this->give(c2_E, gp);

    * ifplas = 0;
    if ( ey == 0. ) {
        s.at(1) = s.at(2) = 0.;
        ep.at(1) = ep.at(2) = 0.;
        return;
    }

    if ( yy >  0.01 ) {
        e0 = ey / ( 1. - yy * yy );

        s.at(1) = e0 * ( e.at(1) + yy * e.at(2) );
        s.at(2) = e0 * ( e.at(2) + yy * e.at(1) );
        so = 0.;
        //
        //     find intersection with the edges of the square plast. co
        //     dition - edge j
        //
        for ( j = 1; j <= 2; j++ ) {
            if ( ( s.at(j) - SCT - so ) > 0. ) {
                jj = j;
                scj = SCT;
                so = s.at(j) - SCT;
                * ifplas = 1;
            } else {
                if ( ( SCC - s.at(j) - so ) > 0. ) {
                    jj = j;
                    scj = SCC;
                    so  = SCC - s.at(j);
                }
            }
        }

        if ( so <= 0. ) {
            //
            // LINEAR ELASTIC - NO INTERSECTION
            //
            ep.at(1) = 0.;
            ep.at(2) = 0.;
            return;
        }

        //
        //     edge jj,j
        //
        k = iperm(jj, 2);
        j = jj;
        s.at(j) = scj;
        s.at(k) = ey * e.at(k) + scj * yy;
        //
        //     any corner?
        //
        if ( ( s.at(k) - SCT ) > 0. ) {
            sck = SCT;
        } else {
            if ( ( SCC - s.at(k) ) > 0. ) {
                sck = SCC;
            } else {
                //
                //    staying on the edge  j
                //
                ep.at(j) = e.at(j) + yy *e.at(k) - scj / e0;
                ep.at(k) = 0.;
                return;
            }
        }

        //
        //     corner  sck  on the edge  j
        //
        s.at(k) = sck;
        ep.at(j) = e.at(j) - ( s.at(j) - yy * s.at(k) ) / ey;
        ep.at(k) = e.at(k) - ( s.at(k) - yy * s.at(j) ) / ey;
        return;
    } else {
        //
        //     poisson ratio =0.
        //
        for ( i = 1; i <= 2; i++ ) {
            sci = 0.;
            ep.at(i) = 0.;
            s.at(i) = e.at(i) * ey;
            if ( s.at(i) > SCT ) {
                sci = SCT;
                * ifplas = 1;
            } else {
                sci = SCC;
            }

            s.at(i) = sci;
            ep.at(i) = e.at(i) - sci / ey;
        }

        return;
    }
}




void
Concrete2 :: strsoft(GaussPoint *gp, double epsult, FloatArray &ep, double &ep1, double &ep2, double &ep3,
                     double SCC, double SCT, int &ifupd)
// material constant of concrete
// stored in this.propertyDictionary
//
// E    - Young modulus
// ny   - Poisson ratio
// Ro   - specific mass
// SCCC<=0 - pressure strength
// SCCT>=0 - tension strength
// EPP>=0  - treshold eff. plastic strain for softening in compress.
// EPU>=epp- ultimate eff. pl. strain
// EOPP>=0 - threshold vlumetric plastic strain for soft. in tension
// EOPU>=EOPP ultimate vol. pl. strain
// SHEARTOL  threshold value of the relative shear deformation
//           (psi**2/eef) at which shear is considered in layers. for
//           lower r.s.d. the transverse shear remains elastic decoupled
//           from bending. default value SHEARTOL = 0.01
// IS_PLASTIC_FLOW   indicates that plastic flow (not deformation theory)
//                   is used in pressure.
// IFAD   IFAD<=0 STATE VARIABLES WILL NOT BE UPDATED
//        >0 UPDATE S.V.
//
//
// state varibles of the layer are stored
// in gp->matInfo->flagDictionary = flags
//
// SCCM      current pressure strenght
// EPM       max. eff. plastic strain
// SCTM      current tension strenght
// E0PM      max. vol. plastic strain
// SRF       current stress in stirrups
// SEZ       current strain in transverse (z) direction.
// if flow teory of compresion then also available are:
// EXX_PLAST, EYY_PLAST, EZZ_PLAST, EXY_PLAST=2.EXY_PLAST, EYZ_PLAST (2*),EXZ_PLAST
//    - components of plastic strain associated with epp in gp->plasticStrainVector.
//
{
    Concrete2MaterialStatus *status = static_cast< Concrete2MaterialStatus * >( this->giveStatus(gp) );
    double eop, eopr, dep, d, eep, eepr;
    //
    // strain softening
    // tension
    //
    if ( this->give(c2_EOPP, gp) != 0. ) {
        eop = max(ep.at(1), 0.) + max(ep.at(2), 0.) +
        max(ep.at(3), 0.);
        if ( eop >  status->giveTempMaxVolPlasticStrain() ) {
            //
            // Current plastic volum. strain is greater than max.
            // reached up to now.
            //
            // reduced plastic strain withe respect to the elastic
            // limit strain of the virgin concrete E0PR
            //
            eopr = eop - ( this->give(c2_SCCT, gp) - SCT ) / this->give(c2_E, gp);
            if ( eopr <= this->give(c2_EOPP, gp) ) {
                status->giveTempMaxVolPlasticStrain() = eop;
                goto label14;
            }

            //
            //   softening
            //
            if ( eop > epsult ) {
                SCT = 0.;
            } else {
                SCT = this->give(c2_SCCT, gp) * ( 1. - ( eopr - this->give(c2_EOPP, gp) ) / ( epsult - this->give(c2_EOPP, gp) ) );
            }

            // When actual strength is lowered the eff. plastic strain
            // increases though no additional strain occurred
            status->giveTempMaxVolPlasticStrain() = eop + ( status->giveTempCurrentTensionStrength() - SCT ) /
                                                    this->give(c2_E, gp);
            status->giveTempCurrentTensionStrength() = SCT;
        }
    }

    //
    // COMPRESSION - EFF. PLASTIC COMPR. STRAIN  EEP
    // IFUPD=1 INDICATES THAT COMPRESSIVE PLASTIC STRAINS ARE UPDATED
    //
    //

label14:

    ifupd = 0;
    if ( !( ( this->give(c2_EPP, gp) == 0. ) && ( !this->give(c2_IS_PLASTIC_FLOW, gp) ) ) ) {
        ep1 = min(ep.at(1), 0.);
        ep2 = min(ep.at(2), 0.);
        ep3 = min(ep.at(3), 0.);
        dep = ep1 + ep2 + ep3;
        //
        if ( dep < 0. ) {
            // There is a compressive plastic deformation.
            // If the plastic
            // flow theory is used then DEP is the
            // increment of the plastic
            // strain.
            if ( this->give(c2_EPP, gp) == 0. ) {
                // No compression softening.
                // For flow theory the plastic strain component
                // must be updated
                if ( this->give(c2_IS_PLASTIC_FLOW, gp) ) {
                    ifupd = 1;
                    return;
                }
            }

            d = 1.5 * ( ep1 * ep1 + ep2 * ep2 + ep3 * ep3 ) - 0.5 * dep * dep;
            dep = sqrt(d);
            if ( this->give(c2_IS_PLASTIC_FLOW, gp) ) {
                // flow theory
                eep = status->giveTempMaxEffPlasticStrain() + dep;
            } else {
                // deformation theory
                eep = dep;
            }

            if ( eep > status->giveTempMaxEffPlasticStrain() ) {
                //
                // current plast.def. EEP is greater than the max. reached
                // upto now EPM - strength must be lowered and
                // plastic strain components updated if flow rule.
                if ( this->give(c2_IS_PLASTIC_FLOW, gp) ) {
                    ifupd = 1;
                }

                //
                //     SOFTENING?
                //
                if ( this->give(c2_EPP, gp) != 0. ) {
                    //
                    // reduced plastic strain with respect to the elastic
                    // limit of the virgin concrete EEPR
                    //
                    eepr = eep + ( this->give(c2_SCCC, gp) - SCC ) / this->give(c2_E, gp);
                    if ( eepr <= this->give(c2_EPP, gp) ) {
                        status->giveTempMaxEffPlasticStrain() = eep;
                    } else {
                        // softening
                        if ( eepr >= this->give(c2_EPU, gp) ) {
                            SCC = 0.;
                        } else {
                            SCC = this->give(c2_SCCC, gp) * ( 1. - ( eepr - this->give(c2_EPP, gp) ) /
                                                             ( this->give(c2_EPU, gp) - this->give(c2_EPP, gp) ) );
                        }

                        //
                        //   When actual strength is lowered the eff. plastic strain
                        //   increases though no additional strain occurred
                        status->giveTempMaxEffPlasticStrain() = eep - ( status->giveTempCurrentPressureStrength() - SCC ) /
                                                                this->give(c2_E, gp);
                        status->giveTempCurrentPressureStrength() = SCC;
                    }
                }
            }
        }
    }
}





void
Concrete2 :: updateStirrups(GaussPoint *gp, FloatArray &strainIncrement, TimeStep *tStep)
// stirr (double dez, double srf)
//
//
//     EVALUATES THE PLASTIC OR VISCOPLSTIC STRESS- STRAIN
//     RELATIONS FOR STIRRUPS AND UPDATES DATA IN MAT STATUS AND GP DATA
//     ACORDINGLY.
//
//     INPUT
//     strainincrement->at(1) = DEZ - INCREMENT OF STRAIN
//     gp->status->giveCurrentStressInStirrups() - CURENT STRESS, ALSO OUTPUT (INCREMENTED)
//
//     this has - MATERIAL CONSTANTS ARRAY:
//     STIRR_E  - YOUNG MOD.
//     STIRR_R  - UNIAX. STRENGTH=ELAST. LIMIT
//     STIRR_A  - STIRRUPS AREA/UNIT LENGHT (BEAM) OR /UNIT AREA (SHELL)
//     STIRR_TOL- TOLERANCE OF THE EQUIL. IN THE Z-DIRECTION (=0 NO ITER.)
//     STIRR_EREF - EREF - REFRENCE STRAIN RATE FOR PERZYNA'S MATER.
//     STIRR_LAMBDA- COEFF. FOR THAT MATER:
//     SHTIRR_H- ISOTROPIC HARDENING FACTOR
//
//     OUTPUT
//         SRF - REAL STIRRUP STRESS IF IFAD=>1, OTHERWISE STRESS
//         BEFORE INCR.
//     AT PRESENT NO HARDENING AVAILABLE
//
//     PLASTIC STRAIN INCREMENT
//
{
    Concrete2MaterialStatus *status = static_cast< Concrete2MaterialStatus * >( this->giveStatus(gp) );

    double dep, srf, dez, ovs, s, dt;

    dez = strainIncrement.at(1);
    srf = status->giveTempCurrentStressInStirrups();

    dep = 0.;
    ovs = fabs(srf) / this->give(stirr_Ft, gp) - 1.;
    if ( ( ovs ) > 0. ) {
        if ( this->give(stirr_EREF, gp) <= 0. ) {
            if ( ( dez * srf ) > 0. ) {
                dep = dez;
            }
        } else {
            dt = tStep->giveTimeIncrement();
            dep = this->give(stirr_EREF, gp) * exp( ovs / this->give(stirr_LAMBDA, gp) ) * dt;
        }
    }

    s = srf + this->give(stirr_E, gp) * ( dez - dep );
    if ( this->give(c2_IFAD, gp) ) {
        srf = s;
    }

    status->giveTempCurrentStressInStirrups() = srf;
}



void
Concrete2 :: givePlateLayerStiffMtrx(FloatMatrix &answer,
                                     MatResponseMode rMode,
                                     GaussPoint *gp,
                                     TimeStep *tStep)
//
// This material is currently unable compute material stiffness
// so it uses slave material (linearElasticMaterial ) to perform this work
// but this material is currently assumed to be linear elastic ->
// plasticity is not taken into account !.
//
{
    // error ("givePlateLayerStiffMtrx: unable to compute");
    linearElasticMaterial->givePlateLayerStiffMtrx(answer, rMode, gp, tStep);
}


MaterialStatus *
Concrete2 :: CreateStatus(GaussPoint *gp) const
/*
 * creates new  material status  corresponding to this class
 */
{
    Concrete2MaterialStatus *status;

    status = new Concrete2MaterialStatus(1, this->giveDomain(), gp);
    return status;
}



Concrete2MaterialStatus :: Concrete2MaterialStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrainVector(), plasticStrainIncrementVector()
    //
    // constructor
    //
{
    SCCM = EPM = E0PM = SRF = SEZ = 0.0;
    SCTM = -1.0;     // init status if SCTM < 0.;
}



Concrete2MaterialStatus :: ~Concrete2MaterialStatus()
//
// DEstructor
//
{ }

contextIOResultType
Concrete2MaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(SCCM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(EPM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(SCTM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(E0PM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(SRF) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(SEZ) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = plasticStrainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // return result back
    return CIO_OK;
}


contextIOResultType
Concrete2MaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restore state variables from stream
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(SCCM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(EPM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(SCTM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(E0PM) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(SRF) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(SEZ) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = plasticStrainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // return result back
    return CIO_OK;
}



void
Concrete2MaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    StructuralMaterialStatus :: initTempStatus();

    tempSCCM = SCCM;
    tempEPM  = EPM;
    tempSCTM = SCTM;
    tempE0PM = E0PM;
    tempSRF  = SRF;
    tempSEZ  = SEZ;

    if ( plasticStrainVector.giveSize() == 0 ) {
        plasticStrainVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrainVector.zero();
    }

    if ( plasticStrainIncrementVector.giveSize() == 0 ) {
        plasticStrainIncrementVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        plasticStrainIncrementVector.zero();
    } else {
        plasticStrainIncrementVector.zero();
    }
}



void
Concrete2MaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    StructuralMaterialStatus :: updateYourself(tStep);

    SCCM = tempSCCM;
    EPM  = tempEPM;
    SCTM = tempSCTM;
    E0PM = tempE0PM;
    SRF  = tempSRF;
    SEZ  = tempSEZ;

    plasticStrainVector = plasticStrainIncrementVector;
}
} // end namespace oofem
