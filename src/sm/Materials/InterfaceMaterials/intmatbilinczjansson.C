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

#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "intmatbilinczjansson.h"

namespace oofem {

REGISTER_Material( IntMatBilinearCZJansson );


IntMatBilinearCZJansson :: IntMatBilinearCZJansson(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{ }


///@todo - need to rearrange traction and stiffness matrix so the first component is normal
FloatArrayF<3>
IntMatBilinearCZJansson :: giveFirstPKTraction_3d(const FloatArrayF<3> &d, const FloatMatrixF<3,3> &F, GaussPoint *gp, TimeStep *tStep) const
{
    // returns real stress vector in 3d stress space of receiver according to
    // previous level of stress and current
    // strain increment, the only way, how to correctly update gp records

    IntMatBilinearCZJanssonStatus *status = static_cast< IntMatBilinearCZJanssonStatus * >( this->giveStatus(gp) );

    auto Finv = inv(F);
    status->letTempInverseDefGradBe(Finv);

    auto dJ = dot(Finv, d);
    status->letTempMaterialJumpBe(dJ);

    double oldDamage = status->giveDamage();
    double dAlpha = 0.0;
    FloatArrayF<3> Qold, Qtemp, Qtemp_comp;

    if ( dJ.at(3) < 0 ) {
        Qtemp_comp.at(3) = this->knc*dJ.at(3);
    }

    auto Qtrial = status->giveEffectiveMandelTraction(); 

    if ( oldDamage < 0.99 ) {
        //FloatArray Qtrial = status->giveEffectiveMandelTraction(); 


        auto Kstiff = diag<3>({this->ks0, this->ks0, this->kn0});
        //} else {
        //    Kstiff.at(3,3) = this->kn0/(1-oldDamage);
        //}

        dJ -= status->giveOldMaterialJump();

        Qtrial += dot(Kstiff, dJ);
        double Qn = Qtrial.at(3);
        FloatArrayF<3> QtrialShear = {Qtrial[0], Qtrial[1], 0.};

        double Qt = norm(QtrialShear);

//        double S = this->GIc/this->sigf;
        double sigf = this->sigf;
        double gamma = this->gamma;
//        double gammaGf = this->GIIc/this->GIc;
        double mu = this->mu;
        double c = mu/gamma; 

        double Qn_M = 0.5*(Qn + fabs(Qn));
//        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;    ///@todo Martin: no use of parameter mu!!!

        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*(1-c)*pow((Qn_M/sigf),2) + sigf*c*(Qn/sigf) - sigf;

        // completely unused?
        //auto Qold = status->giveEffectiveMandelTraction();  
        Qtemp = Qtrial;

        //Qold.rotatedWith(Rot,'n');

        if ( loadFun/sigf < 0.0000001 ) {
            dAlpha = 0.0;   // new_alpha=old_alpha
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp *= 1-oldDamage;
        } else {
            status->letTempDamageDevBe(true);

             // dalpha = datr
            double C1 = (pow((Qt/gamma),2)+(1-c)*pow(Qn_M,2))/(pow(sigf,2));
            double C2 = c*Qn_M/sigf;

            //double xi = (-C2 + sqrt(pow(C2,2)+(1-c*Qn_M/sigf)*4*C1))/(2*C1);
            double xi = 0.0;

            if ( Qn >= 0 ) {
                xi = (-C2 + sqrt(pow(C2,2)+4*C1))/(2*C1);
            } else {
                if (1-c*Qn/sigf>0) {
                    xi = (sigf*gamma/Qt)*sqrt(1-c*Qn/sigf);
                } else {
                    OOFEM_ERROR("Inconsistent cohesive model specification, 1-c*Qn/sigf =  %e", 1 - c*Qn / sigf);
                }
            }

            Qt = xi*Qt;
            Qn = 0.5*((1+xi)-(1-xi)*sgn(Qn))*Qn;
            Qn_M = 0.5*(Qn + fabs(Qn));

            double beta = pow(Qt,2)*Kstiff.at(3,3)/(pow(Qn_M,2)*Kstiff.at(1,1) + pow(Qt,2)*Kstiff.at(3,3));

            double G_beta = beta*(this->GIIc - pow(gamma*sigf,2)/(2*this->ks0)) + (1-beta)*(this->GIc - pow(sigf,2)/(2*this->kn0)); // assuming linear interpolation between mode I and II


            double eta = (pow(Qn_M,2) + pow(Qt,2)*Kstiff.at(3,3)/Kstiff.at(1,1))/(G_beta*sigf);

            dAlpha = (1/xi-1)*sigf/(2*Kstiff.at(3,3))*eta;

            if ( oldDamage + dAlpha > 1 ) {
                dAlpha = 1-oldDamage;
            }

            double Qt_trial = norm(QtrialShear);

            double Qt1, Qt2;

            if ( Qt_trial > 0 ) {
                Qt1 = Qt*QtrialShear.at(1)/Qt_trial;
                Qt2 = Qt*QtrialShear.at(2)/Qt_trial;
            } else {
                Qt1 = 0.0;
                Qt2 = 0.0;
            }


            FloatArrayF<3> Mstar, M;

            Mstar.at(1) = 2*Qt1*this->kn0/Kstiff.at(1,1)/(sigf);    // Qt = sqrt(Qt1^2 + Qt2^2)
            Mstar.at(2) = 2*Qt2*this->kn0/Kstiff.at(2,2)/(sigf);
            Mstar.at(3) = 2*Qn_M/sigf;

            M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);
            M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
            M.at(3) = 2*(1-c)/(sigf)*Qn_M + c;

            FloatMatrix Smat(4,4), Smati(4,4);
            Smat.zero();    // S_mat=0.d0

            Smat.at(1,1) = 1.0 + dAlpha*(1/eta)*2*this->kn0/(sigf);     // S_mat(1:2,1:2) = eye3(1:2,1:2)+ dalpha*S*dMdJtn_e(1:2,1:2)
            Smat.at(2,2) = 1.0 + dAlpha*(1/eta)*2*this->kn0/(sigf);

            if ( Qn_M > 0 ) {
                Smat.at(3,3) = 1.0 + dAlpha*(1/eta)*2*Kstiff.at(3,3)/(sigf);
            } else {
                Smat.at(3,3) = 1.0;
            }

            Smat.at(1,4) = (1/eta)*Mstar.at(1);         // S_mat(1:2,3) = S*M(1:2)
            Smat.at(2,4) = (1/eta)*Mstar.at(2);
            Smat.at(3,4) = (1/eta)*Mstar.at(3);

            Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
            Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
            Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);
            //Smat.printYourself();
            Smati.beInverseOf(Smat);

            Qtemp.at(1) = Qt1;
            Qtemp.at(2) = Qt2;
            Qtemp.at(3) = Qn;

            IntArray Indx(3);
            Indx.at(1) = 1;
            Indx.at(2) = 2;
            Indx.at(3) = 3;

            FloatMatrix Iep(3,3);
            Iep.beSubMatrixOf(Smati,Indx,Indx);

            status->letTempIepBe(Iep);

            FloatArray alpha_v(3);
            alpha_v.at(1) = Smati.at(4,1);              // alpha_v(1:2) = S_mati(3,1:2)
            alpha_v.at(2) = Smati.at(4,2);
            alpha_v.at(3) = Smati.at(4,3);
            status->letTempAlphavBe(alpha_v);

            dJ = status->giveTempJump();
            status->letTempEffectiveMandelTractionBe(Qtemp);

            Qtemp *= 1 - oldDamage - dAlpha;

#if 0
            if (dJ.at(3)>=0) {
                status->letTempEffectiveMandelTractionBe(Qtemp);
                Qtemp.times(1-oldDamage-dAlpha);
            } else {
                if (oldDamage + dAlpha<1) {
                    Qtemp.at(3) = (1-oldDamage)/(1-oldDamage + dAlpha)*Qtemp.at(3);
                    status->letTempEffectiveMandelTractionBe(Qtemp);
                    Qtemp.times(1-oldDamage-dAlpha);
                } else {
                    status->letTempEffectiveMandelTractionBe(Qtemp);        // SHOULD NEVER BE USED
                    Qtemp.times(1-oldDamage-dAlpha);
                    Qtemp.at(3) = (1-oldDamage)*(Qold.at(3) + Kstiff.at(3,3)*dJElastic.at(3)); // dJElastic.at(3) is always equal to dJ.at(3) if dJ.at(3)<0
                }
            }
#endif
        }
    } else {
        dAlpha = 1.0 - oldDamage;
        dJ = status->giveTempJump();
        status->letTempEffectiveMandelTractionBe(Qtemp);    // SHOULD NEVER BE USED!!
        //if (dJ.at(3)<0) {
        //    Qtemp.at(3) = kn0*dJ.at(3);
        //}
    }

    Qtemp += Qtemp_comp;

    auto answer = Tdot(Finv, Qtemp);                    // t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
//    answer *= 1-oldDamage-dAlpha;                 // t1_s = (1-al)*t_1_hat


    status->letTempDamageBe(oldDamage + dAlpha);
//    status->letTempEffectiveMandelTractionBe(Qtemp);  // NEW!

    status->letTempJumpBe(d);
    status->letTempFirstPKTractionBe(answer);
    status->letTempFBe(F);

    if ( mSemiExplicit ) {
        Qtemp = (1-oldDamage) * Qtrial + Qtemp_comp;
        answer = Tdot(Finv, Qtemp);                    // t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
    }

    return answer;
}


FloatMatrixF<3,3>
IntMatBilinearCZJansson :: give3dStiffnessMatrix_dTdj(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    IntMatBilinearCZJanssonStatus *status = static_cast< IntMatBilinearCZJanssonStatus * >( this->giveStatus(gp) );
    //this->give3dStiffnessMatrix_dTdj_num(answer, rMode, gp, tStep);
    //OOFEM_WARNING("numerical tangent");
    //answer.printYourself();

    FloatMatrixF<3,3> answer;

    if ( status->giveOldDamageDev() ) {
        answer = status->giveOlddTdJ();
        //answer.printYourself();
        status->letOldDamageDevBe(false);
    } else {
        double damage = status->giveTempDamage();
        const auto &Finv = status->giveTempInverseDefGrad();
        const auto &J = status->giveTempJump();

        //FloatMatrix Rot = status->giveTempRotationMatrix();
        auto Kstiff = diag<3>({this->ks0, this->ks0, this->kn0});
        //Kstiff.rotatedWith(Rot);

        FloatMatrix help;
        if ( damage >= 1.0 ) {
            if ( J.at(3) < 0 ) {
                Kstiff.at(1,1) = 0.0;
                Kstiff.at(2,2) = 0.0;
                Kstiff.at(3,3) = this->knc;

                answer = rotate(Kstiff, Finv);
            }
        } else {
            if ( status->giveTempDamage() - status->giveDamage()==0.0 ) {
                if ( J.at(3) < 0 ) {
                    Kstiff.at(3,3) += this->knc/(1-damage);
                }

                answer = (1-damage) * rotate(Kstiff, Finv);
            } else {
                const auto &Iep = status->giveTempIep();

                //if (J.at(3)<0) {
                //    Kstiff.at(1,1) = (1-damage)*Kstiff.at(1,1);
                //    Kstiff.at(2,2) = (1-damage)*Kstiff.at(2,2);
                //    Kstiff.at(3,3) = Kstiff.at(3,3);
                //} else {
                //    Kstiff.times((1-damage));
                //}

                answer = (1-damage) * dot(Kstiff, Iep);
                //answer.rotatedWith(Rot);                  // Ea_h = MATMUL(TRANSPOSE(Rot),MATMUL(Keye3,Iep))
                                                            // Ea_h = MATMUL(Ea_h,Rot)

                const auto &alpha_v = status->giveTempAlphav();
                //alpha_v.rotatedWith(Rot, 't');            // alpha_v = MATMUL(TRANSPOSE(Rot),alpha_v)

                const auto &Qtemp = status->giveTempEffectiveMandelTraction();

                auto temp1 = Tdot(Finv, Qtemp);             // CALL gmopen33(MATMUL(TRANSPOSE(Fci),Q),MATMUL(alpha_v,Fci),t1halFci_o)
                auto temp2 = Tdot(Finv, alpha_v);

                auto t1hatFinvOpen = dyad(temp1, temp2);

                if ( J.at(3) < 0 ) {
                    answer.at(3,3) += this->knc;
                }

                // Ea = (1-new_alpha)*MATMUL(TRANSPOSE(Fci),MATMUL(Ea_h,Fci)) -&
                answer = rotate(answer, Finv) - t1hatFinvOpen;     // t1halFci_o
            }
        }
    }
    status->letTempdTdJBe(answer);

    //Finv.printYourself();
    //Kstiff.printYourself();
    //printf("analytical tangent \n");
    //answer.printYourself();
    return answer;

}


//const double tolerance = 1.0e-12; // small number
void
IntMatBilinearCZJansson :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZJansson_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZJansson_knc);

    this->ks0 = kn0;                        // Defaults to kn0
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZJansson_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZJansson_g1c);

    this->GIIc = GIc;                       //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatBilinearCZJansson_g2c);

    IR_GIVE_FIELD(ir, sigf, _IFT_IntMatBilinearCZJansson_sigf);

    IR_GIVE_FIELD(ir, mu, _IFT_IntMatBilinearCZJansson_mu);

    IR_GIVE_FIELD(ir, gamma, _IFT_IntMatBilinearCZJansson_gamma);

    // check validity of the material paramters
    if ( this->kn0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZJansson_kn, "must be positive");
    } else if ( this->ks0 < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZJansson_ks, "must be positive");
    } else if ( this->GIc < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZJansson_g1c, "must be positive");
    } else if ( this->GIIc < 0.0 ) {
        throw ValueInputException(ir, _IFT_IntMatBilinearCZJansson_g2c, "must be positive");
    } else if ( this->gamma < 0.0  ) { 
        throw ValueInputException(ir, _IFT_IntMatBilinearCZJansson_gamma, "must be positive");
    }

    if ( ir.hasField(_IFT_IntMatBilinearCZJansson_semiexplicit) ) {
        mSemiExplicit = true;
        printf("In IntMatBilinearCZJansson::initializeFrom: Semi-explicit time integration activated.\n");
    }
}

int
IntMatBilinearCZJansson :: checkConsistency()
{
    return 1;
}

void
IntMatBilinearCZJansson  :: printYourself()
{
    printf("Paramters for BilinearCZMaterial: \n");

    printf("-Strength paramters \n");
    printf("  sigf  = %e \n", this->sigf);
    printf("  GIc   = %e \n", this->GIc);
    printf("  GIIc  = %e \n", this->GIIc);
    printf("  gamma = sigfs/sigfn = %e \n", this->gamma);
    printf("  mu    = %e \n", this->mu);

    printf("-Stiffness parameters \n");
    printf("  kn0  = %e \n", this->kn0);
    printf("  ks0  = %e \n", this->ks0);
    printf("  knc  = %e \n", this->knc);

}


//IntMatBilinearCZJanssonStatus :: IntMatBilinearCZJanssonStatus(GaussPoint *g) : StructuralMaterialStatus(g)
IntMatBilinearCZJanssonStatus :: IntMatBilinearCZJanssonStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{
    tempFInv = eye<3>();

#if 0
    ///@todo Martin: Very bad implementation of intialisation of Rot
    //*************************************************************
    FloatMatrix Gcov;
    FloatArray N;

    Shell7Base *shell = dynamic_cast< Shell7Base* >(gp->giveElement());
    if ( !shell ) {
        OOFEM_ERROR("BilinearCZMaterialJansson :: giveRealStressVector - oh no wrong element type");
    }
    FloatArray lCoords(3);
    lCoords.at(1) = gp->giveCoordinate(1);
    lCoords.at(2) = gp->giveCoordinate(2);
    shell->evalInitialCovarBaseVectorsAt(lCoords, Gcov);

    FloatArray G1, G2;
    Gcov.copyColumn(G1,1);
    Gcov.copyColumn(G2,2);
    N.beVectorProductOf(G1,G2);
    N.normalize();

    tempRot.resize(3,3);
    G1.normalize();
    G2.beVectorProductOf(N, G1);
    tempRot.setColumn(G1,1);
    tempRot.setColumn(G2,2);
    tempRot.setColumn(N,3);
#endif

    Iep = tempFInv;
}


void
IntMatBilinearCZJanssonStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    ///@todo Martin: check need of this
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    /*
    fprintf(file, "status { ");
    if ( this->damage > 0.0 ) {
    fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    }

    fprintf(file, "}\n");
    */
}


void
IntMatBilinearCZJanssonStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();

    tempMaterialJump = oldMaterialJump;
    tempDamage = damage;
    tempQEffective = QEffective;

    tempFInv = eye<3>();

    Iep = tempFInv;
    alphav = oldMaterialJump;

    tempDamageDev = false;
}

void
IntMatBilinearCZJanssonStatus :: updateYourself(TimeStep *atTime)
{
    StructuralInterfaceMaterialStatus :: updateYourself(atTime);
    damage = tempDamage;
    oldMaterialJump = tempMaterialJump;
    QEffective = tempQEffective;

    old_dTdJ = temp_dTdJ;
    oldDamageDev = tempDamageDev;
}


} // end namespace oofem
