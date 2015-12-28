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


#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "Elements/Shells/shell7base.h"
#include "intmatbilinczfagerstromrate.h"
//#include "vld.h"

namespace oofem {

    REGISTER_Material( IntMatBilinearCZFagerstromRate );

IntMatBilinearCZFagerstromRate :: IntMatBilinearCZFagerstromRate(int n, Domain *d) : IntMatBilinearCZFagerstrom(n, d)
{
    // constructor
}


IntMatBilinearCZFagerstromRate :: ~IntMatBilinearCZFagerstromRate()
{
    // destructor
}

void 
IntMatBilinearCZFagerstromRate :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &d,
                                                     const FloatMatrix &F, TimeStep *tStep)
{
    // returns real stress vector in 3d stress space of receiver according to
    // previous level of stress and current
    // strain increment, the only way, how to correctly update gp records
    //
    
    IntMatBilinearCZFagerstromStatus *status = static_cast< IntMatBilinearCZFagerstromStatus * >( this->giveStatus(gp) );

    FloatMatrix Finv;

    Finv.beInverseOf(F);
    status->letTempInverseDefGradBe(Finv);

    FloatArray dJ;
    dJ.beProductOf(Finv,d);
    status->letTempMaterialJumpBe(dJ);

    double oldDamage = status->giveDamage();
    double dAlpha = 0.0;
    double dt = tStep->giveTimeIncrement();

    FloatArray Qold(3), Qtemp(3), Qtemp_comp(3);
    Qtemp.zero();
    Qtemp_comp.zero();

    if ( dJ.at(3) < 0 ) {
        Qtemp_comp.at(3) = this->knc*dJ.at(3);
    }


    // SUBROUTINE stress_damage_XFEM_direct_Mandel(dJ,N,Qold,old_alpha,Fci,Q,Ea,new_alpha,adtim,dalpha_new,dalpha_old,diss,sig_f,fall);
    if (oldDamage < 0.99) {
        FloatArray Qtrial = status->giveEffectiveMandelTraction(); 

        FloatMatrix Kstiff(3,3);
        FloatArray help;

        Kstiff.zero();
        Kstiff.at(1,1) = this->ks0;
        Kstiff.at(2,2) = this->ks0;
        Kstiff.at(3,3) = this->kn0;

#if 0
        if (dJ.at(3)>=0) {
            Kstiff.at(3,3) = this->kn0;
        } else {
            Kstiff.at(3,3) = this->kn0/(1-oldDamage);
        }
#endif

        dJ.subtract(status->giveOldMaterialJump()); ///@todo Martin: check with Mikael/Jim
        //dJ.rotatedWith(Rot,'n');
        //dJ.printYourself();

        help.beProductOf(Kstiff, dJ);
        Qtrial.add(help);

        double Qn = Qtrial.at(3);
        FloatArray QN(3);
        QN.zero();
        QN.at(3) = Qn;


        FloatArray QtrialShear = Qtrial;
        QtrialShear.subtract(QN);

        double Qt = QtrialShear.computeNorm();


        //double S = this->GIc/this->sigf;
        double sigf = this->sigf;
        double gamma = this->gamma;
        //double gammaGf = this->GIIc/this->GIc;
        double mu = this->mu;
        double c_star = this->c_star;
        double m = this->m;
        //double S = this->GIc/this->sigf;
        double S = (this->GIc - pow(sigf,2)/this->kn0)/sigf;
        double gammaGf = gamma*(this->GIIc - pow(gamma*sigf,2)/this->ks0)/(this->GIc - pow(sigf,2)/this->kn0);



        double Qn_M = 0.5*(Qn + fabs(Qn));


        //double loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;    ///@todo Martin: no use of parameter mu!!!
        double loadFun = sigf*pow(Qt/(gamma*sigf),2) + (sigf/gamma)*(gamma-mu)*pow((Qn_M/sigf),2) - 1/gamma*(gamma*sigf-mu*Qn); // Added support for parameter mu

        Qold = status->giveEffectiveMandelTraction();  ///@todo Martin: check with Mikael/Jim
        Qtemp = Qtrial;

        //Qold.rotatedWith(Rot,'n');

        //double alphaOld = status->giveDamage();

        //if (alphaOld>0.1) {
        //    double bbb=1; ///@todo Should this be used for anything Martin? /JB
        //}

        if (loadFun/sigf < 0.0000001) {
            dAlpha = 0.0;   // new_alpha=old_alpha
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp.times(1-oldDamage);
        } else {
            // dalpha = datr
            double Qt1,Qt2;
            Qt1 = Qtemp.at(1);
            Qt2 = Qtemp.at(2);

            FloatArray M(3); // M = (/2*Q_t/(sig_f*gamma**2), 2*Q_n/sig_f, 0.d0/)
            M.at(1) = 2*Qt1/(pow(gamma,2)*sigf); // Qt = sqrt(Qt1^2 + Qt2^2)
            M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
            M.at(3) = 2*Qn_M/sigf;

            FloatArray dJElastic;
            dJElastic = dJ; // dJtn_e(1:2) = dJtn_v(1:2) - S*dalpha*M(1:2)
            help = M;

            //dAlpha = (dt/(c_star*S))*pow(loadFun_M/sigf,m)*(1/M_norm);
            help.times(S*dAlpha);
            dJElastic.subtract(help);

            FloatMatrix Smat(4,4), Smati(4,4);
            FloatArray R(4), inc(4), dJp(3);
            double loadFunDyn = loadFun;

            const double errorTol = 0.0001;
            for( int iter = 1; fabs(loadFunDyn)/sigf > errorTol; iter++) {
                //printf("loadfun = %e \n",loadFun);
                //double loadFun_M = 0.5*(loadFun + fabs(loadFun));
                //double M_norm = sqrt(M.computeSquaredNorm());

                if (iter>40) {
                    OOFEM_ERROR("BilinearCZMaterialFagerstrom :: giveRealStressVector - no convergence in constitutive driver");
                }
                status->letTempDamageDevBe(true);
                Smat.zero(); // S_mat=0.d0

                R.at(1) = dJElastic.at(1) - (dJ.at(1) - gammaGf*S*dAlpha*M.at(1)); // R(1:2) = dJtn_e(1:2) - (dJtn_v(1:2) - S*dalpha*M(1:2))
                R.at(2) = dJElastic.at(2) - (dJ.at(2) - gammaGf*S*dAlpha*M.at(2));
                R.at(3) = dJElastic.at(3) - (dJ.at(3) - S*dAlpha*M.at(3));
                if (fabs(c_star)>0) {
                    //R.at(4) = loadFun - c_star*pow((dAlpha/dt),m);
                    //R.at(4) = loadFun - c_star*pow((sqrt(dJ.computeSquaredNorm())/dt),m);
                    dJp = dJ;
                    dJp.subtract(dJElastic);
                    R.at(4) = loadFun - c_star*pow((sqrt(dJp.computeSquaredNorm())/dt),m);
                    //R.at(4) = dAlpha - (dt/(c_star*S))*pow(loadFun_M/sigf,m)*(1/M_norm);
                } else {
                    R.at(4) = loadFun;  // R(3) = F/sig_f
                }

                // dMdJtn_e(1,1:2) = (/2/(sig_f*gamma**2),0.d0/)
                // IF (Q_nM>0) THEN
                // dMdJtn_e(2,1:2) = (/0.d0,2/sig_f/)
                // ELSE
                // dMdJtn_e(2,1:2) = (/0.d0,0.d0/)
                // END IF

                // dMdJtn_e = MATMUL(dMdJtn_e,Keye3(1:2,1:2))

                Smat.at(1,1) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(1,1)/(pow(gamma,2)*sigf);  // S_mat(1:2,1:2) = eye3(1:2,1:2)+ dalpha*S*dMdJtn_e(1:2,1:2)
                Smat.at(2,2) = 1.0 + gammaGf*dAlpha*S*2*Kstiff.at(2,2)/(pow(gamma,2)*sigf);
                //dJElastic.printYourself();
                if (Qn_M>0) {
                    Smat.at(3,3) = 1.0 + dAlpha*S*2*Kstiff.at(3,3)/(sigf);
                } else {
                    Smat.at(3,3) = 1.0;
                }

                Smat.at(1,4) = gammaGf*S*M.at(1);   // S_mat(1:2,3) = S*M(1:2)
                Smat.at(2,4) = gammaGf*S*M.at(2);
                Smat.at(3,4) = S*M.at(3);

                double fac1;
                if (sqrt(dJp.computeSquaredNorm())>0) {
                    fac1 = c_star*m/(pow(dt,m))*pow(sqrt(dJp.computeSquaredNorm()),m-2);
                } else {
                    fac1 = 0;
                }


                if ( fabs(c_star) > 0 ) { // Rate-dependent case
                    Smat.at(4,1) = M.at(1)*Kstiff.at(1,1) + fac1*dJp.at(1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
                    Smat.at(4,2) = M.at(2)*Kstiff.at(2,2) + fac1*dJp.at(2);
                    Smat.at(4,3) = (2*(gamma-mu)/(gamma*sigf)*Qn_M + mu/gamma)*Kstiff.at(3,3) + fac1*dJp.at(3); //Martin: included parameter mu
                    //Smat.at(4,4) = -c_star*m*pow((dAlpha/dt),m-1);
                    //Smat.at(4,1) = -fac1*(fac2*M.at(1)*Kstiff.at(1,1) - fac3*M.at(1)*Kstiff.at(1,1));
                    //Smat.at(4,2) = -fac1*(fac2*M.at(2)*Kstiff.at(2,2) - fac3*M.at(2)*Kstiff.at(2,2));
                    //Smat.at(4,1) = -fac1*(fac2*M.at(3)*Kstiff.at(3,3) - fac4*M.at(3)*Kstiff.at(3,3));
                    //Smat.at(4,4) = 1.0;
                } else { // Rate-independent case
                    Smat.at(4,1) = M.at(1)*Kstiff.at(1,1);      // S_mat(3,1:2) = MATMUL(M(1:2),Keye3(1:2,1:2))
                    Smat.at(4,2) = M.at(2)*Kstiff.at(2,2);
                    //Smat.at(4,3) = M.at(3)*Kstiff.at(3,3);
                    Smat.at(4,3) = (2*(gamma-mu)/(gamma*sigf)*Qn_M + mu/gamma)*Kstiff.at(3,3);  //Martin: included parameter mu
                }

                //FloatArray inc;
                //bool transpose = false;
                //Smat.SolveforRhs(R, inc, transpose);

                Smati.beInverseOf(Smat);

                inc.beProductOf(Smati,R);

                dJElastic.at(1) = dJElastic.at(1) - inc.at(1);
                dJElastic.at(2) = dJElastic.at(2) - inc.at(2);
                dJElastic.at(3) = dJElastic.at(3) - inc.at(3);
                //dJElastic.printYourself();
                dAlpha = dAlpha - inc.at(4);

                Qtemp.at(1) = Qold.at(1) + Kstiff.at(1,1)*dJElastic.at(1);
                Qtemp.at(2) = Qold.at(2) + Kstiff.at(2,2)*dJElastic.at(2);
                Qtemp.at(3) = Qold.at(3) + Kstiff.at(3,3)*dJElastic.at(3);

                Qt1 = Qtemp.at(1);
                Qt2 = Qtemp.at(2);
                Qt = sqrt(pow(Qt1,2) + pow(Qt2,2));
                Qn = Qtemp.at(3);                               // Martin: included parameter mu
                Qn_M = 0.5*(Qtemp.at(3)+fabs(Qtemp.at(3)));

                M.at(1) = 2*Qt1/(pow(gamma,2)*sigf);    // Qt = sqrt(Qt1^2 + Qt2^2)
                M.at(2) = 2*Qt2/(pow(gamma,2)*sigf);
                M.at(3) = 2*Qn_M/sigf;

                //loadFun = sigf*pow(Qt/(gamma*sigf),2) + sigf*pow((Qn_M/sigf),2) - sigf;
                loadFun = sigf*pow(Qt/(gamma*sigf),2) + (sigf/gamma)*(gamma-mu)*pow((Qn_M/sigf),2) - 1/gamma*(gamma*sigf-mu*Qn);    //Martin: included parameter mu
                loadFunDyn = loadFun - c_star*pow((sqrt(dJp.computeSquaredNorm())/dt),m);;
            }
            //printf("converged, loadfun = %e, oldDamage = %e \n",loadFun, oldDamage);
            if (oldDamage + dAlpha > 1) {
                dAlpha = 1. - oldDamage;
            }

            FloatMatrix Iep(3,3);
            IntArray Indx(3);
            Indx.at(1) = 1;
            Indx.at(2) = 2;
            Indx.at(3) = 3;

            Iep.beSubMatrixOf(Smati,Indx,Indx);
            status->letTempIepBe(Iep);

            FloatArray alpha_v(3);
            alpha_v.at(1) = Smati.at(4,1);  // alpha_v(1:2) = S_mati(3,1:2)
            alpha_v.at(2) = Smati.at(4,2);
            alpha_v.at(3) = Smati.at(4,3);

            status->letTempAlphavBe(alpha_v);
            status->letTempEffectiveMandelTractionBe(Qtemp);
            Qtemp.times(1-oldDamage-dAlpha);

        }

        //Qtemp.rotatedWith(Rot,'t');           // Q=Qe
        //status->letTempRotationMatrix(Rot);
    } else {
        dAlpha = 1.0 - oldDamage;
        dJ = status->giveTempJump();
        status->letTempEffectiveMandelTractionBe(Qtemp);  // SHOULD NEVER BE USED!!
        //if (dJ.at(3)<0) {
        //    Qtemp.at(3) = kn0*dJ.at(3);
        //}
    }

    Qtemp.at(1) = Qtemp.at(1) + Qtemp_comp.at(1);
    Qtemp.at(2) = Qtemp.at(2) + Qtemp_comp.at(2);
    Qtemp.at(3) = Qtemp.at(3) + Qtemp_comp.at(3);


    answer.beTProductOf(Finv,Qtemp);            // t_1_hat = MATMUL(TRANSPOSE(Fci),Q)
//    answer.times(1-oldDamage-dAlpha);         // t1_s = (1-al)*t_1_hat


    status->letTempDamageBe(oldDamage + dAlpha);
//    status->letTempEffectiveMandelTractionBe(Qtemp);  // NEW!

    status->letTempJumpBe(d);
    status->letTempFirstPKTractionBe(answer);
    status->letTempFBe(F);
}


const double tolerance = 1.0e-12; // small number
IRResultType
IntMatBilinearCZFagerstromRate :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kn0, _IFT_IntMatBilinearCZFagerstrom_kn);
    this->knc = kn0;                        // Defaults to the same stiffness in compression and tension
    IR_GIVE_OPTIONAL_FIELD(ir, this->knc, _IFT_IntMatBilinearCZFagerstrom_knc);

    this->ks0 = kn0;                        // Defaults to kn0
    IR_GIVE_OPTIONAL_FIELD(ir, ks0, _IFT_IntMatBilinearCZFagerstrom_ks);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatBilinearCZFagerstrom_g1c);

    this->GIIc = GIc;                       //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatBilinearCZFagerstrom_g2c);

    IR_GIVE_FIELD(ir, sigf, _IFT_IntMatBilinearCZFagerstrom_sigf);

    IR_GIVE_FIELD(ir, mu, _IFT_IntMatBilinearCZFagerstrom_mu);

    IR_GIVE_FIELD(ir, gamma, _IFT_IntMatBilinearCZFagerstrom_gamma);

    IR_GIVE_FIELD(ir, c_star, _IFT_IntMatBilinearCZFagerstromRate_cstar);

    IR_GIVE_FIELD(ir, m, _IFT_IntMatBilinearCZFagerstromRate_m);

    StructuralInterfaceMaterial ::initializeFrom(ir);

    this->checkConsistency();                                // check validity of the material paramters
    this->printYourself();
    return IRRT_OK;
}

void IntMatBilinearCZFagerstromRate :: giveInputRecord(DynamicInputRecord &input)
{
    IntMatBilinearCZFagerstrom::giveInputRecord(input);

    input.setField(c_star, _IFT_IntMatBilinearCZFagerstromRate_cstar);
    input.setField(knc, _IFT_IntMatBilinearCZFagerstromRate_m);

}

int
IntMatBilinearCZFagerstromRate :: checkConsistency()
{
    if ( this->kn0 < 0.0 ) {
        OOFEM_ERROR("stiffness kn0 is negative (%.2e)", this->kn0);
    } else if ( this->ks0 < 0.0 ) {
        OOFEM_ERROR("stiffness ks0 is negative (%.2e)", this->ks0);
    } else if ( this->GIc < 0.0 ) {
        OOFEM_ERROR("GIc is negative (%.2e)", this->GIc);
    } else if ( this->GIIc < 0.0 ) {
        OOFEM_ERROR("GIIc is negative (%.2e)", this->GIIc);
    } else if ( this->gamma < 0.0  ) { 
        OOFEM_ERROR("gamma (%.2e) is below zero which is unphysical",
            this->gamma);
    } 
    return 1;
}

void
IntMatBilinearCZFagerstromRate  :: printYourself()
{
    IntMatBilinearCZFagerstrom  :: printYourself();
    printf("-Rate parameters \n");
    printf("  c_star  = %e \n", this->c_star);
    printf("  m  = %e \n", this->m);

}


} // end namespace oofem
