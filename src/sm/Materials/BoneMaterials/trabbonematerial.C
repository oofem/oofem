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

#include "trabbonematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "../sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(TrabBoneMaterial);

TrabBoneMaterial :: TrabBoneMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


int
TrabBoneMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _1dMat;
}


void TrabBoneMaterial :: computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep)
{
    TrabBoneMaterialStatus *status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );
    alpha = status->giveTempAlpha();
}


void
TrabBoneMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mode, GaussPoint *gp,
                                          TimeStep *tStep)
{
    TrabBoneMaterialStatus *status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    double epsnew;
    double epsp, depsp;
    double dam, alpha;
    double matconstc;
    FloatArray epsnewArray;

    if ( mode == ElasticStiffness ) {
        answer.resize(1, 1);
        answer.at(1, 1) = E0;
    } else if ( mode == SecantStiffness ) {
        dam = status->giveTempDam();
        matconstc = status->giveMatConstC();

        answer.resize(1, 1);
        answer.at(1, 1) = ( 1.0 - dam ) * E0 + matconstc;
    } else {
        epsnew = status->giveTempStrainVector().at(1);
        epsnewArray = status->giveTempStrainVector();
        epsp = status->giveTempPlasStrainVector().at(1);
        depsp = status->giveTempIncPlasStrainVector().at(1);
        alpha = status->giveTempAlpha();
        dam = status->giveTempDam();

        answer.resize(1, 1);

        computeDensification(gp, epsnewArray);
        matconstc = status->giveMatConstC();

        if ( depsp != 0.0 ) {
            answer.at(1, 1) = ( 1.0 - dam ) * ( E0 * ( Eil + Ek ) ) / ( E0 + Eil + Ek )
                              - E0 * E0 * ( epsnew - epsp ) / ( E0 + Eil + Ek ) * adam * exp(-adam * alpha) * depsp / fabs(depsp) + matconstc;
        } else {
            answer.at(1, 1) = ( 1.0 - dam ) * E0 + matconstc;
        }
    }

    status->setSmtrx( answer.at(1, 1) );
}


void
TrabBoneMaterial :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    double epsnew, epsold;
    double epsp, depsp;
    double alpha;
    double sigp, sigY;
    double gNewton, dgNewton;

    TrabBoneMaterialStatus *status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    epsnew = totalStrain.at(1);
    epsold = status->giveStrainVector().at(1);
    epsp = status->givePlasStrainVector().at(1);
    alpha = status->giveAlpha();
    sigp = E0 * epsnew - ( E0 + Ek ) * epsp;

    if ( sigp < 0.0 ) {
        sigY = SigYn;
    } else {
        sigY = SigYp;
    }

    if ( sigp / fabs(sigp) * sigp > sigY + Eil * alpha + Eie * ( 1 - exp(-kie * alpha) ) ) {
        depsp = epsnew - epsold;
        gNewton = sigp / fabs(sigp) * ( E0 * epsnew - ( E0 + Ek ) * ( epsp + depsp ) ) - ( sigY + Eil * ( alpha +
                                                                                                          sigp / fabs(sigp) * depsp ) + Eie * ( 1 - exp( -kie * ( alpha + sigp / fabs(sigp) * depsp ) ) ) );
        while ( fabs(gNewton) > 10.e-7 ) {
            gNewton = sigp / fabs(sigp) * ( E0 * epsnew - ( E0 + Ek ) * ( epsp + depsp ) ) - ( sigY + Eil * ( alpha +
                                                                                                              sigp / fabs(sigp) * depsp ) + Eie * ( 1 - exp( -kie * ( alpha + sigp / fabs(sigp) * depsp ) ) ) );
            dgNewton = -sigp / fabs(sigp) * ( ( E0 + Ek ) + Eil +
                                             Eie * kie * exp( -kie * ( alpha + sigp / fabs(sigp) * depsp ) ) );
            depsp += -gNewton / dgNewton;
        }

        epsp += depsp;
    } else {
        depsp = 0.0;
    }

    alpha += fabs(depsp);

    status->setTempEpsp(epsp);
    status->setTempDepsp(depsp);
    status->setTempAlpha(alpha);
}


void
TrabBoneMaterial :: computeDensification(GaussPoint *gp, const FloatArray &totalStrain)
{
    double epsnew, sigc;
    double matconstc;

    TrabBoneMaterialStatus *status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    epsnew = totalStrain.at(1);

    if ( epsnew > EpsC ) {
        sigc = 0.0;
        matconstc = 0.0;
    } else {
        sigc = Cc * ( epsnew - EpsC ) + Cc2 * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC );
        matconstc = Cc + 7 * Cc2 * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC ) * ( epsnew - EpsC );
    }

    status->setSigC(sigc);
    status->setMatConstC(matconstc);
}


double
TrabBoneMaterial :: computeDamageParam(double alpha, GaussPoint *gp)
{
    double dam;
    if ( alpha > 0. ) {
        dam = 1.0 - exp(-adam * alpha);
        //    dam = adam*alpha;
    } else {
        dam = 0.0;
    }

    if ( alpha < 0. ) {
        OOFEM_ERROR("Alpha less than zero. Not possible");
    }

    return dam;
}


double
TrabBoneMaterial :: computeDamage(GaussPoint *gp,  TimeStep *tStep)
{
    double alpha;

    computeCumPlastStrain(alpha, gp, tStep);

    double dam = computeDamageParam(alpha, gp);

    return dam;
}


void
TrabBoneMaterial :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp,
                                            const FloatArray &totalStrain,
                                            TimeStep *tStep)
{
    double epsnew, epsp;
    double dam;
    double sig, sigc;

    TrabBoneMaterialStatus *status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    performPlasticityReturn(gp, totalStrain);

    epsnew = totalStrain.at(1);

    computeDensification(gp, totalStrain);

    dam = computeDamage(gp, tStep);

    epsp = status->giveTempPlasStrainVector().at(1);
    sigc = status->giveSigC();
    sig = ( 1 - dam ) * E0 * ( epsnew - epsp ) + sigc;

    answer.resize(1);
    answer.at(1) = sig;
    status->setTempDam(dam);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}


IRResultType
TrabBoneMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Read material properties here

    IR_GIVE_FIELD(ir, E0, _IFT_TrabBoneMaterial_E0);
    IR_GIVE_FIELD(ir, Eil, _IFT_TrabBoneMaterial_Eil);
    IR_GIVE_FIELD(ir, Eie, _IFT_TrabBoneMaterial_Eie);
    IR_GIVE_FIELD(ir, kie, _IFT_TrabBoneMaterial_kie);
    IR_GIVE_FIELD(ir, Ek, _IFT_TrabBoneMaterial_Ek);
    IR_GIVE_FIELD(ir, Cc, _IFT_TrabBoneMaterial_Cc);
    IR_GIVE_FIELD(ir, Cc2, _IFT_TrabBoneMaterial_Cc2);
    IR_GIVE_FIELD(ir, EpsC, _IFT_TrabBoneMaterial_EpsC);
    IR_GIVE_FIELD(ir, SigYp, _IFT_TrabBoneMaterial_SigYp);
    IR_GIVE_FIELD(ir, SigYn, _IFT_TrabBoneMaterial_SigYn);
    IR_GIVE_FIELD(ir, adam, _IFT_TrabBoneMaterial_adam);

    return StructuralMaterial :: initializeFrom(ir);
}


/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////

TrabBoneMaterialStatus :: TrabBoneMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    alpha = 0.0;
    tempAlpha = 0.0;
    dam = 0.0;
    tempDam = 0.0;
    smtrx = 0.0;
    slope = 0.0;
    sigC = 0.0;
    matConstC = 0.0;
    epsp.resize(1);
    tempEpsp.resize(1);
    epsp.at(1) = 0.0;
    tempEpsp.at(1) = 0.0;
    tempDepsp.resize(1);
    tempDepsp.at(1) = 0.0;
}


TrabBoneMaterialStatus :: ~TrabBoneMaterialStatus()
{ }


double
TrabBoneMaterialStatus :: giveAlpha()
{
    return alpha;
}

double
TrabBoneMaterialStatus :: giveTempAlpha()
{
    return tempAlpha;
}

double
TrabBoneMaterialStatus :: giveDam()
{
    return dam;
}

double
TrabBoneMaterialStatus :: giveTempDam()
{
    return tempDam;
}

double
TrabBoneMaterialStatus :: giveSmtrx()
{
    return smtrx;
}

double
TrabBoneMaterialStatus :: giveSlope()
{
    return slope;
}

double
TrabBoneMaterialStatus :: giveSigC()
{
    return sigC;
}

double
TrabBoneMaterialStatus :: giveMatConstC()
{
    return matConstC;
}

const FloatArray &
TrabBoneMaterialStatus :: givePlasStrainVector()
{
    return epsp;
}

const FloatArray &
TrabBoneMaterialStatus :: giveTempPlasStrainVector()
{
    return tempEpsp;
}

const FloatArray &
TrabBoneMaterialStatus :: giveTempIncPlasStrainVector()
{
    return tempDepsp;
}


void
TrabBoneMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "plastrains %f, alpha %f, dam %f, Slope %f, Smtrx %f ", this->tempEpsp.at(1),
            this->tempAlpha, this->tempDam, this->slope, this->smtrx);
    fprintf(file, "}\n");
}


void
TrabBoneMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempAlpha = this->alpha;
    this->tempDam = this->dam;
    this->tempEpsp = this->epsp;
}


void
TrabBoneMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->alpha = this->tempAlpha;
    this->dam = this->tempDam;
    this->epsp = this->tempEpsp;
}


contextIOResultType
TrabBoneMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    //if (fwrite(&kappa,sizeof(double),1,stream) != CIO_OK) THROW_CIOERR(CIO_IOERR);
    //if (fwrite(&damage,sizeof(double),1,stream)!= CIO_OK) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


contextIOResultType
TrabBoneMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if (fread (&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fread (&damage,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


MaterialStatus *TrabBoneMaterial :: CreateStatus(GaussPoint *gp) const
{
    TrabBoneMaterialStatus *status =
        new  TrabBoneMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}
} // end namespace oofem
