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
#include "sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(TrabBoneMaterial);

TrabBoneMaterial :: TrabBoneMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


bool
TrabBoneMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _1dMat;
}


double TrabBoneMaterial :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );
    return status->giveTempAlpha();
}


FloatMatrixF<1,1>
TrabBoneMaterial :: give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    if ( mode == ElasticStiffness ) {
        return {E0};
    } else if ( mode == SecantStiffness ) {
        double dam = status->giveTempDam();
        double matconstc = status->giveMatConstC();

        return {( 1.0 - dam ) * E0 + matconstc};
    } else {
        double epsnew = status->giveTempStrainVector().at(1);
        auto &epsnewArray = status->giveTempStrainVector();
        double epsp = status->giveTempPlasStrainVector().at(1);
        double depsp = status->giveTempIncPlasStrainVector().at(1);
        double alpha = status->giveTempAlpha();
        double dam = status->giveTempDam();
        double matconstc = status->giveMatConstC();

        /// FIXME: this function call modifies the material state (this should never happen in a tangent computation):
        computeDensification(gp, epsnewArray);
        if ( depsp != 0.0 ) {
            return {( 1.0 - dam ) * ( E0 * ( Eil + Ek ) ) / ( E0 + Eil + Ek )
                    - E0 * E0 * ( epsnew - epsp ) / ( E0 + Eil + Ek ) * adam * exp(-adam * alpha) * depsp / fabs(depsp) + matconstc};
        } else {
            return {( 1.0 - dam ) * E0 + matconstc};
        }
    }
}


void
TrabBoneMaterial :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain) const
{
    auto status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    double epsnew = totalStrain.at(1);
    double epsold = status->giveStrainVector().at(1);
    double epsp = status->givePlasStrainVector().at(1);
    double alpha = status->giveAlpha();
    double sigp = E0 * epsnew - ( E0 + Ek ) * epsp;

    double sigY;
    if ( sigp < 0.0 ) {
        sigY = SigYn;
    } else {
        sigY = SigYp;
    }

    double depsp;
    if ( sigp / fabs(sigp) * sigp > sigY + Eil * alpha + Eie * ( 1 - exp(-kie * alpha) ) ) {
        depsp = epsnew - epsold;
        double gNewton = sigp / fabs(sigp) * ( E0 * epsnew - ( E0 + Ek ) * ( epsp + depsp ) ) - ( sigY + Eil * ( alpha +
                                                                                                          sigp / fabs(sigp) * depsp ) + Eie * ( 1 - exp( -kie * ( alpha + sigp / fabs(sigp) * depsp ) ) ) );
        while ( fabs(gNewton) > 10.e-7 ) {
            gNewton = sigp / fabs(sigp) * ( E0 * epsnew - ( E0 + Ek ) * ( epsp + depsp ) ) - ( sigY + Eil * ( alpha +
                                                                                                              sigp / fabs(sigp) * depsp ) + Eie * ( 1 - exp( -kie * ( alpha + sigp / fabs(sigp) * depsp ) ) ) );
            double dgNewton = -sigp / fabs(sigp) * ( ( E0 + Ek ) + Eil +
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
TrabBoneMaterial :: computeDensification(GaussPoint *gp, const FloatArray &totalStrain) const
{
    auto status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );
    double epsnew = totalStrain.at(1);
    double sigc, matconstc;
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
TrabBoneMaterial :: computeDamageParam(double alpha, GaussPoint *gp) const
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
TrabBoneMaterial :: computeDamage(GaussPoint *gp,  TimeStep *tStep) const
{
    double alpha = computeCumPlastStrain(gp, tStep);
    return computeDamageParam(alpha, gp);
}


FloatArrayF<1>
TrabBoneMaterial :: giveRealStressVector_1d(const FloatArrayF<1> &totalStrain,
                                            GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    performPlasticityReturn(gp, totalStrain);

    computeDensification(gp, totalStrain);

    double epsnew = totalStrain.at(1);
    double dam = computeDamage(gp, tStep);
    double epsp = status->giveTempPlasStrainVector().at(1);
    double sigc = status->giveSigC();
    double sig = ( 1 - dam ) * E0 * ( epsnew - epsp ) + sigc;

    FloatArrayF<1> answer = {sig};
    status->setTempDam(dam);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}


void
TrabBoneMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

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
}


/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////

TrabBoneMaterialStatus :: TrabBoneMaterialStatus(GaussPoint *g) : StructuralMaterialStatus(g)
{
}


void
TrabBoneMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "plastrains %f, alpha %f, dam %f, Slope %f ", this->tempEpsp.at(1),
            this->tempAlpha, this->tempDam, this->slope);
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


void
TrabBoneMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    // save parent class status
    StructuralMaterialStatus :: saveContext(stream, mode);

    // write a raw data
    //if (fwrite(&kappa,sizeof(double),1,stream) != CIO_OK) THROW_CIOERR(CIO_IOERR);
    //if (fwrite(&damage,sizeof(double),1,stream)!= CIO_OK) THROW_CIOERR(CIO_IOERR);
}


void
TrabBoneMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    // read parent class status
    StructuralMaterialStatus :: restoreContext(stream, mode);

    // read raw data
    //if (fread (&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fread (&damage,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
}


MaterialStatus *TrabBoneMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new TrabBoneMaterialStatus(gp);
}
} // end namespace oofem
