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

#ifndef idmgrad1_h
#define idmgrad1_h

#include "sm/Materials/ConcreteMaterials/idm1.h"
#include "sm/Materials/graddpmaterialextensioninterface.h"

#define _IFT_IDGMaterial_Name "idmgrad1"

namespace oofem {
/**
 * Gradient-enhanced Isotropic Damage model for concrete in tension,
 */
class IDGMaterial : public IsotropicDamageMaterial1, GradDpMaterialExtensionInterface
{
protected:
    /// Length scale parameter
    double internalLength;
    /// Parameter specifiying averaging type(0 - classical, 1 - distance based, 2 - stress based)
    int averType;
    /// Parameters specifying how the length scale parameter l is adjusted
    double beta, t;

public:
    /// Constructor
    IDGMaterial(int n, Domain * d);
    /// Destructor
    virtual ~IDGMaterial();

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    const char *giveClassName() const override { return "IDGMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_IDGMaterial_Name; }
    IRResultType initializeFrom(InputRecord *ir) override;

    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradDpMaterialExtensionInterfaceType ) {
            return static_cast< GradDpMaterialExtensionInterface * >(this);
        } else {
            return nullptr;
        }
    }
    int hasMaterialModeCapability(MaterialMode mode) override;

    void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    void giveStiffnessMatrix(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode, GaussPoint *gp,  TimeStep *tStep) override;
    void give1dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void give1dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    void giveInternalLengthDerivative(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
};


/**
 * Material status for gradient-enhanced isotropic damage model for concrete in tension.
 */
class IDGMaterialStatus : public IsotropicDamageMaterial1Status, GradDpMaterialStatusExtensionInterface
{
public:
    IDGMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~IDGMaterialStatus();

    const char *giveClassName() const override { return "IDGMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveNonlocalCumulatedStrain() override { return nonlocalCumulatedStrain; }
    void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) override { this->nonlocalCumulatedStrain = nonlocalCumulatedStrain; }
};
} // end namespace oofem
#endif // idmgrad1_h
