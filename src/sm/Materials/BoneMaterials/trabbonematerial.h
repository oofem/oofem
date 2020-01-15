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

#ifndef trabbonematerial_h
#define trabbonematerial_h

#include "sm/Materials/structuralmaterial.h"
#include "floatarrayf.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "sm/Materials/structuralms.h"
#include "cltypes.h"

///@name Input fields for TrabBoneMaterial
//@{
#define _IFT_TrabBoneMaterial_Name "trabbone"
#define _IFT_TrabBoneMaterial_E0 "e0"
#define _IFT_TrabBoneMaterial_Eil "eil"
#define _IFT_TrabBoneMaterial_Eie "eie"
#define _IFT_TrabBoneMaterial_kie "kie"
#define _IFT_TrabBoneMaterial_Ek "ek"
#define _IFT_TrabBoneMaterial_Cc "cc"
#define _IFT_TrabBoneMaterial_Cc2 "cc2"
#define _IFT_TrabBoneMaterial_EpsC "epsc"
#define _IFT_TrabBoneMaterial_SigYp "sigyp"
#define _IFT_TrabBoneMaterial_SigYn "sigyn"
#define _IFT_TrabBoneMaterial_adam "adam"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to TrabBoneMaterial.
 */
class TrabBoneMaterialStatus : public StructuralMaterialStatus
{
protected:
    double tempAlpha = 0., alpha = 0.;
    double tempDam = 0., dam = 0.;
    double slope = 0.;
    double sigC = 0., matConstC = 0.;

    FloatArrayF<1> tempEpsp, epsp, tempDepsp;

public:
    TrabBoneMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    double giveAlpha() const { return alpha; }
    double giveTempAlpha() const { return tempAlpha; }
    double giveDam() const { return dam; }
    double giveTempDam() const { return tempDam; }
    double giveSlope() const { return slope; }
    double giveSigC() const { return sigC; }
    double giveMatConstC() const { return matConstC; }
    const FloatArrayF<1> &givePlasStrainVector() const { return epsp; }
    const FloatArrayF<1> &giveTempPlasStrainVector() const { return tempEpsp; }
    const FloatArrayF<1> &giveTempIncPlasStrainVector() const { return tempDepsp; }

    void setTempAlpha(double al) { tempAlpha = al; }
    void setTempDam(double da) { tempDam = da; }
    void setSlope(double slp) { slope = slp; }
    void setSigC(double sc) { sigC = sc; }
    void setMatConstC(double mcc) { matConstC = mcc; }
    void setTempEpsp(double epsip) { tempEpsp.at(1) = epsip; }
    void setTempDepsp(double depsip) { tempDepsp.at(1) = depsip; }

    const char *giveClassName() const override { return "TrabBoneMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Trabecular bone material model.
 */
class TrabBoneMaterial : public StructuralMaterial
{
protected:
    double E0 = 0., Eil = 0., Eie = 0.;
    double kie = 0., Ek = 0.;
    double Cc = 0., Cc2 = 0., EpsC = 0.;
    double SigYp = 0., SigYn = 0.;
    double adam = 0.;

public:
    TrabBoneMaterial(int n, Domain * d);

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain) const;

    void computeDensification(GaussPoint *gp, const FloatArray &totalStrain) const;

    double computeDamageParam(double alpha, GaussPoint *gp) const;

    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;

    virtual double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp,
                               TimeStep *tStep) const override;

    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode) const override;

    const char *giveInputRecordName() const override { return _IFT_TrabBoneMaterial_Name; }
    const char *giveClassName() const override { return "TrabBoneMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif
