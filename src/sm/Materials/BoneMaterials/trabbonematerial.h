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
    double tempAlpha, alpha;
    double tempDam, dam;
    double smtrx, slope;
    double sigC, matConstC;
    FloatArray tempEpsp, epsp, tempDepsp;

public:
    TrabBoneMaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~TrabBoneMaterialStatus();

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    double giveAlpha();
    double giveTempAlpha();
    double giveDam();
    double giveTempDam();
    double giveSmtrx();
    double giveSlope();
    double giveSigC();
    double giveMatConstC();
    const FloatArray &givePlasStrainVector();
    const FloatArray &giveTempPlasStrainVector();
    const FloatArray &giveTempIncPlasStrainVector();

    void setTempAlpha(double al) { tempAlpha = al; }
    void setTempDam(double da) { tempDam = da; }
    void setSmtrx(double smt) { smtrx = smt; }
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
    double E0, Eil, Eie, kie, Ek, Cc, Cc2, EpsC, SigYp, SigYn, adam;

public:
    TrabBoneMaterial(int n, Domain * d);

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    void computeDensification(GaussPoint *gp, const FloatArray &totalStrain);

    double computeDamageParam(double alpha, GaussPoint *gp);

    double computeDamage(GaussPoint *gp, TimeStep *tStep);

    virtual void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep);

    void give1dStressStiffMtrx(FloatMatrix &answer,
                               MatResponseMode mode, GaussPoint *gp,
                               TimeStep *tStep) override;

    void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp,
                                 const FloatArray &reducedStrain, TimeStep *tStep) override;

    int hasMaterialModeCapability(MaterialMode) override;

    const char *giveInputRecordName() const override { return _IFT_TrabBoneMaterial_Name; }
    const char *giveClassName() const override { return "TrabBoneMaterial"; }

    IRResultType initializeFrom(InputRecord *ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif
