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

#ifndef trabboneembed_h
#define trabboneembed_h

#include "sm/Materials/structuralmaterial.h"
#include "floatarrayf.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "sm/Materials/structuralms.h"
#include "cltypes.h"

///@name Input fields for TrabBoneEmbed
//@{
#define _IFT_TrabBoneEmbed_Name "trabboneembed"
#define _IFT_TrabBoneEmbed_eps0 "eps0"
#define _IFT_TrabBoneEmbed_nu0 "nu0"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to TrabBoneEmbed.
 */
class TrabBoneEmbedStatus : public StructuralMaterialStatus
{
protected:
    double tempAlpha = 0., alpha = 0.;
    double tempDam = 0., dam = 0.;
    double tempPSED = 0., psed = 0.;
    double tempTSED = 0., tsed = 0.;

    FloatMatrixF<6,6> smtrx, matConstD;
    FloatArrayF<6> densStress, tempPlasDef, plasDef, tempIncPlasDef;

public:
    TrabBoneEmbedStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    double giveTempTSED() const { return tempTSED; }
    void setTempDam(double da) { tempDam = da; }
    void setTempTSED(double tse) { tempTSED = tse; }
    void setTempAlpha(double al) { tempAlpha = al; }
    void setTempPlasDef(const FloatArrayF<6> &epsip) { tempPlasDef = epsip; }

    const FloatArrayF<6> &givePlasDef() const { return plasDef; }

    const char *giveClassName() const override { return "TrabBoneEmbedStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Trabecular bone embedding material model.
 */
class TrabBoneEmbed : public StructuralMaterial
{
protected:
    double eps0 = 0., nu0 = 0.;

public:
    TrabBoneEmbed(int n, Domain * d);

    void performPlasticityReturn(GaussPoint *gp, const FloatArrayF<6> &totalStrain) const;

    double computeDamageParam(double alpha, GaussPoint *gp) const;

    double computeDamage(GaussPoint *gp, TimeStep *tStep) const;

    virtual double computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const;

    /// Constructs the anisotropic compliance tensor.
    static FloatMatrixF<6,6> constructIsoComplTensor(double eps0, double nu0);

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                           TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_TrabBoneEmbed_Name; }
    const char *giveClassName() const override { return "TrabBoneEmbed"; }

    void initializeFrom(InputRecord &ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif
