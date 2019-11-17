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

#ifndef m4_h
#define m4_h

#include "microplanematerial_bazant.h"
#include "sm/Materials/structuralms.h"

///@name Input fields for M4Material
//@{
#define _IFT_M4Material_Name "microplane_m4"
#define _IFT_M4Material_c3 "c3"
#define _IFT_M4Material_c4 "c4"
#define _IFT_M4Material_c20 "c20"
#define _IFT_M4Material_k1 "k1"
#define _IFT_M4Material_k2 "k2"
#define _IFT_M4Material_k3 "k3"
#define _IFT_M4Material_k4 "k4"
#define _IFT_M4Material_talpha "talpha"
//@}

namespace oofem {
/**
 * Related material model status to M4Material class
 * for storing history  variables in particular integration point
 * (microplane). Unique copy for each microplane must exist.
 */
class M4MaterialStatus : public StructuralMaterialStatus
{
protected:
    std::vector<MicroplaneState> microplaneStrain, tempMicroplaneStrain;
    std::vector<MicroplaneState> microplaneStress, tempMicroplaneStress;

public:
    M4MaterialStatus(GaussPoint *g, int nplanes);

    const char *giveClassName() const override { return "M4MaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    const MicroplaneState &giveMicroplaneStrain(int mnumber) const { return microplaneStrain[mnumber-1]; }
    const MicroplaneState &giveMicroplaneStress(int mnumber) const { return microplaneStress[mnumber-1]; }

    const MicroplaneState &giveTempMicroplaneStress(int mnumber) const { return tempMicroplaneStress[mnumber-1]; }

    void letTempMicroplaneStrainBe(int mnumber, const MicroplaneState &state) { tempMicroplaneStrain[mnumber-1] = state; }
    void letTempMicroplaneStressBe(int mnumber, const MicroplaneState &state) { tempMicroplaneStress[mnumber-1] = state; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Implementation of microplane material model according to Bazant's boundary curve
 * approach.
 */
class M4Material : public MicroplaneMaterial_Bazant
{
protected:
    double c1 = 0., c2 = 0., c3 = 0., c4 = 0., c5 = 0.;
    double c6 = 0., c7 = 0., c8 = 0., c9 = 0., c10 = 0.;
    double c11 = 0., c12 = 0., c13 = 0., c14 = 0., c15 = 0.;
    double c16 = 0., c17 = 0., c18 = 0., c19 = 0., c20 = 0.; /*c... fixed empirical constants*/
    double k1 = 0., k2 = 0., k3 = 0., k4 = 0., k5 = 0., mu = 0.;
    double EV = 0., ED = 0., ET = 0.;
    double talpha = 0.;

public:
    /**
     * Constructor. Creates  Bazant's Boundary Curve Microplane Material belonging
     * to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M4Material(int n, Domain *d);

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    MicroplaneState giveRealMicroplaneStressVector(GaussPoint *gp, int mnumber, const MicroplaneState &strain, TimeStep *tStep) const override;

    static double FVplus(double ev, double k1, double c13, double c14, double c15, double Ev);
    static double FVminus(double ev, double k1, double k3, double k4, double E);
    static double FDminus(double ed, double k1, double c7, double c8, double c9, double E);
    static double FDplus(double ed, double k1, double c5, double c6, double c7, double c20, double E);
    static double FN(double en, double sv, double k1, double c1, double c2, double c3, double c4,
              double E, double Ev);
    static double FT(double sn, double ev, double k1, double k2, double c10,
              double c11, double c12, double Et);

    void updateVolumetricStressTo(GaussPoint *gp, int mnumber, double sigv) const override;

    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_M4Material_Name; }
    const char *giveClassName() const override { return "M4Material"; }

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new M4MaterialStatus(gp, numberOfMicroplanes); }
};
} // end namespace oofem
#endif // m4_h
