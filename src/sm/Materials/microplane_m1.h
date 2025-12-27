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

#ifndef microplane_m1_h
#define microplane_m1_h

// by commenting this out we can switch back to the old implementation (2D version not inheriting from MicroplaneMaterial)
#define microplane_m1_new_implementation

#ifdef microplane_m1_new_implementation
// ========================= new implementation =========================

 #include "structuralms.h"
 #include "microplanematerial.h"
 #include "gausspoint.h"

///@name Input fields for M1Material
//@{
 #define _IFT_M1Material_Name "microplane_m1"
 #define _IFT_M1Material_s0 "s0"
 #define _IFT_M1Material_hn "hn"
 #define _IFT_M1Material_talpha "talpha"
//@}

namespace oofem {
class M1MaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray sigN, tempSigN, epspN, tempEpspN;
    IntArray plasticState;

public:
    M1MaterialStatus(GaussPoint *g);

    const char *giveClassName() const override { return "M1MaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letTempNormalMplanePlasticStrainsBe(FloatArray epsilonpN) { tempEpspN =  std :: move(epsilonpN); }
    void letPlasticStateIndicatorsBe(IntArray plSt) { plasticState =  std :: move(plSt); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplanePlasticStrains() { return epspN; }
    const FloatArray &giveTempNormalMplanePlasticStrains() { return tempEpspN; }
    const IntArray &givePlasticStateIndicators() { return plasticState; }
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};

/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
class M1Material : public MicroplaneMaterial
{
protected:
    double EN = 0.; // normal microplane elastic modulus
    double ENtan = 0.; // normal microplane tangent (elastoplastic) modulus
    double HN = 0.; // normal microplane hardening/softening modulus
    double s0 = 0.; // normal microplane initial yield stress

public:
    /**
     * Constructor. Creates M1Material belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M1Material(int n, Domain *d);

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "M1Material"; }
    const char *giveInputRecordName() const override { return _IFT_M1Material_Name; }
    void initializeFrom(InputRecord &ir) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<M1MaterialStatus>(gp); }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem

#else
// ========================= old implementation =========================

 #include "structuralms.h"
 #include "structuralmaterial.h"
 #include "gausspoint.h"


///@name Input fields for M1Material
//@{
 #define _IFT_M1Material_Name "microplane_m1"
 #define _IFT_M1Material_e "e"
 #define _IFT_M1Material_s0 "s0"
 #define _IFT_M1Material_hn "hn"
 #define _IFT_M1Material_talpha "talpha"
 #define _IFT_M1Material_nmp "nmp"
//@}

namespace oofem {
/**
 * Status of simple microplane model - version M1, just with normal microplane strains.
 */
class M1MaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray sigN, tempSigN, sigNyield;

public:
    M1MaterialStatus(GaussPoint *g);

    // definition
    const char *giveClassName() const override { return "M1MaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letNormalMplaneYieldStressesBe(FloatArray sigmaNyield) { sigNyield =  std :: move(sigmaNyield); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplaneYieldStresses() { return sigNyield; }
    void initTempStatus() override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
class M1Material : public StructuralMaterial
{
protected:
    double E = 0.; // Young's modulus
    double nu = 0.; // Poisson ratio
    double EN = 0.; // normal microplane elastic modulus
    double HN = 0.; // normal microplane hardening/softening modulus
    double s0 = 0.; // normal microplane initial yield stress
    int nmp = 0; // number of microplanes
    FloatMatrix n; // microplane normals
    FloatMatrix N; // N = n x n in Voigt notation
    FloatMatrix NN; // NN = n x n x n x n in special notation
    FloatArray mw; // microplane weights

public:
    /**
     * Constructor. Microplane Material M1 belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M1Material(int n, Domain *d);

    void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp,
                                          const FloatArray &reducedStrain, TimeStep *tStep) override;
    void giveElasticPlaneStressStiffMtrx(FloatMatrix &answer) override;
    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    const char *giveClassName() const override { return "M1Material"; }
    const char *giveInputRecordName() const override { return _IFT_M1Material_Name; }
    void initializeFrom(InputRecord &ir) override;
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new M1MaterialStatus(gp); }

protected:
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // end of old implementation
#endif // microplane_m1_h
