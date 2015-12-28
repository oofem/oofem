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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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
    M1MaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~M1MaterialStatus();

    // definition
    virtual const char *giveClassName() const { return "M1MaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letTempNormalMplanePlasticStrainsBe(FloatArray epsilonpN) { tempEpspN =  std :: move(epsilonpN); }
    void letPlasticStateIndicatorsBe(IntArray plSt) { plasticState =  std :: move(plSt); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplanePlasticStrains() { return epspN; }
    const FloatArray &giveTempNormalMplanePlasticStrains() { return tempEpspN; }
    const IntArray &givePlasticStateIndicators() { return plasticState; }
    virtual void initTempStatus();
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};

/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
class M1Material : public MicroplaneMaterial
{
protected:

    double EN; // normal microplane elastic modulus
    double ENtan; // normal microplane tangent (elastoplastic) modulus
    double HN; // normal microplane hardening/softening modulus
    double s0; // normal microplane initial yield stress

public:
    /**
     * Constructor. Creates M1Material belonging to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M1Material(int n, Domain *d);
    /// Destructor.
    virtual ~M1Material() { }

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);
    virtual void giveRealMicroplaneStressVector(FloatArray &answer, Microplane *mplane, const FloatArray &strain, TimeStep *tStep) {};

    virtual const char *giveClassName() const { return "M1Material"; }
    virtual const char *giveInputRecordName() const { return _IFT_M1Material_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new M1MaterialStatus(1, domain, gp); }

protected:
    MaterialStatus *CreateMicroplaneStatus(GaussPoint *gp) { return new M1MaterialStatus(1, domain, gp); }
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
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
    M1MaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~M1MaterialStatus();

    // definition
    virtual const char *giveClassName() const { return "M1MaterialStatus"; }
    void letTempNormalMplaneStressesBe(FloatArray sigmaN) { tempSigN =  std :: move(sigmaN); }
    void letNormalMplaneYieldStressesBe(FloatArray sigmaNyield) { sigNyield =  std :: move(sigmaNyield); }
    const FloatArray &giveNormalMplaneStresses() { return sigN; }
    const FloatArray &giveTempNormalMplaneStresses() { return tempSigN; }
    const FloatArray &giveNormalMplaneYieldStresses() { return sigNyield; }
    virtual void initTempStatus();
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * Simple microplane model - version M1, just with normal microplane strains.
 */
class M1Material : public StructuralMaterial
{
protected:

    double E; // Young's modulus
    double nu; // Poisson ratio
    double EN; // normal microplane elastic modulus
    double HN; // normal microplane hardening/softening modulus
    double s0; // normal microplane initial yield stress
    int nmp; // number of microplanes
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
    /// Destructor.
    virtual ~M1Material() { }

    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp,
                                                  const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveElasticPlaneStressStiffMtrx(FloatMatrix &answer);
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual const char *giveClassName() const { return "M1Material"; }
    virtual const char *giveInputRecordName() const { return _IFT_M1Material_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new M1MaterialStatus(1, domain, gp); }

protected:
    MaterialStatus *CreateMicroplaneStatus(GaussPoint *gp) { return new M1MaterialStatus(1, domain, gp); }
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};
} // end namespace oofem
#endif // end of old implementation
#endif // microplane_m1_h
