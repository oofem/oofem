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

#ifndef concrete2_h
#define concrete2_h

#include "femcmpnn.h"
#include "dictionary.h"
#include "material.h"
#include "Materials/deformationtheorymaterial.h"
#include "Materials/isolinearelasticmaterial.h"
#include "floatarray.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for Concrete2
//@{
#define _IFT_Concrete2_Name "concrete2"
#define _IFT_Concrete2_e "e"
#define _IFT_Concrete2_n "n"
#define _IFT_Concrete2_sccc "sccc"
#define _IFT_Concrete2_scct "scct"
#define _IFT_Concrete2_epp "epp"
#define _IFT_Concrete2_epu "epu"
#define _IFT_Concrete2_eopp "eopp"
#define _IFT_Concrete2_eopu "eopu"
#define _IFT_Concrete2_sheartol "sheartol"
#define _IFT_Concrete2_is_plastic_flow "is_plastic_flow"
#define _IFT_Concrete2_ifad "ifad"
#define _IFT_Concrete2_stirr_e "stirr_e"
#define _IFT_Concrete2_stirr_ft "stirr_ft"
#define _IFT_Concrete2_stirr_a "stirr_a"
#define _IFT_Concrete2_stirr_tol "stirr_tol"
#define _IFT_Concrete2_stirr_eref "stirr_eref"
#define _IFT_Concrete2_stirr_lambda "stirr_lambda"
//@}

namespace oofem {
#define c2_SCCC  300
#define c2_SCCT  301
#define c2_EPP   302
#define c2_EPU   303
#define c2_EOPP  304
#define c2_EOPU  305
#define c2_SHEARTOL 306
#define c2_E     307
#define c2_n     308
#define stirr_E  309
#define stirr_Ft 310
#define stirr_A  311
#define stirr_TOL 312
#define stirr_EREF 313
#define stirr_LAMBDA 314
#define c2_IS_PLASTIC_FLOW 315
#define c2_IFAD  316

/**
 * This class implements associated Material Status to Concrete2Material.
 */
class Concrete2MaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray plasticStrainVector; // full form
    FloatArray plasticStrainIncrementVector;

    double tempSCCM, tempEPM, tempSCTM, tempE0PM, tempSRF, tempSEZ;

    double SCCM; ///< Current pressure strength.
    double EPM;  ///< Max. eff. plastic strain.
    double SCTM; ///< Current tension strength.
    double E0PM; ///< Max. vol. plastic strain.
    double SRF;  ///< current stress in stirrups.
    double SEZ;  ///< Current strain in transverse (z) direction.

public:
    Concrete2MaterialStatus(int n, Domain * d, GaussPoint * g);
    virtual ~Concrete2MaterialStatus();
    virtual void printOutputAt(FILE *file, TimeStep *tStep)
    { StructuralMaterialStatus :: printOutputAt(file, tStep); }

    const FloatArray & givePlasticStrainVector() const { return plasticStrainVector; }
    const FloatArray & givePlasticStrainIncrementVector() const
    { return plasticStrainIncrementVector; }
    void letPlasticStrainVectorBe(FloatArray v)
    { plasticStrainVector = std :: move(v); }
    void letPlasticStrainIncrementVectorBe(FloatArray v)
    { plasticStrainIncrementVector = std :: move(v); }

    double &giveTempCurrentPressureStrength() { return tempSCCM; }
    double &giveTempMaxEffPlasticStrain()     { return tempEPM; }
    double &giveTempCurrentTensionStrength()  { return tempSCTM; }
    double &giveTempCurrentStressInStirrups() { return tempSRF; }
    double &giveTempCurrentStrainInZDir()     { return tempSEZ; }
    double &giveTempMaxVolPlasticStrain()     { return tempE0PM; }

    // query for non-tem variables (usefull for postprocessing)
    double &giveCurrentPressureStrength() { return SCCM; }
    double &giveMaxEffPlasticStrain()     { return EPM; }
    double &giveCurrentTensionStrength()  { return SCTM; }
    double &giveCurrentStressInStirrups() { return SRF; }
    double &giveCurrentStrainInZDir()     { return SEZ; }
    double &giveMaxVolPlasticStrain()     { return E0PM; }


    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "Concrete2MaterialStatus"; }
};


/**
 * NonLinear elasto-plastic material model with hardening.
 * Plane stress or uniaxial stress + transverse shear in
 * concrete layers with transverse stirrups.
 */
class Concrete2 : public DeformationTheoryMaterial
{
private:
    double SCCC; ///< Pressure strength.
    double SCCT; ///< Tension strength.
    double EPP;  ///< Threshold eff. plastic strain for softening in compress.
    double EPU;  ///< Ultimate eff. pl. strain.
    double EOPP; ///< Threshold volumetric plastic strain for soft. in tension.
    double EOPU; ///< Ultimate vol. pl. strain.
    /**
     * Threshold value of the relative shear deformation
     * (psi^2/eef) at which shear is considered in layers. for
     * lower r.s.d. the transverse shear remains elastic decoupled
     * from bending. default value SHEARTOL = 0.01
     */
    double SHEARTOL;

    double E, n;
    // stirrups
    double stirrE, stirrFt, stirrA, stirrTOL, stirrEREF, stirrLAMBDA;
    /// Indicates that plastic flow (not deformation theory) is used in pressure.
    int IS_PLASTIC_FLOW;
    /// Determines if state variables should be updated or not (>0 updates).
    int IFAD;

    LinearElasticMaterial *linearElasticMaterial;

public:
    Concrete2(int n, Domain * d);
    virtual ~Concrete2();

    virtual void giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &, TimeStep *tStep);

    virtual void givePlateLayerStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    void dtp3(GaussPoint *gp, FloatArray &e, FloatArray &s, FloatArray &ep,
              double SCC, double SCT, int *ifplas);
    void dtp2(GaussPoint *gp, FloatArray &e, FloatArray &s, FloatArray &ep,
              double SCC, double SCT, int *ifplas);
    void stirr(double dez, double srf);
    void strsoft(GaussPoint *gp, double epsult, FloatArray &ep, double &ep1,
                 double &ep2, double &ep3, double SCC, double SCT, int &ifupd);

    // two functions used to initialize and updating temporary variables in
    // gp's status. These variables are used to control process, when
    // we try to find equilibrium state.
    void updateStirrups(GaussPoint *gp, FloatArray &strainIncrement, TimeStep *tStep);

public:
    virtual double give(int, GaussPoint *gp);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "Concrete2"; }
    virtual const char *giveInputRecordName() const { return _IFT_Concrete2_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // concrete2_h
