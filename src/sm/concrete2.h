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
 *               Copyright (C) 1993 - 2011   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef concrete2_h
#define concrete2_h

#include "femcmpnn.h"
#include "dictionr.h"
#include "material.h"
#include "deformationtheorymaterial.h"
#include "isolinearelasticmaterial.h"
#include "flotarry.h"

#include "structuralms.h"

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
    Concrete2MaterialStatus(int n, Domain *d, GaussPoint *g);
    ~Concrete2MaterialStatus();
    void   printOutputAt(FILE *file, TimeStep *tStep)
    { StructuralMaterialStatus :: printOutputAt(file, tStep); }

    void givePlasticStrainVector(FloatArray &answer) const { answer = plasticStrainVector; }
    void givePlasticStrainIncrementVector(FloatArray &answer) const
    { answer = plasticStrainIncrementVector; }
    void letPlasticStrainVectorBe(FloatArray &v)
    { plasticStrainVector = v; }
    void letPlasticStrainIncrementVectorBe(FloatArray &v)
    { plasticStrainIncrementVector = v; }

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
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    const char *giveClassName() const { return "Concrete2MaterialStatus"; }
    classType giveClassID() const { return Concrete2MaterialStatusClass; }
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
    Concrete2(int n, Domain *d);
    ~Concrete2();

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &, TimeStep *atTime);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime);

    MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    void  giveRealStresses3dShellLayer(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                       const FloatArray &strain, TimeStep *atTime);
    void  dtp3(GaussPoint *gp, FloatArray *e, FloatArray *s, FloatArray *ep,
               double SCC, double SCT, int *ifplas);
    void  dtp2(GaussPoint *gp, FloatArray *e, FloatArray *s, FloatArray *ep,
               double SCC, double SCT, int *ifplas);
    void  stirr(double dez, double srf);
    void  strsoft(GaussPoint *gp, double epsult, FloatArray *ep, double &ep1,
                  double &ep2, double &ep3, double SCC, double SCT, int &ifupd);

    // two functions used to initialize and updating temporary variables in
    // gp's status. These variables are used to control process, when
    // we try to find equilibrium state.
    void updateStirrups(GaussPoint *gp, FloatArray *strainIncrement);

public:
    double give(int, GaussPoint *gp);

    // identification and auxiliary functions
    int hasNonLinearBehaviour() { return 1; }
    const char *giveClassName() const { return "Concrete2"; }
    classType giveClassID() const { return Concrete2Class; }
    IRResultType initializeFrom(InputRecord *ir);
};
} // end namespace oofem
#endif // concrete2_h
