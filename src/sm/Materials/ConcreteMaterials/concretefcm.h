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
 *               Copyright (C) 1993 - 2016   Borek Patzak
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

#ifndef concretefcm_h
#define concretefcm_h

#include "../fcm.h"
#include "randommaterialext.h"

///@name Input fields for ConcreteFCM
//@{
#define _IFT_ConcreteFCM_Name "concretefcm"
#define _IFT_ConcreteFCM_softType "softtype"
#define _IFT_ConcreteFCM_shearType "sheartype"
#define _IFT_ConcreteFCM_shearStrengthType "shearstrengthtype"
#define _IFT_ConcreteFCM_gf "gf"
#define _IFT_ConcreteFCM_ft "ft"
#define _IFT_ConcreteFCM_beta "beta"
#define _IFT_ConcreteFCM_sf "sf"
#define _IFT_ConcreteFCM_fc "fc"
#define _IFT_ConcreteFCM_ag "ag"
#define _IFT_ConcreteFCM_lengthScale "lengthscale"
#define _IFT_ConcreteFCM_soft_w "soft_w"
#define _IFT_ConcreteFCM_soft_function_w "soft(w)"
#define _IFT_ConcreteFCM_soft_eps "soft_eps"
#define _IFT_ConcreteFCM_soft_function_eps "soft(eps)"
#define _IFT_ConcreteFCM_beta_w "beta_w"
#define _IFT_ConcreteFCM_beta_function "beta(w)"
#define _IFT_ConcreteFCM_H "h"
#define _IFT_ConcreteFCM_eps_f "eps_f"
//@}

namespace oofem {
/**
 * This class manages the status of ConcreteFCM
 */
class ConcreteFCMStatus : public FCMMaterialStatus, public RandomMaterialStatusExtensionInterface
{
public:
    ConcreteFCMStatus(int n, Domain *d, GaussPoint *g);
    virtual ~ConcreteFCMStatus();
    /// Writes information into the output file.
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "ConcreteFCMStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual Interface *giveInterface(InterfaceType it);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};



/**
 * This class implements a ConcreteFCM material in a finite element problem.
 *
 * Cracking is described using Fixed Crack Model based on crack band approach.
 * The model offers several predefined functions for softening and for reduction
 * of shear stiffness caused by cracking.
 */
class ConcreteFCM : public FCMMaterial, public RandomMaterialExtensionInterface
{
public:
    ConcreteFCM(int n, Domain *d);
    virtual ~ConcreteFCM() {
        delete linearElasticMaterial;
    }

    // identification and auxiliary functions
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "ConcreteFCM"; }
    virtual const char *giveInputRecordName() const { return _IFT_ConcreteFCM_Name; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new ConcreteFCMStatus(1, domain, gp); }

    virtual double give(int aProperty, GaussPoint *gp);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;


protected:

    /// Fracture energy
    double Gf;
    /// Tensile strenght
    double Ft;
    /// shear retention factor
    double beta;
    /// shear factor
    double sf;

    // 3 parameters for collins aggregate interlock:
    /// Collins' aggregate interlock: compressive strength in MPa
    double fc;
    /// Collins' aggregate interlock: aggregate diameter in appropriate units (same as FE mesh)
    double ag;
    /// Collins' aggregate interlock: 1 for meter, 1000 for analysis in mm
    double lengthScale;

    /// user-defined softening (traction-COD)
    FloatArray soft_w, soft_function_w;
    /// user-defined softening (traction-strain)
    FloatArray soft_eps, soft_function_eps;
    /// user-defined shear retention factor (with respect to crack opening)
    FloatArray beta_w, beta_function;

    /// hardening modulus
    double H;
    /// strain at failure
    double eps_f;

    /// returns tensile strength (can be random)
    virtual double giveTensileStrength(GaussPoint *gp) { return this->give(ft_strength, gp); }

    /// returns fracture energy (can be random)
    virtual double giveFractureEnergy(GaussPoint *gp) { return this->give(gf_ID, gp); }

    /// returns stiffness in the normal direction of the i-th crack
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, int i);

    /// returns Geff which is necessary in the global stiffness matrix
    virtual double computeEffectiveShearModulus(GaussPoint *gp, int i);

    /// shear modulus for a given crack plane (1, 2, 3)
    virtual double computeD2ModulusForCrack(GaussPoint *gp, int icrack);

    /// computes normal stress associated with i-th crack direction
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i);

    /// computes the maximum value of the shear stress; if the shear stress exceeds this value, it is cropped
    virtual double maxShearStress(GaussPoint *gp, int i);

    /// checks possible snap-back
    virtual void checkSnapBack(GaussPoint *gp, int crack);

    /// type of post-peak behavior in the normal direction to the crack plane
    enum SofteningType { ST_NONE, ST_Exponential, ST_Linear, ST_Hordijk, ST_UserDefinedCrack, ST_LinearHardeningStrain, ST_UserDefinedStrain, ST_Unknown };
    SofteningType softType;

    /// type of reduction of the shear stiffness caused by cracking
    enum ShearRetentionType { SHR_NONE, SHR_Const_ShearRetFactor, SHR_Const_ShearFactorCoeff, SHR_UserDefined_ShearRetFactor, SHR_Unknown };
    ShearRetentionType shearType;

    /// defines the maximum value of shear stress
    enum ShearStrengthType { SHS_NONE, SHS_Const_Ft, SHS_Collins_Interlock, SHS_Unknown };
    ShearStrengthType shearStrengthType;
};
} // end namespace oofem
#endif // concretefcm_h
