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

#ifndef steelrelaxmat_h
#define steelrealxmat_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SteelRelaxMat
//@{
#define _IFT_SteelRelaxMat_Name "steelrelaxmat"
#define _IFT_SteelRelaxMat_E "e"
#define _IFT_SteelRelaxMat_reinfClass "reinfclass"
#define _IFT_SteelRelaxMat_k1 "k1"
#define _IFT_SteelRelaxMat_k2 "k2"
#define _IFT_SteelRelaxMat_rho1000 "rho1000"
#define _IFT_SteelRelaxMat_timeFactor "timefactor"
#define _IFT_SteelRelaxMat_charStrength "charstrength"
#define _IFT_SteelRelaxMat_approach "approach"
#define _IFT_SteelRelaxMat_tolerance "tolerance"
#define _IFT_SteelRelaxMat_relRelaxBound "relrelaxbound"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * Implementation of the material model for steel relaxation given in Eurocode 2 
 * (the same as in Model Code 2010) and in Ba\v{z}ant and Yu (J. of Eng. Mech, 2013)
 * which reduces to the first model under constant strain. At variable strain history 
 * the first model uses the approach employing the so-called {\sl{equivalent time}} 
 * approach described in Annex D in the Eurocode 2.
 * The current implementation takes into account only prestress losses
 * due to steel relaxation, other losses (e.g. slip at anchorage,
 * thermal dilation, friction, etc.) need to be treated separately. The same holds
 * for the stress transfer from prestressing reinforcement to concrete in
 * the region called {\sl{transmission length}}. On the other hand,
 * losses due to sequential prestressing, elastic deformation and both
 * short-time and long-time creep and shrinkage are taken into account
 * automatically provided that a suitable material model is chosen for
 * concrete. 
 * See material manual and the above-mentioned documents for details.
 */
class SteelRelaxMat : public StructuralMaterial
{
protected:

    /// Young's modulus
    double E;

    /// constant depending on the reinforcement class
    double k1;

    /// constant depending on the reinforcement class
    double k2;

    /// constant depending on the reinforcement class
    double rho1000;

    /// ratio of prestress vs. characteristic strength
    double mu;

    /**
     * Scaling factor transforming the actual time into
     * appropriate units needed by the formulae of the eurocode. For analysis
     * in days timeFactor = 1, for analysis in seconds timeFactor = 86,400.
     */
    double timeFactor;

    //double stiffnessFactor;

    /// characteristic strength of prestressing steel in appropriate units (not necessarily MPa)
    double charStrength;

    /// tolerance specifying the residual in the stress evaluation algorithm, default value is $10^{-6}$
    double tolerance;

    /**
     * Ratio of stress to characteristic strength
     * under which the relaxation is zero (typically 0.4--0.5); default
     * value is zero.
     */
    double relRelaxBound;

    /**
     * 0 = approach according to Ba\v{z}ant and Yu, 
     * 1 = equivalent time approach according to Eurocode 2 and {\sl{fib}} Model Code 2010
     */
    enum approachType { Bazant_EC2, EquivTime_EC2 } Approach;



public:
    SteelRelaxMat(int n, Domain *d);
    virtual ~SteelRelaxMat();

    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    /**
     * evaluates stress-related strain - subtracts not only temperature strains but also strains caused by steel relaxation
     */
    virtual void giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain,
                                                       TimeStep *tStep, ValueModeType mode);
    /**
     * evaluates eigenstrain due to steel relaxation
     */
    void computeStressRelaxationStrainVector(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain,
                                             TimeStep *tStep, ValueModeType mode);

    /**
     * computes steel relaxation (eurocode formula)
     */
    void evalStressRelaxationAtConstStrain(double &answer, GaussPoint *gp, double dt);

    /**
     * implementation of cumulative time approach according to Eurocode to get prestress loss at variable strain
     */
    void computeIncrOfPrestressLossAtVarStrain(double &answer, GaussPoint *gp, TimeStep *tStep, double stress);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_SteelRelaxMat_Name; }
    virtual const char *giveClassName() const { return "SteelRelaxMat"; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class SteelRelaxMatStatus : public StructuralMaterialStatus
{
protected:

  /**
   * For Bazant's approach this internal variable is a cumulative viscous strain
   * while for Eurocode approach (equivalent time) it is a cumulative prestress loss
   */
    double relaxIntVariable;
    double tempRelaxIntVariable;

    double prestress;

public:
    SteelRelaxMatStatus(int n, Domain *d, GaussPoint *g);
    virtual ~SteelRelaxMatStatus();



    void setTempRelaxIntVariable(double src) { tempRelaxIntVariable = src; }
    double giveTempRelaxIntVariable(void) { return tempRelaxIntVariable; }
    double giveRelaxIntVariable(void) { return relaxIntVariable; }

    void setPrestress(double src) { prestress = src; }
    double givePrestress(void) { return prestress; }


    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "SteelRelaxMatStatus"; }
};
} // end namespace oofem
#endif // steelrelaxmat_h
