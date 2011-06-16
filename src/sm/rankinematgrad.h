/*
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2000   Borek Patzak
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

//   *************************************************
//   *** CLASS GRADIENT RANKINE DAMAGE-PLASTIC MODEL ***
//   *************************************************

#ifndef RankineMatGrad_h

#include "rankinemat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
class GaussPoint;


/////////////////////////////////////////////////////////////////
///////// GRADIENT RANKINE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


class RankineMatGradStatus : public RankineMatStatus
{
protected:
    // STATE VARIABLE DECLARATION

    // for printing only
    double kappa_nl;
    double kappa_hat;

public:
    // CONSTRUCTOR
    RankineMatGradStatus(int n, Domain *d, GaussPoint *g);

    // DESTRUCTOR
    ~RankineMatGradStatus() {; }

    // OUTPUT PRINT
    // Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    // DEFINITION
    const char *giveClassName() const { return "RankineMatGradStatus"; }
    classType   giveClassID() const { return RankineMatClass; }

    // INITIALISATION OF TEMPORARY VARIABLES
    // Initializes the temporary internal variables, describing the current state according to
    // previously reached equilibrium internal variables.
    virtual void initTempStatus();

    // UPDATE VARIABLES
    // Update equilibrium history variables according to temp-variables.
    // Invoked, after new equilibrium state has been reached.
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }
};


/////////////////////////////////////////////////////////////////
//////////// GRADIENT RANKINE MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


class RankineMatGrad : public RankineMat
{
protected:
    // STATE VARIABLE DECLARATION
    // declare material properties here
    double R;
    double mParam;
    double negligible_damage;

public:
    // CONSTRUCTOR
    RankineMatGrad(int n, Domain *d);

    // DESTRUCTOR
    ~RankineMatGrad() {; }

    // INITIALIZATION OF FUNCTION/SUBROUTINE
    const char *giveClassName() const { return "RankineMatGrad"; }
    classType   giveClassID()   const { return RankineMatClass; }
    const char *giveInputRecordName() const { return "RankineMatGrad"; }

    // Initializes the receiver from given record
    IRResultType initializeFrom(InputRecord *ir);
    int  hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    /*******************************************************************/
    //virtual void give1dStressStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    virtual void givePlaneStressStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    //virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    /*********************************************************************/
    //void give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    //void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /*********************************************************************/
    //void give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStressGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    //void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /**********************************************************************/
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /********************************************************************/
    void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);
    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    double giveNonlocalCumPlasticStrain(GaussPoint *gp);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
    int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    InternalStateValueType giveIPValueType(InternalStateType type);
    int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    int giveIPValueSize(InternalStateType type, GaussPoint *gp);

protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RankineMatGradStatus(1, RankineMat :: domain, gp); }
};
} // end namespace oofem
#define RankineMatGrad_h
#endif
