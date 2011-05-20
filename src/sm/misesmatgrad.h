/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.h,v 1.9 2003/04/06 14:08:30 bp Exp $ */
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

//   ********************************************
//   *** CLASS NONLOCAL TRABECULAR BONE MODEL ***
//   ********************************************

#ifndef MisesMatGrad_h

#include "misesmat.h"
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
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


class MisesMatGradStatus : public MisesMatStatus
{
protected:
    // STATE VARIABLE DECLARATION
    // Equivalent strain for avaraging
    double localCumPlastStrainForAverage;

public:
    // CONSTRUCTOR
    MisesMatGradStatus(int n, Domain *d, GaussPoint *g);

    // DESTRUCTOR
    ~MisesMatGradStatus();

    // OUTPUT PRINT
    // Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    // DEFINITION
    const char *giveClassName() const { return "MisesMatGradStatus"; }
    classType   giveClassID() const { return MisesMatClass; }

    // INITIALISATION OF TEMPORARY VARIABLES
    // Initializes the temporary internal variables, describing the current state according to
    // previously reached equilibrium internal variables.
    virtual void initTempStatus();

    // UPDATE VARIABLES
    // Update equilibrium history variables according to temp-variables.
    // Invoked, after new equilibrium state has been reached.
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached
};


/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


class MisesMatGrad : public MisesMat
{
protected:
    // STATE VARIABLE DECLARATION
    // declare material properties here
    double R;
    double mParam;

public:
    // CONSTRUCTOR
    MisesMatGrad(int n, Domain *d);

    // DESTRUCTOR
    ~MisesMatGrad();

    // INITIALIZATION OF FUNCTION/SUBROUTINE
    const char *giveClassName() const { return "MisesMatGrad"; }
    classType   giveClassID()   const { return MisesMatClass; }
    const char *giveInputRecordName() const { return "MisesMatGrad"; }

    // Initializes the receiver from given record
    IRResultType initializeFrom(InputRecord *ir);
    int  hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    /*******************************************************************/
    virtual void give1dStressStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    /*********************************************************************/
    void give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /*********************************************************************/
    void give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void givePlaneStrainGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /**********************************************************************/
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /********************************************************************/
    void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);
    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatGradStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define MisesMatGrad_h
#endif
