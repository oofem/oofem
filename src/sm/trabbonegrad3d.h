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

//   *************************************************
//   *** CLASS GRADIENT BONE DAMAGE-PLASTIC MODEL ***
//   *************************************************

#ifndef TrabBoneGrad3D_h

#include "trabbone3d.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"
#include "linearelasticmaterial.h"
#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
class GaussPoint;


/////////////////////////////////////////////////////////////////
///////// GRADIENT MISES NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneGrad3DStatus : public TrabBone3DStatus
{
protected:
    // STATE VARIABLE DECLARATION
    // Equivalent strain for avaraging
    double nlKappa;
     /// reference to the basic elastic material
    LinearElasticMaterial *linearElasticMaterial;
public:
    // CONSTRUCTOR
    TrabBoneGrad3DStatus(int n, Domain *d, GaussPoint *g);

    // DESTRUCTOR
    ~TrabBoneGrad3DStatus();

    // OUTPUT PRINT
    // Prints the receiver state to stream
    void   printOutputAt(FILE *file, TimeStep *tStep);

    // DEFINITION
    const char *giveClassName() const { return "TrabBoneGrad3DStatus"; }
    classType   giveClassID() const { return TrabBoneGrad3DClass; }
    double giveNlKappa(){return nlKappa;}
    void setNlKappa(double kappa){nlKappa = kappa;}
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
//////////// GRADIENT BONE MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


class TrabBoneGrad3D : public TrabBone3D
{
protected:
    // STATE VARIABLE DECLARATION
    // declare material properties here
    double l;
    double mParam;

public:
    // CONSTRUCTOR
    TrabBoneGrad3D(int n, Domain *d);

    // DESTRUCTOR
    ~TrabBoneGrad3D();

    // INITIALIZATION OF FUNCTION/SUBROUTINE
    const char *giveClassName() const { return "TrabBoneGrad3D"; }
    classType   giveClassID()   const { return TrabBoneGrad3DClass; }
    const char *giveInputRecordName() const { return "TrabBoneGrad3D"; }

    // Initializes the receiver from given record
    IRResultType initializeFrom(InputRecord *ir);
    int  hasMaterialModeCapability(MaterialMode mode);




    void giveCharacteristicMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);
    /*******************************************************************/
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,  MatResponseForm, MatResponseMode, GaussPoint * gp,  TimeStep * atTime);
    /*********************************************************************/
      void give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /*********************************************************************/
     void give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /**********************************************************************/
    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    /********************************************************************/
    void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);
    void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);
    //LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneGrad3DStatus(1, TrabBone3D :: domain, gp); }
};
} // end namespace oofem
#define TrabBoneGrad3D_h
#endif
