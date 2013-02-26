/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.h,v 1.9.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ******************************************************************************
//   *** CLASS NONLOCAL ISOTROPIC DAMAGE MODEL FOR CONCRETE IN TENSION ************
//   ******************************************************************************

#ifndef idmgrad1_h
#define idmgrad1_h

#include "idm1.h"


#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {


class IDGMaterial : public IsotropicDamageMaterial1
{
protected:
  //length scale parameter
  double internalLength;
  ///Parameter specifiying which averaging type(0 - classical, 1 - distance based, 2 - stress based)
  int averType;
  /// Parameters specifying how the length scale parameter l is adjusted
  double beta,t;
    

public:

    /// Constructor
    IDGMaterial(int n, Domain *d);
    /// Destructor
    ~IDGMaterial();

   MaterialStatus *CreateStatus(GaussPoint *gp);
    // identification and auxiliary functions
    const char *giveClassName() const { return "IDGMaterial"; }
    classType giveClassID()         const { return IDGMaterialClass; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "idmgrad1"; }
    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,const FloatArray &totalStrain,TimeStep *atTime);
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    /// Initializes the receiver from given record
    IRResultType initializeFrom(InputRecord *ir);

    int hasMaterialModeCapability(MaterialMode mode);

    void giveCharacteristicMatrix(FloatMatrix &answer,  MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime);

    void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    void givePlaneStressKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);


    void givePlaneStressGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    void giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    void giveInternalLengthDerivative(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

 virtual  void initDamaged(double kappa, FloatArray &strainVector, GaussPoint *gp);
 
  void computeEta(FloatMatrix &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);



};






class IDGMaterialStatus : public IsotropicDamageMaterial1Status
{


public:
    /// Constructor
    IDGMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~IDGMaterialStatus();

   

    // definition


    const char *giveClassName() const { return "IDGMaterialStatus"; }
    classType             giveClassID() const { return IDGMaterialClass; }

    /**
     * Initializes the temporary internal variables, describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

 
    // saves current context(state) into stream
    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
   


};

  
} // end namespace oofem
#endif // idmgrad1_h
