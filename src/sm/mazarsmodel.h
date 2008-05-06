/* $Header: /home/cvs/bp/oofem/sm/src/mazarsmodel.h,v 1.6 2003/04/06 14:08:31 bp Exp $ */
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

//   *******************************************************
//   *** CLASS MAZARS DAMAGE MODEL FOR CONCRETE ************
//   *******************************************************

#ifndef mazarsmodel_h
#define mazarsmodel_h

#include "linearelasticmaterial.h"
#include "idm1.h"
#include "structuralms.h"


/**
 * This class implements associated Material Status to MazarsMaterial.
 */
class MazarsMaterialStatus : public IsotropicDamageMaterial1Status
{
protected:

    /**
     * characteristic element length for compression, fixed as square from element size (for 2d)
     */
    double lec;

public:

    /// Constructor
    MazarsMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~MazarsMaterialStatus() { }

    // void   printOutputAt (FILE *file, TimeStep* tStep) ;

    /// Returns characteristic length stored in receiver
    double giveLec()  { return lec; }
    /// Sets characteristic length to given value
    void   setLec(double ls) { lec = ls; }

    // definition
    const char *giveClassName() const { return "MazarsMaterialStatus"; }
    classType             giveClassID() const { return IsotropicDamageMaterialStatusClass; }

    // saves current context(state) into stream
    /**
     * Stores context of receiver into given stream.
     * Le attribute is stored. Corresponding parent method invoked.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return nonzero if o.k.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * Le attribute is restored. Corresponding parent method invoked.
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return nonzero if o.k.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class implements a Mazars damage model form concrete.
 * It introduces two damage parameters omega_t and omega_c that
 * are computed from the same equivalent strain using two different damage functions
 * g_t and g_c. The g_t is identified from the uniaxial tension tests, while
 * g_c from compressive test. The damage parameter for general stress states
 * omega is obtained as a linear combination of omega_t and omega_c.
 */
class MazarsMaterial : public IsotropicDamageMaterial1
{
    /*
     *
     * DESCRIPTION
     * This class implements a Mazars damage model for concrete
     *
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     * - Storing & restoring Material Staus sored in gp matStatusDictionary.
     */

protected:

    /// Model parametr's related to the shape of uniaxial stress-strain diagrams
    double At, Bt, Ac, Bc;
    /// Strain at the onset of non-linearity
    double eps_0;
    /// Refernce elem-length for objectivity
    double hReft, hRefc;
    /// Beta coefficient reducing the effect of shear; defaul val = 1.06
    double beta;
    /// Model variants
    enum mazarsModelVariant { maz_original, maz_modTension } modelVersion;
public:

    /// Constructor
    MazarsMaterial(int n, Domain *d);
    /// Destructor
    ~MazarsMaterial();

    // identification and auxiliary functions
    const char *giveClassName() const { return "MazarsMaterial"; }
    classType giveClassID()         const { return MazarsMaterialClass; }


    /**
     * Initializes receiver acording to object description stored in input record.
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    /**
     * Computes the equivalent strain measure from given strain vector (full form).
     * @param kappa return param, comtaining the corresponding equivalent strain
     * @param strain total strain vector in full form
     * @param gp integration point
     * @param atTime timestep
     */
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *atTime);
    /**
     * computes the value of damage parameter omega, based on given value of equivalent strain
     * @param omega contains result
     * @param kappa equivalent strain measure
     */
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    /// Creates corresponding status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MazarsMaterialStatus(1, IsotropicDamageMaterial1 :: domain, gp); }

protected:

    /**
     *  Perfoms initialization, when damage first appear. The Le characteristic length is
     *  computed from the direction of largest positive principal strain and stored
     *  in corresponding status.
     *  @param kappa scalar measure of strain level
     *  @param totalStrainVector current total strain vector
     *  @param gp integration point
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);
    /**
     * Computes elastic stiffness for normal stress components
     * @param answer result of size (3,3)
     * @param mode determines the MatResponseMode
     * @param gp integration point
     * @param atTime time step
     */
    void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                          MatResponseMode rMode,
                                          GaussPoint *gp, TimeStep *atTime);
};

#endif // mazarsmodel_h

