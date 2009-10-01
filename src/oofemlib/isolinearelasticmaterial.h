/* $Header: /home/cvs/bp/oofem/oofemlib/src/isolinearelasticmaterial.h,v 1.10 2003/04/06 14:08:24 bp Exp $ */
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


//   ************************************************
//   *** CLASS ISOTROPIC LINEAR ELACSTIC MATERIAL ***
//   ************************************************

#ifndef isolinearelasticmaterial_h
#define isolinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"


class GaussPoint;

/**
 * Class implementing isotropic linear elastic material.
 */
class IsotropicLinearElasticMaterial : public LinearElasticMaterial
{
    /*
     * This class implements an isotropic linear elastic material in a finite
     * element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     *
     * DESCRIPTION
     * ISOTROPIC Linear Elastic Material
     *
     * TASK
     * - Returning standard material stiffness marix for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - methods Give2dPlaneStressMtrx, GivePlaneStrainMtrx, Give1dStressMtrx are
     * introduced since form of this matrices is well known, and for
     * faster response mainly in linear elastic problems.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     */
protected:
    /// Material properties
    double E, nu, G;

public:
    /**
     * Constructor. Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param n material model number in domain
     * @param d domain which receiver belongs to
     */
    IsotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d) { }
    /**
     * Constructor. Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param  n material model number in domain
     * @param d domain which receiver belongs to
     * @param E Young modulus
     * @param nu Poisson ratio
     */
    IsotropicLinearElasticMaterial(int n, Domain *d, double E, double nu);
    /// Destructor
    ~IsotropicLinearElasticMaterial() { }

    /**
     * Computes characteristic matrix of receiver in given integration point.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step
     */
    void  giveCharacteristicMatrix(FloatMatrix &answer,
                                   MatResponseForm form,
                                   MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *atTime);

    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    // identification and auxiliary functions
    /**
     * Test for particular material mode capability.
     * @param mode material mode requested
     * @return nonzero if available
     */
    int hasMaterialModeCapability(MaterialMode mode);
    /// Returns "IsotropicLinearElasticMaterial" - class  name of the receiver.
    const char *giveClassName() const { return "IsotropicLinearElasticMaterial"; }
    /// Returns IsotropicLinearElasticMaterialClass - classType id of receiver.
    classType giveClassID()         const { return IsotropicLinearElasticMaterialClass; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "IsoLE"; }
    /**
     * Initializes receiver acording to object description stored in input record.
     * The E modulus (keyword "E"), Poisson ration ("nu") and coefficient of thermal dilatation
     * alpha ("talpha") are read. The parent class instanciateFrom method is called.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    // non-standart - returns time independent material constant
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id.
     * @param aProperty id of peroperty requested
     * @param gp integration point
     * @return property value
     */
    double   give(int, GaussPoint*);

    /// returns the shear elastic modulus G = E / (2*(1+nu))
    double giveShearModulus()
    {return G;}

    /// returns the bulk elastic modulus K = E / (3*(1-2*nu))
    double giveBulkModulus()
    {return E/(3.*(1.-2.*nu));}


    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer computed results
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);

    /**
     * Creates new copy of associated status (StructuralMaterialStatus class )
     * and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return reference to new status.
     */
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    /**
     * @name Methods for computing material mode contributions
     * These general methods are overloaded, because default implementation computes 3d stiffness
     * matrix using give3dMaterialStiffnessMatrix and
     * reduces it to plane stress stiffness using reduce method described above.
     * Howewer, this reduction is quite time consuming and if it is possible,
     * it is recomended to overload this method and provide direct method for computing
     * particular stiffness matrix.
     */
    //@{
    /**
     * Method for computing plane stress stiffness matrix of receiver.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void givePlaneStressStiffMtrx(FloatMatrix & answer,
                                  MatResponseForm, MatResponseMode, GaussPoint * gp,
                                  TimeStep * atTime);
    /**
     * Method for computing plane strain stiffness matrix of receiver.
     * Note: as already described, if zero strain component is imposed
     * (Plane strain, ..) this condition must be taken into account in geometrical
     * relations, and corresponding component has to be included in reduced vector.
     * (So plane strain conditions are eps_z = gamma_xz = gamma_yz = 0, but relations
     * for eps_z and sigma_z are included).
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void givePlaneStrainStiffMtrx(FloatMatrix & answer,
                                  MatResponseForm, MatResponseMode, GaussPoint * gp,
                                  TimeStep * atTime);

    /**
     * Method for computing 1d  stiffness matrix of receiver.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give1dStressStiffMtrx(FloatMatrix & answer,
                               MatResponseForm, MatResponseMode, GaussPoint * gp,
                               TimeStep * atTime);

    /**
     * Method for computing 2d beam layer stiffness matrix of receiver.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give2dBeamStiffMtrx(FloatMatrix &answer,
                             MatResponseForm form, MatResponseMode rMode,
                             GaussPoint *gp,
                             TimeStep *tStep);

    /**
     * Method for computing 3d beam layer stiffness matrix of receiver.
     * @param answer stiffness matrix
     * @param form material response form
     * @param mode material response mode
     * @param gp integration point, which load history is used
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    void give3dBeamStiffMtrx(FloatMatrix &answer,
                             MatResponseForm form, MatResponseMode rMode,
                             GaussPoint *gp,
                             TimeStep *tStep);
    //@}
    friend class CrossSection;
};


#endif // isolinearelasticmaterial_h
