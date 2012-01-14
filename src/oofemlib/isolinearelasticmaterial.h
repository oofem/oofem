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

#ifndef isolinearelasticmaterial_h
#define isolinearelasticmaterial_h

#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;

/**
 * Class implementing isotropic linear elastic material.
 * This class implements an isotropic linear elastic material in a finite
 * element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 *
 * Tasks:
 * - Returning standard material stiffness matrix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - methods give2dPlaneStressMtrx, givePlaneStrainMtrx, give1dStressMtrx are
 *   introduced since form of this matrices is well known, and for
 *   faster response mainly in linear elastic problems.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class IsotropicLinearElasticMaterial : public LinearElasticMaterial
{
protected:
    /// Young's modulus.
    double E;
    /// Poisson's ratio.
    double nu;
    /// Shear modulus.
    double G;

public:
    /**
     * Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param n material model number in domain
     * @param d domain which receiver belongs to
     */
    IsotropicLinearElasticMaterial(int n, Domain *d) : LinearElasticMaterial(n, d) { }
    /**
     * Creates a new IsotropicLinearElasticMaterial class instance
     * with given number belonging to domain d.
     * @param n Material model number in domain.
     * @param d Domain which receiver belongs to.
     * @param E Young modulus.
     * @param nu Poisson ratio.
     */
    IsotropicLinearElasticMaterial(int n, Domain *d, double E, double nu);
    /// Destructor.
    ~IsotropicLinearElasticMaterial() { }

    void giveCharacteristicMatrix(FloatMatrix &answer,
                                  MatResponseForm form,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *atTime);

    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer Vector of thermal dilatation coefficients.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when atTime is current time step).
     */
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    // identification and auxiliary functions
    int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "IsotropicLinearElasticMaterial"; }
    classType giveClassID() const { return IsotropicLinearElasticMaterialClass; }
    /// Returns input record name of the receiver.
    const char *giveInputRecordName() const { return "IsoLE"; }
    /**
     * Initializes receiver according to object description stored in input record.
     * The E modulus (keyword "E"), Poisson ratio ("nu") and coefficient of thermal dilatation
     * alpha ("talpha") are read. The parent class instanciateFrom method is called.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /**
     * Setups the input record string of receiver
     * @param str String to be filled by input record.
     * @param keyword Print record keyword (default true).
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    // non-standard - returns time independent material constant
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id.
     * @param aProperty ID of property requested.
     * @param gp Integration point.
     * @return Property value.
     */
    double give(int aProperty, GaussPoint *gp);

    /// Returns Young's modulus.
    double giveYoungsModulus() { return E; }

    /// Returns Poisson's ratio.
    double givePoissonsRatio() { return nu; }

    /// Returns the shear elastic modulus G = E / (2*(1+nu)).
    double giveShearModulus() { return G; }

    /// Returns the bulk elastic modulus K = E / (3*(1-2*nu)).
    double giveBulkModulus() { return E / ( 3. * ( 1. - 2. * nu ) ); }

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
    void givePlaneStressStiffMtrx(FloatMatrix & answer,
                                  MatResponseForm, MatResponseMode, GaussPoint * gp,
                                  TimeStep * atTime);

    void givePlaneStrainStiffMtrx(FloatMatrix & answer,
                                  MatResponseForm, MatResponseMode, GaussPoint * gp,
                                  TimeStep * atTime);

    void give1dStressStiffMtrx(FloatMatrix & answer,
                               MatResponseForm, MatResponseMode, GaussPoint * gp,
                               TimeStep * atTime);

    void give2dBeamStiffMtrx(FloatMatrix &answer,
                             MatResponseForm form, MatResponseMode rMode,
                             GaussPoint *gp,
                             TimeStep *tStep);

    void give3dBeamStiffMtrx(FloatMatrix &answer,
                             MatResponseForm form, MatResponseMode rMode,
                             GaussPoint *gp,
                             TimeStep *tStep);

    friend class CrossSection;
};
} // end namespace oofem
#endif // isolinearelasticmaterial_h
