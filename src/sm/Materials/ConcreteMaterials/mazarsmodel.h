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

#ifndef mazarsmodel_h
#define mazarsmodel_h

#include "Materials/linearelasticmaterial.h"
#include "Materials/ConcreteMaterials/idm1.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for MazarsMaterial
//@{
#define _IFT_MazarsMaterial_Name "mazarsmodel"
#define _IFT_MazarsMaterial_version "version"
#define _IFT_MazarsMaterial_e0 "e0"
#define _IFT_MazarsMaterial_ac "ac"
#define _IFT_MazarsMaterial_bc "bc"
#define _IFT_MazarsMaterial_beta "beta"
#define _IFT_MazarsMaterial_at "at"
#define _IFT_MazarsMaterial_bt "bt"
#define _IFT_MazarsMaterial_ef "ef"
#define _IFT_MazarsMaterial_r "r"
#define _IFT_MazarsMaterial_hreft "hreft"
#define _IFT_MazarsMaterial_hrefc "hrefc"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to MazarsMaterial.
 */
class MazarsMaterialStatus : public IsotropicDamageMaterial1Status
{
protected:
    /// Characteristic element length for compression, fixed as square from element size (for 2d).
    double lec;

public:
    /// Constructor.
    MazarsMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor.
    virtual ~MazarsMaterialStatus() { }

    /// Returns characteristic length stored in receiver.
    double giveLec() { return lec; }
    /// Sets characteristic length to given value.
    void setLec(double ls) { lec = ls; }

    // definition
    virtual const char *giveClassName() const { return "MazarsMaterialStatus"; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
protected:
    /// Elastic parameters.
    double E, nu;
    /// Model parameters related to the shape of uniaxial stress-strain diagrams.
    double At, Bt, Ac, Bc;
    /// Reference elem-length for objectivity.
    double hReft, hRefc;
    /// Beta coefficient reducing the effect of shear; default val = 1.06.
    double beta;

    /// Model variants.
    enum mazarsModelVariant { maz_original, maz_modTension } modelVersion;

public:
    /// Constructor
    MazarsMaterial(int n, Domain * d);
    /// Destructor
    virtual ~MazarsMaterial();

    // identification and auxiliary functions
    virtual const char *giveInputRecordName() const { return _IFT_MazarsMaterial_Name; }
    virtual const char *giveClassName() const { return "MazarsMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MazarsMaterialStatus(1, IsotropicDamageMaterial1 :: domain, gp); }

protected:
    /**
     *  Perfoms initialization, when damage first appear. The Le characteristic length is
     *  computed from the direction of largest positive principal strain and stored
     *  in corresponding status.
     *  @param kappa Scalar measure of strain level,
     *  @param totalStrainVector Current total strain vector,
     *  @param gp Integration point,
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp);
    /*
     * Computes elastic stiffness for normal stress components.
     * @param answer Result of size (3,3).
     * @param mode Determines the MatResponseMode.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    /*
     * void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
     *                                    MatResponseMode rMode,
     *                                    GaussPoint *gp, TimeStep *tStep);
     */

    int giveNumberOfSpatialDimensions(GaussPoint *gp);
    void giveNormalBlockOfElasticCompliance(FloatMatrix &answer, GaussPoint *gp);
    double computeGt(double kappa, GaussPoint *gp);
    double computeGc(double kappa, GaussPoint *gp);
};
} // end namespace oofem
#endif // mazarsmodel_h
