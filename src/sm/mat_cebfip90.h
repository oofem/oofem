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

#ifndef mat_cebfip90_h
#define mat_cebfip90_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

namespace oofem {
// material contant's keys for give()
class GaussPoint;

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 * Stores scalar measure of the largest strain level reached.
 */
class CebFipSlip90MaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Scalar measure of the largest slip reached in material.
    double kappa;
    /// Non-equilibrated scalar of the largest slip displacement.
    double tempKappa;

public:
    /// Constructor.
    CebFipSlip90MaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    ~CebFipSlip90MaterialStatus();

    void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    // definition
    const char *giveClassName() const { return "CebFipSlip90MaterialStatus"; }
    classType giveClassID() const { return MaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class CebFipSlip90Material : public StructuralMaterial
{
protected:
    /// Max force (stress).
    double tmax;
    /// Slip valu at begining of yield plateau.
    double s1;
    /// Slip at end of plateau.
    double s2;
    /// Slip when residual force/stress activated.
    double s3;
    /// Residual force/stress.
    double tres;
    /// Alpha coeff.
    double alpha;

public:
    /// Constructor
    CebFipSlip90Material(int n, Domain *d);
    /// Destructor
    ~CebFipSlip90Material();

    int hasNonLinearBehaviour() { return 1; }
    int hasMaterialModeCapability(MaterialMode mode);

    const char *giveClassName() const { return "CebFipSlip90Material"; }
    classType giveClassID() const { return StructuralMaterialClass; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    void giveRealStressVector(FloatArray & answer,  MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual int giveStressStrainComponentIndOf(MatResponseForm form, MaterialMode mmode, int ind);
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode mmode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &strainVector);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the value of bond force/stress, based on given value of slip value.
     * @param kappa Slip value.
     */
    double computeBondForce(double kappa);
    /**
     * Computes the value of bond force/stress stiffness, based on given value of slip value.
     * @param kappa Slip value.
     */
    double computeBondForceStiffness(double kappa);

    IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new CebFipSlip90MaterialStatus(1, domain, gp); }

protected:
    void give1dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // mat_cebfip90_h
