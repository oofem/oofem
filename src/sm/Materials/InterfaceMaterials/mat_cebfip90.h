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

#ifndef mat_cebfip90_h
#define mat_cebfip90_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for CebFipSlip90Material
//@{
#define _IFT_CebFipSlip90Material_Name "cebfipslip90"
#define _IFT_CebFipSlip90Material_tmax "tmax"
#define _IFT_CebFipSlip90Material_tres "tres"
#define _IFT_CebFipSlip90Material_s1 "s1"
#define _IFT_CebFipSlip90Material_s2 "s2"
#define _IFT_CebFipSlip90Material_s3 "s3"
//@}

namespace oofem {

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 * Stores scalar measure of the largest strain level reached.
 */
class CebFipSlip90MaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:
    /// Scalar measure of the largest slip reached in material.
    double kappa;
    /// Non-equilibrated scalar of the largest slip displacement.
    double tempKappa;

public:
    /// Constructor.
    CebFipSlip90MaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor.
    virtual ~CebFipSlip90MaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    // definition
    virtual const char *giveClassName() const { return "CebFipSlip90MaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class CebFipSlip90Material : public StructuralInterfaceMaterial
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
    CebFipSlip90Material(int n, Domain * d);
    /// Destructor
    virtual ~CebFipSlip90Material();

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual const char *giveInputRecordName() const { return _IFT_CebFipSlip90Material_Name; }
    virtual const char *giveClassName() const { return "CebFipSlip90Material"; }

    virtual void giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep);
    virtual void give1dStiffnessMatrix_Eng(FloatMatrix &answer,  MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

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

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new CebFipSlip90MaterialStatus(1, domain, gp); }

};
} // end namespace oofem
#endif // mat_cebfip90_h
