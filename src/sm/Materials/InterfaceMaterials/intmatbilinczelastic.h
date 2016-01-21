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

#ifndef IntMatBilinearCZElastic_h
#define IntMatBilinearCZElastic_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatBilinearCZElastic
//@{
#define _IFT_IntMatBilinearCZElastic_Name "intmatbilinearczelastic"
#define _IFT_IntMatBilinearCZElastic_kn "kn"
#define _IFT_IntMatBilinearCZElastic_ks "ks"
#define _IFT_IntMatBilinearCZElastic_knc "knc"
#define _IFT_IntMatBilinearCZElastic_g1c "g1c"
#define _IFT_IntMatBilinearCZElastic_sigfn "sigfn"
#define _IFT_IntMatBilinearCZElastic_sigfs "sigfs"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for ...
 */
class IntMatBilinearCZElasticStatus : public StructuralInterfaceMaterialStatus
{
protected:

public:
    /// Constructor
    IntMatBilinearCZElasticStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~IntMatBilinearCZElasticStatus();

    double giveDamage() { return 0.0; } // no damage in this model
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "IntMatBilinearCZElasticStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    //virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    //virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * Simple isotropic damage based model for 2d interface elements.
 * In 2d, the interface elements are used to model contact layer between
 * element edges. The generalized strain vector contains two relative displacements
 * (in normal and shear direction). The generalized stress vector contains corresponding
 * tractions in normal and tangent direction.
 *
 * The behaviour of the model is elastic, described by normal and shear stiffness components.
 * Isotropic damage is initiated  when the stress reaches the tensile strength. Damage evolution
 * is governed by normal component of generalized strain vector (normal relative displacement)
 * by an exponential softening law.
 */
class IntMatBilinearCZElastic : public StructuralInterfaceMaterial
{
protected:
    /// Material parameters
    double kn0;
    double ks0;
    double knc;   // stiffness in compression
    double GIc;
    double sigfn;
    double sigfs;

    double gn0;   // normal jump at damage initiation
    double gs0;   // shear jump at damage initiations
    double gnmax; // max normal jump

    double kn1;   // slope during softening part

    virtual int checkConsistency();
    void give3dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);
public:
    /// Constructor
    IntMatBilinearCZElastic(int n, Domain * d);
    /// Destructor
    virtual ~IntMatBilinearCZElastic();

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual const char *giveClassName() const { return "IntMatBilinearCZElastic"; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatBilinearCZElastic_Name; }

    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jumpVector,
                                        const FloatMatrix &F, TimeStep *tStep);

    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new StructuralInterfaceMaterialStatus(1, domain, gp); }
    void printYourself();
    virtual bool hasAnalyticalTangentStiffness() const { return true; }
};
} // end namespace oofem
#endif // isointerfacedamage01_h
