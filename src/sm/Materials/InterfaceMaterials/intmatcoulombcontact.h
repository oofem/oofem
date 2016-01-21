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

#ifndef inmatcoulombcontact_h
#define inmatcoulombcontact_h

#include "material.h"
#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatCoulombContact
//@{
#define _IFT_IntMatCoulombContact_Name "intmatcoulombcontact"
#define _IFT_IntMatCoulombContact_kn "kn"
#define _IFT_IntMatCoulombContact_knt "knt"
#define _IFT_IntMatCoulombContact_frictCoeff "fc"
#define _IFT_IntMatCoulombContact_stiffCoeff "stiffcoeff"
#define _IFT_IntMatCoulombContact_normalClearance "normalclearance"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to IntMatCoulombContact.
 * @author Milan Jirisek
 * @author Vit Smilauer
 * @author Jim Brouzoulis
 */
class IntMatCoulombContactStatus : public StructuralInterfaceMaterialStatus
{
protected:
    FloatArray shearStressShift, tempShearStressShift;

public:
    /// Constructor
    IntMatCoulombContactStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~IntMatCoulombContactStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    FloatArray giveShearStressShift();
    void setTempShearStressShift(FloatArray newShearStressShift) { tempShearStressShift = newShearStressShift; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class represents a "simple" interface material which is linear elastic in 
 * the normal direction. It is similar to a Coulomb contact model in the sense that 
 * the maximum shear stress is bounded by the normal stress times the coefficient 
 * of friction.
 * 
 * @author Milan Jirisek
 * @author Vit Smilauer
 * @author Jim Brouzoulis
 */
class IntMatCoulombContact : public StructuralInterfaceMaterial
{
protected:
    double kn;
    double stiffCoeff;
    double frictCoeff;
    /// Normal distance which needs to be closed when interface element should act in compression (distance is 0 by default).
    double normalClearance;

public:
    /// Constructor
    IntMatCoulombContact( int n, Domain *d );
    /// Destructor
    virtual ~IntMatCoulombContact();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatCoulombContact_Name; }

    virtual void giveEngTraction_3d( FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &jump, TimeStep *tStep);

    virtual void giveEngTraction_2d( FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &jump, TimeStep *tStep);

    virtual void giveEngTraction_1d( FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &jump, TimeStep *tStep);

    virtual void give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                           GaussPoint *gp, TimeStep *tStep)
    { this->giveGeneralStiffnessMatrix(answer, rMode, gp, tStep, 1); }

    virtual void give2dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                           GaussPoint *gp, TimeStep *tStep)
    { this->giveGeneralStiffnessMatrix(answer, rMode, gp, tStep, 2); }

    virtual void give1dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode rMode,
                                           GaussPoint *gp, TimeStep *tStep)
    { this->giveGeneralStiffnessMatrix(answer, rMode, gp, tStep, 3); }
    
    // This method returns the stiffness matrix according to the size of 
    // the spatial jump.
    void giveGeneralStiffnessMatrix(FloatMatrix &answer,
                                      MatResponseMode rMode,
                                      GaussPoint *gp,
                                      TimeStep *tStep, int numSpaceDim);

    void computeEngTraction(double &normalStress, FloatArray &shearStress, 
                             FloatArray &tempShearStressShift,
                             double normalJump, const FloatArray &shearJump );

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual bool hasAnalyticalTangentStiffness( ) const { return true; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatCoulombContactStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
