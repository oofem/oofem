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

#ifndef abaqususermaterial_h
#define abaqususermaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for AbaqusUserMaterial
//@{
#define _IFT_AbaqusUserMaterial_Name "abaqususermaterial"
#define _IFT_AbaqusUserMaterial_numState "numstate"
#define _IFT_AbaqusUserMaterial_properties "properties"
#define _IFT_AbaqusUserMaterial_initialStress "initialstress"
#define _IFT_AbaqusUserMaterial_userMaterial "umat"
#define _IFT_AbaqusUserMaterial_name "name"
#define _IFT_AbaqusUserMaterial_numericalTangent "numericaltangent"
#define _IFT_AbaqusUserMaterial_numericalTangentPerturbation "perturbation"
//@}

namespace oofem {
/**
 * This class allows for custom user materials from Abaqus (UMAT).
 * @note Experimental, subject to change. Many optional arguments haven't been dealt with properly (and few umat functions actually use all of them).
 * @note Nothing for large deformations have been tested.
 * @author Mikael Öhman
 *
 * The umat material should be compiled as a shared library. For example, compile like this;
 * @verbatim
 * gfortran -fPIC -fdefault-real-8 -shared -Wl,-soname,umat.so -o umat.so my_umat_code.f
 * @endverbatim
 *
 * The interface for a user material file for Abaqus is as follows:
 *    SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,
 *       SCD,RPL,DDSDDT,DRPLDE,DRPLDT,
 *       STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
 *       NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
 *       CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
 *
 *    CHARACTER*80 CMNAME
 *    DOUBLE PRECISION STRESS(NTENS),STATEV(NSTATV),
 *       DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
 *       STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
 *       PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
 */
class AbaqusUserMaterial : public StructuralMaterial
{
private:
    /// Gausspoint counter
    static int n;
    /// Dynamically loaded umat.
    void *umatobj;

    /// Pointer to the dynamically loaded umat-function (translated to C)
    void ( * umat )(double *stress, double *statev, double *ddsdde, double *sse, double *spd, // 5
                    double *scd, double *rpl, double *ddsddt, double *drplde, double *drpldt, // 5
                    double *stran, double *dstran, double time [ 2 ], double *dtime, double *temp, // 4
                    double *dtemp, double predef [ 1 ], double dpred [ 1 ], char cmname [ 80 ], int *ndi, // 5
                    int *nshr, int *ntens, int *nstatv, double *props, int *nprops, double coords [ 3 ], // 6
                    double *drot, double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1, // 5
                    int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc); // 6
    /// Name for material routine.
    char cmname [ 80 ];
    /// Size of the state vector.
    int numState;
    /// Material properties.
    FloatArray properties;
    /// Initial stress.
    FloatArray initialStress;
    /**
     * Flag to determine how the stress and Jacobian are interpreted.
     * 0 implies that the P and dPdF are returned from the umat routine.
     */
    int mStressInterpretation;

    /**
     * Flag to determine if numerical tangent should be used.
     */
    bool mUseNumericalTangent;

    /**
     * Size of perturbation if numerical tangent is used.
     */
    double mPerturbation;

    /// Name of the file that contains the umat function
    std :: string filename;

public:
    /// Constructor.
    AbaqusUserMaterial(int n, Domain * d);
    /// Destructor.
    virtual ~AbaqusUserMaterial();

    /**
     * Reads the following values;
     *  - numstate (required, integer): number of state variables.
     *  - properties (required, FloatArray): material property values.
     *  - umat (required, string): Filename of umat file dynamically library.
     *  - name (optional, string, default "umat"): Name of material model (used for input to umat routine).
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp,
                                                    TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx_dPdF(FloatMatrix &answer,
                                               MatResponseMode mmode, GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                            const FloatArray &reducedF, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int hasNonLinearBehaviour() { return true; }
    virtual const char *giveClassName() const { return "AbaqusUserMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_AbaqusUserMaterial_Name; }
};

class AbaqusUserMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Number of state variables.
    int numState;
    /// General state vector.
    FloatArray stateVector;
    /// Temporary state vector.
    FloatArray tempStateVector;
    /// Temporary elastic tangent.
    FloatMatrix tempTangent;

    /// Checker to see if tangent has been computed.
    bool hasTangentFlag;

public:
    /// Constructor.
    AbaqusUserMaterialStatus(int n, Domain * d, GaussPoint * gp, int numState);
    /// Destructor.
    virtual ~AbaqusUserMaterialStatus() { }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    bool hasTangent() const { return hasTangentFlag; }

    const FloatArray &giveStateVector() const { return stateVector; }
    FloatArray &letStateVectorBe(FloatArray &s) { return stateVector = s; }
    const FloatArray &giveTempStateVector() const { return tempStateVector; }
    FloatArray &letTempStateVectorBe(FloatArray &s) { return tempStateVector = s; }
    const FloatMatrix &giveTempTangent() { return tempTangent; }
    void letTempTangentBe(FloatMatrix t) {
        tempTangent = std :: move(t);
        hasTangentFlag = true;
    }

 

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "AbaqusUserMaterialStatus"; }
};
} // end namespace oofem
#endif // abaqususermaterial_h
