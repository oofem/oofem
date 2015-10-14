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

#ifndef structuralpythonmaterial_h
#define structuralpythonmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

///@name Input fields for StructuralPythonMaterial
//@{
#define _IFT_StructuralPythonMaterial_Name "structuralpythonmaterial"
#define _IFT_StructuralPythonMaterial_moduleName "module" /// The name of the module with the supplied functions (i.e. the name of the python script, without file extension)
//@}

namespace oofem {
/**
 * Custom user supplied python scripts for material models.
 * The python module should contain the functions
 * @code{.py}
 * computeStress(oldStrain, oldStress, strain, state, time) # returns stress
 * computePK1Stress(oldF, oldP, F, state, time) # returns P
 * @endcode
 * and optionally
 * @code{.py}
 * computeStressTangent(strain, stress, state, time) # returns ds/de
 * computePK1StressTangent(F, P, state, time) # return dP/dF
 * @endcode
 * else numerical derivatives are used. The state variable should be a dictionary storing either doubles or arrays of doubles.
 * 
 * This code is still experimental, and needs extensive testing.
 * @author Mikael Öhman
 */
class StructuralPythonMaterial : public StructuralMaterial
{
private:
    /// Name of the file that contains the python function
    std :: string moduleName;
    /// The module containing the functions
    PyObject *mpModule;

    /// Compiled function for small deformations
    PyObject *smallDef;
    PyObject *smallDefTangent;

    /// Compiled function for large deformations
    PyObject *largeDef;
    PyObject *largeDefTangent;

    /// Numerical pertubation for numerical tangents
    double pert;
public:
    /// Constructor.
    StructuralPythonMaterial(int n, Domain * d);
    /// Destructor.
    virtual ~StructuralPythonMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    void callStressFunction(PyObject *func, const FloatArray &oldStrain, const FloatArray &oldStress, const FloatArray &strain, FloatArray &stress, PyObject *stateDict, PyObject *tempStateDict, TimeStep *tStep) const;
    void callTangentFunction(FloatMatrix &answer, PyObject *func, const FloatArray &strain, const FloatArray &stress, PyObject *stateDict, PyObject *tempStateDict, TimeStep *tStep) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp,
                                                    TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                            const FloatArray &reducedF, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int hasNonLinearBehaviour() { return true; }
    virtual const char *giveClassName() const { return "StructuralPythonMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_StructuralPythonMaterial_Name; }
};

class StructuralPythonMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Internal state variables
    PyObject *stateDict;
    PyObject *tempStateDict;

public:
    /// Constructor.
    StructuralPythonMaterialStatus(Domain * d, GaussPoint * gp);
    /// Destructor.
    virtual ~StructuralPythonMaterialStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    void reinitTempStateDictionary();

    PyObject *giveStateDictionary() { return stateDict; }
    PyObject *giveTempStateDictionary() { return tempStateDict; }

    virtual const char *giveClassName() const { return "StructuralPythonMaterialStatus"; }
};
} // end namespace oofem
#endif // structuralpythonmaterial_h
