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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //Conversion for lists
#include <pybind11/operators.h>
namespace py = pybind11;

#include <string>

#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "intarray.h"

#include "femcmpnn.h"
#include "timestep.h"
#include "domain.h"
#include "engngm.h"
#include "staggeredproblem.h"
#include "dof.h"
#include "dofmanager.h"
#include "element.h"
#include "Elements/structuralelement.h"
#include "datastream.h"

#include "generalboundarycondition.h"
#include "boundarycondition.h"
#include "initialcondition.h"
#include "function.h"
#include "material.h"
#include "integrationpointstatus.h"
#include "matstatus.h"
#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"

#include "crosssection.h"

#include "field.h"
#include "feinterpol.h"
#include "util.h"
#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "timer.h"

#include "classfactory.h"

#include "chartype.h"
#include "elementgeometrytype.h"
#include "internalstatetype.h"

#include "integrationrule.h"
#include "gausspoint.h"
#include "inputrecord.h"
#include "dynamicinputrecord.h"

#include "classfactory.h"
#include "unknownnumberingscheme.h"
#include "vtkxmlexportmodule.h"
#include "vtkmemoryexportmodule.h"
#include "homexportmodule.h"

#include "uniformgridfield.h"
#include "unstructuredgridfield.h"
#include "dofmanvalfield.h"
#include "pythonfield.h"
#include <iostream>
#include "oofemutil.h"


//test
void test (oofem::Element& e) {
    oofem::IntArray a = {1,2,3};

    e.giveDofManDofIDMask (1, a);
    fprintf (stderr, "test:");
    a.printYourself();
}

/*
    Trampoline classes
*/
template <class ElementBase = oofem::Element> class PyElement : public ElementBase {
        // inherit the constructor
public:
        using ElementBase::ElementBase;
        // trampoline (need one for each virtual method)
        int giveNumberOfDofs() override {
            PYBIND11_OVERLOAD (int, ElementBase, giveNumberOfDofs,);
        }
        int computeNumberOfDofs() override {
            PYBIND11_OVERLOAD (int, ElementBase, computeNumberOfDofs,);
        }

        int giveNumberOfInternalDofManagers() const override {
            PYBIND11_OVERLOAD (int, ElementBase, giveNumberOfInternalDofManagers,);
        }
        oofem::DofManager *giveInternalDofManager(int i) const override {
            PYBIND11_OVERLOAD (oofem::DofManager*, ElementBase, giveInternalDofManager,i);
        }
        void giveCharacteristicMatrix(oofem::FloatMatrix &answer, oofem::CharType type, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, ElementBase, giveCharacteristicMatrix,std::ref(answer),type, tStep);
        }
        void giveCharacteristicVector(oofem::FloatArray &answer, oofem::CharType type, oofem::ValueModeType mode, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, ElementBase, giveCharacteristicVector,std::ref(answer), type, mode, tStep);
        }
        double giveCharacteristicValue(oofem::CharType type, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD (double, ElementBase, giveCharacteristicValue, type, tStep);
        }
        bool computeGtoLRotationMatrix(oofem::FloatMatrix &answer) override {
            PYBIND11_OVERLOAD (bool, ElementBase, computeGtoLRotationMatrix,std::ref(answer));
        }
        void giveDofManDofIDMask(int inode, oofem::IntArray &answer) const override {
            PYBIND11_OVERLOAD (void, ElementBase, giveDofManDofIDMask, inode, std::ref(answer));
        }

        void giveInternalDofManDofIDMask(int inode, oofem::IntArray &answer) const override {
            PYBIND11_OVERLOAD (void, ElementBase, giveInternalDofManDofIDMask,inode, std::ref(answer));
        }
        double computeVolumeAround(oofem::GaussPoint *gp) override {
            PYBIND11_OVERLOAD (double, ElementBase, computeVolumeAround, gp);
        }
        oofem::FEInterpolation *giveInterpolation() const override {
            PYBIND11_OVERLOAD (oofem::FEInterpolation*, ElementBase, giveInterpolation,);
        }
        void printOutputAt(FILE *file, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, ElementBase, printOutputAt, file, tStep);
        }
        oofem::MaterialMode giveMaterialMode() override {
            PYBIND11_OVERLOAD (oofem::MaterialMode, ElementBase, giveMaterialMode,);
        }
        const char *giveClassName() const override {
            PYBIND11_OVERLOAD_PURE (const char*, ElementBase, giveClassName,);
        }
        const char *giveInputRecordName() const override {
            PYBIND11_OVERLOAD_PURE (const char*, ElementBase, giveInputRecordName,);
        }
    };

    template <class StructuralElementBase = oofem::StructuralElement> class PyStructuralElement : public PyElement<StructuralElementBase> {
        // inherit the constructor
    public:
        using PyElement<StructuralElementBase>::PyElement;
        // trampoline (need one for each virtual method)
        void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, StructuralElementBase, computeMassMatrix, std::ref(answer), tStep);
        }
        void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, StructuralElementBase, computeStiffnessMatrix, std::ref(answer), rMode, tStep);
        }
        void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, StructuralElementBase, computeInitialStressMatrix, std::ref(answer), tStep);
        }
        void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override {
            PYBIND11_OVERLOAD_PURE (void, StructuralElementBase, computeStressVector, std::ref(answer), std::ref(strain), gp, tStep);
        }
        void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override {
            PYBIND11_OVERLOAD_PURE (void, StructuralElementBase, computeBmatrixAt, gp, std::ref(answer), lowerIndx, upperIndx);
        }
        void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override {
            PYBIND11_OVERLOAD (void, StructuralElementBase, computeNmatrixAt, iLocCoord, std::ref(answer));
        }
        void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep) override {
            PYBIND11_OVERLOAD_PURE (void, StructuralElementBase, computeConstitutiveMatrixAt, std::ref(answer), rMode, gp, tStep);
        }
      void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override {
        PYBIND11_OVERLOAD (void, StructuralElementBase, giveInternalForcesVector, std::ref(answer), tStep, useUpdatedGpRecord);
      }
    };


    template <class IntegrationPointStatusBase = oofem::IntegrationPointStatus> class PyIntegrationPointStatus : public IntegrationPointStatusBase {
        // inherit the constructor
    public:
        using IntegrationPointStatusBase::IntegrationPointStatusBase;
        // trampoline (need one for each virtual method)
        const char *giveClassName() const override {
            PYBIND11_OVERLOAD_PURE (const char*, IntegrationPointStatusBase, giveClassName,);
        }
    };

    template <class MaterialStatusBase = oofem::MaterialStatus> class PyMaterialStatus : public PyIntegrationPointStatus<MaterialStatusBase> {
        // inherit the constructor
    public:
        using PyIntegrationPointStatus<MaterialStatusBase>::PyIntegrationPointStatus;
        // trampoline (need one for each virtual method)
        void printOutputAt(FILE *file, TimeStep *tStep) const override {
            PYBIND11_OVERLOAD (void, MaterialStatusBase, printOutputAt, file, tStep);
        }
        void initTempStatus() override {
            PYBIND11_OVERLOAD (void, MaterialStatusBase, initTempStatus, );
        }
        void updateYourself(TimeStep *tStep) override {
            PYBIND11_OVERLOAD (void, MaterialStatusBase, updateYourself, tStep);
        }
        bool giveMaterialProperty(int propID, double &value) override {
          PYBIND11_OVERLOAD (bool, MaterialStatusBase, giveMaterialProperty, propID, std::ref(value));
        }
        void setMaterialProperty(int propID, double value) override {
            PYBIND11_OVERLOAD (void, MaterialStatusBase, setMaterialProperty, propID, value);
        }

    };

    template <class StructuralMaterialStatusBase = oofem::StructuralMaterialStatus> class PyStructuralMaterialStatus : public PyMaterialStatus<StructuralMaterialStatusBase> {
        // inherit the constructor
        using PyMaterialStatus<StructuralMaterialStatusBase>::PyMaterialStatus;
    };


    template <class MaterialBase = oofem::Material> class PyMaterial : public MaterialBase {
        // inherit the constructor
        using MaterialBase::MaterialBase;
        // trampoline (need one for each virtual method)
        bool isCharacteristicMtrxSymmetric(oofem::MatResponseMode rMode) const override {
            PYBIND11_OVERLOAD(bool, MaterialBase, isCharacteristicMtrxSymmetric, rMode);
        }
        double give(int aProperty, oofem::GaussPoint *gp) const override {
            PYBIND11_OVERLOAD(bool, MaterialBase, give, aProperty, gp);
        }
        bool hasProperty(int aProperty, oofem::GaussPoint *gp) const override {
            PYBIND11_OVERLOAD(bool, MaterialBase, hasProperty, aProperty, gp);
        }
        void modifyProperty(int aProperty, double value, oofem::GaussPoint *gp) override {
            PYBIND11_OVERLOAD(void, MaterialBase, modifyProperty, aProperty, value, gp);
        }
        bool hasMaterialModeCapability(oofem::MaterialMode mode) const override   {
            PYBIND11_OVERLOAD(bool, MaterialBase, hasMaterialModeCapability, mode);
        }
        bool hasCastingTimeSupport() const override {
            PYBIND11_OVERLOAD(bool, MaterialBase, hasCastingTimeSupport, );
        }
        int setIPValue(const oofem::FloatArray &value, oofem::GaussPoint *gp, oofem::InternalStateType type) override {
            PYBIND11_OVERLOAD(int, MaterialBase, setIPValue, std::ref(value), gp, type);
        }
        int giveIPValue(oofem::FloatArray &answer, oofem::GaussPoint *gp, oofem::InternalStateType type, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(int, MaterialBase, giveIPValue, std::ref(answer), gp, type, tStep);
        }
        void initializeFrom(oofem::InputRecord &ir) override {
            PYBIND11_OVERLOAD(void, MaterialBase, initializeFrom, ir);
        }
        //virtual void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);
        //virtual void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);
        int checkConsistency() override {
            PYBIND11_OVERLOAD(int, MaterialBase, checkConsistency, );
        }
        void restoreConsistency(oofem::GaussPoint *gp) override {
            PYBIND11_OVERLOAD(void, MaterialBase, restoreConsistency, gp);
        }
        int initMaterial(oofem::Element *element) override {
            PYBIND11_OVERLOAD(int, MaterialBase, initMaterial, element);
        }
      /*
        /// Exposing the following won't work (opeened issue https://github.com/pybind/pybind11/issues/1962)
        /// However, giveStatus can be overriden and can create and set status properly. See test4.py.
        std::unique_ptr<oofem::IntegrationPointStatus> CreateStatus(oofem::GaussPoint *gp) const override {
        PYBIND11_OVERLOAD(std::unique_ptr<oofem::IntegrationPointStatus>, MaterialBase, CreateStatus, gp);
        }
      */

      oofem::MaterialStatus* giveStatus (oofem::GaussPoint* gp) const override {
        PYBIND11_OVERLOAD(oofem::MaterialStatus*, MaterialBase, giveStatus, gp);
      }
        void initTempStatus(oofem::GaussPoint *gp) const override {
            PYBIND11_OVERLOAD(void, MaterialBase, initTempStatus, gp);
        }
        const char *giveClassName() const override {
            PYBIND11_OVERLOAD_PURE (const char*, MaterialBase, giveClassName,);
        }
        const char *giveInputRecordName() const override {
            PYBIND11_OVERLOAD_PURE (const char*, MaterialBase, giveInputRecordName,);
        }


    };

    
    class PyField : public oofem::Field {
        // inherit the constructor
        using oofem::Field::Field;
        // trampoline (need one for each virtual method
        int evaluateAt(oofem::FloatArray &answer, const oofem::FloatArray &coords,
        oofem::ValueModeType mode, oofem::TimeStep *tStep) override  {
            PYBIND11_OVERLOAD_PURE(int, oofem::Field, evaluateAt, std::ref(answer), coords, mode, tStep);
        }
        int evaluateAt(oofem::FloatArray &answer, oofem::DofManager *dman,
        oofem::ValueModeType mode, oofem::TimeStep *tStep) override  {
            PYBIND11_OVERLOAD_PURE(int, oofem::Field, evaluateAt, std::ref(answer), dman, mode, tStep);
        }
        const char *giveClassName() const override {
            PYBIND11_OVERLOAD_PURE(const char*, oofem::Field, giveClassName, );
        }
        void saveContext(DataStream &stream) override {
            PYBIND11_OVERLOAD_PURE(void, oofem::Field, saveContext, stream);
        }
        void restoreContext(DataStream &stream) override {
            PYBIND11_OVERLOAD_PURE(void, oofem::Field, restoreContext, stream);
        }
        
        
        
    };
    
    
    template <class StructuralMaterialBase = oofem::StructuralMaterial> class PyStructuralMaterial : public PyMaterial<StructuralMaterialBase> {
        // inherit the constructor
        using PyMaterial<StructuralMaterial>::PyMaterial;
        // trampoline (need one for each virtual method)
        bool hasMaterialModeCapability(oofem::MaterialMode mode) const override {
            PYBIND11_OVERLOAD(bool, StructuralMaterialBase, hasMaterialModeCapability, mode);
        }
        const char *giveClassName() const override {
            PYBIND11_OVERLOAD(const char*, StructuralMaterialBase, giveClassName, );
        }
        void initializeFrom(oofem::InputRecord &ir) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, initializeFrom, ir);
        }
        void giveInputRecord(oofem::DynamicInputRecord &input) override {
           PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveInputRecord, input);
        }
        void giveStiffnessMatrix(oofem::FloatMatrix &answer, oofem::MatResponseMode mode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveStiffnessMatrix, std::ref(answer), mode, gp, tStep);
        }
        void giveRealStressVector (oofem::FloatArray &answer, oofem::GaussPoint *gp, const oofem::FloatArray &reducedStrain, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveRealStressVector, std::ref(answer), gp, reducedStrain, tStep);
        }
        oofem::FloatArrayF<6> giveRealStressVector_3d(const oofem::FloatArrayF<6> &E, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            FloatArray strain = E;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_3d, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_3d
        oofem::FloatArrayF<4> giveRealStressVector_PlaneStrain(const oofem::FloatArrayF<4> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_PlaneStrain, strain, gp, tStep);
        }
        /// Iteratively calls giveRealStressVector_3d to find the stress controlled equal to zero·
        oofem::FloatArray giveRealStressVector_StressControl(const oofem::FloatArray &reducedE, const oofem::IntArray &strainControl, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_StressControl, reducedE, strainControl, gp, tStep);
        }
        oofem::FloatArray giveRealStressVector_ShellStressControl(const oofem::FloatArray &reducedE, const oofem::IntArray &strainControl, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_ShellStressControl, reducedE, strainControl, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<3> giveRealStressVector_PlaneStress(const oofem::FloatArrayF<3> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_PlaneStress, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<1> giveRealStressVector_1d(const oofem::FloatArrayF<1> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_1d, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<2> giveRealStressVector_Warping(const oofem::FloatArrayF<2> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_Warping, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<2> giveRealStressVector_2dBeamLayer(const oofem::FloatArrayF<2> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_2dBeamLayer, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<5> giveRealStressVector_PlateLayer(const oofem::FloatArrayF<5> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_PlateLayer, strain, gp, tStep);
        }
        /// Default implementation relies on giveRealStressVector_StressControl
        oofem::FloatArrayF<3> giveRealStressVector_Fiber(const oofem::FloatArrayF<3> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_Fiber, strain, gp, tStep);
        }
        oofem::FloatArrayF<3> giveRealStressVector_2dPlateSubSoil(const oofem::FloatArrayF<3> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_2dPlateSubSoil, strain, gp, tStep);
        }
        oofem::FloatArrayF<6> giveRealStressVector_3dBeamSubSoil(const oofem::FloatArrayF<6> &reducedE, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            oofem::FloatArray strain = reducedE;
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveRealStressVector_3dBeamSubSoil, strain, gp, tStep);
        }
        oofem::FloatArrayF<9> giveFirstPKStressVector_3d(const oofem::FloatArrayF<9> &reducedF, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveFirstPKStressVector_3d, reducedF, gp, tStep);
        }
        /// Default implementation relies on giveFirstPKStressVector_3d
        oofem::FloatArrayF<5> giveFirstPKStressVector_PlaneStrain(const oofem::FloatArrayF<5> &reducedF, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveFirstPKStressVector_PlaneStrain, reducedF, gp, tStep);
        }
        /// Default implementation relies on giveFirstPKStressVector_3d
        oofem::FloatArrayF<4> giveFirstPKStressVector_PlaneStress(const oofem::FloatArrayF<4> &reducedF, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveFirstPKStressVector_PlaneStress, reducedF, gp, tStep);
        }
        /// Default implementation relies on giveFirstPKStressVector_3d
        oofem::FloatArrayF<1> giveFirstPKStressVector_1d(const oofem::FloatArrayF<1> &reducedF, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveFirstPKStressVector_1d, reducedF, gp, tStep);
        }

        void giveCauchyStressVector_3d(oofem::FloatArray &answer, oofem::GaussPoint *gp, const oofem::FloatArray &reducedF, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveCauchyStressVector_3d, std::ref(answer), gp, reducedF, tStep);
        }
        void giveCauchyStressVector_PlaneStrain(oofem::FloatArray &answer, oofem::GaussPoint *gp, const oofem::FloatArray &reducedF, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveCauchyStressVector_PlaneStrain, std::ref(answer), gp, reducedF, tStep);
        }
        void giveCauchyStressVector_PlaneStress(oofem::FloatArray &answer, oofem::GaussPoint *gp, const oofem::FloatArray &reducedF, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveCauchyStressVector_PlaneStress, std::ref(answer), gp, reducedF, tStep);
        }
        void giveCauchyStressVector_1d(oofem::FloatArray &answer, oofem::GaussPoint *gp, const oofem::FloatArray &reducedF, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, giveCauchyStressVector_1d, std::ref(answer), gp, reducedF, tStep);
        }

        oofem::FloatArrayF<6> giveThermalDilatationVector(oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, giveThermalDilatationVector, gp, tStep);
        }
        oofem::FloatArray computeStressIndependentStrainVector(oofem::GaussPoint *gp, oofem::TimeStep *tStep, oofem::ValueModeType mode) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, computeStressIndependentStrainVector, gp, tStep, mode);
        }
#if 0
        // I don't think this method ought to be overloaded / Mikael
        oofem::FloatArrayF<6> computeStressIndependentStrainVector_3d(oofem::GaussPoint *gp, oofem::TimeStep *tStep, oofem::ValueModeType mode) const override {
            PYBIND11_OVERLOAD(oofem::FloatArray, StructuralMaterialBase, computeStressIndependentStrainVector_3d, gp, tStep, mode);
        }
#endif

        oofem::FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(oofem::MatResponseMode mode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give3dMaterialStiffnessMatrix, mode, gp, tStep);
        }
        oofem::FloatMatrixF<9,9> give3dMaterialStiffnessMatrix_dPdF(oofem::MatResponseMode mode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give3dMaterialStiffnessMatrix_dPdF, mode, gp, tStep);
        }
        void give3dMaterialStiffnessMatrix_dCde(oofem::FloatMatrix &answer, oofem::MatResponseMode mode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, give3dMaterialStiffnessMatrix_dCde, std::ref(answer), mode, gp, tStep);
        }

        FloatMatrixF<3,3> givePlaneStressStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, givePlaneStressStiffMtrx, mmode, gp, tStep);
        }
        oofem::FloatMatrixF<4,4> givePlaneStressStiffMtrx_dPdF(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, givePlaneStressStiffMtrx_dPdF, mmode, gp, tStep);
        }
        void givePlaneStressStiffMtrx_dCde(oofem::FloatMatrix &answer, oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, givePlaneStressStiffMtrx_dCde, std::ref(answer), mmode, gp, tStep);
        }

        oofem::FloatMatrixF<4,4> givePlaneStrainStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, givePlaneStrainStiffMtrx, mmode, gp, tStep);
        }
        oofem::FloatMatrixF<5,5> givePlaneStrainStiffMtrx_dPdF(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, givePlaneStrainStiffMtrx_dPdF, mmode, gp, tStep);
        }
        void givePlaneStrainStiffMtrx_dCde(oofem::FloatMatrix &answer, oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, givePlaneStrainStiffMtrx_dCde, std::ref(answer), mmode, gp, tStep);
        }

        oofem::FloatMatrixF<1,1> give1dStressStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give1dStressStiffMtrx, mmode, gp, tStep);
        }

        oofem::FloatMatrixF<1,1> give1dStressStiffMtrx_dPdF(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give1dStressStiffMtrx_dPdF, mmode, gp, tStep);
        }
        void give1dStressStiffMtrx_dCde(oofem::FloatMatrix &answer, oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) override {
            PYBIND11_OVERLOAD(void, StructuralMaterialBase, give1dStressStiffMtrx_dCde, std::ref(answer), mmode, gp, tStep);
        }

        oofem::FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give2dBeamLayerStiffMtrx, mmode, gp, tStep);
        }
        oofem::FloatMatrixF<3,3> give2dPlateSubSoilStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give2dPlateSubSoilStiffMtrx, mmode, gp, tStep);
        }
        oofem::FloatMatrixF<6,6> give3dBeamSubSoilStiffMtrx(oofem::MatResponseMode mmode, oofem::GaussPoint *gp, oofem::TimeStep *tStep) const override {
            PYBIND11_OVERLOAD(oofem::FloatMatrix, StructuralMaterialBase, give3dBeamSubSoilStiffMtrx, mmode, gp, tStep);
        }
     };

PYBIND11_MODULE(oofempy, m) {
    m.doc() = "oofem python bindings module"; // optional module docstring

    py::class_<oofem::FloatArray>(m, "FloatArray")
        .def(py::init<int>(), py::arg("n")=0)
        .def(py::init([](py::sequence s){
            oofem::FloatArray* ans = new oofem::FloatArray((int) py::len(s));
            for (unsigned int i=0; i<py::len(s); i++) {
                (*ans)[i]=s[i].cast<double>();
            }
            return ans;
        }
        ))
        .def("printYourself", (void (oofem::FloatArray::*)() const) &oofem::FloatArray::printYourself, "Prints receiver")
        .def("printYourself", (void (oofem::FloatArray::*)(const std::string &) const) &oofem::FloatArray::printYourself, "Prints receiver")
        .def("pY", &oofem::FloatArray::pY)
        .def("__setitem__", [](oofem::FloatArray &s, size_t i, double v) {
            s[i] = v;
        })
        .def("__getitem__", [](const oofem::FloatArray &s, size_t i) {
            if (i >= (size_t) s.giveSize()) throw py::index_error();
            return s[i];
        })
        .def("__repr__",
            [](const oofem::FloatArray &s) {
                std::ostringstream streamObj;
                std::string strObj;
                std::string a = "<oofempy.FloatArray: {";
                for ( int i = 0; i < s.giveSize(); ++i ) {
                    if ( i > 40 ) {
                        a.append("...");
                        break;
                    } else {
                        streamObj.str("");
                        streamObj.clear();
                        streamObj << s[i];//convert to scientific notation if necessary
                        strObj = streamObj.str();
                        a.append(strObj);
                        //a.append(std::to_string(s[i]));
                        a.append(", ");
                    }
                }
                a.append("}>");
                return a;
        })
        .def("resize", &oofem::FloatArray::resize)
        .def("assemble", &oofem::FloatArray::assemble)
        .def("distance", (double (oofem::FloatArray::*)(const oofem::FloatArray &) const) &oofem::FloatArray::distance)
        .def("normalize", &oofem::FloatArray::normalize)
        .def("computeNorm", &oofem::FloatArray::computeNorm)
        .def("product", &oofem::FloatArray::product)
        .def("zero", &oofem::FloatArray::zero)
        .def("beProductOf", &oofem::FloatArray::beProductOf)

        // expose FloatArray operators
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * float())
        .def(float() * py::self)
        .def(py::self += py::self)
        .def(py::self -= py::self)
        .def(py::self *= float())
        ;
     py::implicitly_convertible<py::sequence, oofem::FloatArray>();

     py::class_<oofem::FloatMatrix>(m, "FloatMatrix")
        .def(py::init<>())
        .def(py::init<int,int>())
        .def("printYourself", (void (oofem::FloatMatrix::*)() const) &oofem::FloatMatrix::printYourself, "Prints receiver")
        .def("printYourself", (void (oofem::FloatMatrix::*)(const std::string &) const) &oofem::FloatMatrix::printYourself, "Prints receiver")
        .def("pY", &oofem::FloatMatrix::pY)
        .def("__setitem__", [](oofem::FloatMatrix &s, py::tuple indx, double v) {
            s(((py::handle)(indx[0])).cast<int>(), ((py::handle)(indx[1])).cast<int>()) = v;
        })
        .def("__getitem__", [](const oofem::FloatMatrix &s, py::tuple indx) {
            if (((py::handle)(indx[0])).cast<int>() >= s.giveNumberOfRows()) throw py::index_error();
            if (((py::handle)(indx[1])).cast<int>() >= s.giveNumberOfColumns()) throw py::index_error();
            return s(((py::handle)(indx[0])).cast<int>(), ((py::handle)(indx[1])).cast<int>());
        })
        .def("__repr__",
            [](const oofem::FloatArray &s) {
                std::string a = "<oofempy.FloatArray: {";
                for ( int i = 0; i < s.giveSize(); ++i ) {
                    if ( i > 40 ) {
                        a.append("...");
                        break;
                    } else {
                        a.append(std::to_string(s[i]));
                        a.append(", ");
                    }
                }
                a.append("}>");
                return a;
            })
        .def("resize", &oofem::FloatMatrix::resize)
        .def("isSquare", &oofem::FloatMatrix::isSquare)
        .def("assemble", (void (oofem::FloatMatrix::*)(const FloatMatrix &, const IntArray &)) &oofem::FloatMatrix::assemble, "Assembles the contribution to the receiver")
        .def("assemble", (void (oofem::FloatMatrix::*)(const FloatMatrix &, const IntArray &, const IntArray&)) &oofem::FloatMatrix::assemble, "Assembles the contribution to the receiver")
        .def("computeFrobeniusNorm", &oofem::FloatMatrix::computeFrobeniusNorm)
        .def("computeNorm", &oofem::FloatMatrix::computeNorm)
        .def("beDiagonal", &oofem::FloatMatrix::beDiagonal)
        .def("giveTrace", &oofem::FloatMatrix::giveTrace)
        .def("giveDeterminant", &oofem::FloatMatrix::giveDeterminant)
        .def("zero", &oofem::FloatMatrix::zero)
        .def("beUnitMatrix", &oofem::FloatMatrix::beUnitMatrix)
        .def("beTranspositionOf", &oofem::FloatMatrix::beTranspositionOf)
        .def("beProductOf", &oofem::FloatMatrix::beProductOf)
        .def("addProductOf", &oofem::FloatMatrix::addProductOf)
        .def("addTProductOf", &oofem::FloatMatrix::addTProductOf)
        .def("beTProductOf", &oofem::FloatMatrix::beTProductOf)
        .def("beProductTOf", &oofem::FloatMatrix::beProductTOf)
        .def("beDyadicProductOf", &oofem::FloatMatrix::beDyadicProductOf)
        .def("beInverseOf", &oofem::FloatMatrix::beInverseOf)
        .def("solveForRhs", (bool (oofem::FloatMatrix::*)(const FloatArray &, FloatArray &, bool)) &oofem::FloatMatrix::solveForRhs)
        .def("solveForRhs", (void (oofem::FloatMatrix::*)(const FloatMatrix &, FloatMatrix &, bool)) &oofem::FloatMatrix::solveForRhs)
        .def("plusProductSymmUpper", &oofem::FloatMatrix::plusProductSymmUpper)
        .def("plusDyadSymmUpper", &oofem::FloatMatrix::plusDyadSymmUpper)
        .def("plusProductUnsym", &oofem::FloatMatrix::plusProductUnsym)
        .def("plusDyadUnsym", &oofem::FloatMatrix::plusDyadUnsym)
        // expose FloatArray operators
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * py::self)
        .def(py::self * FloatArray())
        .def(py::self += py::self)
        .def(py::self -= py::self)
        ;

    py::class_<oofem::IntArray>(m, "IntArray")

        .def(py::init<int>(), py::arg("n")=0)
        .def(py::init<const oofem::IntArray&>())
        .def(py::init([](py::sequence s){
            oofem::IntArray* ans = new oofem::IntArray((int) py::len(s));
            for (unsigned int i=0; i<py::len(s); i++) {
                (*ans)[i]=s[i].cast<int>();
            }
            return ans;
        }
        ))
        .def("resize", &oofem::IntArray::resize)
        .def("clear", &oofem::IntArray::clear)
        .def("preallocate", &oofem::IntArray::preallocate)
        .def("followedBy", (void (oofem::IntArray::*)(int , int)) &oofem::IntArray::followedBy, "Appends number to receiver")
        .def("followedBy", (void (oofem::IntArray::*)(const oofem::IntArray& , int)) &oofem::IntArray::followedBy, "Appends array to receiver")
        .def("isEmpty", &oofem::IntArray::isEmpty)
        .def("containsOnlyZeroes", &oofem::IntArray::containsOnlyZeroes)
        .def("containsSorted", &oofem::IntArray::containsSorted)
        .def("insertSorted", &oofem::IntArray::insertSorted)
        .def("eraseSorted", &oofem::IntArray::eraseSorted)
        .def("findFirstIndexOf", &oofem::IntArray::findFirstIndexOf, "Finds index of first occurrence of given value in array")
        .def("contains", &oofem::IntArray::contains)
        .def("sort", &oofem::IntArray::sort)
        .def("zero", &oofem::IntArray::zero)
        .def("pY", &oofem::IntArray::pY)
        .def("__repr__",
            [](const oofem::IntArray &s) {
                std::string a = "<oofempy.IntArray: {";
                for ( int i = 0; i < s.giveSize(); ++i ) {
                    if ( i > 40 ) {
                        a.append("...");
                        break;
                    } else {
                        a.append(std::to_string(s[i]));
                        a.append(", ");
                    }
                }
                a.append("}>");
                return a;
            })
        .def("__setitem__", [](oofem::IntArray &s, size_t i, int v) {
            s[i] = v;
        })
        .def("__getitem__", [](const oofem::IntArray &s, size_t i) {
            if (i >= (size_t) s.giveSize()) throw py::index_error();
            return s[i];
        })

    ;
    py::implicitly_convertible<py::sequence, oofem::IntArray>();


    py::class_<oofem::DataReader>(m, "DataReader")
    ;

    py::class_<oofem::OOFEMTXTDataReader, oofem::DataReader>(m, "OOFEMTXTDataReader")
        .def(py::init<std::string>())
    ;


    py::class_<oofem::InputRecord>(m, "InputRecord")
    ;

    typedef const char *InputFieldType;
    py::class_<oofem::DynamicInputRecord, oofem::InputRecord>(m, "DynamicInputRecord")
        .def(py::init<std::string, int>(), py::arg("answer") = "", py::arg("value")=0)
        .def("finish", &oofem::DynamicInputRecord::finish, py::arg("wrn")=true)
        .def("setRecordKeywordField", &oofem::DynamicInputRecord::setRecordKeywordField)
        .def("setRecordKeywordNumber", &oofem::DynamicInputRecord::setRecordKeywordNumber)
        .def("setField", (void (oofem::DynamicInputRecord::*)(int, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(double, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(bool, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(std::string, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatArray, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatMatrix, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::FloatMatrix, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(oofem::IntArray, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(std::vector<std::string>, InputFieldType)) &oofem::DynamicInputRecord::setField)
        .def("setField", (void (oofem::DynamicInputRecord::*)(InputFieldType)) &oofem::DynamicInputRecord::setField)
    ;

    py::class_<oofem::OOFEMTXTInputRecord, oofem::InputRecord>(m, "OOFEMTXTInputRecord")
        .def(py::init<>())
        .def(py::init<int, std::string>())
        .def("finish", &oofem::OOFEMTXTInputRecord::finish, py::arg("wrn")=true)
        .def("setRecordString", &oofem::OOFEMTXTInputRecord::setRecordString)
    ;


    py::class_<oofem::FEMComponent>(m, "FEMComponent")
        .def("giveClassName", &oofem::FEMComponent::giveClassName)
        .def("giveInputRecordName", &oofem::FEMComponent::giveInputRecordName)
        .def("giveDomain", &oofem::FEMComponent::giveDomain, py::return_value_policy::reference)
        .def("setDomain", &oofem::FEMComponent::setDomain)
        .def("giveNumber", &oofem::FEMComponent::giveNumber)
        .def("setNumber", &oofem::FEMComponent::setNumber)
        .def("checkConsistency", &oofem::FEMComponent::checkConsistency)
        .def("printYourself", &oofem::FEMComponent::printYourself)
        .def_property_readonly("number", &oofem::FEMComponent::giveNumber)

        ;

    py::class_<oofem::MetaStep>(m, "MetaStep")
        .def("setNumberOfSteps", &oofem::MetaStep::setNumberOfSteps)
        .def("giveNumberOfSteps", &oofem::MetaStep::giveNumberOfSteps)
    ;
        
    py::class_<oofem::TimeStep>(m, "TimeStep")
        .def("giveNumber", &oofem::TimeStep::giveNumber)
        .def("giveTargetTime", &oofem::TimeStep::giveTargetTime)
        .def("giveIntrinsicTime", &oofem::TimeStep::giveIntrinsicTime)
        .def("giveTimeIncrement", &oofem::TimeStep::giveTimeIncrement)
        .def("setTimeIncrement", &oofem::TimeStep::setTimeIncrement)
        .def("setTime", &oofem::TimeStep::setTime)
        .def("setTargetTime", &oofem::TimeStep::setTargetTime)
        .def("setIntrinsicTime", &oofem::TimeStep::setIntrinsicTime)
        .def("isNotTheLastStep", &oofem::TimeStep::isNotTheLastStep)
        .def("isTheFirstStep", &oofem::TimeStep::isTheFirstStep)
        .def("isTheCurrentTimeStep", &oofem::TimeStep::isTheCurrentTimeStep)
        ;

    py::class_<oofem::FieldManager>(m, "FieldManager")
        .def("registerField", &oofem::FieldManager::registerField, py::keep_alive<1, 2>())
        ;

    py::class_<oofem::EngngModelContext>(m, "EngngModelContext")
        .def("giveFieldManager", &oofem::EngngModelContext::giveFieldManager, py::return_value_policy::reference)

        ;

    py::class_<oofem::EngngModel>(m, "EngngModel")
        .def("giveDomain", &oofem::EngngModel::giveDomain, py::return_value_policy::reference)
        .def("setDomain", &oofem::EngngModel::setDomain, py::keep_alive<1, 3>())
        .def("giveNumberOfDomains", &oofem::EngngModel::giveNumberOfDomains)
        .def("giveMetaStep", &oofem::EngngModel::giveMetaStep, py::return_value_policy::reference)
        .def("terminateAnalysis", &oofem::EngngModel::terminateAnalysis)
        .def("terminate", &oofem::EngngModel::terminate)
        .def("solveYourself", &oofem::EngngModel::solveYourself)
        .def("solveYourselfAt", &oofem::EngngModel::solveYourselfAt)
        .def("terminate",&oofem::EngngModel::terminate)
        .def("giveField", &oofem::EngngModel::giveField)
        //.def("giveCurrentStep", &oofem::EngngModel::giveCurrentStep, py::return_value_policy::reference)
        .def("giveCurrentStep", &oofem::EngngModel::giveCurrentStep, py::return_value_policy::reference, py::arg("force") = false)
        .def("givePreviousStep", &oofem::EngngModel::givePreviousStep, py::return_value_policy::reference)
        .def("giveNextStep", &oofem::EngngModel::giveNextStep, py::return_value_policy::reference)
        .def("giveNumberOfSteps", &oofem::EngngModel::giveNumberOfSteps, py::arg("force")=false)
        .def("giveUnknownComponent", &oofem::EngngModel::giveUnknownComponent)
        .def("checkProblemConsistency", &oofem::EngngModel::checkProblemConsistency)
        .def("giveTimer", &oofem::EngngModel::giveTimer, py::return_value_policy::reference)
        .def("init", &oofem::EngngModel::init)
        .def("initializeYourself", &oofem::EngngModel::initializeYourself)
        .def("initMetaStepAttributes", &oofem::EngngModel::initMetaStepAttributes)
        .def("preInitializeNextStep", &oofem::EngngModel::preInitializeNextStep)
        .def("updateYourself", &oofem::EngngModel::updateYourself)
        .def("postInitialize", &oofem::EngngModel::postInitialize)
        .def("setRenumberFlag", &oofem::EngngModel::setRenumberFlag)
        .def("giveContext", &oofem::EngngModel::giveContext, py::return_value_policy::reference)
        .def("forceEquationNumbering", py::overload_cast<int>(&oofem::EngngModel::forceEquationNumbering))
        .def("forceEquationNumbering", py::overload_cast<>(&oofem::EngngModel::forceEquationNumbering))
        ;
    
    py::class_<oofem::StaggeredProblem, oofem::EngngModel>(m, "StaggeredProblem")
        .def("giveSlaveProblem", &oofem::StaggeredProblem::giveSlaveProblem, py::return_value_policy::reference)
        .def("giveTimeControl", &oofem::StaggeredProblem::giveTimeControl, py::return_value_policy::reference)
    ;
        
    py::class_<oofem::Domain>(m, "Domain")
        .def(py::init<int, int, oofem::EngngModel*>())
        .def("giveNumber", &oofem::Domain::giveNumber)
        .def("setNumber", &oofem::Domain::setNumber)
        .def("giveDofManager", &oofem::Domain::giveDofManager, py::return_value_policy::reference)
        //.def("giveDofManagers", &oofem::Domain::giveDofManagers)
        .def("giveElement", &oofem::Domain::giveElement, py::return_value_policy::reference)
        //.def("giveElements", &oofem::Domain::giveElements, py::return_value_policy::reference)
        .def("giveNumberOfElements", &oofem::Domain::giveNumberOfElements)
        .def("giveNumberOfDofManagers", &oofem::Domain::giveNumberOfDofManagers)
        .def("giveNumberOfBoundaryConditions", &oofem::Domain::giveNumberOfBoundaryConditions)
        .def("giveNumberOfFunctions", &oofem::Domain::giveNumberOfFunctions)
        .def("giveNumberOfMaterialModels", &oofem::Domain::giveNumberOfMaterialModels)
        .def("giveNumberOfCrossSectionModels", &oofem::Domain::giveNumberOfCrossSectionModels)
        .def("giveNumberOfInitialConditions", &oofem::Domain::giveNumberOfInitialConditions)
        .def("giveNumberOfRegions", &oofem::Domain::giveNumberOfRegions)
        .def("giveNumberOfNonlocalBarriers", &oofem::Domain::giveNumberOfNonlocalBarriers)
        .def("giveNumberOfSets", &oofem::Domain::giveNumberOfSets)
        .def("giveBc", &oofem::Domain::giveBc, py::return_value_policy::reference)
        .def("giveIc", &oofem::Domain::giveIc, py::return_value_policy::reference)
        .def("giveFunction", &oofem::Domain::giveFunction, py::return_value_policy::reference)
        .def("giveMaterial", &oofem::Domain::giveMaterial, py::return_value_policy::reference)
        .def("giveCrossSection", &oofem::Domain::giveCrossSection, py::return_value_policy::reference)
        .def("resizeDofManagers", &oofem::Domain::resizeDofManagers)
        .def("setDofManager", &oofem::Domain::py_setDofManager, py::keep_alive<0, 2>())
        .def("resizeElements", &oofem::Domain::resizeElements)
        .def("setElement", &oofem::Domain::py_setElement, py::keep_alive<0, 2>())
        .def("resizeMaterials", &oofem::Domain::resizeMaterials)
        .def("setMaterial", &oofem::Domain::py_setMaterial, py::keep_alive<0, 2>())
        .def("resizeCrossSectionModels", &oofem::Domain::resizeCrossSectionModels)
        .def("setCrossSection", &oofem::Domain::py_setCrossSection, py::keep_alive<0, 2>())
        .def("resizeBoundaryConditions", &oofem::Domain::resizeBoundaryConditions)
        .def("setBoundaryCondition", &oofem::Domain::py_setBoundaryCondition, py::keep_alive<0, 2>())
        .def("resizeInitialConditions", &oofem::Domain::resizeInitialConditions)
        .def("setInitialCondition", &oofem::Domain::py_setInitialCondition, py::keep_alive<0, 2>())
        .def("resizeFunctions", &oofem::Domain::resizeFunctions)
        .def("setFunction", &oofem::Domain::py_setFunction, py::keep_alive<0, 2>())
        .def("resizeSets", &oofem::Domain::resizeSets)
        .def("setSet", &oofem::Domain::py_setSet, py::keep_alive<0, 2>())
    ;

    py::class_<oofem::Dof>(m, "Dof")
        .def("giveDofManNumber", &oofem::Dof::giveDofManNumber)
        .def("giveDofManager", &oofem::Dof::giveDofManager, py::return_value_policy::reference)
        .def("giveDofManGlobalNumber", &oofem::Dof::giveDofManGlobalNumber)
        .def("giveUnknown", (double (oofem::Dof::*)(oofem::ValueModeType, oofem::TimeStep*)) &oofem::Dof::giveUnknown)
        .def("giveDofID", &oofem::Dof::giveDofID)
        .def("printYourself", &oofem::Dof::printYourself)
    ;

    py::class_<oofem::DofManager, oofem::FEMComponent>(m, "DofManager")
        .def("giveDofWithID", &oofem::DofManager::giveDofWithID, py::return_value_policy::reference)
        .def("giveNumberOfDofs", &oofem::DofManager::giveNumberOfDofs)
        .def("giveUnknownVector", (void (oofem::DofManager::*)(oofem::FloatArray&, const oofem::IntArray&, oofem::ValueModeType, oofem::TimeStep*, bool)) &oofem::DofManager::giveUnknownVector)
        .def("givePrescribedUnknownVector", &oofem::DofManager::givePrescribedUnknownVector)
        .def("giveCoordinates", &oofem::DofManager::giveCoordinates, py::return_value_policy::reference)
        .def("appendDof", &oofem::DofManager::appendDof, py::keep_alive<1, 2>())
        .def("removeDof", &oofem::DofManager::removeDof)
        .def("hasDofID", &oofem::DofManager::hasDofID)
        .def("giveGlobalNumber", &oofem::DofManager::giveGlobalNumber)
    ;

    /*
        Element
    */

    py::class_<oofem::Element, oofem::FEMComponent, PyElement<>>(m, "Element")
        .def(py::init<int, oofem::Domain*>())
        .def("giveLocationArray", (void (oofem::Element::*)(oofem::IntArray &, const oofem::UnknownNumberingScheme &, oofem::IntArray *dofIds) const) &oofem::Element::giveLocationArray)
        .def("giveLocationArray", (void (oofem::Element::*)(oofem::IntArray &, const oofem::IntArray &, const oofem::UnknownNumberingScheme &, oofem::IntArray *) const) &oofem::Element::giveLocationArray)
        //.def("giveCharacteristicMatrix", &oofem::Element::giveCharacteristicMatrix)
        //.def("giveCharacteristicVector", &oofem::Element::giveCharacteristicVector)
        //.def("giveCharacteristicValue", &oofem::Element::giveCharacteristicValue)
        .def("computeVectorOf", (void (oofem::Element::*)(oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray &)) &oofem::Element::computeVectorOf)
        .def("computeVectorOf", (void (oofem::Element::*)(const oofem::IntArray &, oofem::ValueModeType , oofem::TimeStep*, oofem::FloatArray &, bool)) &oofem::Element::computeVectorOf)
        .def("computeVectorOfPrescribed", (void (oofem::Element::*)(oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray &)) &oofem::Element::computeVectorOfPrescribed)
        .def("computeVectorOfPrescribed", (void (oofem::Element::*)(const oofem::IntArray &, oofem::ValueModeType , oofem::TimeStep *, oofem::FloatArray & )) &oofem::Element::computeVectorOfPrescribed)
        .def("giveDofManagerNumber", &oofem::Element::giveDofManagerNumber)
        .def("giveDofManArray", &oofem::Element::giveDofManArray, py::return_value_policy::reference)
        .def("giveDofManager", &oofem::Element::giveDofManager, py::return_value_policy::reference)
        .def("giveMaterial", &oofem::Element::giveMaterial, py::return_value_policy::reference)
        .def("giveCrossSection", &oofem::Element::giveCrossSection, py::return_value_policy::reference)
        .def("giveMaterialNumber", &oofem::Element::giveMaterialNumber)
        .def("giveNumberOfDofManagers", &oofem::Element::giveNumberOfDofManagers)
        .def("giveNumberOfNodes", &oofem::Element::giveNumberOfNodes)
        .def("giveRegionNumber", &oofem::Element::giveRegionNumber)
        .def("giveNode", &oofem::Element::giveNode)
        .def("updateYourself", &oofem::Element::updateYourself)
        .def("isActivated", &oofem::Element::isActivated)
        .def("isCast", &oofem::Element::isCast)
        .def("giveGeometryType", &oofem::Element::giveGeometryType)
        .def("giveIntegrationRule", &oofem::Element::giveIntegrationRule, py::return_value_policy::reference)
        .def("giveDefaultIntegrationRulePtr", &oofem::Element::giveDefaultIntegrationRulePtr, py::return_value_policy::reference)
        .def("giveIPValue", &oofem::Element::giveIPValue)
        .def("giveLabel", &oofem::Element::giveLabel)
        .def("initializeFrom", &oofem::Element::initializeFrom)
        .def("setDofManagers", &oofem::Element::setDofManagers)
        .def("setNumberOfDofManagers", &oofem::Element::setNumberOfDofManagers)
        .def("giveDofManDofIDMask", &oofem::Element::giveDofManDofIDMask)
        .def("getActivityTimeFunctionNumber", &oofem::Element::getActivityTimeFunctionNumber)
        .def("setActivityTimeFunctionNumber", &oofem::Element::setActivityTimeFunctionNumber)
    ;

    py::class_<oofem::StructuralElement, oofem::Element, PyStructuralElement<>>(m, "StructuralElement")
        .def(py::init<int, oofem::Domain*>())
        .def("giveDofManDofIDMask", &oofem::StructuralElement::giveDofManDofIDMask)
    ;


    py::class_<oofem::GeneralBoundaryCondition, oofem::FEMComponent>(m, "GeneralBoundaryCondition")
        .def("getIsImposedTimeFunctionNumber", &oofem::GeneralBoundaryCondition::getIsImposedTimeFunctionNumber)
        .def("setIsImposedTimeFunctionNumber", &oofem::GeneralBoundaryCondition::setIsImposedTimeFunctionNumber)
    ;


    py::class_<oofem::BoundaryCondition, oofem::GeneralBoundaryCondition>(m, "BoundaryCondition")
        .def("setPrescribedValue", &oofem::BoundaryCondition::setPrescribedValue)
    ;

    py::class_<oofem::InitialCondition, oofem::FEMComponent>(m, "InitialCondition")
    ;

//     py::class_<oofem::Timer>(m, "Timer")
//     ;
    
    py::class_<oofem::EngngModelTimer> (m, "EngngModelTimer")
        .def("startTimer", &oofem::EngngModelTimer::startTimer)
        .def("stopTimer", &oofem::EngngModelTimer::stopTimer)
        .def("getUtime", &oofem::EngngModelTimer::getUtime)
    ;
    
    py::enum_<oofem::EngngModelTimer::EngngModelTimerType>(m, "EngngModelTimerType")
        .value("EMTT_AnalysisTimer", oofem::EngngModelTimer::EMTT_AnalysisTimer)
        .value("EMTT_SolutionStepTimer", oofem::EngngModelTimer::EMTT_SolutionStepTimer)
        .value("EMTT_NetComputationalStepTimer", oofem::EngngModelTimer::EMTT_NetComputationalStepTimer)
        .value("EMTT_LoadBalancingTimer", oofem::EngngModelTimer::EMTT_LoadBalancingTimer)
        .value("EMTT_DataTransferTimer", oofem::EngngModelTimer::EMTT_DataTransferTimer)
        .value("EMTT_LastTimer", oofem::EngngModelTimer::EMTT_LastTimer)
        .export_values()
    ;
    
    
    /*
        Function class
    */

    /* Trampoline classes allowing to define custom derived types from within Python */
    class PyFunction : public oofem::Function {
        // inherit the constructor
        using oofem::Function::Function;
        // trampoline (need one for each virtual method)
        double evaluateAtTime(double t) override {
            PYBIND11_OVERLOAD (
                double, // return type
                oofem::Function, // parent class
                evaluateAtTime, // name of method in c++ (must match python name)
                t // argument(s)
            );
        }
        double evaluateVelocityAtTime(double t) override {
            PYBIND11_OVERLOAD_PURE(
                double,
                oofem::Function,
                evaluateVelocityAtTime,
                t
            )
        }
        double evaluateAccelerationAtTime(double t) override {
            PYBIND11_OVERLOAD_PURE(
                double,
                oofem::Function,
                evaluateAccelerationAtTime,
                t
            )
        }
        const char* giveClassName() const override {
           PYBIND11_OVERLOAD_PURE(
                const char*,
                oofem::Function,
                giveClassName,
            )
        }
        const char* giveInputRecordName() const override {
           PYBIND11_OVERLOAD_PURE(
                const char*,
                oofem::Function,
                giveInputRecordName,
            )
        }


    };

    py::class_<oofem::Function, oofem::FEMComponent, PyFunction>(m, "Function")
    .def(py::init<int, oofem::Domain*>())
    .def("evaluate", (double (oofem::Function::*)(TimeStep *, ValueModeType)) &oofem::Function::evaluate)
    //abstract services
    .def ("evaluateAtTime", &oofem::Function::evaluateAtTime)
    .def ("evaluateVelocityAtTime", &oofem::Function::evaluateVelocityAtTime)
    .def ("evaluateAccelerationAtTime", &oofem::Function::evaluateVelocityAtTime)
    ;

    py::class_<oofem::IntegrationPointStatus, PyIntegrationPointStatus<>>(m, "IntegrationPointStatus")
      .def(py::init<oofem::GaussPoint*>())
        .def("giveClassName", &oofem::IntegrationPointStatus::giveClassName)
    ;
    py::class_<oofem::MaterialStatus, oofem::IntegrationPointStatus, PyMaterialStatus<>>(m, "MaterialStatus")
        .def(py::init<oofem::GaussPoint *>())
    ;

    py::class_<oofem::StructuralMaterialStatus, oofem::MaterialStatus, PyStructuralMaterialStatus<>>(m, "StructuralMaterialStatus")
        .def(py::init<oofem::GaussPoint *>())
        .def("giveStrainVector", &oofem::StructuralMaterialStatus::giveStrainVector)
        .def("giveStressVector", &oofem::StructuralMaterialStatus::giveStressVector)
        .def("letTempStressVectorBe", &oofem::StructuralMaterialStatus::letTempStressVectorBe)
        .def("letTempStrainVectorBe", &oofem::StructuralMaterialStatus::letTempStrainVectorBe)
    ;

    py::class_<oofem::Material, oofem::FEMComponent>(m, "Material")
        .def("giveStatus", &oofem::Material::giveStatus, py::return_value_policy::reference)
        .def("CreateStatus", &oofem::Material::CreateStatus, py::return_value_policy::reference)
        .def("giveIPValue", &oofem::Material::giveIPValue)
    ;

    py::class_<oofem::StructuralMaterial, oofem::Material, PyStructuralMaterial<>>(m, "StructuralMaterial")
        .def(py::init<int, oofem::Domain*>())
        .def("giveStatus", &oofem::StructuralMaterial::giveStatus, py::return_value_policy::reference)
        .def("CreateStatus", &oofem::StructuralMaterial::CreateStatus, py::return_value_policy::reference)
    ;

    py::class_<oofem::CrossSection, oofem::FEMComponent>(m, "CrossSection")
    ;

    py::class_<oofem::Set, oofem::FEMComponent>(m, "Set")
    ;

    
    py::class_<oofem::UnknownNumberingScheme>(m, "UnknownNumberingScheme")
    ;

    py::class_<oofem::IntegrationRule>(m, "IntegrationRule")
        .def("giveNumberOfIntegrationPoints", &oofem::IntegrationRule::giveNumberOfIntegrationPoints)
        .def("getIntegrationPoint", &oofem::IntegrationRule::getIntegrationPoint, py::return_value_policy::reference)
    ;
    py::class_<oofem::GaussPoint>(m, "GaussPoint")
      .def ("giveMaterialStatus", (const oofem::IntegrationPointStatus* (oofem::GaussPoint::*)() const) &oofem::GaussPoint::giveMaterialStatus, py::return_value_policy::reference)
      .def ("setMaterialStatus", (oofem::IntegrationPointStatus* (oofem::GaussPoint::*)(oofem::IntegrationPointStatus*)) &oofem::GaussPoint::setMaterialStatus, py::keep_alive<0, 2>())
    ;

    py::class_<oofem::ExportModule>(m, "ExportModule")
      .def("initialize", &oofem::ExportModule::initialize)
      .def("doOutput", &oofem::ExportModule::doOutput)
      .def("terminate", &oofem::ExportModule::terminate)
      ;
    
    py::class_<oofem::HOMExportModule, oofem::ExportModule>(m, "HOMExportModule")
      ;

    py::class_<oofem::VTKPiece>(m, "VTKPiece")
      .def("giveNumberOfNodes", &oofem::VTKPiece::giveNumberOfNodes)
      .def("giveNodeCoords", &oofem::VTKPiece::giveNodeCoords)
      .def("giveNumberOfCells", &oofem::VTKPiece::giveNumberOfCells)
      .def("giveCellConnectivity", &oofem::VTKPiece::giveCellConnectivity)
      .def("giveCellType", &oofem::VTKPiece::giveCellType)
      .def("giveCellOffset", &oofem::VTKPiece::giveCellOffset)
      .def("givePrimaryVarInNode", &oofem::VTKPiece::givePrimaryVarInNode)
      .def("giveInternalVarInNode", &oofem::VTKPiece::giveInternalVarInNode)
      .def("giveCellVar", &oofem::VTKPiece::giveCellVar)

      .def("getVertices", &oofem::VTKPiece::getVertices, py::return_value_policy::move)
      .def("getCellConnectivity", &oofem::VTKPiece::getCellConnectivity, py::return_value_policy::move)
      .def("getCellTypes", &oofem::VTKPiece::getCellTypes, py::return_value_policy::move)
      .def("getPrimaryVertexValues", &oofem::VTKPiece::getPrimaryVertexValues, py::return_value_policy::move)
      .def("getInternalVertexValues", &oofem::VTKPiece::getInternalVertexValues, py::return_value_policy::move)
      .def("getCellValues", &oofem::VTKPiece::getCellValues, py::return_value_policy::move)
      ;

    py::class_<oofem::VTKBaseExportModule, oofem::ExportModule>(m, "VTKBaseExportModule")
      ;
    
    py::class_<oofem::VTKXMLExportModule, oofem::VTKBaseExportModule>(m, "VTKXMLExportModule")
      .def("getVTKPieces", &oofem::VTKXMLExportModule::getVTKPieces,  py::return_value_policy::reference)
      ;

    py::class_<oofem::VTKMemoryExportModule, oofem::VTKBaseExportModule>(m, "VTKMemoryExportModule")
      .def("getVTKPieces", &oofem::VTKMemoryExportModule::getVTKPieces,  py::return_value_policy::reference)
      ;


    py::class_<oofem::ClassFactory>(m, "ClassFactory")
        .def("createElement", &oofem::ClassFactory::createElement)
        .def("createEngngModel", &oofem::ClassFactory::createEngngModel)
    ;

    m.def("getClassFactory", &oofem::GiveClassFactory, py::return_value_policy::reference);
    m.def("InstanciateProblem", &oofem::InstanciateProblem);
    //std::unique_ptr<EngngModel> InstanciateProblem(DataReader &dr, problemMode mode, int contextFlag, EngngModel *master = 0, bool parallelFlag = false);

    py::enum_<oofem::FieldType>(m, "FieldType")
        .value("FT_Unknown", oofem::FieldType::FT_Unknown)
        .value("FT_Velocity", oofem::FieldType::FT_Velocity)
        .value("FT_Displacements", oofem::FieldType::FT_Displacements)
        .value("FT_VelocityPressure", oofem::FieldType::FT_VelocityPressure)
        .value("FT_Pressure", oofem::FieldType::FT_Pressure)
        .value("FT_Temperature", oofem::FieldType::FT_Temperature)
        .value("FT_HumidityConcentration", oofem::FieldType::FT_HumidityConcentration)
        .value("FT_TransportProblemUnknowns", oofem::FieldType::FT_TransportProblemUnknowns)
        .value("FT_TemperatureAmbient", oofem::FieldType::FT_TemperatureAmbient)
    ;


    py::enum_<oofem::UnknownType>(m, "UnknownType")
      .value("DisplacementVector", oofem::UnknownType::DisplacementVector)
      .value("Temperature", oofem::UnknownType::Temperature)
      ;
    
    py::enum_<oofem::ValueModeType>(m, "ValueModeType")
        .value("VM_Unknown", oofem::ValueModeType::VM_Unknown)
        .value("VM_Total", oofem::ValueModeType::VM_Total)
        .value("VM_Velocity", oofem::ValueModeType::VM_Velocity)
        .value("VM_Acceleration", oofem::ValueModeType::VM_Acceleration)
        .value("VM_Incremental", oofem::ValueModeType::VM_Incremental)
        .value("VM_RhsTotal", oofem::ValueModeType::VM_RhsTotal)
        .value("VM_RhsIncremental", oofem::ValueModeType::VM_RhsIncremental)
        .value("VM_RhsInitial", oofem::ValueModeType::VM_RhsInitial)
        .value("VM_Intermediate", oofem::ValueModeType::VM_Intermediate)
        .value("VM_TotalIntrinsic", oofem::ValueModeType::VM_TotalIntrinsic)
    ;

    py::enum_<oofem::DofIDItem>(m, "DofIDItem")
      .value("Undef", oofem::DofIDItem::Undef)
      .value("D_u", oofem::DofIDItem::D_u)
      .value("D_v", oofem::DofIDItem::D_v)
      .value("D_w", oofem::DofIDItem::D_w)
      .value("R_u", oofem::DofIDItem::R_u)
      .value("R_v", oofem::DofIDItem::R_v)
      .value("R_w", oofem::DofIDItem::R_w)
      .value("V_u", oofem::DofIDItem::V_u)
      .value("V_v", oofem::DofIDItem::V_v)
      .value("V_w", oofem::DofIDItem::V_w)
      .value("T_f", oofem::DofIDItem::T_f)
      .value("P_f", oofem::DofIDItem::P_f)
      .value("G_0", oofem::DofIDItem::G_0)
      .value("G_1", oofem::DofIDItem::G_1)
      .value("C_1", oofem::DofIDItem::C_1)
      .value("W_u", oofem::DofIDItem::W_u)
      .value("W_v", oofem::DofIDItem::W_v)
      .value("W_w", oofem::DofIDItem::W_w)
      .value("Gamma", oofem::DofIDItem::Gamma)
      .value("D_u_edge_const", oofem::DofIDItem::D_u_edge_const)
      .value("D_u_edge_lin", oofem::DofIDItem::D_u_edge_lin)
      .value("D_v_edge_const", oofem::DofIDItem::D_v_edge_const)
      .value("D_v_edge_lin", oofem::DofIDItem::D_v_edge_lin)
      .value("Warp_PsiTheta", oofem::DofIDItem::Warp_PsiTheta)
      .value("Warp_Theta", oofem::DofIDItem::Warp_Theta)
      .value("LMP_u", oofem::DofIDItem::LMP_u)
      .value("LMP_v", oofem::DofIDItem::LMP_v)
      .value("LMP_w", oofem::DofIDItem::LMP_w)
      .value("Trac_u", oofem::DofIDItem::Trac_u)
      .value("Trac_v", oofem::DofIDItem::Trac_v)
      .value("Trac_w", oofem::DofIDItem::Trac_w)
      ;

    py::enum_<oofem::problemMode>(m, "problemMode")
        .value("processor", oofem::problemMode::_processor)
        .value("postProcessor", oofem::problemMode::_postProcessor)
    ;

    py::enum_<oofem::domainType>(m, "domainType")
        .value("_unknownMode", oofem::domainType::_unknownMode)
        .value("_2dPlaneStressMode", oofem::domainType::_2dPlaneStressMode)
        .value("_PlaneStrainMode", oofem::domainType::_PlaneStrainMode)
        .value("_2dPlaneStressRotMode", oofem::domainType::_2dPlaneStressRotMode)
        .value("_3dMode", oofem::domainType::_3dMode)
        .value("_3dAxisymmMode", oofem::domainType::_3dAxisymmMode)
        .value("_2dMindlinPlateMode", oofem::domainType::_2dMindlinPlateMode)
        .value("_3dDegeneratedShellMode", oofem::domainType::_3dDegeneratedShellMode)
        .value("_3dShellMode", oofem::domainType::_3dShellMode)
        .value("_2dTrussMode", oofem::domainType::_2dTrussMode)
        .value("_1dTrussMode", oofem::domainType::_1dTrussMode)
        .value("_2dBeamMode", oofem::domainType::_2dBeamMode)
        .value("_HeatTransferMode", oofem::domainType::_HeatTransferMode)
        .value("_Mass1TransferMode", oofem::domainType::_Mass1TransferMode)
        .value("_HeatMass1Mode", oofem::domainType::_HeatMass1Mode)
        .value("_2dIncompressibleFlow", oofem::domainType::_2dIncompressibleFlow)
        .value("_3dIncompressibleFlow", oofem::domainType::_3dIncompressibleFlow)
        .value("_2dLatticeMode", oofem::domainType::_2dLatticeMode)
        .value("_2dLatticeMassTransportMode", oofem::domainType::_2dLatticeMassTransportMode)
        .value("_3dLatticeMode", oofem::domainType::_3dLatticeMode)
        .value("_3dLatticeMassTransportMode", oofem::domainType::_3dLatticeMassTransportMode)
        .value("_2dLatticeHeatTransferMode", oofem::domainType::_2dLatticeHeatTransferMode)
        .value("_3dLatticeHeatTransferMode", oofem::domainType::_3dLatticeHeatTransferMode)
        .value("_3dDirShellMode", oofem::domainType::_3dDirShellMode)
        .value("_WarpingMode", oofem::domainType::_WarpingMode)
    ;



    py::enum_<oofem::CharType>(m, "CharType")
      .value("UnknownCharType", oofem::CharType::UnknownCharType)
      .value("StiffnessMatrix", oofem::CharType::StiffnessMatrix)
      .value("TangentStiffnessMatrix", oofem::CharType::TangentStiffnessMatrix)
      .value("SecantStiffnessMatrix", oofem::CharType::SecantStiffnessMatrix)
      .value("ElasticStiffnessMatrix", oofem::CharType::ElasticStiffnessMatrix)
      .value("MassMatrix", oofem::CharType::MassMatrix)
      .value("LumpedMassMatrix", oofem::CharType::LumpedMassMatrix)
      .value("ConductivityMatrix", oofem::CharType::ConductivityMatrix)
      .value("CapacityMatrix", oofem::CharType::CapacityMatrix)
      .value("InitialStressMatrix", oofem::CharType::InitialStressMatrix)
      .value("ExternalForcesVector", oofem::CharType::ExternalForcesVector)
      .value("InternalForcesVector", oofem::CharType::InternalForcesVector)
      .value("LastEquilibratedInternalForcesVector", oofem::CharType::LastEquilibratedInternalForcesVector)
      .value("InertiaForcesVector", oofem::CharType::InertiaForcesVector)
      .value("AuxVelocityLhs", oofem::CharType::AuxVelocityLhs)
      .value("VelocityLhs", oofem::CharType::VelocityLhs)
      .value("PressureGradientMatrix", oofem::CharType::PressureGradientMatrix)
      .value("DivergenceMatrix", oofem::CharType::DivergenceMatrix)
      .value("VelocityLaplacianMatrix", oofem::CharType::VelocityLaplacianMatrix)
      .value("PressureLaplacianMatrix", oofem::CharType::PressureLaplacianMatrix)
      .value("StabilizationMassMatrix", oofem::CharType::StabilizationMassMatrix)
      .value("PressureGradientVector", oofem::CharType::PressureGradientVector)
      .value("MassVelocityVector", oofem::CharType::MassVelocityVector)
      .value("MassAuxVelocityVector", oofem::CharType::MassAuxVelocityVector)
      .value("LaplacePressureVector", oofem::CharType::LaplacePressureVector)
      .value("LaplaceVelocityVector", oofem::CharType::LaplaceVelocityVector)
      .value("DivergenceAuxVelocityVector", oofem::CharType::DivergenceAuxVelocityVector)
      .value("DivergenceVelocityVector", oofem::CharType::DivergenceVelocityVector)
      ;

    py::enum_<oofem::Element_Geometry_Type>(m, "Element_Geometry_Type")
      .value("EGT_point", oofem::Element_Geometry_Type::EGT_point)
      .value("EGT_line_1", oofem::Element_Geometry_Type::EGT_line_1)
      .value("EGT_line_2", oofem::Element_Geometry_Type::EGT_line_2)
      .value("EGT_triangle_1", oofem::Element_Geometry_Type::EGT_triangle_1)
      .value("EGT_triangle_2", oofem::Element_Geometry_Type::EGT_triangle_2)
      .value("EGT_quad_1", oofem::Element_Geometry_Type::EGT_quad_1)
      .value("EGT_quad_1_interface", oofem::Element_Geometry_Type::EGT_quad_1_interface)
      .value("EGT_quad_21_interface", oofem::Element_Geometry_Type::EGT_quad_21_interface)
      .value("EGT_quad_2", oofem::Element_Geometry_Type::EGT_quad_2)
      .value("EGT_quad9_2", oofem::Element_Geometry_Type::EGT_quad9_2)
      .value("EGT_tetra_1", oofem::Element_Geometry_Type::EGT_tetra_1)
      .value("EGT_tetra_2", oofem::Element_Geometry_Type::EGT_tetra_2)
      .value("EGT_hexa_1", oofem::Element_Geometry_Type::EGT_hexa_1)
      .value("EGT_hexa_2", oofem::Element_Geometry_Type::EGT_hexa_2)
      .value("EGT_hexa_27", oofem::Element_Geometry_Type::EGT_hexa_27)
      .value("EGT_wedge_1", oofem::Element_Geometry_Type::EGT_wedge_1)
      .value("EGT_wedge_2", oofem::Element_Geometry_Type::EGT_wedge_2)
      .value("EGT_Composite", oofem::Element_Geometry_Type::EGT_Composite)
      .value("EGT_unknown", oofem::Element_Geometry_Type::EGT_unknown)
      ;

    py::enum_<oofem::InternalStateType>(m, "InternalStateType")
      .value("IST_Undefined", oofem::InternalStateType::IST_Undefined)
      .value("IST_StressTensor", oofem::InternalStateType::IST_StressTensor)
      .value("IST_PrincipalStressTensor", oofem::InternalStateType::IST_PrincipalStressTensor)
      .value("IST_PrincipalStressTempTensor", oofem::InternalStateType::IST_PrincipalStressTempTensor)
      .value("IST_StrainTensor", oofem::InternalStateType::IST_StrainTensor)
      .value("IST_PrincipalStrainTensor", oofem::InternalStateType::IST_PrincipalStrainTensor)
      .value("IST_PrincipalStrainTempTensor", oofem::InternalStateType::IST_PrincipalStrainTempTensor)
      .value("IST_BeamForceMomentTensor", oofem::InternalStateType::IST_BeamForceMomentTensor)
      .value("IST_BeamStrainCurvatureTensor", oofem::InternalStateType::IST_BeamStrainCurvatureTensor)
      .value("IST_ShellMomentTensor", oofem::InternalStateType::IST_ShellMomentTensor)
      .value("IST_ShellForceTensor", oofem::InternalStateType::IST_ShellForceTensor)
      .value("IST_CurvatureTensor", oofem::InternalStateType::IST_CurvatureTensor)
      .value("IST_DisplacementVector", oofem::InternalStateType::IST_DisplacementVector)
      .value("IST_DamageTensor", oofem::InternalStateType::IST_DamageTensor)
      .value("IST_DamageInvTensor", oofem::InternalStateType::IST_DamageInvTensor)
      .value("IST_PrincipalDamageTensor", oofem::InternalStateType::IST_PrincipalDamageTensor)
      .value("IST_PrincipalDamageTempTensor", oofem::InternalStateType::IST_PrincipalDamageTempTensor)
      .value("IST_CrackState", oofem::InternalStateType::IST_CrackState)
      .value("IST_StressTensorTemp", oofem::InternalStateType::IST_StressTensorTemp)
      .value("IST_StrainTensorTemp", oofem::InternalStateType::IST_StrainTensorTemp)
      .value("IST_ShellForceTensorTemp", oofem::InternalStateType::IST_ShellForceTensorTemp)
      .value("IST_ShellMomentTensorTemp", oofem::InternalStateType::IST_ShellMomentTensorTemp)
      .value("IST_CurvatureTensorTemp", oofem::InternalStateType::IST_CurvatureTensorTemp)
      .value("IST_DisplacementVectorTemp", oofem::InternalStateType::IST_DisplacementVectorTemp)
      .value("IST_DamageTensorTemp", oofem::InternalStateType::IST_DamageTensorTemp)
      .value("IST_DamageInvTensorTemp", oofem::InternalStateType::IST_DamageInvTensorTemp)
      .value("IST_CrackStateTemp", oofem::InternalStateType::IST_CrackStateTemp)
      .value("IST_PlasticStrainTensor", oofem::InternalStateType::IST_PlasticStrainTensor)
      .value("IST_PrincipalPlasticStrainTensor", oofem::InternalStateType::IST_PrincipalPlasticStrainTensor)
      .value("IST_CylindricalStressTensor", oofem::InternalStateType::IST_CylindricalStressTensor)
      .value("IST_CylindricalStrainTensor", oofem::InternalStateType::IST_CylindricalStrainTensor)
      .value("IST_MaxEquivalentStrainLevel", oofem::InternalStateType::IST_MaxEquivalentStrainLevel)
      .value("IST_ErrorIndicatorLevel", oofem::InternalStateType::IST_ErrorIndicatorLevel)
      .value("IST_InternalStressError", oofem::InternalStateType::IST_InternalStressError)
      .value("IST_PrimaryUnknownError", oofem::InternalStateType::IST_PrimaryUnknownError)
      .value("IST_RelMeshDensity", oofem::InternalStateType::IST_RelMeshDensity)
      .value("IST_MicroplaneDamageValues", oofem::InternalStateType::IST_MicroplaneDamageValues)
      .value("IST_Temperature", oofem::InternalStateType::IST_Temperature)
      .value("IST_MassConcentration_1", oofem::InternalStateType::IST_MassConcentration_1)
      .value("IST_HydrationDegree", oofem::InternalStateType::IST_HydrationDegree)
      .value("IST_Humidity", oofem::InternalStateType::IST_Humidity)
      .value("IST_Velocity", oofem::InternalStateType::IST_Velocity)
      .value("IST_Pressure", oofem::InternalStateType::IST_Pressure)
      .value("IST_VOFFraction", oofem::InternalStateType::IST_VOFFraction)
      .value("IST_Density", oofem::InternalStateType::IST_Density)
      .value("IST_MaterialInterfaceVal", oofem::InternalStateType::IST_MaterialInterfaceVal)
      .value("IST_MaterialNumber", oofem::InternalStateType::IST_MaterialNumber)
      .value("IST_ElementNumber", oofem::InternalStateType::IST_ElementNumber)
      .value("IST_BoneVolumeFraction", oofem::InternalStateType::IST_BoneVolumeFraction)
      .value("IST_PlasStrainEnerDens", oofem::InternalStateType::IST_PlasStrainEnerDens)
      .value("IST_ElasStrainEnerDens", oofem::InternalStateType::IST_ElasStrainEnerDens)
      .value("IST_TotalStrainEnerDens", oofem::InternalStateType::IST_TotalStrainEnerDens)
      .value("IST_DamageScalar", oofem::InternalStateType::IST_DamageScalar)
      .value("IST_MaterialOrientation_x", oofem::InternalStateType::IST_MaterialOrientation_x)
      .value("IST_MaterialOrientation_y", oofem::InternalStateType::IST_MaterialOrientation_y)
      .value("IST_MaterialOrientation_z", oofem::InternalStateType::IST_MaterialOrientation_z)
      .value("IST_TemperatureFlow", oofem::InternalStateType::IST_TemperatureFlow)
      .value("IST_MassConcentrationFlow_1", oofem::InternalStateType::IST_MassConcentrationFlow_1)
      .value("IST_HumidityFlow", oofem::InternalStateType::IST_HumidityFlow)
      .value("IST_CrackStatuses", oofem::InternalStateType::IST_CrackStatuses)
      .value("IST_CrackedFlag", oofem::InternalStateType::IST_CrackedFlag)
      .value("IST_CrackDirs", oofem::InternalStateType::IST_CrackDirs)
      .value("IST_CumPlasticStrain", oofem::InternalStateType::IST_CumPlasticStrain)
      .value("IST_CumPlasticStrain_2", oofem::InternalStateType::IST_CumPlasticStrain_2)
      .value("IST_StressWorkDensity", oofem::InternalStateType::IST_StressWorkDensity)
      .value("IST_DissWorkDensity", oofem::InternalStateType::IST_DissWorkDensity)
      .value("IST_FreeEnergyDensity", oofem::InternalStateType::IST_FreeEnergyDensity)
      .value("IST_ThermalConductivityIsotropic", oofem::InternalStateType::IST_ThermalConductivityIsotropic)
      .value("IST_HeatCapacity", oofem::InternalStateType::IST_HeatCapacity)
      .value("IST_AverageTemperature", oofem::InternalStateType::IST_AverageTemperature)
      .value("IST_YoungModulusVirginPaste", oofem::InternalStateType::IST_YoungModulusVirginPaste)
      .value("IST_PoissonRatioVirginPaste", oofem::InternalStateType::IST_PoissonRatioVirginPaste)
      .value("IST_YoungModulusConcrete", oofem::InternalStateType::IST_YoungModulusConcrete)
      .value("IST_PoissonRatioConcrete", oofem::InternalStateType::IST_PoissonRatioConcrete)
      .value("IST_VolumetricPlasticStrain", oofem::InternalStateType::IST_VolumetricPlasticStrain)
      .value("IST_DeviatoricStrain", oofem::InternalStateType::IST_DeviatoricStrain)
      .value("IST_DeviatoricStress", oofem::InternalStateType::IST_DeviatoricStress)
      .value("IST_Viscosity", oofem::InternalStateType::IST_Viscosity)
      .value("IST_CharacteristicLength", oofem::InternalStateType::IST_CharacteristicLength)
      .value("IST_DeviatoricStrainMeasure", oofem::InternalStateType::IST_DeviatoricStrainMeasure)
      .value("IST_DeviatoricStressMeasure", oofem::InternalStateType::IST_DeviatoricStressMeasure)
      .value("IST_vonMisesStress", oofem::InternalStateType::IST_vonMisesStress)
      .value("IST_CrackVector", oofem::InternalStateType::IST_CrackVector)
      .value("IST_PressureGradient", oofem::InternalStateType::IST_PressureGradient)
      .value("IST_DissWork", oofem::InternalStateType::IST_DissWork)
      .value("IST_DeltaDissWork", oofem::InternalStateType::IST_DeltaDissWork)
      .value("IST_StressCapPos", oofem::InternalStateType::IST_StressCapPos)
      .value("IST_TangentNorm", oofem::InternalStateType::IST_TangentNorm)
      .value("IST_Tangent", oofem::InternalStateType::IST_Tangent)
      .value("IST_DirectorField", oofem::InternalStateType::IST_DirectorField)
      .value("IST_CrackWidth", oofem::InternalStateType::IST_CrackWidth)
      .value("IST_DeformationGradientTensor", oofem::InternalStateType::IST_DeformationGradientTensor)
      .value("IST_FirstPKStressTensor", oofem::InternalStateType::IST_FirstPKStressTensor)
      .value("IST_XFEMEnrichment", oofem::InternalStateType::IST_XFEMEnrichment)
      .value("IST_XFEMNumIntersecPoints", oofem::InternalStateType::IST_XFEMNumIntersecPoints)
      .value("IST_XFEMLevelSetPhi", oofem::InternalStateType::IST_XFEMLevelSetPhi)
      .value("IST_Maturity", oofem::InternalStateType::IST_Maturity)
      .value("IST_CauchyStressTensor", oofem::InternalStateType::IST_CauchyStressTensor)
      .value("IST_InterfaceJump", oofem::InternalStateType::IST_InterfaceJump)
      .value("IST_InterfaceTraction", oofem::InternalStateType::IST_InterfaceTraction)
      .value("IST_InterfaceFirstPKTraction", oofem::InternalStateType::IST_InterfaceFirstPKTraction)
      .value("IST_StressTensor_Reduced", oofem::InternalStateType::IST_StressTensor_Reduced)
      .value("IST_StrainTensor_Reduced", oofem::InternalStateType::IST_StrainTensor_Reduced)
      .value("IST_CrossSectionNumber", oofem::InternalStateType::IST_CrossSectionNumber)
      .value("IST_ShellStrainTensor", oofem::InternalStateType::IST_ShellStrainTensor)
      .value("IST_AbaqusStateVector", oofem::InternalStateType::IST_AbaqusStateVector)
      .value("IST_AutogenousShrinkageTensor", oofem::InternalStateType::IST_AutogenousShrinkageTensor)
      .value("IST_DryingShrinkageTensor", oofem::InternalStateType::IST_DryingShrinkageTensor)
      .value("IST_TotalShrinkageTensor", oofem::InternalStateType::IST_TotalShrinkageTensor)
      .value("IST_ThermalStrainTensor", oofem::InternalStateType::IST_ThermalStrainTensor)
      .value("IST_CreepStrainTensor", oofem::InternalStateType::IST_CreepStrainTensor)
      .value("IST_TensileStrength", oofem::InternalStateType::IST_TensileStrength)
      .value("IST_ResidualTensileStrength", oofem::InternalStateType::IST_ResidualTensileStrength)
      .value("IST_LocalEquivalentStrain", oofem::InternalStateType::IST_LocalEquivalentStrain)
      .value("IST_CrackIndex", oofem::InternalStateType::IST_CrackIndex)
      .value("IST_EigenStrainTensor", oofem::InternalStateType::IST_EigenStrainTensor)
      .value("IST_CrackStrainTensor", oofem::InternalStateType::IST_CrackStrainTensor)
      .value("IST_2ndCrackWidth", oofem::InternalStateType::IST_2ndCrackWidth)
      .value("IST_2ndCrackVector", oofem::InternalStateType::IST_2ndCrackVector)
      .value("IST_3rdCrackWidth", oofem::InternalStateType::IST_3rdCrackWidth)
      .value("IST_3rdCrackVector", oofem::InternalStateType::IST_3rdCrackVector)
      .value("IST_FiberStressLocal", oofem::InternalStateType::IST_FiberStressLocal)
      .value("IST_FiberStressNL", oofem::InternalStateType::IST_FiberStressNL)
      .value("IST_EnergyMassCapacity", oofem::InternalStateType::IST_EnergyMassCapacity)
      .value("IST_PrincStressVector1", oofem::InternalStateType::IST_PrincStressVector1)
      .value("IST_PrincStressVector2", oofem::InternalStateType::IST_PrincStressVector2)
      .value("IST_PrincStressVector3", oofem::InternalStateType::IST_PrincStressVector3)
      .value("IST_InterfaceNormal", oofem::InternalStateType::IST_InterfaceNormal)
      .value("IST_MomentTensor", oofem::InternalStateType::IST_MomentTensor)
      .value("IST_MomentTensorTemp", oofem::InternalStateType::IST_MomentTensorTemp)
      .value("IST_YieldStrength", oofem::InternalStateType::IST_YieldStrength)
      .value("IST_ElasticStrainTensor", oofem::InternalStateType::IST_ElasticStrainTensor)
      .value("IST_MoistureContent", oofem::InternalStateType::IST_MoistureContent)
      .value("IST_CrackStatusesTemp", oofem::InternalStateType::IST_CrackStatusesTemp)
      .value("IST_CrackSlip", oofem::InternalStateType::IST_CrackSlip)
      .value("IST_EquivalentTime", oofem::InternalStateType::IST_EquivalentTime)
      .value("IST_IncrementCreepModulus", oofem::InternalStateType::IST_IncrementCreepModulus)
      ;

    py::register_exception<oofem::InputException>(m, "InputException");
    py::register_exception<oofem::MissingKeywordInputException>(m, "MissingKeywordInputException");
    py::register_exception<oofem::BadFormatInputException>(m, "BadFormatInputException");
    py::register_exception<oofem::ValueInputException>(m, "ValueInputException");

   
    
    py::enum_<oofem::MatResponseMode>(m, "MatResponseMode")
        .value("TangentStiffness", oofem::MatResponseMode::TangentStiffness)
        .value("SecantStiffness", oofem::MatResponseMode::SecantStiffness)
        .value("ElasticStiffness", oofem::MatResponseMode::ElasticStiffness)
        .value("Conductivity", oofem::MatResponseMode::Conductivity)
        .value("Conductivity_ww", oofem::MatResponseMode::Conductivity_ww)
        .value("Conductivity_hh", oofem::MatResponseMode::Conductivity_hh)
        .value("Conductivity_hw", oofem::MatResponseMode::Conductivity_hw)
        .value("Conductivity_wh", oofem::MatResponseMode::Conductivity_wh)
        .value("Capacity", oofem::MatResponseMode::Capacity)
        .value("Capacity_ww", oofem::MatResponseMode::Capacity_ww)
        .value("Capacity_hh", oofem::MatResponseMode::Capacity_hh)
        .value("Capacity_hw", oofem::MatResponseMode::Capacity_hw)
        .value("Capacity_wh", oofem::MatResponseMode::Capacity_wh)
        .value("IntSource", oofem::MatResponseMode::IntSource)
        .value("IntSource_ww", oofem::MatResponseMode::IntSource_ww)
        .value("IntSource_hh", oofem::MatResponseMode::IntSource_hh)
        .value("IntSource_hw", oofem::MatResponseMode::IntSource_hw)
        .value("IntSource_wh", oofem::MatResponseMode::IntSource_wh)
    ;


    m.def("linearStatic", &linearStatic, py::return_value_policy::move);
    m.def("staticStructural", &staticStructural, py::return_value_policy::move);
    m.def("domain", &domain, py::return_value_policy::move);
    m.def("truss1d", &truss1d, py::return_value_policy::move);
    m.def("beam2d", &beam2d, py::return_value_policy::move);
    m.def("trPlaneStress2d", &trPlaneStress2d, py::return_value_policy::move);
    m.def("planeStress2d", &planeStress2d, py::return_value_policy::move);
    m.def("transientTransport", &transientTransport, py::return_value_policy::move);
    m.def("qBrick1ht", &qBrick1ht, py::return_value_policy::move);

    m.def("node", &node, py::return_value_policy::move);
    m.def("boundaryCondition", &boundaryCondition, py::return_value_policy::move);
    m.def("initialCondition", &initialCondition, py::return_value_policy::move);
    m.def("constantEdgeLoad", &constantEdgeLoad, py::return_value_policy::move);
    m.def("constantSurfaceLoad", &constantSurfaceLoad, py::return_value_policy::move);
    m.def("nodalLoad", &nodalLoad, py::return_value_policy::move);
    m.def("structTemperatureLoad", &structTemperatureLoad, py::return_value_policy::move);
    m.def("structEigenstrainLoad", &structEigenstrainLoad, py::return_value_policy::move);
    m.def("isoLE", &isoLE, py::return_value_policy::move);
    m.def("idm1", &idm1, py::return_value_policy::move);
    m.def("isoHeat", &isoHeat, py::return_value_policy::move);

    m.def("simpleCS", &simpleCS, py::return_value_policy::move);
    m.def("simpleTransportCS", &simpleTransportCS, py::return_value_policy::move);
    m.def("peakFunction", &peakFunction, py::return_value_policy::move);
    m.def("constantFunction", &constantFunction, py::return_value_policy::move);
    m.def("piecewiseLinFunction", &piecewiseLinFunction, py::return_value_policy::move);
    m.def("vtkxml", &vtkxml, py::return_value_policy::move);
    m.def("vtkmemory", &vtkmemory, py::return_value_policy::move);
    m.def("homExport", &homExport, py::return_value_policy::move);
    m.def("createSet", &createSet, py::return_value_policy::move);


//std::shared_ptr<oofem::Field>
    py::class_<oofem::Field, PyField, std::shared_ptr<oofem::Field>>(m, "Field")
        .def(py::init<oofem::FieldType>())  
//         .def("evaluateAt", (int (oofem::Field::*)(oofem::FloatArray &answer, const oofem::FloatArray &coords, oofem::ValueModeType mode, oofem::TimeStep *tStep)) &oofem::Field::evaluateAt)
        
        .def("evaluateAt", [](oofem::FloatArray &answer, const oofem::FloatArray &coords, oofem::ValueModeType mode, oofem::TimeStep *tStep){
            return std::make_tuple(answer,coords);//TODO-how to invoke the function if it does not exist yet?
        })
        .def("giveType", &oofem::Field::giveType)
        .def("setType", &oofem::Field::setType)
        ;

    py::class_<oofem::UniformGridField, oofem::Field, std::shared_ptr<oofem::UniformGridField>>(m, "UniformGridField")
        .def(py::init<>())
        .def("setGeometry", &oofem::UniformGridField::setGeometry)
        .def("setValues", &oofem::UniformGridField::setValues)
        ;

    py::class_<oofem::UnstructuredGridField, oofem::Field, std::shared_ptr<oofem::UnstructuredGridField>>(m, "UnstructuredGridField")
        .def(py::init<int, int, double>(), py::arg().noconvert(), py::arg().noconvert(), py::arg("octreeOriginShift") = 0.0)
        .def("addVertex", &oofem::UnstructuredGridField::addVertex)
        .def("setVertexValue", &oofem::UnstructuredGridField::setVertexValue)
        ;
    
    py::class_<oofem::DofManValueField, oofem::Field,  std::shared_ptr<oofem::DofManValueField>>(m, "DofManValueField")
        .def(py::init<oofem::FieldType,oofem::Domain*>())
        .def(py::init<oofem::FieldType,int,int, const std::string, const std::string>(), py::arg().noconvert(), py::arg().noconvert(), py::arg().noconvert(), py::arg("engngModel")="transienttransport", py::arg("domainDofsDefaults")="heattransfer")
        .def("addNode", &oofem::DofManValueField::addNode )
        .def("addElement", &oofem::DofManValueField::addElement )
        .def("setDofManValue", &oofem::DofManValueField::setDofManValue )
        .def("getNodeCoordinates", &oofem::DofManValueField::getNodeCoordinates )
        ;
    
        
//depends on Python.h
#ifdef _PYBIND_BINDINGS       
    py::class_<oofem::PythonField, oofem::Field, std::shared_ptr<oofem::PythonField>>(m, "PythonField")
        .def(py::init<>())
        .def("setModuleName", &oofem::PythonField::setModuleName)
        .def("setFunctionName", &oofem::PythonField::setFunctionName)
        ;   
#endif


    m.def("test", &test);
 }
