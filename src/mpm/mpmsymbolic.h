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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#ifndef mpmsymbolic_h
#define mpmsymbolic_h

/**
 * Multiphysics module 
 * Classes:
 * - ElementBase(Element) defining geometry
 * - Variable class representing unknown field (or test feld) in a weak psolution. The variable has its interpolation, type (scalar, vector), size.
    When test field, it keeps reference to its primary (unknown) variable. The history parameter dermines how many time steps to remember. 
 * - Term class represnting a term to evaluate on element. Paramaters element(geometry), variables
 * - Element - responsible for defining and performing integration (of terms), assembly of term contributions. 
 */

#include "classfactory.h"
#include "mpm.h"
#include "mpmevaluator2.h"
#include "logger.h"
#include "feinterpol.h"
#include "CrossSections/structuralcrosssection.h"
#include "matresponsemode.h"
#include "engngm.h"



namespace oofem {

    /* Define custom functors for evaluator */
    auto MPMfunctor_Grad_s = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the symmetric gradient of the first argument (assumed to be a vector field) 
        // ARGS: args[0] - pointer to VarSlot containing the vector field (as a user pointer)
        //       args[1] - pointer to GaussPoint (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting symmetric gradient matrix
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Grad_s functor with %ld arguments\n", args.size());
        if (args.size() != 2) {
            OOFEM_ERROR("MPMfunctor_Grad_s functor expects exactly 2 arguments: vector field (Variable class) and GaussPoint.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        // functor logic
        const MPElement* cell = static_cast<const MPElement*>(gp->giveElement());
        FloatMatrix answer;
        cell->computeGradSymMatrixAt(answer, v, gp);

        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };

    // define gradient of the unknown field variable (displacement) functor
    auto MPMfunctor_Grad = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the gradient of the first argument (assumed to be a scalar field) 
        // ARGS: args[0] - pointer to VarSlot containing the scalar field (as a user pointer)
        //       args[1] - pointer to GaussPoint (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting symmetric gradient matrix
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Grad functor with %ld arguments\n", args.size());
        if (args.size() != 2) {
            OOFEM_ERROR("MPMfunctor_Grad functor expects exactly 2 arguments: scalar field (Variable class) and GaussPoint.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        const MPElement* cell = static_cast<const MPElement*>(gp->giveElement());

        // functor logic
        FloatMatrix answer;
        cell->computeGradMatrixAt(answer, v, gp);

        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };

    // define divergence op functor on vector field variable (e.g. velocity) 
    auto MPMfunctor_Div = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the divergence of the first argument (assumed to be a vector field) 
        // ARGS: args[0] - pointer to VarSlot containing the vector field (as a user pointer)
        //       args[1] - pointer to GaussPoint (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting divergence (scalar)
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Div functor with %ld arguments\n", args.size());
        if (args.size() != 2) {
            OOFEM_ERROR("MPMfunctor_Div functor expects exactly 2 arguments: vector field (Variable class) and GaussPoint.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        const MPElement* cell = static_cast<const MPElement*>(gp->giveElement());
        const MaterialMode mmode = gp->giveMaterialMode();
        // functor logic
        FloatMatrix answer;
        const FEInterpolation* interpol = v->interpolation;

        FloatMatrix dndx;
        interpol->evaldNdx(dndx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(cell));
        int nnodes = interpol->giveNumberOfNodes(cell->giveGeometryType());
        answer.resize(1, nnodes*v->size);
        int nsd = 1*mmodeIs1D(mmode) + 2*mmodeIs2D(mmode) + 3*mmodeIs3D(mmode);
        for (int i = 0; i< nnodes; i++) {
            for (int j = 0; j< nsd; j++) {
                answer(0, i*v->size+j) = dndx(i, j);
            }
        }
        out.value = answer;
        out.type = VarSlot::Type::MATRIX;   
    };

    // define Interpolation op (interpolation matrix) of the unknown field variable functor
    auto MPMfunctor_N= [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the gradient of the first argument (assumed to be a scalar field) 
        // ARGS: args[0] - pointer to VarSlot containing the scalar field (as a user pointer)
        //       args[1] - pointer to GaussPoint (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting interpolation matrix
        OOFEM_LOG_DEBUG("    [C++ Callback] Called N functor with %ld arguments\n", args.size());
        if (args.size() != 2) {
            OOFEM_ERROR("MPMfunctor_N functor expects exactly 2 arguments: scalar field (Variable class) and GaussPoint.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        const MPElement* cell = static_cast<const MPElement*>(gp->giveElement());
        //const MaterialMode mmode = gp->giveMaterialMode();
        

        // functor logic
        FloatMatrix N;
        cell->computeNMatrixAt(N, v, gp);

        out.value = N;
        out.type = VarSlot::Type::MATRIX;
    };

auto MPMfunctor_MProp = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the constitutive property  
        // ARGS: args[0] - pointer to GaussPoint (as a user pointer)
        //       args[1] - pointer to TimeStep (as a user pointer)
        //       args[2] - property ID (as a double, to be casted to MaterialResponseMode enum)
        // OUTPUT: out - VarSlot to store the resulting property (e.g. stress, etc depending on property ID)
        OOFEM_LOG_DEBUG("    [C++ Callback] Called MProp functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_MProp functor expects exactly 3 arguments: GaussPoint, TimeStep and PropertyID.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        double raw_val2 = std::get<double>(args[2]->value);
        // 2. Cast back to your specific application type (Variable class)
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr0);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr1);
        MatResponseMode propertyID = static_cast<MatResponseMode>(raw_val2);

        // functor logic
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());
        StructuralCrossSection* cs = static_cast<StructuralCrossSection*>(cell->giveCrossSection());

        FloatArray charVec;

        double property = cs->giveMaterial(gp)->giveCharacteristicValue(propertyID, gp, tstep);
        
        out.value = property;
        out.type = VarSlot::Type::SCALAR;
        std::ostringstream oss;
        //oss << "MProp result: " << std::get<double>(out.value)  << "\n\n";
        OOFEM_LOG_DEBUG("%s", oss.str().c_str());
    };

    auto MPMfunctor_MVec = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the constitutive variable/property  
        // ARGS: args[0] - pointer to GaussPoint (as a user pointer)
        //       args[1] - pointer to TimeStep (as a user pointer)
        //       args[2] - property ID (as a double, to be casted to MaterialResponseMode enum)
        // OUTPUT: out - VarSlot to store the resulting characteristic vector (e.g. stress, etc depending on property ID)
        // The generalized state is no longer passed in: it is pushed to the material once per
        // iteration before any term is evaluated, so this is a cache read.
        OOFEM_LOG_DEBUG("    [C++ Callback] Called MVec functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_MVec functor expects exactly 3 arguments: GaussPoint, TimeStep and PropertyID.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        double raw_val2 = std::get<double>(args[2]->value);
        // 2. Cast back to your specific application type (Variable class)
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr0);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr1);
        MatResponseMode propertyID = static_cast<MatResponseMode>(raw_val2);

        // functor logic
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());
        StructuralCrossSection* cs = static_cast<StructuralCrossSection*>(cell->giveCrossSection());

        FloatArray charVec;

        cs->giveMaterial(gp)->giveCharacteristicVector(charVec, propertyID, gp, tstep);
        
        out.value = FloatMatrix::fromArray(charVec);
        out.type = VarSlot::Type::MATRIX;
        std::ostringstream oss;
        oss << "MVec result: " << std::get<FloatMatrix>(out.value)  << "\n\n";
        OOFEM_LOG_DEBUG("%s", oss.str().c_str());
    };

    auto MPMfunctor_MDer = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the constitutive derivative of the given property 
        // ARGS: args[0] - pointer to GaussPoint (as a user pointer)
        //       args[1] - pointer to TimeStep (as a user pointer)
        //       args[2] - property ID (as a double, to be casted to MaterialResponseMode enum)
        // OUTPUT: out - VarSlot to store the resulting symmetric gradient matrix
        OOFEM_LOG_DEBUG("    [C++ Callback] Called MDer functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_MDer functor expects exactly 3 arguments: GaussPoint, TimeStep and PropertyID.");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        double raw_val2 = std::get<double>(args[2]->value);
        // 2. Cast back to your specific application type (Variable class)
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr0);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr1);
        MatResponseMode propertyID = static_cast<MatResponseMode>(raw_val2);

        // functor logic
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());
        StructuralCrossSection* cs = static_cast<StructuralCrossSection*>(cell->giveCrossSection());

        FloatMatrix D;
        if (propertyID == MatResponseMode::DeviatoricStiffness) {
            cs->giveMaterial(gp)->giveCharacteristicMatrix(D, propertyID, gp, tstep);
        } else {
            cs->giveMaterial(gp)->giveCharacteristicMatrix(D, propertyID, gp, tstep);
        }
        out.value = D;
        out.type = VarSlot::Type::MATRIX;
    };

    auto MPMfunctor_Sig = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the Stress vector of the first argument (assumed to be a vector field) 
        // ARGS: args[0] - pointer to GaussPoint (as a user pointer)
        //       args[1] - pointer to TimeStep (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting stress vector
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Sig functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_Sig functor expects exactly 3 arguments: Field Variable, GaussPoint (gp), TimeStep (ts).");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        void* raw_ptr2 = std::get<void*>(args[2]->value);

        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr2);
        // functor logic
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());
        StructuralCrossSection* cs = static_cast<StructuralCrossSection*>(cell->giveCrossSection());

        // The strain is no longer derived here: it was pushed to the material as part of the
        // generalized state once per iteration, so this is a cache read. The field argument is
        // kept for backward compatibility of the input syntax and is validated against the field
        // that actually supplies the strain on this cell, so that a deck naming the wrong one
        // fails instead of silently reading another field's state.
        const Variable* strainSource = cell->giveStateVariableSource(FT_Displacements);
        if (strainSource != nullptr && strainSource != v) {
            OOFEM_ERROR("Sig(%s, ...) does not match the field supplying the strain on element %d ('%s')",
                        v->name.c_str(), cell->giveNumber(), strainSource->name.c_str());
        }

        FloatMatrix answer;
        FloatArray sig;
        cs->giveMaterial(gp)->giveCharacteristicVector(sig, MatResponseMode::Stress, gp, tstep);
        answer = FloatMatrix::fromArray(sig);
        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };

    auto MPMfunctor_Sig_dev = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the deviatoric stress vector of the first argument (assumed to be a vector field) 
        // ARGS: args[0] - pointer to GaussPoint (as a user pointer)
        //       args[1] - pointer to TimeStep (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting deviatoric stress vector
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Sig_dev functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_Sig_dev functor expects exactly 3 arguments: Field Variable, GaussPoint (gp), TimeStep (ts).");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        void* raw_ptr2 = std::get<void*>(args[2]->value);

        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr2);
        // functor logic
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());
        StructuralCrossSection* cs = static_cast<StructuralCrossSection*>(cell->giveCrossSection());

        // The strain is no longer derived here: it was pushed to the material as part of the
        // generalized state once per iteration, so this is a cache read. The field argument is
        // kept for backward compatibility of the input syntax and is validated against the field
        // that actually supplies the strain on this cell, so that a deck naming the wrong one
        // fails instead of silently reading another field's state.
        const Variable* strainSource = cell->giveStateVariableSource(FT_Displacements);
        if (strainSource != nullptr && strainSource != v) {
            OOFEM_ERROR("Sig_dev(%s, ...) does not match the field supplying the strain on element %d ('%s')",
                        v->name.c_str(), cell->giveNumber(), strainSource->name.c_str());
        }

        FloatMatrix answer;
        FloatArray sig;
        cs->giveMaterial(gp)->giveCharacteristicVector(sig, MatResponseMode::DeviatoricStress, gp, tstep);
        answer = FloatMatrix::fromArray(sig);
        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };

    // define functor to return element field nodal values (e.g. temperature at nodes) as a column matrix
    auto MPMfunctor_FieldNodalValues = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the nodal values of the first argument (assumed to be a field) at intrinsic time of time step.
        // ARGS: args[0] - pointer to Variable (as a user pointer)
        //       args[1] - element (cell) (as a user pointer)
        //       args[2] - timep step (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting nodal values column matrix (dof ordering determined by field dof ordering)
        OOFEM_LOG_DEBUG("    [C++ Callback] Called FieldNodalValues functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_FieldNodalValues functor expects exactly 3 arguments: Variable, Element(cell), and TimeStep(ts)");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        void* raw_ptr2 = std::get<void*>(args[2]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        MPElement* cell = static_cast<MPElement*>(raw_ptr1);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr2);

        // functor logic
        FloatArray u;
        cell->getUnknownVector(u, v, VM_TotalIntrinsic, tstep); // get nodal values of the variable at current time step
        FloatMatrix answer = FloatMatrix::fromArray(u);

        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    }; 
    // define functor to return element field nodal velocities (e.g. temperature rates at nodes) as a column matrix
    auto MPMfunctor_FieldNodalVelocities = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the nodal values of the first argument (assumed to be a field) at intrinsic time of time step.
        // ARGS: args[0] - pointer to Variable (as a user pointer)
        //       args[1] - element (cell) (as a user pointer)
        //       args[2] - timep step (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting nodal values column matrix (dof ordering determined by field dof ordering)
        OOFEM_LOG_DEBUG("    [C++ Callback] Called FieldNodalVelocities functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_FieldNodalVelocities functor expects exactly 3 arguments: Variable, Element(cell), and TimeStep(ts)");
        }
        // 1. Retrieve the generic pointers to arguments
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        void* raw_ptr2 = std::get<void*>(args[2]->value);
        // 2. Cast back to your specific application type (Variable class)
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        MPElement* cell = static_cast<MPElement*>(raw_ptr1);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr2);

        // functor logic
        FloatArray u;
        cell->getUnknownVector(u, v, VM_Velocity, tstep); // get nodal velocities of the variable at current time step
        FloatMatrix answer = FloatMatrix::fromArray(u);

        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    }; 

    // define functor to evaluate field variable at given integration point
    auto MPMfunctor_Eval = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        // Compute the field value at Gauss point.
        // ARGS: args[0] - pointer to Variable (as a user pointer)
        //       args[1] - GaussPoint (as a user pointer)
        //       args[2] - time step (as a user pointer)
        // OUTPUT: out - VarSlot to store the resulting field value as a column matrix
        OOFEM_LOG_DEBUG("    [C++ Callback] Called EvalField functor with %ld arguments\n", args.size());
        if (args.size() != 3) {
            OOFEM_ERROR("MPMfunctor_Eval functor expects exactly 3 arguments: Variable, GaussPoint, and TimeStep.");
        }
        
        void* raw_ptr0 = std::get<void*>(args[0]->value);
        void* raw_ptr1 = std::get<void*>(args[1]->value);
        void* raw_ptr2 = std::get<void*>(args[2]->value);
        
        const Variable* v = static_cast<const Variable*>(raw_ptr0);
        GaussPoint* gp = static_cast<GaussPoint*>(raw_ptr1);
        TimeStep* tstep = static_cast<TimeStep*>(raw_ptr2);
        MPElement* cell = static_cast<MPElement*>(gp->giveElement());

        FloatArray u, nvec;
        FloatMatrix N, uMat, answer;
        
        cell->getUnknownVector(u, v, VM_TotalIntrinsic, tstep);
        v->interpolation->evalN(nvec, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(cell));
        N.beNMatrixOf(nvec, v->size);
        
        uMat = FloatMatrix::fromArray(u);
        answer.beProductOf(N, uMat);

        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };

    // define functor to concatenate vectors and scalars vertically
    auto MPMfunctor_vcat = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        OOFEM_LOG_DEBUG("    [C++ Callback] Called Concat functor with %ld arguments\n", args.size());
        if (args.empty()) {
            OOFEM_ERROR("MPMfunctor_vcat functor expects at least 1 argument.");
        }
        
        int total_rows = 0;
        int num_cols = -1;
        
        for (const auto& arg : args) {
            if (arg->type == VarSlot::Type::MATRIX) {
                const auto& m = std::get<FloatMatrix>(arg->value);
                total_rows += m.rows();
                if (num_cols == -1) num_cols = m.cols();
                else if (num_cols != m.cols()) OOFEM_ERROR("MPMfunctor_vcat: Matrix column count mismatch.");
            } else if (arg->type == VarSlot::Type::SCALAR) {
                total_rows += 1;
                if (num_cols == -1) num_cols = 1;
                else if (num_cols != 1) OOFEM_ERROR("MPMfunctor_vcat: Scalar cannot be concatenated with matrix having cols != 1.");
            } else {
                OOFEM_ERROR("MPMfunctor_vcat: Unsupported argument type.");
            }
        }
        
        if (num_cols == -1) num_cols = 1;
        FloatMatrix answer(total_rows, num_cols);
        int current_row = 0;
        
        for (const auto& arg : args) {
            if (arg->type == VarSlot::Type::MATRIX) {
                const auto& m = std::get<FloatMatrix>(arg->value);
                for (int r = 0; r < m.rows(); ++r) {
                    for (int c = 0; c < num_cols; ++c) {
                        answer(current_row + r, c) = m(r, c);
                    }
                }
                current_row += m.rows();
            } else if (arg->type == VarSlot::Type::SCALAR) {
                answer(current_row++, 0) = std::get<double>(arg->value);
            }
        }
        
        out.value = answer;
        out.type = VarSlot::Type::MATRIX;
    };


    auto MPMfunctor_LumpMatrix = [](const std::vector<const VarSlot*>& args, VarSlot& out) {
        /*
        HRZ Lumping (Hintz-Rock-Zienkiewicz)
        Instead of summing the rows, HRZ lumping preserves only the diagonal terms of the consistent matrix 
        and scales them so that their total sum equals the total mass of the element.
        */

        // 1. Validation
        if (args.size() != 2) {
            throw std::runtime_error("LumpMatrix requires 2 arguments: (ConsistentMatrix, TargetTotalMass)");
        }
        if (args[0]->type != VarSlot::Type::MATRIX || args[1]->type != VarSlot::Type::SCALAR) {
            throw std::runtime_error("LumpMatrix inputs must be (Matrix, Scalar).");
        }

        const FloatMatrix& M_c = std::get<FloatMatrix>(args[0]->value);
        double target_mass = std::get<double>(args[1]->value);

        if (M_c.rows() != M_c.cols()) {
            throw std::runtime_error("LumpMatrix requires a square consistent mass matrix.");
        }

        // 2. Calculate the sum of the diagonal
        double diag_sum = 0.0;
        for (int i = 0; i < M_c.rows(); ++i) {
            diag_sum += M_c(i, i);
        }

        // Prevent division by zero if the matrix is empty/invalid
        if (std::abs(diag_sum) < 1e-14) {
            throw std::runtime_error("LumpMatrix failed: Sum of diagonal is practically zero.");
        }

        // 3. Calculate the HRZ scaling factor
        double scale_factor = target_mass / diag_sum;

        // 4. Create the lumped matrix by scaling the diagonal
        FloatMatrix M_lumped(M_c.rows(), M_c.cols());
        M_lumped.zero();
        for (int i = 0; i < M_c.rows(); ++i) {
            M_lumped(i, i) = M_c(i, i) * scale_factor;
        }

        // 5. Output to VM
        out.value = M_lumped;
        out.type = VarSlot::Type::MATRIX;
    };


    auto MPMfunctor_print =[](const std::vector<const VarSlot*>& args, VarSlot& out) {
        
        std::cout << "[Debug Output: ";
        
        for (size_t i = 0; i < args.size(); ++i) {
            const VarSlot& slot = *args[i];
            
            if (slot.type == VarSlot::Type::SCALAR) {
                std::cout << "Scalar: " << std::get<double>(slot.value);
            } 
            else if (slot.type == VarSlot::Type::MATRIX) {
                // Add a newline for matrices so they format nicely
                std::cout << "Matrix:" << std::get<FloatMatrix>(slot.value);
            } 
            else if (slot.type == VarSlot::Type::USER_PTR) {
                std::cout << "<C++ Object at " << std::get<void*>(slot.value) << ">";
            }
            
            if (i < args.size() - 1) std::cout << ", ";
        }
        
        std::cout << "]" << std::endl;

        // --- CRITICAL: Satisfy the VM by returning a dummy scalar ---
        out.value = 0.0;
        out.type = VarSlot::Type::SCALAR;
    };
/**
 * @brief Symbolic term allowing to parse and evaluate user defined expressions
 * 
 */
class SymbolicTerm : public GenericCellTerm {
    protected:
        std::string lhsExpression, rhsExpression;
        mutable int pool_ptr=0;
        struct VMContext {
            mutable std::vector<Instruction> program;
            mutable std::map<std::string, int> symbols;
            mutable std::map<int, VarData> constants;
            /**
             * Immutable execution environment for this expression, established at the end of
             * initializeFrom and read-only thereafter; see buildEnvironment and _evaluateVM.
             */
            std::unique_ptr<MPMEnvironment> env;
        };
        mutable VMContext lhsExpressionContext, rhsExpressionContext;
        
        EngngModel *problem;
 
        struct TestField {
            FloatMatrix values; // Nodal values (e.g., Temperature at 3 nodes)
            TestField(double start_val) { 
            values=FloatMatrix::fromIniList({{start_val}, {start_val + 10.0}, {start_val + 20.0}}); 
        }
    };

    public:
    SymbolicTerm() : GenericCellTerm() {}
    SymbolicTerm (const Variable *testField, const Variable* unknownField, const std::string &lexpr, const std::string& rexpr, MaterialMode m=MaterialMode::_Unknown)  : GenericCellTerm(testField, unknownField, m), lhsExpression(lexpr), rhsExpression(rexpr) {}

    void initializeFrom(const std::shared_ptr<InputRecord> &ir, EngngModel* problem) override {
        GenericCellTerm::initializeFrom(ir, problem);
        IR_GIVE_FIELD(ir, lhsExpression, "lexpression");
        IR_GIVE_FIELD(ir, rhsExpression, "rexpression");

        MPMCompiler compiler;

        // Declare functions (functors)
        compiler.register_function("Grad_s");
        compiler.register_function("Grad");
        compiler.register_function("Div");
        compiler.register_function("N");
        compiler.register_function("Sig");
        compiler.register_function("Sig_dev");
        compiler.register_function("MDer");
        compiler.register_function("MVec"); // characteristic vector (e.g. stress) from material response
        compiler.register_function("MProp"); // material property
        compiler.register_function("vcat"); // matrix/vector vertical concatenation
        compiler.register_function("LumpMatrix"); // HRZ lumping of consistent mass matrix
        compiler.register_function("eval"); 
        compiler.register_function("print");


        compiler.register_function("ru");
        compiler.register_function("rv");

        try {
            compiler.compile_script(lhsExpression, lhsExpressionContext.program, lhsExpressionContext.symbols, lhsExpressionContext.constants, pool_ptr);
        } catch (const std::exception& e) {
            std::string msg = "SymbolicTerm: Compilation error in expression '" + lhsExpression + "': " + e.what();
            OOFEM_ERROR("%s", msg.c_str());
        }
        try {
            compiler.compile_script(rhsExpression, rhsExpressionContext.program, rhsExpressionContext.symbols, rhsExpressionContext.constants, pool_ptr);
        } catch (const std::exception& e) {
            std::string msg = "SymbolicTerm: Compilation error in expression '" + rhsExpression + "': " + e.what();
            OOFEM_ERROR("%s", msg.c_str());
        }
        this->problem = problem;

        // Both expressions are compiled and the problem is known, so the invariant part of their
        // execution environment can be established now, once, instead of at every evaluation.
        this->buildEnvironment(lhsExpressionContext);
        this->buildEnvironment(rhsExpressionContext);
    }

    /**
     * Establishes the immutable execution environment of one compiled expression.
     *
     * Everything except the point being evaluated is invariant: the symbol table, the functor
     * table, the compiled constants, the problem's variables and the response-mode literals.
     * Setting all of that up per Gauss point, per term, per sweep dominated the cost of evaluating
     * cheap expressions.
     *
     * Called at the end of initializeFrom rather than lazily on first evaluation, so that the
     * environment is fully built before any assembly starts and evaluation needs no locking.
     */
    void buildEnvironment(VMContext& context) const {
        context.env = std::make_unique<MPMEnvironment>();
        MPMEnvironment& env = *context.env;

        env.symbols = context.symbols;

        env.functors["Grad_s"] = MPMfunctor_Grad_s;
        env.functors["Grad"] = MPMfunctor_Grad;
        env.functors["Div"] = MPMfunctor_Div;
        env.functors["N"] = MPMfunctor_N;
        env.functors["Sig"] = MPMfunctor_Sig;
        env.functors["Sig_dev"] = MPMfunctor_Sig_dev;
        env.functors["MDer"] = MPMfunctor_MDer;
        env.functors["MVec"] = MPMfunctor_MVec;
        env.functors["MProp"] = MPMfunctor_MProp;
        env.functors["vcat"] = MPMfunctor_vcat;
        env.functors["eval"] = MPMfunctor_Eval;
        env.functors["LumpMatrix"] = MPMfunctor_LumpMatrix;
        env.functors["print"] = MPMfunctor_print;
        env.functors["ru"] = MPMfunctor_FieldNodalValues;
        env.functors["rv"] = MPMfunctor_FieldNodalVelocities;

        // Seed the template slot pool with the invariant bindings, by doing them once on a
        // throwaway evaluator and keeping its pool.
        MPMEvaluator seed(pool_ptr, context.symbols);
        for (auto const& [idx, val] : context.constants) {
            seed.init_slot(idx, val);
        }
        // all problem variables, as user pointers
        for (auto &i : problem->giveVariables()) {
            seed.set_variable(i.first, (void*)i.second.get());
        }
        // Every response mode by name, so that a deck can say
        // MDer(gp, ts, MatResponseMode::Permeability) rather than MDer(gp, ts, 19) and stop
        // encoding enum values. Driven off the {value, name} table that enum-impl.h already
        // generates for its ToString helper, so new modes need no work here; set_variable ignores
        // names the script does not use.
        for (auto &item : EnumData<MatResponseMode>::value_to_name) {
            seed.set_variable(std::string("MatResponseMode::") + item.name, (double)item.value);
        }

        env.pool = seed.givePool();
        env.is_set = seed.giveIsSet();
    }

    void _evaluateVM (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep, VMContext& context) const {
        try {
            // Private scratch over the shared, read-only environment: only the slot pool is
            // copied, so concurrent evaluations of the same expression do not interfere.
            MPMEvaluator vm(*context.env);

            // The only genuinely per-evaluation bindings.
            vm.set_variable("gp", (void*)gp);
            vm.set_variable("ts", (void*)tStep);
            vm.set_variable("cell", (void*)&cell);

            vm.execute(context.program);
            if (vm.get_result().type == VarSlot::Type::MATRIX) {
                answer = std::get<FloatMatrix>(vm.get_result().value);
                std::ostringstream oss;
                oss << "Result: " << answer << "\n\n";
                OOFEM_LOG_DEBUG("%s", oss.str().c_str());
            }

        } catch (const std::exception& e) {
            OOFEM_ERROR("VM ERROR: %s", e.what());
        }
    }
    void evaluate_lin (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override {
        _evaluateVM(answer, cell, gp, tStep, lhsExpressionContext);
    }
    void evaluate (FloatArray&answer, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const override {
        FloatMatrix help;
        _evaluateVM(help, cell, gp, tStep, rhsExpressionContext);
        // convert result to array
        if (help.isNotEmpty()) {
            help.copyColumn(answer, 1);
        }
    }
    void getDimensions(Element& cell) const override {}

}; // end class SymbolicTerm
}  // end namespace oofem
#endif // mpmsymbolic_h