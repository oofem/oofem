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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#ifndef structuralslipfe2material_h
#define structuralslipfe2material_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralfe2material.h"

#include <memory>

///@name Input fields for StructuralSlipFE2Material
//@{
#define _IFT_StructuralSlipFE2Material_Name "structslipfe2material"
#define _IFT_StructuralSlipFE2Material_useExternalStiffness "use_ext_stiffness"
#define _IFT_StructuralSlipFE2Material_allGPResults "export_all_gps"
#define _IFT_StructuralSlipFE2Material_outputSelectedResults "output_selected_el_gps"
#define _IFT_StructuralSlipFE2Material_dStressdEps "dsde"

//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;

class StructuralSlipFE2MaterialStatus : public StructuralFE2MaterialStatus
{
protected:
    int nel=1; //macroscopic element number
    int gpn=1; //gauss point number

public:
    StructuralSlipFE2MaterialStatus(int rank, GaussPoint * g,  const std :: string & inputfile);

    void setTangent(const FloatMatrix &iTangent) {tangent = iTangent; oldTangent = false;}
    const char *giveClassName() const override { return "StructuralSlipFE2MaterialStatus"; }
};


/**
 * Multiscale constitutive model for structural problems on reinforced concrete structures.
 * At the macroscale, both the displacement and reinforcement slip can be treated as variables.
 * Currently, only plane stress mode is supported.
 *
 * This material uses the PrescribedGradient boundary conditions to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedGradient boundary condition.
 * - It must be the first boundary condition
 *
 * For the macroscopic displacement field, the following boundary conditions on the RVE are supported (must use PrescribedGradientHomogenization class):
 * - PrescribedGradientBCDirichletRC (macroscopic displacement gradient prescribed with Dirichlet BCs on concrete and optionally on reinforcement)
 * - PrescribedGradientBCNeumannRC (macroscopic displacement gradient prescribed with Neumann BCs on concrete)
 * - PrescribedGradientMultiple (for combining more than one PrescribedGradient boundary condition on the RVE)
 * ref: Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2018). Two‐scale finite element modelling of reinforced concrete structures: Effective response and subscale fracture development. International Journal for Numerical Methods in Engineering, 114(10), 1074–1102. https://doi.org/10.1002/nme.5776
 * *
 * @author Adam Sciegaj
 */
class StructuralSlipFE2Material : public StructuralFE2Material
{
protected:
    bool useExtStiff = false;
    bool allGPRes = false;
    /**
     * FloatArray that contains element and Gauss point number for the desired results to be exported.
     * The syntax in the input file should follow:
     * output_selected_el_gp size elNumber gpNumber
     * e.g. output_selected_el_gp 6 22 4 21 2 20 3 will export Gauss point 4 from el 22, Gauss point 2 from el 21 etc.
     */
    IntArray outputSelected;

    FloatMatrix givendStressdEpsTangent; //3x3 stiffness matrix (plane stress)

public:
    StructuralSlipFE2Material(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveInputRecordName() const override { return _IFT_StructuralSlipFE2Material_Name; }
    const char *giveClassName() const override { return "StructuralSlipFE2Material"; }

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};

} // end namespace oofem
#endif // structuralslipfe2material_h
