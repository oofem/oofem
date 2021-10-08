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
#include "prescribeddispsliphomogenization.h"

#include <memory>

///@name Input fields for StructuralSlipFE2Material
//@{
#define _IFT_StructuralSlipFE2Material_Name "structslipfe2material"
#define _IFT_StructuralSlipFE2Material_useExternalStiffness "use_ext_stiffness"
#define _IFT_StructuralSlipFE2Material_allGPResults "export_all_gps"
#define _IFT_StructuralSlipFE2Material_outputSelectedResults "output_selected_el_gps"
#define _IFT_StructuralSlipFE2Material_dStressdEps "dsde"
#define _IFT_StructuralSlipFE2Material_dBStressdEps "dbsde"
#define _IFT_StructuralSlipFE2Material_dRStressdEps "drsde"
#define _IFT_StructuralSlipFE2Material_dStressdS "dsds"
#define _IFT_StructuralSlipFE2Material_dBStressdS "dbsds"
#define _IFT_StructuralSlipFE2Material_dRStressdS "drsds"
#define _IFT_StructuralSlipFE2Material_dStressdG "dsdg"
#define _IFT_StructuralSlipFE2Material_dBStressdG "dbsdg"
#define _IFT_StructuralSlipFE2Material_dRStressdG "drsdg"

//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;
class PrescribedDispSlipHomogenization;

class StructuralSlipFE2MaterialStatus : public StructuralMaterialStatus
{
protected:
    /// The RVE
    std :: unique_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    PrescribedDispSlipHomogenization *bc = nullptr;

    std :: string mInputFile;

    int nel=1; //macroscopic element number
    int gpn=1; //gauss point number

    //Tangents wrt to strain, slip and slip gradient
    FloatMatrix dStressdEpsTangent, dBStressdEpsTangent, dRStressdEpsTangent;
    bool olddSdETangent = true, olddBSdETangent = true, olddRSdETangent =true;

    FloatMatrix dStressdSTangent, dBStressdSTangent, dRStressdSTangent;
    bool olddSdSTangent = true, olddBSdSTangent = true, olddRSdSTangent = true;

    FloatMatrix dStressdGTangent, dBStressdGTangent, dRStressdGTangent;
    bool olddSdGTangent = true, olddBSdGTangent = true, olddRSdGTangent = true;

    FloatArray slipVector;
    FloatArray bStressVector;
    FloatArray tempSlipVector;
    FloatArray tempBStressVector;

    FloatArray slipGradVector;
    FloatArray rStressVector;
    FloatArray tempSlipGradVector;
    FloatArray tempRStressVector;

public:
    StructuralSlipFE2MaterialStatus(int rank, GaussPoint * g,  const std :: string & inputfile, int el, int gp);

    EngngModel *giveRVE() const { return this->rve.get(); }
    PrescribedDispSlipHomogenization *giveBC();

    bool createRVE(const std :: string &inputfile, int rank, int el, int gp);
    void setTimeStep(TimeStep *tStep);

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    const char *giveClassName() const override { return "StructuralSlipFE2MaterialStatus"; }

    /// Setters and getters
    FloatMatrix &givedStressdEpsTangent() { return dStressdEpsTangent; }
    void setdStressdEpsTangent(const FloatMatrix &iTangent) {dStressdEpsTangent = iTangent; olddSdETangent = false;}
    void setdBStressdEpsTangent(const FloatMatrix &iTangent) {dBStressdEpsTangent = iTangent; olddBSdETangent = false;}
    void setdRStressdEpsTangent(const FloatMatrix &iTangent) {dRStressdEpsTangent = iTangent; olddRSdETangent = false;}

    void setdStressdSTangent(const FloatMatrix &iTangent) {dStressdSTangent = iTangent; olddSdSTangent = false;}
    void setdBStressdSTangent(const FloatMatrix &iTangent) {dBStressdSTangent = iTangent; olddBSdSTangent = false;}
    void setdRStressdSTangent(const FloatMatrix &iTangent) {dRStressdSTangent = iTangent; olddRSdSTangent = false;}

    void setdStressdGTangent(const FloatMatrix &iTangent) {dStressdGTangent = iTangent; olddSdGTangent = false;}
    void setdBStressdGTangent(const FloatMatrix &iTangent) {dBStressdGTangent = iTangent; olddBSdGTangent = false;}
    void setdRStressdGTangent(const FloatMatrix &iTangent) {dRStressdGTangent = iTangent; olddRSdGTangent = false;}

    const FloatArray &giveSlipVector() const { return slipVector; }
    const FloatArray &giveTransferStressVector() const { return bStressVector; }
    const FloatArray &giveTempSlipVector() const { return tempSlipVector; }
    void letTempTransferStressVectorBe(const FloatArray &v) { tempBStressVector = v; }
    void letTempSlipVectorBe(const FloatArray &v) { tempSlipVector = v; }

    const FloatArray &giveSlipGradVector() const { return slipGradVector; }
    const FloatArray &giveReinfStressVector() const { return rStressVector; }
    const FloatArray &giveTempSlipGradVector() const { return tempSlipGradVector; }
    void letTempReinfStressVectorBe(const FloatArray &v) { tempRStressVector = v; }
    void letTempSlipGradVectorBe(const FloatArray &v) { tempSlipGradVector = v; }
};


/**
 * Multiscale constitutive model for structural problems on reinforced concrete structures.
 * At the macroscale, both the displacement and reinforcement slip can be treated as variables.
 * Currently, only plane stress mode is supported.
 *
 * This material uses the PrescribedDispSlip boundary conditions to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedDispSlip boundary condition.
 * - It must be the first boundary condition
 *
 * The following boundary conditions on the RVE are supported (must use PrescribedDispSlipHomogenization class):
 * - PrescribedDispSlipBCDirichletRC (macroscopic displacement gradient, slip and slip gradient prescribed with Dirichlet BCs on concrete and optionally on reinforcement)
 * - PrescribedDispSlipBCNeumannRC (macroscopic displacement gradient, slip and slip gradient prescribed with Neumann BCs on concrete and optionally on reinforcement)
 * - PrescribedDispSlipMultiple (for combining more than one PrescribedDispSlip boundary condition on the RVE)
 * Input fields slip and slipGrad are optional. If left unspecified, the BCs work then as usual PrescribedGradient BCs (PrescribedGradient, PrescribedGradientBCNeumann)
 *
 * refs: Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2018). Two‐scale finite element modelling of reinforced concrete structures: Effective response and subscale fracture development. International Journal for Numerical Methods in Engineering, 114(10), 1074–1102. https://doi.org/10.1002/nme.5776
 * Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2019). A multiscale model for reinforced concrete with macroscopic variation of reinforcement slip. Computational Mechanics, 63(2), 139–158. https://doi.org/10.1007/s00466-018-1588-3
 * Sciegaj, A., Larsson, F., Lundgren, K., & Runesson, K. (2020). On a volume averaged measure of macroscopic reinforcement slip in two-scale modeling of reinforced concrete. International Journal for Numerical Methods in Engineering, 121(8), 1822–1846. https://doi.org/10.1002/nme.6288
 *
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

    //sensitivities
    FloatMatrix givendStressdEpsTangent; //3x3 stiffness matrix (plane stress)
    FloatMatrix givendBStressdEpsTangent, givendRStressdEpsTangent;
    FloatMatrix givendStressdSTangent, givendBStressdSTangent, givendRStressdSTangent;
    FloatMatrix givendStressdGTangent, givendBStressdGTangent, givendRStressdGTangent;

public:
    StructuralSlipFE2Material(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveInputRecordName() const override { return _IFT_StructuralSlipFE2Material_Name; }
    const char *giveClassName() const override { return "StructuralSlipFE2Material"; }

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    /**
    * Computes the homogenized stress, homogenized bond stress and reinforcement stress for the createRVE
    */
    virtual void giveHomogenizedFields(FloatArray &stress, FloatArray &bStress, FloatArray &rStress, const FloatArray &strain, const FloatArray &slip,
                                       const FloatArray &slipGradient, GaussPoint *gp, TimeStep *tStep);
    /**
     * Computes the sensitivity matrices for the RVE. Necessary for building the stiffness matrix.
     */
    virtual void giveSensitivities(FloatMatrix &dStressdEps, FloatMatrix &dStressdS, FloatMatrix &dStressdG, FloatMatrix &dBStressdEps, FloatMatrix &dBStressdS,
                           FloatMatrix &dBStressdG, FloatMatrix &dRStressdEps, FloatMatrix &dRStressdS, FloatMatrix &dRStressdG, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

};

} // end namespace oofem
#endif // structuralslipfe2material_h
