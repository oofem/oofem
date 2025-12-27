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

// Use the following flags and specify location of the cmake files.
// cmake -DUSE_MFRONT=ON -DMFrontGenericInterface_DIR=/usr/local/share/mgis/cmake/ ..

#include "mfrontusermaterial.h"
#include "gausspoint.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif

#include <cstring>

namespace oofem {
REGISTER_Material(MFrontUserMaterial);

// eventual reindexing
int const MFrontUserMaterial :: mfront2oo9[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
int const MFrontUserMaterial :: mfront2oo6[6] = {0, 1, 2, 3, 4, 5};

MFrontUserMaterial :: MFrontUserMaterial(int n, Domain *d) : StructuralMaterial(n, d) {}

MFrontUserMaterial :: ~MFrontUserMaterial()
{

}

void MFrontUserMaterial :: initializeFrom(InputRecord &ir)
{
    using namespace mgis::behaviour;

    StructuralMaterial :: initializeFrom(ir);


    FloatArray s(6);
    this->initialStress = s;

    std :: string libname;
    std :: string modelname;
    IR_GIVE_FIELD(ir, libname, _IFT_MFrontUserMaterial_libpath);
    IR_GIVE_FIELD(ir, modelname, _IFT_MFrontUserMaterial_modelname);
    strncpy(this->libname, libname.c_str(), 199);
    strncpy(this->modelname, modelname.c_str(), 199);
    // loading the MFront library using MGIS
    this->behaviour = std::make_unique<Behaviour>(load(this->libname, this->modelname, Hypothesis::TRIDIMENSIONAL));
    // initialize material properties from user data
    IR_GIVE_FIELD(ir, this->properties, _IFT_MFrontUserMaterial_properties);
    // check that the number of material properties given by the user is correct
    const auto nprops = mgis::behaviour::getArraySize(this->behaviour->mps, this->behaviour->hypothesis);
    if (this->properties.giveSize() != nprops) {
        throw(std::runtime_error("wrong number of material properties"));
    }
    // treating the case for external state variables
    const auto nesvs = mgis::behaviour::getArraySize(this->behaviour->esvs, this->behaviour->hypothesis);
    if (nesvs != 1) {
        throw(std::runtime_error("wrong number of external state variables"));
    }
    if (this->behaviour->esvs[0].name != "Temperature") {
        throw(std::runtime_error("the temperature is the only external state variable supported so far"));
    }

}

void MFrontUserMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);
}

MaterialStatus *MFrontUserMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MFrontUserMaterialStatus(gp, *(this->behaviour));
}


FloatMatrixF<6,6>
MFrontUserMaterial :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = dynamic_cast< MFrontUserMaterialStatus * >( this->giveStatus(gp) );
    if (!ms->hasTangent()) {
        // Evaluating the function once, so that the tangent can be obtained.
        const_cast<MFrontUserMaterial*>(this)->giveRealStressVector_3d(zeros<6>(), gp, tStep);
    }

    double h = mPerturbation;
    FloatArray strain, strainh, stress, stressh;
    strain = static_cast<StructuralMaterialStatus *>(gp->giveMaterialStatus())->giveTempStrainVector();
    stress = static_cast<StructuralMaterialStatus *>(gp->giveMaterialStatus())->giveTempStressVector();
    FloatMatrix En(strain.giveSize(), strain.giveSize());
    for (int i = 1; i <= strain.giveSize(); ++i) {
        strainh = strain;
        strainh.at(i) += h;
        stressh = this->giveRealStressVector_3d(strainh, gp, tStep);
        stressh.subtract(stress);
        stressh.times(1.0 / h);
        En.setColumn(stressh, i);
    }
    stress = this->giveRealStressVector_3d(strain, gp, tStep);

    return En;
}


FloatMatrixF<9,9>
MFrontUserMaterial :: give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = dynamic_cast<MFrontUserMaterialStatus * >(this->giveStatus(gp));
    if (!ms->hasTangent()) {
        // Evaluating the function once, so that the tangent can be obtained.
        const auto &vF = ms->giveTempFVector();
        this->giveFirstPKStressVector_3d(vF, gp, tStep);
    }

    double h = mPerturbation;
    auto const &vF = ms->giveTempFVector();
    auto const &stress = ms->giveTempPVector();
    FloatMatrixF<9,9> En;
    for ( int i = 1; i <= 9; ++i ) {
        auto vF_h = vF;
        vF_h.at(i) += h;
        auto stressh = this->giveFirstPKStressVector_3d(vF_h, gp, tStep);
        auto ds = (stressh - stress) * (1. / h);
        En.setColumn(ds, i);
    }

    // Reset
    this->giveFirstPKStressVector_3d(vF, gp, tStep);

    return En;
}

FloatMatrixF<5,5>
MFrontUserMaterial :: givePlaneStrainStiffMtrx_dPdF(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const
{
    auto dPdF3D = this->give3dMaterialStiffnessMatrix_dPdF(mmode, gp, tStep);
    return dPdF3D({0, 1, 2, 5, 8}, {0, 1, 2, 5, 8});
}


FloatArrayF<6>
MFrontUserMaterial :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< MFrontUserMaterialStatus * >( this->giveStatus(gp) );

    FloatArrayF<6> old_strain = ms->giveStrainVector();
    FloatArrayF<6> old_stress = initialStress + ms->giveStressVector();

    auto old_state = ms->giveStateVector();
    auto new_state = old_state;
    FloatMatrixF<6,6> mfront_jacobian;

    // Times and increment
    double dtime = tStep->giveTimeIncrement();

    // Change the component order
    auto mfront_stress = old_stress[mfront2oo6];
    auto mfront_old_strain = old_strain[mfront2oo6];
    auto mfront_new_strain = strain[mfront2oo6];

    // Beware of MFront conventions here !!!!
    mfront_stress[3] *= 2.;
    mfront_stress[4] *= 2.;
    mfront_stress[5] *= 2.;
    mfront_old_strain[3] /= 2.;
    mfront_old_strain[4] /= 2.;
    mfront_old_strain[5] /= 2.;
    mfront_new_strain[3] /= 2.;
    mfront_new_strain[4] /= 2.;
    mfront_new_strain[5] /= 2.;

    //
    double externalStateVariables[1] = {293.15};
    double stored_energy = 0;
    double dissipated_energy = 0;

    mgis::behaviour::BehaviourDataView v;
    v.dt = dtime;

// modifications du to API changes in MGIS
#ifdef MGIS_BEHAVIOUR_API_VERSION
#if MGIS_BEHAVIOUR_API_VERSION == 1
    auto rdt = double{};
    v.rdt = &rdt;
    v.speed_of_sound = nullptr;
#else /* MGIS_BEHAVIOUR_API_VERSION == 1 */
#error "Unsupported MGIS' behaviour API"
#endif  /* MGIS_BEHAVIOUR_API_VERSION == 1 */
#else /* MGIS_BEHAVIOUR_API_VERSION */
    v.rdt = 1;
#endif /* MGIS_BEHAVIOUR_API_VERSION */

    v.K = &mfront_jacobian(0, 0);
    // perform an integration with computation of the consistent tangent operator
    v.K[0] = 4;
    // setting the initial state
    v.s0.gradients                = &mfront_old_strain[0];
    v.s0.thermodynamic_forces     = &old_stress[0];
    v.s0.material_properties      = nullptr; // &(this->properties[0])
    v.s0.internal_state_variables = &old_state[0]; // shall be null if numState==0
    v.s0.external_state_variables = externalStateVariables;
    v.s0.stored_energy            = &stored_energy;
    v.s0.dissipated_energy        = &dissipated_energy;
    // setting the state at the end of the time step
    v.s1.gradients                = &mfront_new_strain[0];
    v.s1.thermodynamic_forces     = &mfront_stress[0];
    v.s1.material_properties      = nullptr; // &(this->properties[0])
    v.s1.internal_state_variables = &new_state[0]; // shall be null if numState==0
    v.s1.external_state_variables = externalStateVariables;
    v.s1.stored_energy            = &stored_energy;
    v.s1.dissipated_energy        = &dissipated_energy;

    integrate(v, *(this->behaviour));

    // Change to OOFEM's component order
    auto jacobian = mfront_jacobian(mfront2oo6, mfront2oo6);
    // subtracking the initial stress
    auto stress = mfront_stress[mfront2oo6] - initialStress;
    // Here MFront->OOFEM conventions for the strain **and** the jacobian
    //

    ms->letTempStrainVectorBe(strain);
    ms->letTempStressVectorBe(stress);
    ms->letTempStateVectorBe(new_state);
    ms->letTempTangentBe(jacobian);
    return stress;
}


int MFrontUserMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto ms = static_cast<MFrontUserMaterialStatus *>(this->giveStatus(gp));
    if ( type == IST_Undefined || type == IST_AbaqusStateVector ) {
        // The undefined value is used to just dump the entire state vector.
        answer = ms->giveStateVector();
        return 1;
    } else if (type == IST_StressTensor) {
        answer = ms->giveStressVector();
        answer.add(initialStress);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

void MFrontUserMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    tempStateVector = stateVector;
}

MFrontUserMaterialStatus :: MFrontUserMaterialStatus(GaussPoint *gp, const mgis::behaviour::Behaviour &b) :
    StructuralMaterialStatus(gp), hasTangentFlag(false)
{
    strainVector.resize(6);
    strainVector.zero();
    // getting number of state variables
    const auto numState = mgis::behaviour::getArraySize(b.isvs, b.hypothesis);
    // allocating  number of state variables
    this->stateVector.resize(numState);
    this->tempStateVector.resize(numState);
    this->stateVector.zero(); // ??
    this->tempStateVector.zero(); // ??
}

void MFrontUserMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    stateVector = tempStateVector;
}

void MFrontUserMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(File, tStep);

    fprintf(File, "  stateVector ");
    FloatArray state = this->giveStateVector();
    for ( auto &var : state ) {
        fprintf( File, " % .4e", var );
    }

    fprintf(File, "\n");
}

} // end namespace oofem
