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

#include "abaqususermaterial.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#ifdef _WIN32 //_MSC_VER and __MINGW32__ included
 #include <Windows.h>
#else
 #include <dlfcn.h>
#endif

#include <cstring>

namespace oofem {
REGISTER_Material(AbaqusUserMaterial);

int AbaqusUserMaterial :: n = 1;

AbaqusUserMaterial :: AbaqusUserMaterial(int n, Domain *d) :
    StructuralMaterial(n, d),
    umatobj(NULL), umat(NULL),
    mStressInterpretation(0),
    mUseNumericalTangent(false),
    mPerturbation(1.0e-7)
{ }

AbaqusUserMaterial :: ~AbaqusUserMaterial()
{
#ifdef _WIN32
    if ( this->umatobj ) {
        FreeLibrary( ( HMODULE ) this->umatobj );
    }
#else
    if ( this->umatobj ) {
        dlclose(this->umatobj);
    }

#endif
}

IRResultType AbaqusUserMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    std :: string umatname;

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, this->numState, _IFT_AbaqusUserMaterial_numState);
    IR_GIVE_FIELD(ir, this->properties, _IFT_AbaqusUserMaterial_properties);
    IR_GIVE_FIELD(ir, this->filename, _IFT_AbaqusUserMaterial_userMaterial);
    umatname = "umat";
    IR_GIVE_OPTIONAL_FIELD(ir, umatname, _IFT_AbaqusUserMaterial_name);
    strncpy(this->cmname, umatname.c_str(), 80);
    IR_GIVE_OPTIONAL_FIELD(ir, this->initialStress, _IFT_AbaqusUserMaterial_initialStress);

#ifdef _WIN32
    ///@todo Check all the windows support.
    this->umatobj = ( void * ) LoadLibrary( filename.c_str() );
    if ( !this->umatobj ) {
        DWORD dlresult = GetLastError(); //works for MinGW 32bit and MSVC
        OOFEM_ERROR("Couldn't load \"%s\",\nerror code = %d", filename.c_str(), dlresult);
    }

    //     * ( void ** )( & this->umat ) = GetProcAddress( ( HMODULE ) this->umatobj, "umat_" );
    * ( FARPROC * ) ( & this->umat ) = GetProcAddress( ( HMODULE ) this->umatobj, "umat_" ); //works for MinGW 32bit
    if ( !this->umat ) {
        //         char *dlresult = GetLastError();
        DWORD dlresult = GetLastError(); //works for MinGW 32bit
        OOFEM_ERROR("Couldn't load symbol umat,\nerror code: %d\n", dlresult);
    }

#else
    this->umatobj = dlopen(filename.c_str(), RTLD_NOW);
    if ( !this->umatobj ) {
        OOFEM_ERROR("couldn't load \"%s\",\ndlerror: %s", filename.c_str(), dlerror() );
    }

    * ( void ** ) ( & this->umat ) = dlsym(this->umatobj, "umat_");

    char *dlresult = dlerror();
    if ( dlresult ) {
        OOFEM_ERROR("couldn't load symbol umat,\ndlerror: %s\n", dlresult);
    }

#endif

    if ( ir->hasField(_IFT_AbaqusUserMaterial_numericalTangent) ) {
        mUseNumericalTangent = true;
    }

    if ( ir->hasField(_IFT_AbaqusUserMaterial_numericalTangentPerturbation) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, mPerturbation, _IFT_AbaqusUserMaterial_numericalTangentPerturbation);
        printf("mPerturbation: %e\n", mPerturbation);
    }

    return IRRT_OK;
}

void AbaqusUserMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralMaterial :: giveInputRecord(input);

    input.setField(this->numState, _IFT_AbaqusUserMaterial_numState);
    input.setField(this->properties, _IFT_AbaqusUserMaterial_properties);
    input.setField(this->filename, _IFT_AbaqusUserMaterial_userMaterial);
    input.setField(std::string(this->cmname), _IFT_AbaqusUserMaterial_name);
}

MaterialStatus *AbaqusUserMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new AbaqusUserMaterialStatus(n++, this->giveDomain(), gp, this->numState);
}


void AbaqusUserMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = dynamic_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    if ( !ms->hasTangent() ) { ///@todo Make this hack fit more nicely into OOFEM in general;
        // Evaluating the function once, so that the tangent can be obtained.
        FloatArray stress(6), strain(6);
        strain.zero();
        this->giveRealStressVector_3d(stress, gp, strain, tStep);
    }

    answer = ms->giveTempTangent();

#if 0
    double h = 1e-7;
    FloatArray strain, strainh, stress, stressh;
    strain = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStrainVector();
    stress = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveTempStressVector();
    FloatMatrix En( strain.giveSize(), strain.giveSize() );
    for ( int i = 1; i <= strain.giveSize(); ++i ) {
        strainh = strain;
        strainh.at(i) += h;
        this->giveRealStressVector_3d(stressh, form, gp, strainh, tStep);
        stressh.subtract(stress);
        stressh.times(1.0 / h);
        En.setColumn(stressh, i);
    }
    this->giveRealStressVector_3d(stress, form, gp, strain, tStep);

    printf("En = ");
    En.printYourself();
    printf("Tangent = ");
    answer.printYourself();
#endif
}


void AbaqusUserMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                              MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = dynamic_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    if ( !ms->hasTangent() ) { ///@todo Make this hack fit more nicely into OOFEM in general;
        // Evaluating the function once, so that the tangent can be obtained.
        FloatArray stress, vF;
        vF = ms->giveTempFVector();
        this->giveFirstPKStressVector_3d(stress, gp, vF, tStep);
    }

    if ( !mUseNumericalTangent ) {
        //    if(mStressInterpretation == 0) {
        answer = ms->giveTempTangent();
        /*
         *  }
         *  else {
         *              // The Abaqus Documentation of User Subroutines for UMAT Section 1.1.31 says that DDSDDE is defined as
         *              // partial(Delta(sigma))/partial(Delta(epsilon)).
         *              FloatMatrix dSdE;
         *              dSdE = ms->giveTempTangent();
         *              this->give_dPdF_from(dSdE, answer, gp);
         *  }
         */
    } else   {
        double h = mPerturbation;
        FloatArray vF, vF_h, stress, stressh;
        vF = ms->giveTempFVector();
        stress = ( ( StructuralMaterialStatus * ) gp->giveMaterialStatus() )->giveTempPVector();
        FloatMatrix En(9, 9);
        for ( int i = 1; i <= 9; ++i ) {
            vF_h = vF;
            vF_h.at(i) += h;
            this->giveFirstPKStressVector_3d(stressh, gp, vF_h, tStep);
            stressh.subtract(stress);
            stressh.times(1.0 / h);
            En.setColumn(stressh, i);
        }

        // Reset
        this->giveFirstPKStressVector_3d(stressh, gp, vF, tStep);

        /*
         *      printf("En = ");
         *      En.printYourself();
         *
         *      answer = ms->giveTempTangent();
         *      printf("Tangent = ");
         *      answer.printYourself();
         *
         *      FloatMatrix diff = En;
         *      diff.subtract(answer);
         *      printf("diff: "); diff.printYourself();
         */

        answer = En;
    }
}

void AbaqusUserMaterial :: givePlaneStrainStiffMtrx_dPdF(FloatMatrix &answer,
                                                         MatResponseMode mmode, GaussPoint *gp,
                                                         TimeStep *tStep)
{
    FloatMatrix dPdF3d;
    this->give3dMaterialStiffnessMatrix_dPdF(dPdF3d, mmode, gp, tStep);
    StructuralMaterial :: giveReducedMatrixForm(answer, dPdF3d, _PlaneStrain);
}

void AbaqusUserMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &reducedStrain, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = static_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    /* User-defined material name, left justified. Some internal material models are given names starting with
     * the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for
     * CMNAME. */
    //char cmname[80];

    // Sizes of the tensors.
    int ndi = 3;
    int nshr = 3;

    int ntens = ndi + nshr;
    FloatArray strain = ms->giveStrainVector();
    FloatArray stress = ms->giveStressVector();
    // adding the initial stress
    stress.add(initialStress);
    FloatArray strainIncrement;
    strainIncrement.beDifferenceOf(reducedStrain, strain);
    FloatArray state = ms->giveStateVector();
    FloatMatrix jacobian(ntens, ntens);
    int numProperties = this->properties.giveSize();

    // Temperature and increment
    double temp = 0.0, dtemp = 0.0;

    // Times and increment
    double dtime = tStep->giveTimeIncrement();
    ///@todo Check this. I'm just guessing. Maybe intrinsic time instead?
    double time [ 2 ] = {
        tStep->giveTargetTime() - dtime, tStep->giveTargetTime()
    };
    double pnewdt = 1.0; ///@todo Right default value? umat routines may change this (although we ignore it)

    /* Specific elastic strain energy, plastic dissipation, and “creep” dissipation, respectively. These are passed
     * in as the values at the start of the increment and should be updated to the corresponding specific energy
     * values at the end of the increment. They have no effect on the solution, except that they are used for
     * energy output. */
    double sse, spd, scd;

    // Outputs only in a fully coupled thermal-stress analysis:
    double rpl = 0.0; // Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material.
    FloatArray ddsddt(ntens); // Variation of the stress increments with respect to the temperature.
    FloatArray drplde(ntens); // Variation of RPL with respect to the strain increments.
    double drpldt = 0.0; // Variation of RPL with respect to the temperature.

    /* An array containing the coordinates of this point. These are the current coordinates if geometric
     * nonlinearity is accounted for during the step (see “Procedures: overview,” Section 6.1.1); otherwise,
     * the array contains the original coordinates of the point */
    FloatArray coords;
    gp->giveElement()->computeGlobalCoordinates( coords, gp->giveNaturalCoordinates() );

    /* Rotation increment matrix. This matrix represents the increment of rigid body rotation of the basis
     * system in which the components of stress (STRESS) and strain (STRAN) are stored. It is provided so
     * that vector- or tensor-valued state variables can be rotated appropriately in this subroutine: stress and
     * strain components are already rotated by this amount before UMAT is called. This matrix is passed in
     * as a unit matrix for small-displacement analysis and for large-displacement analysis if the basis system
     * for the material point rotates with the material (as in a shell element or when a local orientation is used).*/
    FloatMatrix drot(3, 3);
    drot.beUnitMatrix();

    /* Characteristic element length, which is a typical length of a line across an element for a first-order
     * element; it is half of the same typical length for a second-order element. For beams and trusses it is a
     * characteristic length along the element axis. For membranes and shells it is a characteristic length in
     * the reference surface. For axisymmetric elements it is a characteristic length in the
     * plane only.
     * For cohesive elements it is equal to the constitutive thickness.*/
    double celent = 0.0; /// @todo Include the characteristic element length

    /* Array containing the deformation gradient at the beginning of the increment. See the discussion
     * regarding the availability of the deformation gradient for various element types. */
    FloatMatrix dfgrd0(3, 3);
    dfgrd0.beUnitMatrix();
    /* Array containing the deformation gradient at the end of the increment. The components of this array
     * are set to zero if nonlinear geometric effects are not included in the step definition associated with
     * this increment. See the discussion regarding the availability of the deformation gradient for various
     * element types. */
    FloatMatrix dfgrd1(3, 3);
    dfgrd1.beUnitMatrix();

    int noel = gp->giveElement()->giveNumber(); // Element number.
    int npt = 0; // Integration point number.

    int layer = 0; // Layer number (for composite shells and layered solids)..
    int kspt = 0; // Section point number within the current layer.
    int kstep = 0; // Step number.
    int kinc = 0; // Increment number.

    ///@todo No idea about these parameters
    double predef;
    double dpred;

    // Change to Abaqus's component order
    stress.changeComponentOrder();
    strain.changeComponentOrder();
    strainIncrement.changeComponentOrder();

    //    printf("stress oofem: "); stress.printYourself();

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector_3d - Calling subroutine");
    this->umat(stress.givePointer(), // STRESS
               state.givePointer(), // STATEV
               jacobian.givePointer(), // DDSDDE
               & sse, // SSE
               & spd, // SPD
               & scd, // SCD
               & rpl, // RPL
               ddsddt.givePointer(), // DDSDDT
               drplde.givePointer(), // DRPLDE
               & drpldt, // DRPLDT
               strain.givePointer(), // STRAN
               strainIncrement.givePointer(), // DSTRAN
               time, // TIME
               & dtime, // DTIME
               & temp, // TEMP
               & dtemp, // DTEMP
               & predef, // PREDEF
               & dpred, // DPRED
               this->cmname, // CMNAME
               & ndi, // NDI
               & nshr, // NSHR
               & ntens, // NTENS
               & numState, // NSTATV
               properties.givePointer(), // PROPS
               & numProperties, // NPROPS
               coords.givePointer(), // COORDS
               drot.givePointer(), // DROT
               & pnewdt, // PNEWDT
               & celent, // CELENT
               dfgrd0.givePointer(), // DFGRD0
               dfgrd1.givePointer(), // DFGRD1
               & noel, // NOEL
               & npt, // NPT
               & layer, // LAYER
               & kspt, // KSPT
               & kstep, // KSTEP
               & kinc); // KINC

    // Change to OOFEM's component order
    stress.changeComponentOrder();
    // subtracking the initial stress
    stress.subtract(initialStress);
    strain.changeComponentOrder();
    strainIncrement.changeComponentOrder();
    jacobian.changeComponentOrder();

    ms->letTempStrainVectorBe(reducedStrain);
    ms->letTempStressVectorBe(stress);
    ms->letTempStateVectorBe(state);
    ms->letTempTangentBe(jacobian);
    answer = stress;

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector - Calling subroutine was successful");
}


void AbaqusUserMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &vF, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = static_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    /* User-defined material name, left justified. Some internal material models are given names starting with
     * the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for
     * CMNAME. */
    //char cmname[80];

    // Sizes of the tensors.
    int ndi = 3;
    int nshr = 3;

    int ntens = ndi + nshr;
    FloatArray stress(9); // PK1
    FloatArray strainIncrement;

    // compute Green-Lagrange strain
    FloatArray strain;
    FloatArray vS;
    FloatMatrix F, E;
    F.beMatrixForm(vF);
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    strain.beSymVectorFormOfStrain(E);

    strainIncrement.beDifferenceOf(strain, ms->giveStrainVector());
    FloatArray state = ms->giveStateVector();
    FloatMatrix jacobian(9, 9); // dPdF
    int numProperties = this->properties.giveSize();

    // Temperature and increment
    double temp = 0.0, dtemp = 0.0;

    // Times and increment
    double dtime = tStep->giveTimeIncrement();
    ///@todo Check this. I'm just guessing. Maybe intrinsic time instead?
    double time [ 2 ] = {
        tStep->giveTargetTime() - dtime, tStep->giveTargetTime()
    };
    double pnewdt = 1.0; ///@todo Right default value? umat routines may change this (although we ignore it)

    /* Specific elastic strain energy, plastic dissipation, and “creep” dissipation, respectively. These are passed
     * in as the values at the start of the increment and should be updated to the corresponding specific energy
     * values at the end of the increment. They have no effect on the solution, except that they are used for
     * energy output. */
    double sse, spd, scd;

    // Outputs only in a fully coupled thermal-stress analysis:
    double rpl = 0.0; // Volumetric heat generation per unit time at the end of the increment caused by mechanical working of the material.
    FloatArray ddsddt(ntens); // Variation of the stress increments with respect to the temperature.
    FloatArray drplde(ntens); // Variation of RPL with respect to the strain increments.
    double drpldt = 0.0; // Variation of RPL with respect to the temperature.

    /* An array containing the coordinates of this point. These are the current coordinates if geometric
     * nonlinearity is accounted for during the step (see “Procedures: overview,” Section 6.1.1); otherwise,
     * the array contains the original coordinates of the point */
    FloatArray coords;
    gp->giveElement()->computeGlobalCoordinates( coords, gp->giveNaturalCoordinates() ); ///@todo Large deformations?

    /* Rotation increment matrix. This matrix represents the increment of rigid body rotation of the basis
     * system in which the components of stress (STRESS) and strain (STRAN) are stored. It is provided so
     * that vector- or tensor-valued state variables can be rotated appropriately in this subroutine: stress and
     * strain components are already rotated by this amount before UMAT is called. This matrix is passed in
     * as a unit matrix for small-displacement analysis and for large-displacement analysis if the basis system
     * for the material point rotates with the material (as in a shell element or when a local orientation is used).*/
    FloatMatrix drot(3, 3);
    drot.beUnitMatrix();

    /* Characteristic element length, which is a typical length of a line across an element for a first-order
     * element; it is half of the same typical length for a second-order element. For beams and trusses it is a
     * characteristic length along the element axis. For membranes and shells it is a characteristic length in
     * the reference surface. For axisymmetric elements it is a characteristic length in the
     * plane only.
     * For cohesive elements it is equal to the constitutive thickness.*/
    double celent = 0.0; /// @todo Include the characteristic element length

    /* Array containing the deformation gradient at the beginning of the increment. See the discussion
     * regarding the availability of the deformation gradient for various element types. */
    FloatMatrix dfgrd0(3, 3);
    /* Array containing the deformation gradient at the end of the increment. The components of this array
     * are set to zero if nonlinear geometric effects are not included in the step definition associated with
     * this increment. See the discussion regarding the availability of the deformation gradient for various
     * element types. */
    FloatMatrix dfgrd1(3, 3);

    dfgrd0.beMatrixForm( ms->giveFVector() );
    dfgrd1.beMatrixForm(vF);

    int noel = gp->giveElement()->giveNumber(); // Element number.
    int npt = 0; // Integration point number.

    // We intentionally ignore the layer number since that is handled by the layered cross-section in OOFEM.
    int layer = 0; // Layer number (for composite shells and layered solids)..
    int kspt = 0; // Section point number within the current layer.
    int kstep = tStep->giveMetaStepNumber(); // Step number.
    int kinc = 0; // Increment number.

    ///@todo No idea about these parameters
    double predef;
    double dpred;

    // Change to Abaqus's component order
    stress.changeComponentOrder();
    strain.changeComponentOrder();
    strainIncrement.changeComponentOrder();

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector - Calling subroutine");
    this->umat(stress.givePointer(), // STRESS
               state.givePointer(), // STATEV
               jacobian.givePointer(), // DDSDDE
               & sse, // SSE
               & spd, // SPD
               & scd, // SCD
               & rpl, // RPL
               ddsddt.givePointer(), // DDSDDT
               drplde.givePointer(), // DRPLDE
               & drpldt, // DRPLDT
               strain.givePointer(), // STRAN
               strainIncrement.givePointer(), // DSTRAN
               time, // TIME
               & dtime, // DTIME
               & temp, // TEMP
               & dtemp, // DTEMP
               & predef, // PREDEF
               & dpred, // DPRED
               this->cmname, // CMNAME
               & ndi, // NDI
               & nshr, // NSHR
               & ntens, // NTENS
               & numState, // NSTATV
               properties.givePointer(), // PROPS
               & numProperties, // NPROPS
               coords.givePointer(), // COORDS
               drot.givePointer(), // DROT
               & pnewdt, // PNEWDT
               & celent, // CELENT
               dfgrd0.givePointer(), // DFGRD0
               dfgrd1.givePointer(), // DFGRD1
               & noel, // NOEL
               & npt, // NPT
               & layer, // LAYER
               & kspt, // KSPT
               & kstep, // KSTEP
               & kinc); // KINC


    // Change to OOFEM's component order
    stress.changeComponentOrder();
    strain.changeComponentOrder();
    strainIncrement.changeComponentOrder();
    jacobian.changeComponentOrder();


    if ( mStressInterpretation == 0 ) {
        FloatMatrix P, cauchyStress;
        P.beMatrixForm(stress);

        FloatArray vP;
        vP.beVectorForm(P);

        cauchyStress.beProductTOf(P, F);
        cauchyStress.times( 1.0 / F.giveDeterminant() );

        FloatArray vCauchyStress;
        vCauchyStress.beSymVectorForm(cauchyStress);

        ms->letTempStrainVectorBe(strain);
        ms->letTempStressVectorBe(vCauchyStress);
        ms->letTempStateVectorBe(state);
        ms->letTempTangentBe(jacobian);
        ms->letTempPVectorBe(vP);
        ms->letTempFVectorBe(vF);

        answer = vP;
    } else   {
        FloatArray vP;
        FloatMatrix P, sigma, F_inv;
        F_inv.beInverseOf(F);

        sigma.beMatrixForm(stress);
        P.beProductOf(F, sigma);
        vP.beVectorForm(P);
        answer = vP;

        // Convert from sigma to S
        FloatMatrix S;
        S.beProductOf(F_inv, P);
        vS.beSymVectorForm(S);

        // @todo Should convert from dsigmadE to dSdE here
        // L2=detF*matmul( matmul ( inv(op_a_V9(F,F), cm-op_a_V9(ident,Tau)-op_b_V9(Tau,ident) ), inv(op_a_V9(Ft,Ft)))

        ms->letTempStrainVectorBe(strain);
        ms->letTempStressVectorBe(vS);
        ms->letTempStateVectorBe(state);
        ms->letTempTangentBe(jacobian);
        ms->letTempPVectorBe(vP);
        ms->letTempFVectorBe(vF);
    }

    OOFEM_LOG_DEBUG("AbaqusUserMaterial :: giveRealStressVector_3d - Calling subroutine was successful");
}


int AbaqusUserMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    AbaqusUserMaterialStatus *ms = static_cast< AbaqusUserMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_Undefined || type == IST_AbaqusStateVector ) {
        // The undefined value is used to just dump the entire state vector.
        answer = ms->giveStateVector();
        return 1;
    } else if ( type == IST_StressTensor ) {
        answer = ms->giveStressVector();
        answer.add(initialStress);
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}



void AbaqusUserMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    tempStateVector = stateVector;
}

AbaqusUserMaterialStatus :: AbaqusUserMaterialStatus(int n, Domain *d, GaussPoint *gp, int numState) :
    StructuralMaterialStatus(n, d, gp),
    numState(numState), stateVector(numState), tempStateVector(numState), hasTangentFlag(false)
{
    strainVector.resize(6);
    strainVector.zero();
}

void AbaqusUserMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    stateVector = tempStateVector;
}

void AbaqusUserMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
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
