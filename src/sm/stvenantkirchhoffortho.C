/*
 * stvenantkirchhoffortho.C
 *
 *  Created on: May 28, 2013
 *      Author: svennine
 */

#include "stvenantkirchhoffortho.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( StVenantKirchhoffOrtho );

StVenantKirchhoffOrtho::StVenantKirchhoffOrtho(int n, Domain *d) : StructuralMaterial(n, d)
{

}

StVenantKirchhoffOrtho::~StVenantKirchhoffOrtho() {

}


int
StVenantKirchhoffOrtho :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports the given mode
//
{
    if ( mode == _3dMat ) {
        return 1;
    }

    return 0;
}

void
StVenantKirchhoffOrtho :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode, GaussPoint *gp, TimeStep *atTime)

// returns the 6x6 tangent stiffness matrix

{

	// Isotropic part
	double lambda_i = E_i*nu_i/((1.0+nu_i)*(1.0-2.0*nu_i));
	double G = E_i/(2.0*(1.0+nu_i));
	double c1 = lambda_i + 2.0*G;

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = c1;
    answer.at(2, 2) = c1;
    answer.at(3, 3) = c1;
    answer.at(4, 4) = 1.0*G;
    answer.at(5, 5) = 1.0*G;
    answer.at(6, 6) = 1.0*G;
    answer.at(1, 2) = answer.at(2, 1) = lambda_i;
    answer.at(1, 3) = answer.at(3, 1) = lambda_i;
    answer.at(1, 4) = answer.at(4, 1) = 0.0;
    answer.at(1, 5) = answer.at(5, 1) = 0.0;
    answer.at(1, 6) = answer.at(6, 1) = 0.0;
    answer.at(2, 3) = answer.at(3, 2) = lambda_i;
    answer.at(2, 4) = answer.at(4, 2) = 0.0;
    answer.at(2, 5) = answer.at(5, 2) = 0.0;
    answer.at(2, 6) = answer.at(6, 2) = 0.0;
    answer.at(3, 4) = answer.at(4, 3) = 0.0;
    answer.at(3, 5) = answer.at(5, 3) = 0.0;
    answer.at(3, 6) = answer.at(6, 3) = 0.0;
    answer.at(4, 5) = answer.at(5, 4) = 0.0;
    answer.at(4, 6) = answer.at(6, 4) = 0.0;
    answer.at(5, 6) = answer.at(6, 5) = 0.0;



    // Fiber part

    // Help variables following eqns 55a-55g in Bonet & Burton (1998)
//    double n = EA/E_i;
    double n 		= E_i/EA;
    double m 		= 1.0 - nu_i - 2.0*n*nu_i*nu_i;
    double lambda 	= E_i*(nu_i + n*nu_i*nu_i)/( m*(1.0+nu_i) );
    double mu 		= E_i/(2.0*(1.0+nu_i));
    double alpha 	= (mu - GA);
    double beta 	= ( E_i*nu_i*nu_i*(1.0-n)/( 4.0*m*(1.0+nu_i) ) );
    double gamma 	= ( EA*(1.0-nu_i)/(8.0*m) - (lambda+2.0*mu)/8.0 + 0.5*alpha - beta );

    double A1_2 = VecA.at(1)*VecA.at(1);
    double A2_2 = VecA.at(2)*VecA.at(2);
    double A3_2 = VecA.at(3)*VecA.at(3);



    answer.at(1, 1) += 8.0*gamma*A1_2*A1_2 + (8.0*beta-4.0*alpha)*A1_2;

    double tmp = 8.0*gamma*A1_2*A2_2 + 4.0*beta*A1_2 + 4.0*beta*A2_2;
    answer.at(1, 2) += tmp;
    answer.at(2, 1) += tmp;

    tmp = 8.0*gamma*A1_2*A3_2 + 4.0*beta*A1_2 + 4.0*beta*A3_2;
    answer.at(1, 3) += tmp;
    answer.at(3, 1) += tmp;

    tmp = 8.0*gamma*A1_2*VecA.at(2)*VecA.at(3) + 4.0*beta*VecA.at(2)*VecA.at(3);
    answer.at(1, 4) += tmp;
    answer.at(4, 1) += tmp;

    tmp = 8.0*gamma*A1_2*VecA.at(1)*VecA.at(3) + (4.0*beta-2.0*alpha)*VecA.at(1)*VecA.at(3);
    answer.at(1, 5) += tmp;
    answer.at(5, 1) += tmp;

    tmp = 8.0*gamma*A1_2*VecA.at(1)*VecA.at(2) + (4.0*beta-2.0*alpha)*VecA.at(1)*VecA.at(2);
    answer.at(1, 6) += tmp;
    answer.at(6, 1) += tmp;



    answer.at(2, 2) += 8.0*gamma*A2_2*A2_2 + (8.0*beta-4.0*alpha)*A2_2;

    tmp = 8.0*gamma*A2_2*A3_2 + 4.0*beta*A2_2 + 4.0*beta*A3_2;
    answer.at(2, 3) += tmp;
    answer.at(3, 2) += tmp;

    tmp = 8.0*gamma*A2_2*VecA.at(2)*VecA.at(3) + (4.0*beta-2.0*alpha)*VecA.at(2)*VecA.at(3);
    answer.at(2, 4) += tmp;
    answer.at(4, 2) += tmp;

    tmp = 8.0*gamma*A2_2*VecA.at(1)*VecA.at(3) + 4.0*beta*VecA.at(1)*VecA.at(3);
    answer.at(2, 5) += tmp;
    answer.at(5, 2) += tmp;

    tmp = 8.0*gamma*A2_2*VecA.at(1)*VecA.at(2) + 4.0*beta*VecA.at(1)*VecA.at(2) -2.0*alpha*( VecA.at(1)*VecA.at(2) );
    answer.at(2, 6) += tmp;
    answer.at(6, 2) += tmp;

    answer.at(3, 3) += 8.0*gamma*A3_2*A3_2 + (8.0*beta-4.0*alpha)*A3_2;

    tmp = 8.0*gamma*A3_2*VecA.at(2)*VecA.at(3) + 4.0*beta*VecA.at(2)*VecA.at(3) -2.0*alpha*( VecA.at(2)*VecA.at(3) );
    answer.at(3, 4) += tmp;
    answer.at(4, 3) += tmp;

    tmp = 8.0*gamma*A3_2*VecA.at(1)*VecA.at(3) + 4.0*beta*VecA.at(1)*VecA.at(3) -2.0*alpha*( VecA.at(1)*VecA.at(3) );
    answer.at(3, 5) += tmp;
    answer.at(5, 3) += tmp;

    tmp = 8.0*gamma*A3_2*VecA.at(1)*VecA.at(2) + 4.0*beta*VecA.at(1)*VecA.at(2);
    answer.at(3, 6) += tmp;
    answer.at(6, 3) += tmp;



    answer.at(4, 4) += 8.0*gamma*VecA.at(2)*VecA.at(3)*VecA.at(2)*VecA.at(3) - alpha*(VecA.at(3)*VecA.at(3) + VecA.at(2)*VecA.at(2) );

    tmp = 8.0*gamma*VecA.at(2)*VecA.at(3)*VecA.at(1)*VecA.at(3) - alpha*VecA.at(1)*VecA.at(2);
    answer.at(4, 5) += tmp;
    answer.at(5, 4) += tmp;

    tmp = 8.0*gamma*VecA.at(2)*VecA.at(3)*VecA.at(1)*VecA.at(2) - alpha*VecA.at(1)*VecA.at(3);
    answer.at(4, 6) += tmp;
    answer.at(6, 4) += tmp;


    answer.at(5, 5) += 8.0*gamma*VecA.at(1)*VecA.at(3)*VecA.at(1)*VecA.at(3) - alpha*(VecA.at(3)*VecA.at(3) + VecA.at(1)*VecA.at(1));

    tmp = 8.0*gamma*VecA.at(1)*VecA.at(3)*VecA.at(1)*VecA.at(2) - alpha*VecA.at(3)*VecA.at(2);
    answer.at(5, 6) += tmp;
    answer.at(6, 5) += tmp;

    answer.at(6, 6) += 8.0*gamma*VecA.at(1)*VecA.at(2)*VecA.at(1)*VecA.at(2) - alpha*(VecA.at(2)*VecA.at(2) + VecA.at(1)*VecA.at(1));


}

void
StVenantKirchhoffOrtho :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime)

// returns 6 components of the stress corresponding to the given total strain

{
	// Isotropic part

	double lambda_i = E_i*nu_i/((1.0+nu_i)*(1.0-2.0*nu_i));
	double G = E_i/(2.0*(1.0+nu_i));

    FloatArray strainVector;

    StVenantKirchhoffOrthoMaterialStatus *status = static_cast< StVenantKirchhoffOrthoMaterialStatus * >( this->giveStatus(gp) );
    this->giveStressDependentPartOfStrainVector(strainVector, gp,
                                                totalStrain,
                                                atTime, VM_Total);

    answer.resize(6);
    double c1 = lambda_i*( strainVector.at(1) + strainVector.at(2) + strainVector.at(3) );

    answer.at(1) = c1 + 2.0*G*strainVector.at(1);
    answer.at(2) = c1 + 2.0*G*strainVector.at(2);
    answer.at(3) = c1 + 2.0*G*strainVector.at(3);

    answer.at(4) = 		1.0*G*strainVector.at(4);
    answer.at(5) = 		1.0*G*strainVector.at(5);
    answer.at(6) = 		1.0*G*strainVector.at(6);



    // Fiber part

    // Help variables following eqns 55a-55g in Bonet & Burton (1998)
//    double n = EA/E_i;
    double n 		= E_i/EA;
    double m 		= 1.0 - nu_i - 2.0*n*nu_i*nu_i;
    double lambda 	= E_i*(nu_i + n*nu_i*nu_i)/( m*(1.0+nu_i) );
    double mu 		= E_i/(2.0*(1.0+nu_i));
    double alpha 	= (mu - GA);
    double beta 	= (E_i*nu_i*nu_i*(1.0-n)/( 4.0*m*(1.0+nu_i) ) );
    double gamma 	= ( EA*(1.0-nu_i)/(8.0*m) - (lambda+2.0*mu)/8.0 + 0.5*alpha - beta );



    FloatMatrix C(3, 3);
    C.at(1, 1) = 1. + 2. * strainVector.at(1);
    C.at(2, 2) = 1. + 2. * strainVector.at(2);
    C.at(3, 3) = 1. + 2. * strainVector.at(3);
    C.at(1, 2) = C.at(2, 1) = strainVector.at(6);
    C.at(1, 3) = C.at(3, 1) = strainVector.at(5);
    C.at(2, 3) = C.at(3, 2) = strainVector.at(4);


    double I1 = C.at(1,1) + C.at(2,2) + C.at(3,3);

    FloatArray CA(3);
    CA.at(1) = C.at(1,1)*VecA.at(1) + C.at(1,2)*VecA.at(2) + C.at(1,3)*VecA.at(3);
    CA.at(2) = C.at(2,1)*VecA.at(1) + C.at(2,2)*VecA.at(2) + C.at(2,3)*VecA.at(3);
    CA.at(3) = C.at(3,1)*VecA.at(1) + C.at(3,2)*VecA.at(2) + C.at(3,3)*VecA.at(3);

    double I4 = VecA.at(1)*CA.at(1) + VecA.at(2)*CA.at(2) + VecA.at(3)*CA.at(3);

    FloatMatrix S_trn(3,3);
    S_trn.zero();

    double c2 = 2.0*beta*(I4-1.0);
    S_trn.at(1,1) += c2;
    S_trn.at(2,2) += c2;
    S_trn.at(3,3) += c2;

    double c3 = 2.0*( 1.0*alpha + beta*(I1-3.0) + 2.0*gamma*(I4-1.0) );
    S_trn.at(1,1) += c3*VecA.at(1)*VecA.at(1);
    S_trn.at(1,2) += c3*VecA.at(1)*VecA.at(2);
    S_trn.at(1,3) += c3*VecA.at(1)*VecA.at(3);

    S_trn.at(2,1) += c3*VecA.at(2)*VecA.at(1);
    S_trn.at(2,2) += c3*VecA.at(2)*VecA.at(2);
    S_trn.at(2,3) += c3*VecA.at(2)*VecA.at(3);

    S_trn.at(3,1) += c3*VecA.at(3)*VecA.at(1);
    S_trn.at(3,2) += c3*VecA.at(3)*VecA.at(2);
    S_trn.at(3,3) += c3*VecA.at(3)*VecA.at(3);






    S_trn.at(1,1) += -alpha*CA.at(1)*VecA.at(1);
    S_trn.at(1,2) += -alpha*CA.at(1)*VecA.at(2);
    S_trn.at(1,3) += -alpha*CA.at(1)*VecA.at(3);

    S_trn.at(2,1) += -alpha*CA.at(2)*VecA.at(1);
    S_trn.at(2,2) += -alpha*CA.at(2)*VecA.at(2);
    S_trn.at(2,3) += -alpha*CA.at(2)*VecA.at(3);

    S_trn.at(3,1) += -alpha*CA.at(3)*VecA.at(1);
    S_trn.at(3,2) += -alpha*CA.at(3)*VecA.at(2);
    S_trn.at(3,3) += -alpha*CA.at(3)*VecA.at(3);




    S_trn.at(1,1) += -alpha*VecA.at(1)*CA.at(1);
    S_trn.at(1,2) += -alpha*VecA.at(1)*CA.at(2);
    S_trn.at(1,3) += -alpha*VecA.at(1)*CA.at(3);

    S_trn.at(2,1) += -alpha*VecA.at(2)*CA.at(1);
    S_trn.at(2,2) += -alpha*VecA.at(2)*CA.at(2);
    S_trn.at(2,3) += -alpha*VecA.at(2)*CA.at(3);

    S_trn.at(3,1) += -alpha*VecA.at(3)*CA.at(1);
    S_trn.at(3,2) += -alpha*VecA.at(3)*CA.at(2);
    S_trn.at(3,3) += -alpha*VecA.at(3)*CA.at(3);


    answer.at(1) += 1.0*S_trn.at(1,1);
    answer.at(2) += 1.0*S_trn.at(2,2);
    answer.at(3) += 1.0*S_trn.at(3,3);

    answer.at(4) += 0.5*(S_trn.at(2,3) + S_trn.at(3,2) );
    answer.at(5) += 0.5*(S_trn.at(1,3) + S_trn.at(3,1) );
    answer.at(6) += 0.5*(S_trn.at(1,2) + S_trn.at(2,1) );







    // update gp
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

}

MaterialStatus *
StVenantKirchhoffOrtho :: CreateStatus(GaussPoint *gp) const
{
    StructuralMaterialStatus *status;

    status = new StructuralMaterialStatus(1, this->giveDomain(), gp);
    return status;
}


IRResultType
StVenantKirchhoffOrtho :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    StructuralMaterial :: initializeFrom(ir);

    // Read material properties here
    VecA.resize(3);

    // Fiber direction vector
    IR_GIVE_FIELD(ir, VecA.at(1)	, _IFT_StVenantKirchhoffOrtho_ex);
    IR_GIVE_FIELD(ir, VecA.at(2)	, _IFT_StVenantKirchhoffOrtho_ey);
    IR_GIVE_FIELD(ir, VecA.at(3)	, _IFT_StVenantKirchhoffOrtho_ez);
    VecA.normalize();
//    printf("A: (%e, %e, %e)\n", VecA.at(1), VecA.at(2), VecA.at(3) );

    IR_GIVE_FIELD(ir, E_i		, _IFT_StVenantKirchhoffOrtho_E_i);
    IR_GIVE_FIELD(ir, nu_i		, _IFT_StVenantKirchhoffOrtho_nu_i);
    IR_GIVE_FIELD(ir, EA		, _IFT_StVenantKirchhoffOrtho_EA);
    IR_GIVE_FIELD(ir, GA		, _IFT_StVenantKirchhoffOrtho_GA);
//    printf("E_i: %e nu_i: %e EA: %e GA: %e\n", E_i, nu_i, EA, GA);

    return IRRT_OK;
}


StVenantKirchhoffOrthoMaterialStatus :: StVenantKirchhoffOrthoMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    // init state variables
}


StVenantKirchhoffOrthoMaterialStatus :: ~StVenantKirchhoffOrthoMaterialStatus()
{}


void
StVenantKirchhoffOrthoMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // print state to output stream

    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf(file, "}\n");
}

// initialize temporary state variables according to equilibrated state vars
void
StVenantKirchhoffOrthoMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
}


// Called when equilibrium reached, set equilibrated vars according to temporary (working) ones.
void
StVenantKirchhoffOrthoMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
}







} // end namespace oofem
