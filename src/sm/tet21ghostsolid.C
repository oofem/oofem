#include "tet21ghostsolid.h"
#include "nlstructuralelement.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "structuralcrosssection.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "load.h"
#include "element.h"
#include "dofmanager.h"

namespace oofem {
REGISTER_Element(tet21ghostsolid);

FEI3dTetQuad tet21ghostsolid :: interpolation;
FEI3dTetLin tet21ghostsolid :: interpolation_lin;

IntArray tet21ghostsolid :: momentum_ordering(30);
IntArray tet21ghostsolid :: conservation_ordering(4);
IntArray tet21ghostsolid :: ghostdisplacement_ordering(30);

tet21ghostsolid::tet21ghostsolid(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{

    numberOfGaussPoints = 5;
    numberOfDofMans = 10;

    double nu=.25, E=1;
    Dghost.resize(6,6);
    Dghost.zero();
    Dghost.at(1,1) = 1-nu; Dghost.at(2,2) = Dghost.at(3,2) = nu;
    Dghost.at(2,2) = 1-nu; Dghost.at(2,1) = Dghost.at(2,3) = nu;
    Dghost.at(3,3) = 1-nu; Dghost.at(3,1) = Dghost.at(3,2) = nu;
    Dghost.at(4,4) = Dghost.at(5,5) = Dghost.at(6,6) = .5*(1-2*nu);
    Dghost.times(E/(1+nu)/(1-2*nu));

    conservation_ordering.setValues(4, 7, 14, 21, 28);

    for ( int i = 0, j = 1; i < 10; ++i ) {
        momentum_ordering(i * 3 + 0) = j++;
        momentum_ordering(i * 3 + 1) = j++;
        momentum_ordering(i * 3 + 2) = j++;
        if ( i <= 3 ) { j++; }
        j=j+3;
    }

    for ( int i = 0, j = 1; i < 10; ++i ) {
        j=j+3;
        ghostdisplacement_ordering(i * 3 + 0) = j++;
        ghostdisplacement_ordering(i * 3 + 1) = j++;
        ghostdisplacement_ordering(i * 3 + 2) = j++;
        if ( i <= 3 ) { j++; }
    }

}

FEInterpolation *
tet21ghostsolid :: giveInterpolation() const
{
    return & interpolation;
}

void
tet21ghostsolid :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
tet21ghostsolid :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    answer.resize(64, 64);
    answer.zero();

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();

    FloatMatrix Kf, G, Kx, D, B, Ed, EdB, dNx;
    FloatArray Nlin, dNv;

    for (int j = 0; j<iRule->giveNumberOfIntegrationPoints(); j++) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);

        double detJ = fabs( ( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) ) );
        double weight = gp->giveWeight();

        this->interpolation.evaldNdx( dNx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( Nlin, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        dNv.resize(30);
        for (int k = 0; k<dNx.giveNumberOfColumns(); k++) {
            dNv.at(k*3+1) = dNx.at(1,k+1);
            dNv.at(k*3+2) = dNx.at(2,k+1);
            dNv.at(k*3+3) = dNx.at(3,k+1);
        }

        if (nlGeometry == 0) {

            this->computeBmatrixAt(gp, B);

            // Fluid part
            gp->setMaterialMode(_3dFlow);
            fluidMaterial->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep);
            gp->setMaterialMode(_3dMat);

            EdB.beProductOf(Ed, B);
            Kf.plusProductSymmUpper(B, EdB, detJ*weight);

            // Ghost solid part
            EdB.beProductOf(Dghost, B);
            Kx.plusProductSymmUpper(B, EdB, detJ*weight);

            // Incompressibility part
            G.plusDyadUnsym(dNv, Nlin, -detJ*weight);

        } else {
            OOFEM_CLASS_ERROR("No support for large deformations yet!");
        }

    }

    double deltat = tStep->giveTimeIncrement();

    FloatMatrix GT, GTdeltat;

    GT.beTranspositionOf(G);
    GTdeltat.beTranspositionOf(G);
    GTdeltat.times(deltat);
    Kf.symmetrized();
    Kx.symmetrized();

    answer.assemble(Kf, ghostdisplacement_ordering, momentum_ordering);
    answer.assemble(G, ghostdisplacement_ordering, conservation_ordering);
    // answer.assemble(GT, conservation_ordering, ghostdisplacement_ordering);
    answer.assemble(GT, conservation_ordering, momentum_ordering);
    answer.assemble(Kx, momentum_ordering, ghostdisplacement_ordering);

    //answer.printYourself();

}


void
tet21ghostsolid :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.resize(64);
    answer.zero();

    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatArray N, gVector, temparray(30), dNv;
    FloatMatrix dNx, G;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interpolation.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) ) );
        double dA = detJ * gp->giveWeight();

        // Body load
        if ( gVector.giveSize() ) {

            double rho = mat->give('d', gp);
            this->interpolation.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );

            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dA;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dA;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dA;
            }
        }

        // "load" from previous step
        this->interpolation.evaldNdx( dNx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        dNv.resize(30);
        for (int k = 0; k<dNx.giveNumberOfColumns(); k++) {
            dNv.at(k*3+1) = dNx.at(1,k+1);
            dNv.at(k*3+2) = dNx.at(2,k+1);
            dNv.at(k*3+3) = dNx.at(3,k+1);
        }

        G.plusDyadUnsym(N, dNv, -dA);
    }

    answer.assemble(temparray, this->momentum_ordering);

    if (!tStep->isTheFirstStep()) {
        G.printYourself();
        FloatArray u_n;
        IntArray id;
        id.setValues(3, 1, 2, 3);

        for (int i = 0; i<this->giveNumberOfDofManagers(); i++) {
            this->giveDofManager(i+1)->giveUnknownVector(u_n, id, VM_Velocity, tStep );
            u_n.printYourself();
        }
    }

    // answer.printYourself();
}

void
tet21ghostsolid :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    answer.resize(64);
    answer.zero();
}

void
tet21ghostsolid :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the mask for node number inode of this element.

    if ( inode <= 4 ) {
        if ( ut == EID_MomentumBalance ) {
            answer.setValues(7, V_u, V_v, V_w, D_u, D_u, D_w, P_f);
        } else {
            OOFEM_ERROR("tet21ghostsolid :: giveDofManDofIDMask: Unknown equation id encountered");
        }
    } else {
        if ( ut == EID_MomentumBalance || ut == EID_MomentumBalance_ConservationEquation ) {
            answer.setValues(6, V_u, V_v, V_w, D_u, D_u, D_w);
        } else {
            OOFEM_ERROR("tet21ghostsolid :: giveDofManDofIDMask: Unknown equation id encountered");
        }
    }
}

void
tet21ghostsolid :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x30] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(6, 30);
    answer.zero();

    for ( int i = 1; i <= 10; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(1, i);
        answer.at(2, 3 * i - 1) = dnx.at(2, i);
        answer.at(3, 3 * i - 0) = dnx.at(3, i);

        answer.at(4, 3 * i - 1) = dnx.at(3, i);
        answer.at(4, 3 * i - 0) = dnx.at(2, i);

        answer.at(5, 3 * i - 2) = dnx.at(3, i);
        answer.at(5, 3 * i - 0) = dnx.at(1, i);

        answer.at(6, 3 * i - 2) = dnx.at(2, i);
        answer.at(6, 3 * i - 1) = dnx.at(1, i);
    }
}

}
