#include "FsbInterpolatorPlaneStress.h"

namespace oofem {

double FsbInterpolatorPlaneStress :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2), inv, dn;

    this->giveDerivatives(dn, lcoords);
    double x, y, i;
    for (i = 1; i <= 4; i++ ) {
        x = cellgeo.giveVertexCoordinates(i).at(xind); // This used to be x = cellgeo.giveVertexCoordinates(i)->at(xind); but could not build. Check if it works properly with dot operator.
        y = cellgeo.giveVertexCoordinates(i).at(yind); // This used to be y = cellgeo.giveVertexCoordinates(i)->at(yind); but could not build. Check if it works properly with dot operator.

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }

    for (i = 5; i <= 7; i++) {
        x = 0.5*(cellgeo.giveVertexCoordinates(i - 3).at(xind) + cellgeo.giveVertexCoordinates(i - 4).at(xind)); // This used to be x = 0.5*(cellgeo.giveVertexCoordinates(i - 3)->at(xind) + cellgeo.giveVertexCoordinates(i - 4)->at(xind)); but could not build. Check if it works properly with dot operator.
        y = 0.5*(cellgeo.giveVertexCoordinates(i - 3).at(yind) + cellgeo.giveVertexCoordinates(i - 4).at(yind)); // This used to be y = 0.5*(cellgeo.giveVertexCoordinates(i - 3)->at(yind) + cellgeo.giveVertexCoordinates(i - 4)->at(yind)); but could not build. Check if it works properly with dot operator.

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }

    // i = 8
        x = 0.5*(cellgeo.giveVertexCoordinates(1).at(xind) + cellgeo.giveVertexCoordinates(i - 4).at(xind)); // This used to be x = 0.5*(cellgeo.giveVertexCoordinates(1).at(xind) + cellgeo.giveVertexCoordinates(i - 4).at(xind)); but could not build. Check if it works properly with dot operator.
        y = 0.5*(cellgeo.giveVertexCoordinates(1).at(yind) + cellgeo.giveVertexCoordinates(i - 4).at(yind)); // This used to be y = 0.5*(cellgeo.giveVertexCoordinates(1).at(yind) + cellgeo.giveVertexCoordinates(i - 4).at(yind)); but could not build. Check if it works properly with dot operator.

        jacobianMatrix.at(1, 1) += dn.at(8, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(8, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(8, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(8, 2) * y;

    inv.beInverseOf(jacobianMatrix);

    answer.beProductTOf(dn, inv);
    return jacobianMatrix.giveDeterminant();
}

void FsbInterpolatorPlaneStress::evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    const double &ksi = lcoords.at(1);
    const double &eta = lcoords.at(2);

    answer = {
        ( 1. + ksi ) * ( 1. + eta ) * 0.25,
        ( 1. - ksi ) * ( 1. + eta ) * 0.25,
        ( 1. - ksi ) * ( 1. - eta ) * 0.25,
        ( 1. + ksi ) * ( 1. - eta ) * 0.25,

        ( 1. - ksi*ksi) * ( 1. + eta) * 0.5,
        ( 1. - ksi) * ( 1. - eta*eta) * 0.5,
        ( 1. - ksi*ksi) * ( 1. - eta) * 0.5,
        ( 1. + ksi) * ( 1. - eta*eta) * 0.5
    };
}

void FsbInterpolatorPlaneStress :: giveDerivatives(FloatMatrix &dn, const FloatArray &lc)
{
    const double &ksi = lc[0];
    const double &eta = lc[1];

    dn.resize(8, 2);

    // dn/dxi
    dn.at(1, 1) =  0.25 * ( 1. + eta );
    dn.at(2, 1) = -0.25 * ( 1. + eta );
    dn.at(3, 1) = -0.25 * ( 1. - eta );
    dn.at(4, 1) =  0.25 * ( 1. - eta );

    dn.at(5, 1) = -1.0 * ksi * ( 1. + eta);
    dn.at(6, 1) = -0.5 * ( 1. - eta*eta);
    dn.at(7, 1) = -1.0 * ksi * ( 1. - eta);
    dn.at(8, 1) =  0.5 * ( 1. - eta*eta);

    // dn/deta
    dn.at(1, 2) =  0.25 * ( 1. + ksi );
    dn.at(2, 2) =  0.25 * ( 1. - ksi );
    dn.at(3, 2) = -0.25 * ( 1. - ksi );
    dn.at(4, 2) = -0.25 * ( 1. + ksi );

    dn.at(5, 2) =  0.5 * (1. - ksi*ksi);
    dn.at(6, 2) = -1.0 * eta * ( 1. - ksi);
    dn.at(7, 2) = -0.5 * (1. - ksi*ksi);
    dn.at(8, 2) = -1.0 * eta * ( 1. + ksi);
}

} // end namespace oofem
