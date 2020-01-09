#include "FsbInterpolatorPlate.h"

namespace oofem {

void FsbInterpolatorPlate::evalN(oofem::FloatArray &answer, const oofem::FloatArray &lcoords, const oofem::FEICellGeometry &cellgeo)
{
    const double &ksi = lcoords.at(1);
    const double &eta = lcoords.at(2);

    double n5 = 0.5*(1 - ksi*ksi)*(1 - eta);
    double n7 = 0.5*(1 - ksi*ksi)*(1 + eta);
    double n6 = 0.5*(1 - eta*eta)*(1 + ksi);
    double n8 = 0.5*(1 - eta*eta)*(1 - ksi);

    answer = {
        0.25*(1 - ksi)*(1 - eta) - 0.5*(n8 + n5),
        0.25*(1 + ksi)*(1 - eta) - 0.5*(n5 + n6),
        0.25*(1 + ksi)*(1 + eta) - 0.5*(n6 + n7),
        0.25*(1 - ksi)*(1 + eta) - 0.5*(n7 + n8),

        n5,
        n6,
        n7,
        n8
    };
}

double FsbInterpolatorPlate::evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(2, 2), inv, dn;

    this->giveDerivatives(dn, lcoords);
    double x, y, i;
    for (i = 1; i <= 4; i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(xind);
        y = cellgeo.giveVertexCoordinates(i)->at(yind);

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }

    for (i = 5; i <= 7; i++) {
        x = 0.5*(cellgeo.giveVertexCoordinates(i - 3)->at(xind) + cellgeo.giveVertexCoordinates(i - 4)->at(xind));
        y = 0.5*(cellgeo.giveVertexCoordinates(i - 3)->at(yind) + cellgeo.giveVertexCoordinates(i - 4)->at(yind));

        jacobianMatrix.at(1, 1) += dn.at(i, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(i, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(i, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(i, 2) * y;
    }

    // i = 8
        x = 0.5*(cellgeo.giveVertexCoordinates(1)->at(xind) + cellgeo.giveVertexCoordinates(i - 4)->at(xind));
        y = 0.5*(cellgeo.giveVertexCoordinates(1)->at(yind) + cellgeo.giveVertexCoordinates(i - 4)->at(yind));

        jacobianMatrix.at(1, 1) += dn.at(8, 1) * x;
        jacobianMatrix.at(1, 2) += dn.at(8, 1) * y;
        jacobianMatrix.at(2, 1) += dn.at(8, 2) * x;
        jacobianMatrix.at(2, 2) += dn.at(8, 2) * y;

    inv.beInverseOf(jacobianMatrix);

    answer.beProductTOf(dn, inv);
    return jacobianMatrix.giveDeterminant();
}

void FsbInterpolatorPlate::giveDerivatives(FloatMatrix &dn, const FloatArray &lc)
{
    const double &ksi = lc[0];
    const double &eta = lc[1];

    dn.resize(8, 2);

    // dn/dxi
    dn.at(1, 1) = -0.25*(1 - eta) + 0.25*(1 - eta*eta) + 0.5*ksi*(1 - eta);
    dn.at(2, 1) =  0.25*(1 - eta) + 0.5*ksi*(1 - eta) - 0.25*(1 - eta*eta);
    dn.at(3, 1) =  0.25*(1 + eta) - 0.25*(1 - eta*eta) + 0.5*ksi*(1 + eta);
    dn.at(4, 1) = -0.25*(1 + eta) + 0.5*ksi*(1 + eta) + 0.25*(1 - eta*eta);

    dn.at(5, 1) = -ksi*(1 - eta);
    dn.at(6, 1) =  0.5*(1 - eta*eta);
    dn.at(7, 1) = -ksi*(1 + eta);
    dn.at(8, 1) = -0.5*(1 - eta*eta);

    // dn/deta
    dn.at(1, 2) = -0.25*(1 - ksi) + 0.5*eta*(1 - ksi) + 0.25*(1 - ksi*ksi);
    dn.at(2, 2) = -0.25*(1 + ksi) + 0.25*(1 - ksi*ksi) - 0.25*(1 - ksi*ksi);
    dn.at(3, 2) =  0.25*(1 + ksi) - 0.25*(1 - ksi*ksi) + 0.5*eta*(1 + ksi);
    dn.at(4, 2) =  0.25*(1 + ksi) + 0.5*eta*(1 + ksi) + 0.5*eta*(1 - ksi);

    dn.at(5, 2) = -0.5*(1 - ksi*ksi);
    dn.at(6, 2) =  0.5*(1 - ksi*ksi);
    dn.at(7, 2) = -eta*(1 + ksi);
    dn.at(8, 2) = -eta*(1 - ksi);
}

} // end namespace oofem
