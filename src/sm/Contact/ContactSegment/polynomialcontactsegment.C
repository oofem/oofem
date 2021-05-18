
#include "polynomialcontactsegment.h"
#include "mathfem.h"

namespace oofem {
REGISTER_ContactSegment(PolynomialContactSegment)

void PolynomialContactSegment::initializeFrom(InputRecord &ir)
{
    FunctionContactSegment::initializeFrom(ir);
    IR_GIVE_FIELD(ir, coeffs, _IFT_PolynomialContactSegment_coeffs);
    order = coeffs.giveSize() - 1;
}

double PolynomialContactSegment::functionValue(const double x) const
{
    double answer = 0.0;
    for ( int i = order; i >= 0; i-- ) {
        answer += coeffs(order - i) * pow(x, i);
    }
    return answer;
}

double PolynomialContactSegment::derivativeValue(const double x) const
{
    double answer = 0.0;
    for ( int i = order; i >= 1; i-- ) {
        answer += i * coeffs(order - i) * pow(x, i - 1);
    }
    return answer;
}

double PolynomialContactSegment::doubleDerivativeValue(const double x) const
{
    double answer = 0.0;
    for ( int i = order; i >= 2; i-- ) {
        answer += ( i - 1 ) * i * coeffs(order - i) * pow(x, i - 2);
    }
    return answer;
}
}
