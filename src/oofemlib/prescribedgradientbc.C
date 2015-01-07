/*
 * prescribedgradientbc.C
 *
 *  Created on: Mar 5, 2014
 *      Author: svennine
 */

#include "prescribedgradientbc.h"
#include "dynamicinputrecord.h"
#include "set.h"
#include "feinterpol.h"
#include "element.h"

#include <cmath>

namespace oofem {
PrescribedGradientBC :: PrescribedGradientBC(int n, Domain *d) : ActiveBoundaryCondition(n, d) {}

PrescribedGradientBC :: ~PrescribedGradientBC() {}

IRResultType PrescribedGradientBC :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    ActiveBoundaryCondition :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mGradient, _IFT_PrescribedGradientBC_gradient);

    mCenterCoord.resize( mGradient.giveNumberOfColumns() );
    mCenterCoord.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, mCenterCoord, _IFT_PrescribedGradientBC_centercoords)

    return IRRT_OK;
}

void PrescribedGradientBC :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    input.setField(mGradient, _IFT_PrescribedGradientBC_gradient);
    input.setField(mCenterCoord, _IFT_PrescribedGradientBC_centercoords);
}


void PrescribedGradientBC :: setPrescribedGradientVoigt(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 3 ) { // Then 2D
        this->mGradient.resize(2, 2);
        this->mGradient.at(1, 1) = t.at(1);
        this->mGradient.at(2, 2) = t.at(2);
        this->mGradient.at(1, 2) = this->mGradient.at(2, 1) = t.at(3);
    } else if ( n == 6 ) { // Then 3D
        this->mGradient.resize(3, 3);
        this->mGradient.at(1, 1) = t.at(1);
        this->mGradient.at(2, 2) = t.at(2);
        this->mGradient.at(3, 3) = t.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->mGradient.at(1, 2) = this->mGradient.at(2, 1) = t.at(6) * 0.5;
        this->mGradient.at(1, 3) = this->mGradient.at(3, 1) = t.at(5) * 0.5;
        this->mGradient.at(2, 3) = this->mGradient.at(3, 2) = t.at(4) * 0.5;
    } else {
        OOFEM_ERROR("Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedTensor directly if needed.");
    }
}


void PrescribedGradientBC :: giveGradientVoigt(FloatArray &oGradient) const
{
    int numRows = mGradient.giveNumberOfRows();
    switch ( numRows ) {
    case 1:
        oGradient = FloatArray {
            mGradient.at(1, 1)
        };
        break;
    case 2:
        // Do not assume symmetry
        oGradient = {
            mGradient.at(1, 1), mGradient.at(2, 2), mGradient.at(1, 2), mGradient.at(2, 1)
        };
        break;
    case 3:
        OOFEM_ERROR("PrescribedGradientBC :: giveGradientVoigt() not implemented for 3 rows.\n")
        break;
    }
}

double PrescribedGradientBC :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}
} /* namespace oofem */
