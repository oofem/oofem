#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"
#include "dofdistributedprimaryfield.h"
#include "unknownnumberingscheme.h"
#include "integral.h"
#include "fei2dquadlin.h"
#include "fei2dlinelin.h"
#include "fei2dlinequad.h"
#include "termlibrary.h"
#include "crosssection.h"
#include "material.h"
#include "masterdof.h"
#include "gaussintegrationrule.h"
#include "prototype2.h"

namespace oofem {
    REGISTER_EngngModel(TestProblem);    
    REGISTER_Element(Q1Element)
    REGISTER_Element(L1Element)
    REGISTER_Term(BTSigmaTerm2)
    REGISTER_Term(NTfTerm)
    const FEInterpolation & Q1Element::gInterpol = FEI2dQuadLin(1,2);
    const FEInterpolation & L1Element::gInterpol = FEI2dLineLin(1,2);
    const FEInterpolation & L1Element::gQuadInterpol = FEI2dLineQuad(1,2);

    void L1Element::initializeFinish()
    {
        MPElement::initializeFinish();
        int n = this->dofManArray.giveSize();
        if ( n != 2 && n != 3 ) {
            OOFEM_ERROR("element %d: 2 nodes (ends) or 3 nodes (ends and inner node) expected, %d given", this->giveNumber(), n);
        }
        // a given 3rd node is the inner node of a quadratic interpolation (which generates it at
        // the midpoint when not given, or given as 0), and makes the geometry quadratic
        if ( !this->dofManArray.contains(0) ) {
            this->numberOfDofMans = n;
            this->quadraticGeometry = ( n == 3 );
        }
    }

    void Q1Element::initializeFinish()
    {
        MPElement::initializeFinish();
        int n = this->dofManArray.giveSize();
        if ( n != 4 && n != 8 ) {
            OOFEM_ERROR("element %d: 4 nodes (corners) or 8 nodes (corners and edge inner nodes) expected, %d given", this->giveNumber(), n);
        }
        // nodes 1-4 define the geometry; given edge nodes 5-8 are the inner nodes of a quadratic
        // interpolation, which generates those not given (or given as 0)
        if ( !this->dofManArray.contains(0) ) {
            this->numberOfDofMans = n;
        }
    }

} // end namespace oofem
