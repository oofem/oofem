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

} // end namespace oofem
