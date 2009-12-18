#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "geometry.h"

namespace oofem {

PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, Patch *patch) : GaussIntegrationRule(n, e) {
    this->patch = patch;
}

PatchIntegrationRule :: ~PatchIntegrationRule() {
    delete patch;
}


int
PatchIntegrationRule :: SetUpPointsOnTriagle(int nPoints, MaterialMode mode, GaussPoint ***arry)
{
    numberOfIntegrationPoints = GaussIntegrationRule :: SetUpPointsOnTriagle(nPoints, mode, arry);
    firstLocalStrainIndx = 1;
    lastLocalStrainIndx = 3;
    // convert patch coordinates into element based, update weights accordingly
    for ( int j = 0; j <  numberOfIntegrationPoints; j++ ) {
        GaussPoint *gp = ( * arry ) [ j ];
        patch->convertGPIntoParental(gp); // convert coordinates into parental
        ElementGeometry *elg = ( ElementGeometry * ) patch->giveParent();
        double parentArea = elg->giveArea();
        Triangle *tr = ( Triangle * ) patch;
        gp->setWeight(8.0 * gp->giveWeight() * tr->getArea() / parentArea); // update integration weight
    }

    return numberOfIntegrationPoints;
}

} // end namespace oofem
