#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "geometry.h"

PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, Patch *patch) : GaussIntegrationRule(n, e) {
    this->patch = patch;
}

PatchIntegrationRule :: ~PatchIntegrationRule() {
    delete patch;
}


int
PatchIntegrationRule :: SetUpPointsOnTriagle(int nPoints, MaterialMode mode, GaussPoint ***arry)
{
  numberOfIntegrationPoints = GaussIntegrationRule::SetUpPointsOnTriagle (nPoints, mode, arry);
  firstLocalStrainIndx = 1;
  lastLocalStrainIndx = 3;
  // convert patch coordinates into element based, update weights accordingly
  for ( int j = 0; j <  numberOfIntegrationPoints; j++ ) {
      GaussPoint *gp = (*arry)[j];
      patch->convertGPIntoParental(gp); // convert coordinates into parental
      Triangle * tr = (Triangle*) patch;
      // here to get the area of the whole element
      gp->setWeight(gp->giveWeight()/tr->getArea()); // update integration weight 
  } 
  return numberOfIntegrationPoints;
}
