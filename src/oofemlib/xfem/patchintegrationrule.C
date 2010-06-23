#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "geometry.h"
#include "usrdefsub.h"
#include "contextioerr.h"
#include "datastream.h"

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
        Element *elg = ( Element * ) patch->giveParent();
        double parentArea = elg->computeVolume();
        Triangle *tr = ( Triangle * ) patch;
        gp->setWeight(8.0 * gp->giveWeight() * tr->getArea() / parentArea); // update integration weight
    }

    return numberOfIntegrationPoints;
}

contextIOResultType
PatchIntegrationRule :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
  //
  // saves full  context (saves state variables, that completely describe
  // current state)
  //
  
  // save parent data
  contextIOResultType iores;
  
  if ( ( iores = IntegrationRule :: saveContext(stream, mode, obj) ) != CIO_OK ) {
    THROW_CIOERR(iores);
  }
  
  // save patch data
  if (this->patch) {
    // store patch type
    int _type = this->patch->givePatchType();
    if ( !stream->write(& _type, 1) ) {
      THROW_CIOERR(CIO_IOERR);
    }
    patch->saveContext (stream, mode, obj);
  } else {
    OOFEM_ERROR("saveContex : can't store NULL patch");
  }

  return CIO_OK;  
}

contextIOResultType
PatchIntegrationRule :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    //
    // restores full element context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;

    if ( stream == NULL ) {
        OOFEM_ERROR("restoreContex : can't write into NULL stream");
    }

    if ( ( iores = IntegrationRule :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
      THROW_CIOERR(iores);
    }
    
    // restore patch data
    if (this->patch) {
      delete this->patch;
    }
    int _ptype;
    if ( !stream->read(& _ptype, 1) ) {
      THROW_CIOERR(CIO_IOERR);
    }
    
    // create new patch
    this->patch = CreateUsrDefPatch ((Patch::PatchType)_ptype, this->giveElement());
    this->patch->restoreContext (stream, mode, obj);
    
    return CIO_OK;  
}



} // end namespace oofem

