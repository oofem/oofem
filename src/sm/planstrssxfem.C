#include "planstrssxfem.h"
#include "engngm.h"
#include "dof.h"
#include "masterdof.h"
#include "planstrss.h"
#include "structuralmaterial.h"
#include "patchintegrationrule.h"
#include <math.h>
#include "interfacetype.h"
#include "xfemelementinterface.h"
#include "structuralcrosssection.h"

Interface *
PlaneStress2dXfem::giveInterface(InterfaceType interface) {
    if (interface != XfemElementInterfaceType)  return PlaneStress2d::giveInterface(interface);
    else if (interface == XfemElementInterfaceType) {
        return (XfemElementInterface *) this;
    } else return NULL;
}


void PlaneStress2dXfem::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                   int li, int ui){
    
    FloatMatrix* simple = new FloatMatrix();
    PlaneStress2d::computeBmatrixAt(gp, *simple);    
    int start = simple->giveNumberOfColumns();
    // evaluation of N,dNdx
    FloatMatrix dNdx;
    interpolation.evaldNdx(dNdx, domain, dofManArray, * gp->giveCoordinates(), 0.0);
    FloatArray N;
    interpolation.evalN(N, * gp->giveCoordinates(), 0.0);
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    int counter = 0;
    AList<FloatMatrix> additionals;
    additionals.put(1,simple);
    counter += simple->giveNumberOfColumns();
    for (int i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++) {
        EnrichmentItem *er = xf->giveEnrichmentItem(i);
        int erndofs = er->giveNumberOfDofs();
        // enrichment function at the gauss point
        FloatArray efgp;
        er->giveEnrichmentFunction()->evaluateFunctionAt(efgp, gp);
        // derivative of enrichment function at the gauss point
        FloatMatrix efgpD;
        er->giveEnrichmentFunction()->evaluateDerivativeAt(efgpD, gp);
        // adds up the number of the dofs from an enrichment item
        // for each node
        for (int j = 1; j <= this->giveNumberOfDofManagers(); j++) {
            if (er->isDofManEnriched(dofManArray.at(j))) {   
                FloatArray *nodecoords = domain->giveDofManager(dofManArray.at(j))->giveCoordinates();
                // ef is a FloatArray containing the value of EnrichmentFunction in a specific for all enriched dofs
                FloatArray efnode;
                er->giveEnrichmentFunction()->evaluateFunctionAt(efnode, nodecoords);
                // matrix to be added anytime a node is enriched
                FloatMatrix* toAdd = new FloatMatrix(3, erndofs);
                toAdd->zero();
                FloatArray help;
                help.resize(2);
                for (int p = 1; p <= 2; p++) {
                    help.at(p) = dNdx.at(j, p) * (efgp.at(p) - efnode.at(p)) + N.at(j) * efgpD.at(p, p);
                }
                for (int k = 1; k <= erndofs; k++) {
                    toAdd->at(k, k) = help.at(k);
                    if (k == 1) {
                        toAdd->at(3, k) = help.at(2);
                    }

                    if (k == 2) {
                        toAdd->at(3, k) = help.at(1);
                    }
                }
                int sz = additionals.giveSize();
                additionals.put(sz+1,toAdd);
                counter += toAdd->giveNumberOfColumns();
            }
        }
    }
    answer.resize(3, counter);
    answer.zero();
    int columns = 1;
    for(int i = 1; i <= additionals.giveSize(); i++){    
        for(int j = 1; j <= additionals.at(i)->giveNumberOfColumns(); j++){
            for(int k = 1;  k <= 3; k++){
                answer.at(k, columns) = additionals.at(i)->at(k,j);
            }
            columns++;
        }
    }
} 

void PlaneStress2dXfem::giveLocationArray(IntArray &locationArray, EquationID, const UnknownNumberingScheme& s) const {
    IntArray interactedEI;
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    xf->getInteractedEI(interactedEI, const_cast<PlaneStress2dXfem *> (this));
    int count = 0;
    // for all enrichment items which element interacts
    for (int i = 1; i <= (const_cast<PlaneStress2dXfem *> (this))->giveNumberOfDofManagers(); i++) {
        DofManager *dm = this->giveDomain()->giveDofManager(dofManArray.at(i));
        for(int j = 1; j <= xf->giveNumberOfEnrichmentItems(); j++){
           EnrichmentItem *er = xf->giveEnrichmentItem(j);
           if(er->isDofManEnriched(dofManArray.at(i))){
              IntArray *dofIdAr = er->getDofIdArray();
              for (int k = 1; k <= dofIdAr->giveSize(); k++) {
                if(dm->hasDofID(dofIdAr->at(k)) == false){
                   int sz = dm->giveNumberOfDofs();
                   Dof *df = new MasterDof( sz + 1, dm, 0, 0, dofIdAr->at(k));
                   int eqN = xf->giveFictPosition(dofManArray.at(i))->at(k);
                   df->setEquationNumber(eqN);  
                   dm->addDof(sz+1,df);                  
                }
              }
           }
        }
    }
    locationArray.resize(0);
    IntArray enriched;
    enriched.resize(0);
    for(int i = 1; i <= (const_cast<PlaneStress2dXfem *> (this))->giveNumberOfDofManagers(); i++){
       DofManager *dm = (const_cast<PlaneStress2dXfem *> (this))->giveDomain()->giveDofManager(dofManArray.at(i));
       for(int j = 1; j <= dm->giveNumberOfDofs(); j++){
		int eqN = dm->giveDof(j)->giveEquationNumber(s);
       		if(j <= 2) locationArray.followedBy(eqN);
       		else enriched.followedBy(eqN);
       }
    }
    locationArray.followedBy(enriched);
} 


int PlaneStress2dXfem::computeNumberOfDofs(EquationID ut) {
   int ret = 0;
   for(int i = 1; i <= giveNumberOfDofManagers(); i++){
       DofManager *dm = this->giveDomain()->giveDofManager(dofManArray.at(i));
       for(int j = 1; j <= dm->giveNumberOfDofs(); j++){
           ret++;
       }
   }
   return ret;
} 


void
PlaneStress2dXfem::giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    PlaneStress2d::giveDofManDofIDMask(inode, EID_MomentumBalance, answer);
    XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    for (int i = 1; i <= xf->giveNumberOfEnrichmentItems(); i++) {
        EnrichmentItem *er = xf->giveEnrichmentItem(i);
        if (er->isDofManEnriched(dofManArray.at(inode))) {
                IntArray *dofIdAr = er->getDofIdArray();
                for (int j = 1; j <= dofIdAr->giveSize(); j++) {
                    answer.followedBy(dofIdAr->at(j));
                }
            }
    }
    return;
} 

void PlaneStress2dXfem::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep){
     XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
     if(xf->isInteracted(this)){
       PatchIntegrationRule *pir = (PatchIntegrationRule*) gp->giveIntegrationRule();
       StructuralMaterial * sm = (StructuralMaterial*) this->giveDomain()->giveMaterial(pir->giveMaterial());
       sm->giveCharacteristicMatrix(answer, ReducedForm, rMode, gp, tStep);
     }
     else PlaneStress2d::computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}


double
PlaneStress2dXfem :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double volume = 0;
    volume = PlaneStress2d::computeVolumeAround(aGaussPoint);
    return volume;
}

void
PlaneStress2dXfem :: computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the receiver's nodes (in nodal cs).
// Dofs cointaining expected unknowns (of expected type) are determined
// using this->GiveNodeDofIDMask function
{
    // FloatArray *answer ;
    int i, j, k, nDofs, size;
    IntArray elementNodeMask;
    FloatArray vec;
    answer.resize( size = this->computeGlobalNumberOfDofs(type) );
    k      = 0;
    int m = 0;
    FloatArray p1 (size - 8);
    
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, type, elementNodeMask);
        this->giveDofManager(i)->giveUnknownVector(vec, elementNodeMask, type, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            if(j <= 2){
               answer.at(++k) = vec.at(j);
            }
            else {
            p1.at(++m) = vec.at(j);
            }
        }
    }
    for(int i = 1; i <= m; i++){
       answer.at(k+i) = p1.at(i);
    }
    return;
}


void
PlaneStress2dXfem :: computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
{
     FloatArray Epsilon;
     this->computeStrainVector(Epsilon, gp, stepN);
     XfemManager *xf = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
     if(xf->isInteracted(this)){
           PatchIntegrationRule *pir = (PatchIntegrationRule*) gp->giveIntegrationRule();
           StructuralMaterial * sm = (StructuralMaterial*) this->giveDomain()->giveMaterial(pir->giveMaterial());
           sm->giveRealStressVector(answer, ReducedForm, gp, Epsilon, stepN);
     }
     else {
          StructuralCrossSection *cs = ( StructuralCrossSection * ) this->giveCrossSection();
          cs->giveRealStresses(answer, ReducedForm, gp, Epsilon, stepN);
     }
}

