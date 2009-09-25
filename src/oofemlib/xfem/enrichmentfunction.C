#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "domain.h"
#include <math.h>
#include "xfemmanager.h"
#include "mathfem.h"

/*
void EnrichmentFunction :: insertEnrichmentItem(EnrichmentItem *er) {
    int sz = assocEnrItemArray.giveSize();
    assocEnrItemArray.resize(sz + 1);
    this->assocEnrItemArray.at(sz + 1) = er->giveNumber();
}

void EnrichmentFunction :: setActive(EnrichmentItem *er) {
    this->activeEnrItem = er;
}
*/
IRResultType EnrichmentFunction :: initializeFrom(InputRecord *ir) {
    return IRRT_OK;
}

double EnrichmentFunction::evaluateFunctionAt(GaussPoint *gp, EnrichmentItem* ei) {
  FloatArray gcoords;
  gp->giveElement()->computeGlobalCoordinates (gcoords, *gp->giveCoordinates());
  return this->evaluateFunctionAt(&gcoords, ei);
}

void EnrichmentFunction::evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem* ei) {
  FloatArray gc;
  gp->giveElement()->computeGlobalCoordinates (gc, *gp->giveCoordinates());
  this->evaluateDerivativeAt(answer, &gc, ei);
}  


double DiscontinuousFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem* ei) {
    double dist = ei->giveGeometry()->computeDistanceTo(point);
    return sgn(dist);
}

void DiscontinuousFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem* ei) {
    answer.resize(2);
    answer.zero();
}

double RampFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem* ei) {
    double dist = ei->giveGeometry()->computeDistanceTo(point);
    double absDist;
    if (dist < 0.0000) absDist = (-1)*dist;
    else absDist = dist;
    return absDist;
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem* ei) {
    double dist = ei->giveGeometry()->computeDistanceTo(point);
    double absDist;
    if (dist < 0.0000) absDist = (-1)*dist;
    else absDist = dist;
    answer.resize(2);
    answer.zero();
    answer.at(1) = answer.at(2) = absDist / dist;
}

double RampFunction :: evaluateFunctionAt(GaussPoint *gp, EnrichmentItem* ei){
  FloatArray N;
  Element *el = gp->giveElement();
  el->giveInterpolation()->evalN(N, * gp->giveCoordinates(), 0.0);
  double dist = 0;
  double absMember = 0;
  double member = 0;
  for(int i = 1; i <= el->giveNumberOfDofManagers(); i++){
      dist = ei->giveGeometry()->computeDistanceTo(el->giveDofManager(i)->giveCoordinates());
      member += N.at(i)*dist;
  }
  if (member < 0.0000) absMember = (-1)*member;
  else absMember = member;
  return absMember;
}

void RampFunction :: evaluateDerivativeAt(FloatArray &answer, GaussPoint *gp, EnrichmentItem* ei){
  FloatArray N;
  Element *el = gp->giveElement();
  el->giveInterpolation()->evalN(N, * gp->giveCoordinates(), 0.0);
  IntArray dofManArray(el->giveNumberOfDofManagers());
  for(int i = 1; i <= el->giveNumberOfDofManagers(); i++){
      dofManArray.at(i) = el->giveDofManagerNumber(i);
  }
  FloatMatrix dNdx;
  el->giveInterpolation()->evaldNdx(dNdx, el->giveDomain(), dofManArray, * gp->giveCoordinates(), 0.0);
  double dist = 0;
  double dfdx = 0;
  double dfdy = 0;
  double phi = 0;
  for(int i = 1; i <= el->giveNumberOfDofManagers(); i++){
      dist = ei->giveGeometry()->computeDistanceTo(el->giveDofManager(i)->giveCoordinates());
      phi += N.at(i)*dist;
      dfdx += dNdx.at(i,1)*dist;
      dfdy += dNdx.at(i,2)*dist;
  }
  answer.resize(2);
  answer.zero();
  answer.at(1) = dfdx*sgn(phi);
  answer.at(2) = dfdy*sgn(phi);
}

// to change
double BranchFunction :: evaluateFunctionAt(FloatArray *point, EnrichmentItem* ei) {
    /*
     *  double ret = 0;
     *  CrackTip *cr = (CrackTip*) activeEnrItem;
     *  FloatArray transfCoor;
     *  Line *l = (Line*) cr->giveGeometry();
     *  l->transformIntoPolar(point, transfCoor);
     *  double rS = sqrt(transfCoor.at(1));
     *  double theta = transfCoor.at(2);
     *  switch (number) {
     *      case 1: ret = rS * sin(0.5 * theta);
     *      case 2: ret = rS * cos(0.5 * theta);
     *      case 3: ret = rS * sin(0.5 * theta) * cos(theta);
     *      case 4: ret = rS * cos(0.5 * theta) * cos(theta);
     *  }
     *  return ret;
     */
  return 0.0;
}

// to change
void BranchFunction :: evaluateDerivativeAt(FloatArray &answer, FloatArray *point, EnrichmentItem* ei) {
    /*
     *  answer.resize(2);
     *  CrackTip *cr = (CrackTip*) activeEnrItem;
     *  Line *l = (Line*) cr->giveGeometry();
     *  FloatArray transfCoor;
     *  l->transformIntoPolar(point, transfCoor);
     *  double alpha = l->computeInclinationAngle();
     *  double fac = 0.5 / sqrt(transfCoor.at(1));
     *  double theta = transfCoor.at(2);
     *  double dPhidr = 0;
     *  double dPhidt = 0;
     *  double drdx = cos(alpha);
     *  double drdy = sin(alpha);
     *  double dtdx = (-1) * sin(alpha);
     *  double dtdy = cos(alpha);
     *
     *  switch (number) {
     *      case 1:
     *          dPhidr = (-1) * fac * sin(0.5 * theta);
     *          dPhidt = fac * cos(0.5 * theta);
     *          break;
     *      case 2:
     *          dPhidr = fac * cos(0.5 * theta);
     *          dPhidt = fac * sin(0.5 * theta);
     *          break;
     *      case 3:
     *          dPhidr = fac * sin(0.5 * theta)*(2 * sin(theta) * sin(theta) - cos(theta));
     *          dPhidt = fac * cos(theta) * cos(1.5 * theta);
     *          break;
     *      case 4:
     *          dPhidr = fac * cos(0.5 * theta)*(cos(theta) + 2 * sin(theta) * sin(theta));
     *          dPhidt = (-1) * fac * cos(theta) * sin(1.5 * theta);
     *          break;
     *  }
     *  answer.at(1) = dPhidr * drdx + dPhidt*dtdx;
     *  answer.at(2) = dPhidr * drdy + dPhidt*dtdy;
     */
}

