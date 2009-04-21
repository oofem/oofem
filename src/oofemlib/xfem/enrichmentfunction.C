#include "enrichmentfunction.h"
#include "enrichmentitem.h"
#include "domain.h"
#include <math.h>
#include "xfemmanager.h"
#include "mathfem.h"

void EnrichmentFunction :: insertEnrichmentItem(EnrichmentItem *er) {
    int sz = assocEnrItemArray.giveSize();
    assocEnrItemArray.resize(sz + 1);
    this->assocEnrItemArray.at(sz + 1) = er->giveNumber();
}

void EnrichmentFunction :: setActive(EnrichmentItem *er) {
    this->activeEnrItem = er;
}

IRResultType EnrichmentFunction :: initializeFrom(InputRecord *ir) {
    return IRRT_OK;
}

void DiscontinuousFunction :: evaluateFunctionAt(FloatArray &answer, FloatArray *point) {
    answer.resize(2);
    double dist = activeEnrItem->giveGeometry()->computeDistanceTo(point);
    answer.at(1) = answer.at(2) = sgn(dist);
}

void DiscontinuousFunction :: evaluateDerivativeAt(FloatMatrix &answer, FloatArray *point) {
    answer.resize(3, 2);
    answer.zero();
}

// to change
void BranchFunction :: evaluateFunctionAt(FloatArray &answer, FloatArray *point) {
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
}

// to change
void BranchFunction :: evaluateDerivativeAt(FloatMatrix &answer, FloatArray *point) {
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

