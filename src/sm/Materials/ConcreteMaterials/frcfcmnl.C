/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "frcfcmnl.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "datastream.h"
#include "node.h"

#include "nonlocalmaterialext.h"

#include "sm/Elements/structuralelement.h"

#include "floatmatrix.h"
#include "floatarray.h"

#include "dynamicinputrecord.h"
#include "unknownnumberingscheme.h"
#include "stressvector.h"
#include "strainvector.h"


namespace oofem {
REGISTER_Material(FRCFCMNL);

FRCFCMNL :: FRCFCMNL(int n, Domain *d) : FRCFCM(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{}


void
FRCFCMNL :: initializeFrom(InputRecord &ir)
{
    FRCFCM :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    //    IR_GIVE_FIELD(ir, participAngle, _IFT_FRCFCMNL_participAngle);
    participAngle = 90.;
}


void
FRCFCMNL :: giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *tStep) const
{
    // first try conventional way as if no sourrounding cracks existed
    FCMMaterial :: giveRealStressVector(answer, gp, totalStrain, tStep);

    FRCFCMNLStatus *status = static_cast< FRCFCMNLStatus * >( this->giveStatus(gp) );
    int numberOfActiveCracks = status->giveNumberOfTempCracks();

    double sigma_f_nonlocal = 0.;

    FloatArray sigma, crackVec;
    FloatMatrix sigmaL2G, sigmaG2L, crackDirs;
    double sigma_f_local = 0.;
    double crackStrain;

    if ( numberOfActiveCracks == 0 ) {
        sigma_f_nonlocal = this->computeNonlocalStressInFibersInUncracked(gp, tStep);
        sigma_f_nonlocal = max( sigma_f_nonlocal, status->giveFiberStressNL(1) );
        status->setTempFiberStressNL(1, sigma_f_nonlocal);

        return;
    } else {
        sigma = answer; // global stress


        sigmaL2G = status->giveL2GStressVectorTransformationMtrx();


        sigmaG2L = status->giveG2LStressVectorTransformationMtrx();

        sigma.rotatedWith(sigmaG2L, 'n'); // global stress transformed to local stress

        sigma_f_local = 0.;

        for ( int iCrack = 1; iCrack <= ( status->giveNumberOfTempCracks() ); iCrack++ ) {
            // get cracking strain for i-th crack
            crackStrain = status->giveTempCrackStrain(iCrack);

            if ( crackStrain > 0. ) {
                // get local fiber stress
	            sigma_f_local = this->computeStressInFibersInCracked(gp, tStep, crackStrain, iCrack);
                // set local fiber stress to status
                status->setTempFiberStressLoc(iCrack, sigma_f_local);
            } else {
                // set zero fiber stress
                status->setTempFiberStressLoc(iCrack, 0.);
            }

            // and compare it to stress from surroundings sigma_f_nonlocal

            crackDirs = status->giveCrackDirs();

            crackVec.resize( status->giveMaxNumberOfCracks(gp) );
            crackVec.zero();

            for ( int i = 1; i <= status->giveMaxNumberOfCracks(gp); i++ ) {
                crackVec.at(i) = crackDirs.at(i, iCrack);
            }

            // compute nonlocal stress
            sigma_f_nonlocal = this->computeNonlocalStressInFibers(crackVec, gp, tStep);

            //test !!!!
            sigma_f_nonlocal = max( sigma_f_nonlocal, status->giveFiberStressNL(iCrack) );


            // set nonlocal fiber stress to status
            status->setTempFiberStressNL(iCrack, sigma_f_nonlocal);

            sigma.at(iCrack) += sigma_f_nonlocal;
        } // loop over crack directions

        sigma.rotatedWith(sigmaL2G, 'n'); // modified local stress transformed to global stress

        answer = sigma;
        status->letTempStressVectorBe(sigma);
    }

    return;
}




/// computes debonded length of the fibers from the crack opening. Delta is here one half of the crack opening.
double
FRCFCMNL :: computeDebondedLength(double delta) const {
    double a = 0.;

    if ( this->fiberType == FT_CAF ) { // continuous aligned fibers
        a = sqrt( ( this->Ef * this->Df * delta ) / ( 2. * this->tau_0 * ( 1. + this->eta ) ) );
    } else if (  this->fiberType == FT_SAF ) { // short aligned fibers
        // no snubbing here - debonding length is simply related to pull-out displacement and not to force...

        a = sqrt( ( this->Ef * this->Df * delta ) / ( 2. * this->tau_0 * ( 1. + this->eta ) ) );
        // threshold is related to the fibre length
        a = min( a, this->Lf / ( 2. * ( 1. + this->eta ) ) );
    } else if (  this->fiberType == FT_SRF ) { // short random fibers
        a = sqrt( ( this->Ef * this->Df * delta ) / ( 2. * this->tau_0 * ( 1. + this->eta ) ) );
        // threshold is related to the fibre length
        a = min( a, this->Lf / ( 2. * ( 1. + this->eta ) ) );
    } else {
        OOFEM_ERROR("Unknown fiber type");
    }

    return a;
}


double
FRCFCMNL :: computeDecreaseInFibreStress(double distance, double delta, double targetDebondedLength) const {
    // compute decrease in fiber stress due to stress transfer to matrix
    // delta_sigma_f = 4 * tau_eff * Vf * Xca'/ ( Df * ( 1 - Vf))
    // Xca' = ( 4 * x^2 - 4 * x * Lf ) / ( -2 * pi * Lf * lambda )
    // lambda = 4 / ( pi * g )
    // where x is the distance from the crack, in this case the distance between the gauss-points


    double delta_sigma_f = 0.;
    double tau, x;

    if ( this->fiberType == FT_CAF ) { // continuous aligned fibers
        delta_sigma_f = this->Vf * 4. * this->tau_0 * distance / this->Df;
    } else {
        if ( delta < this->w_star / 2. ) {
            tau = this->tau_0;
        } else {
            // maybe delta should be changed to delta_max?!?
            tau = this->computeFiberBond(2. * delta);
        }

        x = min(distance, targetDebondedLength);

        if (  this->fiberType == FT_SAF ) { // short aligned fibers
            // delta_sigma_f is computed here in direction of the fibers, therefore no snubbing and no reduction in fibre volume
            // this expression is derived for x < a (we are in the debonding zone) and for fibers perpendicular to the crack
            delta_sigma_f = this->Vf * 4. * tau * ( this->Lf * x - x * x ) / ( this->Df * this->Lf );
        } else {
            // this expression uses the same hypothesis as in Lu & Leung 2016 (taken from Aveston 1974) that the stress transfer by friction
            // does not continue after x (section where we want to now the stress) = xd (length of the debonded zone)
            delta_sigma_f = this->Vf * 4. * tau * ( this->Lf * x - x * x ) / ( 3 * this->Df * this->Lf );
        }
    }

    return delta_sigma_f;
}


void
FRCFCMNL :: computeElementCentroid(FloatArray &answer, GaussPoint *gp) const {
    double A, cx, cy, x_i, x_ii, y_i, y_ii;
    Element *elem;
    elem = gp->giveElement();
    int nnode = elem->giveNumberOfNodes();

    answer.resize(2);
    answer.zero();

    A = 0.;
    cx = cy = 0.;

    if ( ( nnode == 3 ) || ( nnode == 4 ) ) {
        for ( int inod = 1; inod < nnode; inod++ ) {
            x_i = elem->giveNode(inod)->giveCoordinate(1);
            x_ii = elem->giveNode(inod + 1)->giveCoordinate(1);

            y_i = elem->giveNode(inod)->giveCoordinate(2);
            y_ii = elem->giveNode(inod + 1)->giveCoordinate(2);

            A += 0.5 * ( x_i * y_ii - x_ii * y_i );

            cx += ( x_i + x_ii ) * ( x_i * y_ii - x_ii * y_i );
            cy += ( y_i + y_ii ) * ( x_i * y_ii - x_ii * y_i );
        }

        x_i = elem->giveNode(nnode)->giveCoordinate(1);
        x_ii = elem->giveNode(1)->giveCoordinate(1);

        y_i = elem->giveNode(nnode)->giveCoordinate(2);
        y_ii = elem->giveNode(1)->giveCoordinate(2);

        A += 0.5 * ( x_i * y_ii - x_ii * y_i );

        cx += ( x_i + x_ii ) * ( x_i * y_ii - x_ii * y_i );
        cy += ( y_i + y_ii ) * ( x_i * y_ii - x_ii * y_i );
    }

    cx /= ( 6. * A );
    cy /= ( 6. * A );

    answer.at(1) = cx;
    answer.at(2) = cy;
}


bool
FRCFCMNL :: isInElementProjection(GaussPoint *homeGp, GaussPoint *nearGp, int iNlCrack) const
{
    FRCFCMNLStatus *nonlocStatus;
    nonlocStatus = static_cast< FRCFCMNLStatus * >( nearGp->giveMaterialStatus() );

    Element *nearElem;
    FloatArray examinedCoords;
    int nnode;
    double x, y, b_min, b_max, b_home, slope, x_min, x_max;

    nearElem = nearGp->giveElement();

    if ( ( fiberType == FT_CAF ) || ( fiberType == FT_SAF ) ) {
        slope = this->orientationVector.at(2) / this->orientationVector.at(1);
    } else {
        slope = nonlocStatus->giveCrackDirs().at(2, iNlCrack) / nonlocStatus->giveCrackDirs().at(1, iNlCrack);
    }

    this->computeElementCentroid(examinedCoords, homeGp);

    nnode = nearElem->giveNumberOfNodes();

    if ( fabs(slope) < 1.e6 ) {
        b_min = b_max = b_home = 0.;

        for ( int inod = 1; inod <= nnode; inod++ ) {
            x = nearElem->giveNode(inod)->giveCoordinate(1);
            y = nearElem->giveNode(inod)->giveCoordinate(2);

            // y = slope * x + b

            if ( inod == 1 ) {
                b_min = y - slope * x;
                b_max = y - slope * x;
            } else {
                b_min = min(b_min, y - slope * x);
                b_max = max(b_max, y - slope * x);
            }
        }

        b_home = examinedCoords.at(2) - slope *examinedCoords.at(1);

        if ( ( b_min < b_home ) && ( b_home < b_max ) ) {
            return true;
        } else {
            return false;
        }
    } else {
        x_min = x_max = 0.;

        for ( int inod = 1; inod < nnode; inod++ ) {
            x = nearElem->giveNode(inod)->giveCoordinate(1);

            if ( inod == 1 ) {
                x_min = x;
                x_max = x;
            } else {
                x_min = min(x_min, x);
                x_max = max(x_max, x);
            }
        }

        if ( ( x_min < examinedCoords.at(1) ) && ( examinedCoords.at(1) < x_max ) ) {
            return true;
        } else {
            return false;
        }
    }

    // happy compiler
    return true;
}



double
FRCFCMNL :: computeNonlocalStressInFibers(const FloatArray &crackVectorHome, GaussPoint *gp, TimeStep *tStep) const {
    // using "non-local" approach

    // computes stress in fibers in the direction of "crackVectorHome"
    // this can be used either to check the stress-strength criterion or when computing the stress from strain

    // they key is to find the maximum stress in fibers which
    // - diminishes with the distance from crack
    // - diminishes with increasing angle wrt to its normal


    FRCFCMNLStatus *nonlocStatus;

    this->buildNonlocalPointTable(gp);

    auto *list = this->giveIPIntegrationList(gp); // !

    // two fiber stresses - from left and right part of the crack - rewrite if it turns out that one stress is fully sufficientaveraged
    double sigma_f_left = 0.;
    double sigma_f_right = 0.;

    FloatArray coordsHome, coordsTarget;

    double alpha, k_alpha_Home;

    FloatArray crackVectorTarget;
    FloatArray pointVector;

    bool side;

    double distance, targetCrackOpening, delta, a, targetDebondedLength, theta;
    double sigma_f_iCr, delta_sigma_f;


    this->computeElementCentroid(coordsHome, gp);

    crackVectorTarget.resize( gp->giveMaterial()->giveDomain()->giveNumberOfSpatialDimensions() );
    pointVector = crackVectorTarget;

    for ( auto &lir : *list ) {
        GaussPoint *neargp = lir.nearGp;

        if ( neargp->giveElement()->giveNumber() != gp->giveElement()->giveNumber() ) { // it must not be the same element
            if ( neargp->giveMaterial()->giveNumber() == gp->giveMaterial()->giveNumber() ) { // it must be the same material
                nonlocStatus = static_cast< FRCFCMNLStatus * >( neargp->giveMaterialStatus() );

                // crack(s) exist in the near gp
                if ( nonlocStatus->giveNumberOfTempCracks() > 0 ) {
                    // do just once for all cracks
                    // get coordinates

                    this->computeElementCentroid(coordsTarget, neargp);

                    // compute the distance between evaluated element centers
                    distance = 0.;
                    for ( int i = 1; i <= coordsHome.giveSize(); i++ ) {
                        distance += ( coordsHome.at(i) - coordsTarget.at(i) ) * ( coordsHome.at(i) - coordsTarget.at(i) );
                    }

                    distance = sqrt(distance);

                    // zero the orientation vector between two points
                    pointVector.zero();


                    for ( int iNlCrack = 1; iNlCrack <= nonlocStatus->giveNumberOfTempCracks(); iNlCrack++ ) {
                        // length of the tunnel (zone where the fibers are debonded from the matrix and where friction acts) at the target gauss-point
                        targetCrackOpening = this->computeNormalCrackOpening(neargp, iNlCrack);


                        // pull-out displacement
                        delta = ( targetCrackOpening - fibreActivationOpening ) / 2.;

                        if ( delta > 0. ) {
                            a = this->computeDebondedLength(delta);

                            targetDebondedLength = a;

                            if ( distance < targetDebondedLength ) {
                                if ( this->isInElementProjection(gp, neargp, iNlCrack) ) {
                                    // get stress in fibers bridging iNlCrack - local stress
                                    sigma_f_iCr = nonlocStatus->giveTempFiberStressLoc(iNlCrack);


                                    // get vector from points' coordinates - if has not been done before for the previous crack
                                    if ( pointVector.containsOnlyZeroes() ) {
                                        for ( int i = 1; i <= coordsHome.giveSize(); i++ ) {
                                            pointVector.at(i) = ( coordsHome.at(i) - coordsTarget.at(i) ) / distance;
                                        }
                                    }

                                    // the change in the bridging stress due to pulley and snubbing effects
                                    if ( ( this->fiberType == FT_CAF ) ||  ( this->fiberType == FT_SAF ) ) {
                                        theta = this->computeCrackFibreAngle(neargp, iNlCrack);
                                        sigma_f_iCr /= ( fabs( cos(theta) ) * exp(fabs(theta) * this->f) );
                                    } else if ( this->fiberType == FT_SRF ) {
                                        sigma_f_iCr *= 2. / ( 3. * this->g );
                                    } else {
                                        OOFEM_ERROR("Unknown fiber type");
                                    }

                                    // change in fiber stress due to bond friction
                                    delta_sigma_f = this->computeDecreaseInFibreStress(distance, delta, targetDebondedLength);

                                    sigma_f_iCr -= delta_sigma_f;

                                    // only positive values - can even happen to be negative?
                                    sigma_f_iCr = max(sigma_f_iCr, 0.);

                                    // decide the side
                                    alpha = this->computeAngleBetweenVectors(pointVector, crackVectorHome);

                                    if ( alpha < M_PI / 2. ) {
                                        side = true;
                                    } else {
                                        side = false;
                                    }

                                    if ( this->fiberType == FT_CAF ) {
                                        alpha = this->computeAngleBetweenVectors(crackVectorHome, this->orientationVector);
                                    } else if ( this->fiberType == FT_SAF ) {
                                        alpha = this->computeAngleBetweenVectors(crackVectorHome, this->orientationVector);
                                    } else if ( this->fiberType == FT_SRF ) {
                                        // angle factor - point vector vs. target crack
                                        crackVectorTarget.zero();
                                        for ( int i = 1; i <= coordsHome.giveSize(); i++ ) {
                                            crackVectorTarget.at(i) = nonlocStatus->giveCrackDirs().at(i, iNlCrack);
                                        }

                                        alpha = this->computeAngleBetweenVectors(crackVectorHome, crackVectorTarget);
                                    } else {
                                        OOFEM_ERROR("Unknown fiber type");
                                    }

                                    // necessary correction
                                    if ( alpha > M_PI / 2. ) {
                                        alpha = M_PI - alpha;
                                    }

                                    // add some snubbing effect?
                                    if ( alpha >  participAngle * M_PI / 180. ) {
                                        k_alpha_Home = 0.;
                                    } else {
                                        k_alpha_Home = cos( alpha * (  90. / participAngle ) );
                                    }


                                    // get fiber stress
                                    if ( side ) {
                                        sigma_f_left = max(k_alpha_Home * sigma_f_iCr, sigma_f_left);
                                    } else {
                                        sigma_f_right = max(k_alpha_Home * sigma_f_iCr, sigma_f_right);
                                    }
                                } // is in element projection
                            } // is in the debonded zone
                        } // positive delta
                    } // loop over target crack directions
                } // condition for at least crack
            } // the same material material
        } // not the same element
    } // loop over all gauss-points


    return ( max(sigma_f_left, sigma_f_right) );
    //  return ( sigma_f_left + sigma_f_right );
}


double
FRCFCMNL :: computeNonlocalStressInFibersInUncracked(GaussPoint *gp, TimeStep *tStep) const
{
    // using "non-local" approach

    FRCFCMNLStatus *nonlocStatus;

    this->buildNonlocalPointTable(gp);

    auto *list = this->giveIPIntegrationList(gp); // !

    FloatArray coordsHome, coordsTarget;

    double distance, targetCrackOpening, delta, a, targetDebondedLength, theta, sigma_f_iCr, delta_sigma_f;
    double sigma_f = 0.;

    this->computeElementCentroid(coordsHome, gp);


    for ( auto &lir : *list ) {
        GaussPoint *neargp = lir.nearGp;

        if ( neargp->giveElement()->giveNumber() != gp->giveElement()->giveNumber() ) { // it must not be the same element
            if ( neargp->giveMaterial()->giveNumber() == gp->giveMaterial()->giveNumber() ) { // it must be the same material
                nonlocStatus = static_cast< FRCFCMNLStatus * >( neargp->giveMaterialStatus() );

                // crack(s) exist in the near gp
                if ( nonlocStatus->giveNumberOfTempCracks() > 0 ) {
                    // do just once for all cracks
                    // get coordinates

                    this->computeElementCentroid(coordsTarget, neargp);

                    // compute the distance between evaluated element centers
                    distance = 0.;
                    for ( int i = 1; i <= coordsHome.giveSize(); i++ ) {
                        distance += ( coordsHome.at(i) - coordsTarget.at(i) ) * ( coordsHome.at(i) - coordsTarget.at(i) );
                    }

                    distance = sqrt(distance);


                    for ( int iNlCrack = 1; iNlCrack <= nonlocStatus->giveNumberOfTempCracks(); iNlCrack++ ) {
                        // length of the tunnel (zone where the fibers are debonded from the matrix and where friction acts) at the target gauss-point
                        targetCrackOpening = this->computeNormalCrackOpening(neargp, iNlCrack);


                        // pull-out displacement
                        delta = ( targetCrackOpening - fibreActivationOpening ) / 2.;

                        if ( delta > 0. ) {
                            a = this->computeDebondedLength(delta);

                            targetDebondedLength = a;

                            if ( distance < targetDebondedLength ) {
                                if ( this->isInElementProjection(gp, neargp, iNlCrack) ) {
                                    // get stress in fibers bridging iNlCrack - local stress
                                    sigma_f_iCr = nonlocStatus->giveTempFiberStressLoc(iNlCrack);

                                    // the change in the bridging stress due to pulley and snubbing effects
                                    if ( ( this->fiberType == FT_CAF ) ||  ( this->fiberType == FT_SAF ) ) {
                                        theta = this->computeCrackFibreAngle(neargp, iNlCrack);
                                        sigma_f_iCr /= ( fabs( cos(theta) ) * exp(fabs(theta) * this->f) );
                                    } else if ( this->fiberType == FT_SRF ) {
                                        sigma_f_iCr /= this->g;
                                    } else {
                                        OOFEM_ERROR("Unknown fiber type");
                                    }

                                    // change in fiber stress due to bond friction
                                    delta_sigma_f = this->computeDecreaseInFibreStress(distance, delta, targetDebondedLength);

                                    sigma_f_iCr -= delta_sigma_f;

                                    // only positive values - can even happen to be negative?
                                    sigma_f_iCr = max(sigma_f_iCr, 0.);

                                    // get fiber stress
                                    sigma_f = max(sigma_f_iCr, sigma_f);
                                } // is in element projection
                            } // is in the debonded zone
                        } // positive delta
                    } // loop over target crack directions
                } // condition for at least crack
            } // the same material material
        } // not the same element
    } // loop over all gauss-points


    return sigma_f;
}





bool
FRCFCMNL :: isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress) const {
    // evaluates cheaper function than is the nonlocal approach
    if ( !FRCFCM :: isStrengthExceeded(base, gp, tStep, iCrack, trialStress) ) {
        return false;
    }

    FRCFCMNLStatus *status = static_cast< FRCFCMNLStatus * >( this->giveStatus(gp) );

    double sigma_f, sigma_m;

    FloatArray crackVec;

    crackVec.resize( status->giveMaxNumberOfCracks(gp) );
    crackVec.zero();

    for ( int i = 1; i <= status->giveMaxNumberOfCracks(gp); i++ ) {
        crackVec.at(i) = base.at(i, iCrack);
    }

    sigma_f = this->computeNonlocalStressInFibers(crackVec, gp, tStep);

    //test !!!!
    sigma_f = max( sigma_f, status->giveFiberStressNL(iCrack) );

    status->setTempFiberStressNL(iCrack, sigma_f);

    // get neat stress in matrix
    sigma_m = ( trialStress - sigma_f ) / ( 1. - this->Vf );

    // compare matrix strength with computed stress in the matrix

    if ( FCMMaterial :: isStrengthExceeded(base, gp, tStep, iCrack, sigma_m) ) {
        return true;
    } else {
        return false;
    }
}



void
FRCFCMNL :: giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode,
                                        GaussPoint *gp,
                                        TimeStep *tStep) const
//
// returns effective material stiffness matrix in full form
// for gp stress strain mode
//
{
    if (  rMode == ElasticStiffness ) {
        FCMMaterial :: giveMaterialStiffnessMatrix(answer, rMode, gp, tStep);
    } else if ( rMode == SecantStiffness ) {
        /*FRCFCMNLStatus *status = static_cast< FRCFCMNLStatus * >( this->giveStatus(gp) );
         * int numberOfActiveCracks = status->giveNumberOfTempCracks();
         * double localSigmaF;
         * double NLsigmaF;
         *
         * for (int iCrack = 1; iCrack <= numberOfActiveCracks; iCrack++) {
         *
         * localSigmaF = status->giveTempFiberStressLoc(iCrack);
         * NLsigmaF =  status->giveTempFiberStressNL(iCrack);
         *
         * //      if ( NLsigmaF > localSigmaF ) {
         * if ( NLsigmaF > 0. ) {
         *  FCMMaterial :: giveMaterialStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
         *  return;
         * }
         * }*/

        FCMMaterial :: giveMaterialStiffnessMatrix(answer, rMode, gp, tStep);
    } else { // tangent stiffness
        FCMMaterial :: giveMaterialStiffnessMatrix(answer, rMode, gp, tStep);
    }

    return;
}


double
FRCFCMNL :: computeAngleBetweenVectors(const FloatArray &vec1, const FloatArray &vec2) const
{
    // compute angle between two vectors

    if ( vec1.giveSize() != vec2.giveSize() ) {
        OOFEM_ERROR("the length of vec1 and vec2 is not of the same");
    }


    //return 0.;


    double length1 = 0., length2 = 0.;
    double theta = 0.;

    for ( int i = 1; i <= min( vec1.giveSize(), vec2.giveSize() ); i++ ) {
        length1 += pow(vec1.at(i), 2);
        length2 += pow(vec2.at(i), 2);
        theta += vec1.at(i) * vec2.at(i);
    }
    length1 = sqrt(length1);
    length2 = sqrt(length2);

    // normalize
    theta /= ( length1 * length2 );


    // can be exceeded due to truncation error
    if ( theta > 1 ) {
        return 0.;
    } else if ( theta < -1. ) {
        return M_PI;
    }


    return acos(theta);
}


Interface *
FRCFCMNL :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return NULL;
    }
}



int
FRCFCMNL :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FRCFCMNLStatus *status = static_cast< FRCFCMNLStatus * >( this->giveStatus(gp) );

    double local = 0.;
    double nl = 0.;

    for ( int i = 1; i <= status->giveMaxNumberOfCracks(gp); i++ ) {
        local = max(status->giveFiberStressLoc(i), local);
        nl = max(status->giveFiberStressNL(i), nl);
    }


    // stress in fibers (per continuum)
    if  ( type == IST_FiberStressLocal ) {
        answer.resize(1);
        answer.at(1) = local;
        return 1;

        // stress in fibers (per continuum)
    } else if  ( type == IST_FiberStressNL ) {
        answer.resize(1);
        answer.at(1) = nl;
        return 1;
    } else {
        return FRCFCM :: giveIPValue(answer, gp, type, tStep);
    }
}



///////////////////////////////////////////////////////////////////
//                      FRC FCM NL STATUS                       ///
///////////////////////////////////////////////////////////////////


FRCFCMNLStatus :: FRCFCMNLStatus(GaussPoint *gp) :
    FRCFCMStatus(gp), StructuralNonlocalMaterialStatusExtensionInterface(),
    fiberStressLoc(this->nMaxCracks), tempFiberStressLoc(), fiberStressNL(), tempFiberStressNL()
{
    tempFiberStressLoc = fiberStressLoc;
    fiberStressNL = fiberStressLoc;
    tempFiberStressNL = fiberStressLoc;
}


void
FRCFCMNLStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    FRCFCMStatus :: printOutputAt(file, tStep);

    fprintf(file, "maxFiberStressLocal: {");
    for ( double s: fiberStressLoc ) {
        fprintf( file, " %f", s );
    }
    fprintf(file, "}\n");

    fprintf(file, "maxFiberStressNL: {");
    for ( double s: fiberStressLoc ) {
        fprintf( file, " %f", s );
    }
    fprintf(file, "}\n");
}


void
FRCFCMNLStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    FRCFCMStatus :: initTempStatus();

    this->tempFiberStressLoc = this->fiberStressLoc;
    this->tempFiberStressNL = this->fiberStressNL;
}



void
FRCFCMNLStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    FRCFCMStatus :: updateYourself(tStep);

    this->fiberStressLoc = this->tempFiberStressLoc;
    this->fiberStressNL = this->tempFiberStressNL;
}



void
FRCFCMNLStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    FRCFCMStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = fiberStressLoc.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = fiberStressNL.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}

void
FRCFCMNLStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    FRCFCMStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = fiberStressLoc.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = fiberStressNL.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


Interface *
FRCFCMNLStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return ConcreteFCMStatus :: giveInterface(type);
    }
}
} // end namespace oofem
