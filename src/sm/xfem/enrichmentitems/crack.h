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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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
#ifndef CRACK_H_
#define CRACK_H_

#include "xfem/enrichmentitem.h"
#include "xfem/hybridei.h"

#define _IFT_Crack_Name "crack"

namespace oofem {
class XfemManager;
class Domain;
class InputRecord;
class GaussPoint;
class GnuplotExportModule;

/**
 * Crack.
 * @author Erik Svenning
 * @date Feb 14, 2014
 */
class OOFEM_EXPORT Crack : public HybridEI
{
public:
    Crack(int n, XfemManager *xm, Domain *aDomain);

    virtual const char *giveClassName() const { return "Crack"; }
    virtual const char *giveInputRecordName() const { return _IFT_Crack_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    void AppendCohesiveZoneGaussPoint(GaussPoint *ipGP);

    virtual void callGnuplotExportModule(GnuplotExportModule &iExpMod, TimeStep *tStep);

    const std :: vector< GaussPoint * > &giveCohesiveZoneGaussPoints() const { return mCohesiveZoneGaussPoints; }
    const std :: vector< double > &giveCohesiveZoneArcPositions() const { return mCohesiveZoneArcPositions; }

    void computeCrackIntersectionPoints(Crack &iCrack, std :: vector< FloatArray > &oIntersectionPoints, std :: vector< double > &oArcPositions);
    void computeArcPoints(const std :: vector< FloatArray > &iIntersectionPoints, std :: vector< double > &oArcPositions);
    double computeLength();
    virtual int giveDofPoolSize() const;

protected:
    /**
     * Array of pointers to the Gauss points related to the
     * cohesive zone. The array is used for data extraction
     * and visualization only. The reason for keeping an array
     * of pointers here is as follows: the cohesive zone Gauss
     * points are created in the XFEMElementInterface, that
     * (of course) only keeps track of GPs in that element.
     * However, for visualization it is very valuable to be able
     * to plot cohesive zone data (e.g. damage or crack opening)
     * vs the arc length coordinate of the crack. This must be
     * accomplished at the level of the EnrichmentItem, because
     * here we know about the geometry of the crack.
     */
    std :: vector< GaussPoint * >mCohesiveZoneGaussPoints;
    std :: vector< double >mCohesiveZoneArcPositions;
};
} // end namespace oofem

#endif /* CRACK_H_ */
