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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef WEAKPERIODICBC_H_
#define WEAKPERIODICBC_H_

#include <vector>
#include <iostream>

#include "activebc.h"
#include "inputrecord.h"

namespace oofem {
enum basisType { monomial=0, trigonometric=1 };
/**
 * Imposes weak periodicity on the doftype of choice. 2D. It is required that the two edges are parallel and either horizontal or vertical.
 *
 * @author Carl Sandström
 */
class WeakPeriodicbc : public ActiveBoundaryCondition
{
private:

    basisType useBasisType;
    int bcID;
    int orderOfPolygon;

    int direction;

    int normalDirection;

    double smax, smin;
    bool doUpdateSminmax;
    /** Number of Gausspoints used when integrating along the element edges */
    int ngp;


    /** ID of dofs on which weak periodicity is imposed */
    int dofid;

    /** sideSign is the sign of the normal for each side */
    signed int sideSign [ 2 ];

    /** side[] keeps track of which side of the triangle is located along the boundary. element[] keeps track of what element is located along the boundary */
    std :: vector< int >side [ 2 ], element [ 2 ];

    void giveEdgeNormal(FloatArray &answer, int element, int side);

    void updateSminmax();

    void updateDirection();

    double computeBaseFunctionValue(int baseID, double coordinate);

    Node *gammaDman;
    IntArray DofIDList;

public:
    WeakPeriodicbc(int n, Domain *d);
    ~WeakPeriodicbc() { };

    IRResultType initializeFrom(InputRecord *ir);

    basisType giveBasisType() {return useBasisType; };

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain);

    virtual double assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain);

    virtual int giveNumberOfInternalDofManagers();

    virtual DofManager *giveInternalDofManager(int i);

    void addElementSide(int elem, int side);
};
}
#endif /* WEAKPERIODICBC_H_ */
