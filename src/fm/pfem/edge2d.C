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

#include "edge2d.h"
#include "delaunaytriangle.h"

namespace oofem {

	Edge2D :: Edge2D (int node1, int node2)
		: nodes(2)
	{
		nodes.at(1) = node1;
		nodes.at(2) = node2;
	}

	Edge2D :: ~Edge2D()
	{}

	bool
	Edge2D::operator== (const Edge2D& right)
	{
		if ((this->nodes.at(1) == right.nodes.at(1) && this->nodes.at(2) == right.nodes.at(2)) || 
			(this->nodes.at(1) == right.nodes.at(2) && this->nodes.at(2) == right.nodes.at(1)))
		{
			return true;
		}
		else
		{
			return false;
		}
	}



	AlphaEdge2D::AlphaEdge2D(int node1, int node2, double _length)
		: Edge2D(node1, node2)
		, isOnConvexHull(true)
		, outerAlphaBound(0.0)
		, innerAlphaBound(0.0)
		, length(_length)
	{
		sharedByTriangles[0] = NULL;
		sharedByTriangles[1] = NULL;
	}

	AlphaEdge2D::~AlphaEdge2D()
	{}

	void
	AlphaEdge2D :: setSharing(int n, DelaunayTriangle *pTE)
	{
		sharedByTriangles[n-1] = pTE;
	}
}