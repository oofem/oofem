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

#include "aabb.h"
#include <stdexcept>

namespace oofem {

Vector :: Vector() : Vector(0.0, 0.0, 0.0)
{
}

Vector :: Vector(double x, double y, double z) : x(x), y(y), z(z)
{
}

double
Vector :: operator[](int index) const {
	if (index == 0) return this->x;
	if (index == 1) return this->y;
	if (index == 2) return this->z;
	throw std::out_of_range("Index out of range for Vector");
}

AABB :: AABB(const Vector& min, const Vector& max)
{
	this->min = min;
	this->max = max;
}

bool
AABB :: contains(const Vector& v)
{
	return (
		this->min.x <= v.x
		&&
		v.x <= this->max.x
		&&
		this->min.y <= v.y
		&&
		v.y <= this->max.y
		&&
		this->min.z <= v.z
		&&
		v.z <= this->max.z
	);
}

bool
AABB :: contains(double x, double y, double z)
{
	return this->contains(Vector(x, y, z));
}

void
AABB :: merge(double x, double y, double z)
{
	this->merge(Vector(x, y, z));
}

void
AABB :: merge(const Vector& v)
{
	if (v.x < this->min.x) this->min.x = v.x;
	if (v.y < this->min.y) this->min.y = v.y;
	if (v.z < this->min.z) this->min.z = v.z;
	if (v.x > this->max.x) this->max.x = v.x;
	if (v.y > this->max.y) this->max.y = v.y;
	if (v.z > this->max.z) this->max.z = v.z;
}

} // end namespace oofem
