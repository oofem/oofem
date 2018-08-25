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

#ifndef eigenvectorprimaryfield_h
#define eigenvectorprimaryfield_h

#include "dofdistributedprimaryfield.h"

namespace oofem {
/**
 * Class representing the mode shapes of eigen vectors.
 * The values are stored as VM_Total only.
 * Active vector is determined by the time step number.
 */
class OOFEM_EXPORT EigenVectorPrimaryField : public DofDistributedPrimaryField
{
public:
    /**
     * Constructor. Creates a field of given type associated to given domain.
     * @param a Engineering model which field belongs to.
     * @param idomain Index of domain for field.
     * @param ft Type of stored field.
     * @param nHist Number of old time steps to store (minimum 1), i.e. the number of eigen vectors.
     */
    EigenVectorPrimaryField(EngngModel * a, int idomain, FieldType ft, int nHist);
    virtual ~EigenVectorPrimaryField();

    double giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep) override;

    /**
     * Stores all the eigenvectors in one call.
     * @param eigenVectors Matrix with all eigen vectors (stored as columns)
     * @param s Equation numbering for the rows of the vectors.
     */
    void updateAll(const FloatMatrix &eigenVectors, const UnknownNumberingScheme &s);

    void applyDefaultInitialCondition() override;
    void advanceSolution(TimeStep *tStep) override;
};
} // end namespace oofem
#endif // eigenvectorprimaryfield_h
