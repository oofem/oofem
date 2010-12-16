/* $Header: /home/cvs/bp/oofem/oofemlib/src/primaryfield.h,v 1.1 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#ifndef dofdistributedprimaryfield_h
#define dofdistributedprimaryfield_h

#include "primaryfield.h"

#include "flotarry.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {
class PrimaryField;
class Dof;

/**
 * Class representing field of primary variables, which are typically allocated on nodes.
 * The field is determined by DOF values stored in DOF repositories (unknown dictionary).
 * These repositories are maintained and updated by engng models since the algorithms are very specific to each model.
 * The class can return several variables stored in DOF. The purpose of this class is to provide
 * a shell that allows to access these repositories using field services.
 * The class contains also a solution vector for temporal storage of unknowns. The vector needs to be projected back to DOFs.
 */
class DofDistributedPrimaryField : public PrimaryField
{
public:


protected:
public:
    /** Constructor. Creates a field of given type associated to given domain.
     * Not using pointer to domain, because this will prevent the use of PrimaryField as an
     * EngngModel attribute. This is because the domain does not exists when
     * PrimaryField is created (this is when EngngModel is created).
     */
    DofDistributedPrimaryField(EngngModel *a, int idomain, FieldType ft, EquationID ut, int nHist);
    ~DofDistributedPrimaryField();

    /** Copy unknowns from previous solution or DOF's dictionary to the solution vector
     * @param mode what the unknown desribes (increment, total value etc.)
     * @param atTime time of interest
     * @param answer the resulting vector
     */
    virtual void initialize(ValueModeType mode, TimeStep *atTime, FloatArray &answer);

    /** Return value of interest at given DOF
     * @param dof pointer to DOF
     * @param mode what the unknown desribes (increment, total value etc.)
     * @param atTime time of interest
     */
    virtual double giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *atTime);

    /** Project @param vectorToStore back to DOF's dictionary
     * @param mode what the unknown desribes (increment, total value etc.)
     * @param atTime time
     * @param vectorToStore vector with the size of number of equations
     */
    virtual void update(ValueModeType mode, TimeStep *atTime, FloatArray &vectorToStore);

    virtual FloatArray *giveSolutionVector(TimeStep *atTime);
    /**
     */
    virtual void advanceSolution(TimeStep *atTime);

    /** Stores receiver state to output stream.
     * Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream output stream
     * @param mode determines ammount of info in stream
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    saveContext(DataStream *stream, ContextMode mode);
    /** Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @see saveContext member function.
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    restoreContext(DataStream *stream, ContextMode mode);

protected:
};
} // end namespace oofem
#endif // dofdistributedprimaryfield_h
