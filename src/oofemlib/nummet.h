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

#ifndef nummet_h
#define nummet_h

#include "oofemcfg.h"
#include "inputrecord.h"
#include "contextioresulttype.h"
#include "contextmode.h"

namespace oofem {
class EngngModel;
class Domain;
class DataStream;

/**
 * This base class is an abstraction for numerical algorithm.
 *
 * Generally, the particular instances (instances of derived
 * classes) perform some sequence of numerical operations on given data
 * to obtain the solution.
 * The derived class should declare the interface for specific problem type
 * (like solution of linear system). The interface usually consist in
 * declaring virtual abstract function solve, with parameters corresponding
 * to problem under consideration. The solve method should return
 * value of NM_Status type. (Other parameters can be provided via instanciateFrom
 * service, which receives the init record of Engng method).
 * The data are specified using parameters passed to solve method (so called
 * mapping). Typically, each particular Engng model instance is responsible
 * for mapping of its governing equation components to corresponding
 * numerical components. Such mapping allows the numerical method implementation
 * to be independent of a particular physical problem by strictly dealing
 * with numerical components, which are mapped to corresponding physical
 * components of governing equation, that are  hidden to numerical method.
 *
 * It should be pointed out, that all numerical methods solving the same
 * numerical problem use the same compulsory names for the corresponding numerical components - this is
 * enforced by using the same base problem-specific class.
 * It is therefore possible to use any suitable
 * instance of the NumericalMethod class to the solve problem, and leave the whole engineering model code,
 * including mapping, unchanged, because all instances of the Numerical
 * method class provide the common interface.
 * Similarly, a high-level numerical method instance may use services of another
 * low-level numerical method instance. The numerical method
 * instance may also represent interface to an existing procedure
 * written in C or Fortran.
 */
class OOFEM_EXPORT NumericalMethod
{
protected:
    /// Pointer to domain
    Domain *domain;
    /// Pointer to engineering model.
    EngngModel *engngModel;

public:
    /**
     * Constructor.
     * @param d Domain which the receiver belongs to.
     * @param m Engineering model which the receiver belongs to.
     */
    NumericalMethod(Domain * d, EngngModel * m) : domain(d), engngModel(m) { }
    /// Destructor
    virtual ~NumericalMethod() { }

    /// @return Engineering model receiver is connected to.
    EngngModel *giveEngngModel() { return engngModel; }

    virtual IRResultType initializeFrom(InputRecord *ir)  { return IRRT_OK; }

    /**
     * Reinitializes the receiver. This is used, when topology of problem has changed
     * (for example after adaptive refinement or load transfer in parallel applications).
     * This is necessary for numerical methods, that cache some data between solution
     * steps and that may depend on domain or problem topology. The caching of data by
     * receiver is intended only for speeding up the calculation, but numerical method
     * must be always able to generate this data again.
     * This method clears receiver cached data dependent on topology, when it changes.
     */
    virtual void reinitialize() { }

    virtual void setDomain(Domain *d) { domain = d; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL) { return CIO_OK; }
};
} // end namespace oofem
#endif // nummet_h
