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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef Shell7BaseXFEM_h
#define Shell7BaseXFEM_h


//#include "eleminterpmapperinterface.h"
//#include "nodalaveragingrecoverymodel.h"
//#include "layeredcrosssection.h"

//#include "nlstructuralelement.h"
#include "shell7base.h"
#include "xfemelementinterface.h"
namespace oofem {

class FEI3dTrQuad;
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */

class Shell7BaseXFEM : public Shell7Base, public XfemElementInterface
{
protected:
    int numberOfGaussPoints;	

    
    

    
public:
    //Shell7BaseXFEM(int n, Domain *d);	// constructor
    Shell7BaseXFEM(int n, Domain *d) : Shell7Base(n, d), XfemElementInterface(this){};
    virtual ~Shell7BaseXFEM() { }		// destructor -> declaring as virtual will make each subclass call their respective destr.
    // definition & identification
    virtual const char *giveClassName()                const { return "Shell7BaseXFEM"; }
    virtual classType giveClassID()                    const { return Shell7BaseXFEMClass; }
    
	


};



} // end namespace oofem
#endif 
