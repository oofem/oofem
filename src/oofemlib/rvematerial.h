/*
 * RVEMaterial.h
 *
 *  Created on: Mar 26, 2010
 *      Author: carl
 */

#ifndef RVEMATERIAL_H_
#define RVEMATERIAL_H_

#include <stdio.h>
#include <stdlib.h>

#include "material.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "oofem_limits.h"


namespace oofem {

class RVEMaterial //: virtual public Material
{
private:
	int stdoutFID;
	fpos_t stdoutPos;

protected:
	/** name of .in file containing the RVE
	 */
	std::string rveFilename; //char rveFilename [ MAX_FILENAME_LENGTH ];
    std::string rveLogFilename;
	/** Type of boundary condition.
	 */
	int BCType;

public:

	EngngModel *rve;

	// Constructor
    RVEMaterial(int n, Domain *d) { };// : Material(n, d) { };

	// Destructor
	~RVEMaterial(){ free (rve); };

	int SupressRVEoutput;

	IRResultType initializeFrom(InputRecord *ir);

	//MaterialStatus *CreateStatus(GaussPoint *gp) const;

    	/// Returns class name of the receiver.
	void suppressStdout();
	void enableStdout();

   	const char *giveClassName() const { return "RVEMaterial"; };
    	/// Returns classType id of receiver.
    classType giveClassID()         const { return MaterialClass; };

};

}

#endif /* RVEMATERIAL_H_ */
