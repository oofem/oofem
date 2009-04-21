#ifndef planstrssxfem_h
#define planstrssxfem_h

#include "planstrss.h"
#include "xfemelementinterface.h"
#include "usrdefsub.h"

/** temporary class for testing
 * in the usual case instead of PlaneStress2dXfem
 * there will be the standard PlaneStress2d
 */
class PlaneStress2dXfem : public PlaneStress2d, public XfemElementInterface {
public:
    /// Constructor
    PlaneStress2dXfem(int n, Domain * d) : PlaneStress2d(n, d), XfemElementInterface(this) {}
    /// Destructor
    ~PlaneStress2dXfem() {};
    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    /** function implementing the interface service
     * computes enriched part of the strain displacement matrix
     */
    void XfemElementInterface_createEnrBmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    /// computes the enriched part of the location array
    void giveLocationArray(IntArray &locationArray, EquationID);
    const char *giveClassName() const { return "PlaneStress2dXfem"; }
    classType giveClassID() const { return PlaneStress2dXfemClass; }
};

#endif
