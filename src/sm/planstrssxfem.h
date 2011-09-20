#ifndef planstrssxfem_h
#define planstrssxfem_h

#include "planstrss.h"
#include "xfemelementinterface.h"

namespace oofem {
/** temporary class for testing
 * in the usual case instead of PlaneStress2dXfem
 * there will be the standard PlaneStress2d
 */
class PlaneStress2dXfem : public PlaneStress2d, public XfemElementInterface
{
public:
    /// Constructor
    PlaneStress2dXfem(int n, Domain *d) : PlaneStress2d(n, d), XfemElementInterface(this) { }
    /// Destructor
    ~PlaneStress2dXfem() { };
    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    /// computes the enriched part of the location array
    void giveLocationArray(IntArray & locationArray, EquationID, const UnknownNumberingScheme & s) const;
    const char *giveClassName() const { return "PlaneStress2dXfem"; }
    classType giveClassID() const { return PlaneStress2dXfemClass; }
    int computeNumberOfDofs(EquationID ut);
    void computeBmatrixAt(GaussPoint *, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    void giveDofManDofIDMask(int inode, EquationID, IntArray & answer) const;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
    void computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    double computeArea() const;

#ifdef __OOFEG
    void          drawRawGeometry(oofegGraphicContext &);
    //void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    virtual void  drawScalar(oofegGraphicContext &context);
    //virtual void  drawSpecial(oofegGraphicContext &);
    //     void          drawInternalState (oofegGraphicContext&);
#endif
};
} // end namespace oofem
#endif
