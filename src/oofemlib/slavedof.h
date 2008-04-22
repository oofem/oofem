//
//  termitovo - zaloha
//

//   ************************
//   *** CLASS SLAVE DOF ***
//   ************************

#ifndef slavedof_h
#define slavedof_h

#include "dof.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif


/**
 * Class representing "slave" degree of freedom. This dof is generally linked to some master dofs
 * with various weights - contributions - combination (link slave-slave is allowed).
 * The slave dof is fully dependent upon master dofs, therefore it's equation number is undefined.
 * It can have own boundary conditions - it is not yet implemented.
 */
class SlaveDof : public Dof
{
protected:
    /// count of Master DofManagers
    int countOfMasterDofs;
    /// count of Master DofManagers
    int countOfPrimaryMasterDofs;
    /// array of pointers on Master Dofs
    Dof **masterDof;
    /// vector of master contribution coefficients
    FloatArray masterContribution;

public:
    /**
     * Constructor. Creates slave dof with number n, belonging to aNode dof manager.
     * @param i dof number
     * @param aNode receiver will belong to aNode dof manager
     * @param id DofID of slave dof
     */
    SlaveDof(int n, DofManager *aNode, DofID id);
    /**
     * Destructor.
     */
    ~SlaveDof(void) { }

    void initialize(int cntOfMstrDfMngr, Node **mstrNode, const IntArray *mstrDofID, const FloatArray *mstrContribution);
    int giveNumberOfPrimaryMasterDofs(void);
    void giveUnknowns(FloatArray &masterUnknowns, EquationID type, ValueModeType mode, TimeStep *stepN);
    void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    void giveBcValues(FloatArray &masterBcValues, ValueModeType mode, TimeStep *stepN);
    void computeDofTransformation(FloatArray &masterContribs);
    void giveEquationNumbers(IntArray &masterEqNumbers);
    void givePrescribedEquationNumbers(IntArray &masterEqNumbers);

    double giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN);
    double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Slave simply asks necessary master dofs and
     * computes the results.
     * @see MasterDof::giveUnknown function
     */
    double giveLocalUnknown(EquationID, ValueModeType, TimeStep *) { _error("HangingDof :: giveLocalUnknown: local coordinate system doesn't exist");
                                                                     return 0.0; }

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int giveEquationNumber(void) { _error("giveEquationNumber: undefined");
                                   return 0; }

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int givePrescribedEquationNumber(void) { _error("givePrescribedEquationNumber: undefined");
                                             return 0; }
    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    int askNewEquationNumber(TimeStep *tStep) { return 1; }

    /**
     * Returns boundary condition of dof if it is precsribed.
     * HangingDof can not be subjected to bc, it is only mapping to master
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    int hasBc(TimeStep *tStep) { return 0; }

    /**
     * Returns initial condition of dof if it is precsribed.
     * HangingDof can not be subjected to ic, it is only mapping to master
     * @see MasterDof::hasIc
     */
    int hasIc() { return 0; }

    /**
     * RigidArmSlaveDof can not be subjected to ic - it is only mapping to master.
     * @see MasterDof::hasIc
     */
    int hasIcOn(ValueModeType) { return 0; }

    /**
     * Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id culd not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boubdaray condition, zero otherwise
     */
    int giveBcIdValue() { return 0; }

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType saveContext(FILE *stream, void *obj = NULL) { return CIO_OK; }

    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType restoreContext(FILE *stream, void *obj = NULL) { return CIO_OK; }

    /**
     * Returns class name of the receiver.
     */
    const char *giveClassName() const { return "SlaveDof"; }

    /**
     * Returns classType id of receiver.
     */
    classType giveClassID() const { return SlaveDofClass; }
};


#endif // slavedof_h
