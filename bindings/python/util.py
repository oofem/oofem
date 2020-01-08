import sys
sys.path.insert(0, '../../default')
import oofempy


#
# Initializes domain with given sequences of nodes, elements, .....
#
def setupDomain(domain, nodes, elems, mats, css, bcs, ltfs, ics):
    domain.resizeDofManagers(len(nodes))  # process node list
    for n in nodes:
        domain.setDofManager(n.number, n)
    domain.resizeElements(len(elems))  # process elements
    for e in elems:
        domain.setElement(e.number, e)
    domain.resizeMaterials(len(mats))  # process materials
    for m in mats:
        domain.setMaterial(m.number, m)
    domain.resizeCrossSectionModels(len(css))  # process cross sections
    for c in css:
        domain.setCrossSection(c.number, c)
    domain.resizeBoundaryConditions(len(bcs))  # process boundary conditions
    for bc in bcs:
        domain.setBoundaryCondition(bc.number, bc)
    domain.resizeFunctions(len(ltfs))  # process (load time) functions
    for ltf in ltfs:
        domain.setFunction(ltf.number, ltf)
    domain.resizeInitialConditions(len(ics))  # process initial conditions
    for ic in ics:
        domain.setInitialCondition(ic.number, ic)
    return

