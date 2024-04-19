import oofempy

#
# Initializes domain with given sequences of nodes, elements, .....
#
def setupDomain (domain, nodes, elems, css, mats, bcs, ics, ltfs, sets):
    
    domain.resizeDofManagers(len(nodes)) # process nodes
    for n in nodes:
        domain.setDofManager(n.number, n)
    
    domain.resizeElements(len(elems)) # process elements
    for e in elems:
        domain.setElement(e.number, e)
    
    domain.resizeCrossSectionModels(len(css)) #process cross sections
    for c in css:
        domain.setCrossSection(c.number, c)
    
    domain.resizeMaterials(len(mats)) # process materials
    for m in mats:
        domain.setMaterial(m.number, m)
    
    domain.resizeBoundaryConditions(len(bcs)) #process boundary conditions
    for bc in bcs:
        domain.setBoundaryCondition(bc.number, bc)
    
    domain.resizeInitialConditions(len(ics)) #process initial conditions
    for ic in ics:
        domain.setInitialCondition(ic.number, ic)
    
    domain.resizeFunctions(len(ltfs)) # process (load time) functions
    for ltf in ltfs:
        domain.setFunction(ltf.number, ltf)
    
    domain.resizeSets(len(sets)) # process sets
    for s in sets:
        domain.setSet(s.number, s)
       
    return

