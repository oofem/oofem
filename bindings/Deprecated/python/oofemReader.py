# this will raise ImportError right away, if oofem wrapper is not found
import liboofem

from mupif import Cell, Field, FieldID, Mesh

class OofemReader(object):
    """
    Reads OOFEM problem through the Python interface, and return mirror of the data as mupif objects. OOFEM is general object oriented FEM solver (www.oofem.org). 
    Use this class like this::

    # create and solve a problem in oofem
    import liboofem
    dr=liboofem.OOFEMTXTDataReader("tmpatch42.in")
    pb=liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
    pb.checkProblemConsistency()
    pb.setRenumberFlag()
    pb.solveYourself()
    pb.terminateAnalysis()

    # instantiate the reader
    reader=mupif.OofemReader.OofemReader(model=pb)
    # return several fields; they will all share the same mesh
    f1=reader.makeField(fieldID=mupif.FieldID.FID_Displacement)
    f2=reader.makeField(fieldID=mupif.FieldID.FID_Strain)

    .. note:: This class has not been tested yet, as liboofem wrapper is not yet fully functional.
    """
    # shorthands
    _EGT=liboofem.Element_Geometry_Type 
    _FT=liboofem.FieldType
    # mapping from oofem element type enumeration to corresponding mupif cell class
    elementTypeMap={
        _EGT.EGT_triangle_1: Cell.Triangle_2d_lin,
        _EGT.EGT_quad_1:     Cell.Quad_2d_lin,
        _EGT.EGT_tetra_1:    Cell.Tetrahedron_3d_lin,
        _EGT.EGT_hexa_1:     Cell.Brick_3d_lin,
    }
    # mapping from mupif field ID to oofem field type
    fieldTypeMap={
        FieldID.FID_Displacement:  (_FT.FT_Displacements,0),
        FieldID.FID_Stress: None,
        FieldID.FID_Strain: None,
        FieldID.FID_Temperature:   (_FT.FT_Temperature,0),
        FieldID.FID_Humidity:      (_FT.FT_HumidityConcentration,0),
        FieldID.FID_Concentration: (_FT.FT_HumidityConcentration,1),
    }
    
    def __init__(self,model,domain=1):
        '''Construct the reader for later re-use; 
        
        :param liboofem.EngngModel model: oofem problem class representation 
        :param int domain: allows to select problem domain (default is 1, TODO: more user-friendly interface than number?)
        '''
        self.model=model
        if domain<1 or domain>model.giveNumberOfDomains(): raise ValueError('Invalid domain value (must be 1..%d, for this model instance)'%mode.giveNumberOfDomains())
        self.domain=mode.giveDomain(domain)
        self.mesh=None

    def _getMesh(self):
        '''Return mesh object from the model. Called internally from makeField, and the result will be re-used for multiple fields.

        :return: Mesh
        :rtype: mupif.Mesh.UnstructuredMesh
        '''
        if self.mesh: return self.mesh

        dom=self.domain
        self.mesh=Mesh.UnstructuredMesh()
        verts,cells=[],[]

        # vertices
        for dn in range(1,dom.giveNumberOfDofManagers()+1):
            dof=dom.giveDofManager(dn)
            cc=dof.giveCoordinates()
            assert len(cc) in (1,2,3) # 1D, 2D, 3D
            verts.append(tuple([cc[i] for i in range(len(cc))]))
            # note: seems that the previous line should read
            # verts.append(Vertex.Vertex(dn-1, dof.giveLabel(), tuple([cc[i] for i in range(len(cc))])))
        # cells
        for en in range(1,dom.giveNumberOfElements()+1):
            elt=dom.giveElement(en)
            gt=elt.giveGeometryType()
            eltType=elementTypeMap.get(gt) # mupif class which will be instantiated
            if not eltType: raise ValueError('Element type %s not mapped to any mupif cell type.'%(str(gt)))
            # convert to 0-based indices
            vv=tuple([elt.giveDofManager(n)-1 for n in range(1,elt.numberOfDofManagers+1)])
            cells.append(eltType(mesh=self.mesh,number=en-1,label=en-1,vertices=vv))
        # build mesh, cache it and return
        self.mesh.setup(verts,cells)
        return self.mesh

    def makeField(fieldID,timestep=None,valueModeType=liboofem.ValueModeType.VM_Total,units=None):
        '''Return field object for the model.

        TODO: pass timestep (uses current step at the moment), valueModeType and units as arguments?

        :param FieldID: type of field to be returned
        :param timestep: time step representation
        :param liboofem.ValueModeType valueModeType: determines the mode of requested value (total or incremental)
        :return: Field
        :rtype: mupif.Field.Field
        '''
        if timestep is None: timestep=model.giveCurrentStep()
        elif not isinstance(timestep,liboofem.TimeStep): raise ValueError("timestep must be a liboofem.TimeStep.")

        # user can give units in which oofem field is represented
        if units: raise RuntimeError("units not yet implemented.")

        dom=self.domain
        oofemFieldInfo=fieldTypeMap.get(fieldID)
        if not oofemFieldInfo: raise ValueError('Field type %d (%s) not found in this oofem model.'%(fieldID, fieldID.value ))
        field=self.model.giveField(oofemFieldInfo[0],timestep)

        # determine valuye type from the first DOF
        if dom.giveNumberOfDofManagers()<1: raise ValueError("Field has zero DofManagers, unable to determine value type.")
        val=liboofem.FloatArray()
        field.evaluateAtDman(val,dom.giveDofManager(1),valueModeType,timestep)
        valueType={1:ValueType.Scalar,3:ValueType.Vector,9:ValueType.Tensor}.get(len(val))
        if not valueType: raise ValueError("Unhandled value length in field: %d (should be 1, 3 or 9)"%(len(val)))
        checkLen=len(val)

        # create empty field, values set in the loop below
        ret=Field.Field(self._getMesh(),fieldID=fieldID,valueType=valueType,units=units,time=timestep.getTargetTime(),values=None,fieldType=FieldType.FT_vertexBased)

        # assign values at DOFs
        for dn in range(1,dom.giveNumberOfDofManagers()+1):
            dof=dom.giveDofManager(dn)
            field.evaluateAtDman(val,dof,valueModeType,timestep)
            if len(val)!=checkLen: raise ValueError("Invalid length at DofManager #%d: %d, should be %d."%(dn,len(val),checkLen))
            ret.setValue(dn-1,val)

        return ret


