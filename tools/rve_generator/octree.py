import numpy

# Class is named Oct but support arbitrary number of dimensions.
class Octtree:
    def __init__(self, ndim, depthlimit = 100):
        self.ndim = ndim
        self.midpoint = []
        # Subcells are stored in binary encoding
        # The positive side of the midpoint is 1, the negative 0, e.g. [0,1,1] -> 110b = 5
        self.subcells = []
        self.objects = list()
        self.depthlimit = depthlimit

    def __str__(self):
        if len(self.subcells) > 0:
            result += 'Midpoint: %s. '%self.midpoint
            result += '{\n'
            for s in self.subcells:
                result += s.__str__() + ';\n'
            result += '}'
        else:
            result += 'Leaf with %d objects ('%len(self.objects)
            for o in self.objects:
                result += '%f : %s,'%(o[0],o[1])
            result += ')'
        return result

    def insert(self, data, coord):
        if len(self.subcells) > 0:
            # This reinterprets the array of bools as an integer: e.g. [False, True, True] -> 110b = 5
            index = 0
            for k in range(self.ndim): 
                if coord[k] > self.midpoint[k]: index += 1 << k
            self.subcells[index].insert(data, coord)
        else:
            self.objects.append( (data,coord) )
            if self.depthlimit > 0 and len(self.objects) > 10:
                self.subdivide()

    def subdivide(self):
        self.subcells = [ Octtree(self.ndim, self.depthlimit - 1) for _ in range(1 << self.ndim) ]
        # Compute the midpoint:
        #self.midpoint = numpy.mean(self.coords, axis=0)
        self.midpoint = numpy.zeros(self.ndim)
        for o in self.objects:
            self.midpoint += o[1]
        self.midpoint /= len(self.objects)
        # Insert objects into subcells
        for o in self.objects:
            index = 0
            for k in range(self.ndim):
                if o[1][k] > self.midpoint[k]: index += 1 << k
            self.subcells[ index ].insert(o[0], o[1])
        self.objects = list()

    def getObjectsWithin(self, bbox0, bbox1):
        if len(self.subcells) > 0:
            result = list()
            check = [bbox0 > self.midpoint, bbox1 < self.midpoint] 
            # Loop through the subcells and check to see if the bbox overlaps
            for i in range(1 << self.ndim):
                # Easiest to check if any axis does *not* overlap (not overlap = outside)
                outside = False
                for j in range(self.ndim):
                    # Convert to the binary representation to see if which side we are on w.r.t the midpoint
                    if check[i >> j & 1][j]: # We thus shift to obtain the j'th component of i
                        outside = True
                        break
                if not outside:
                    result.extend(self.subcells[i].getObjectsWithin(bbox0, bbox1))
            return result
        else:
            return self.objects # Just take all values. Shouldnt be many false positives. 
