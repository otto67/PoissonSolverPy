import numpy as np
import matplotlib.pyplot as plt

# Helper classes for discretization 
# General grid parameters for a rectangular box
class DiscPrms:

    def __init__(self, nnx=41, nny=41, dt=0.05, x_max=1, y_max=1, t_max=1, x_min=0, y_min=0):
        self.nno_x = nnx
        self.nno_y = nny
        self.dt = dt
        self.x_max = x_max
        self.y_max = y_max
        self.t_max = t_max
        self.x_min = x_min
        self.y_min = y_min
        self.dx = x_max / (nnx - 1)
        self.dy = y_max / (nny - 1)


# Two-dimensional elements for FEM
class Elm2d:

    def __init__(self, elem_type="4bf4gf2d"):
        self.elm_type = elem_type
        self.curr_itg_pt = [0.0, 0.0]
        self.coor_at_itg_pt = [0.0, 0.0]
        self.detJ = 0.0
        self.detJxW = 0.0

    def setPrms(self, x0, lenx, leny, elm_no, lb):
        raise NotImplementedError()

    def N(self, i):
        raise NotImplementedError()

    def dN(self, i, dir_):
        raise NotImplementedError()

    def initAtItgPt(self, coor):
        raise NotImplementedError()

    def initAtItgPtBound(self, pt):
        raise NotImplementedError()

# Bilinear element 
class Elm4bn4gf(Elm2d): 

    def __init__(self):
        super().__init__("4bf4gf2d")
        self.nbf = 4
        self.ngf = 4
        self.l_l_c = [0,0] # Global coords for lower left corner of element
        self.dx = 0
        self.dy = 0
        self.lb = 0


    def setPrms(self, x0, lenx, leny, elm_no, lb):
        self.l_l_c = x0
        self.dx = lenx
        self.dy = leny
        self.detJ = (lenx*leny)/4
        self.detJxW = (lenx*leny)/4
        self.elem_no = elm_no
        self.lb = lb

# Basis function in integration domain coordinates
    def N(self, i):
        a, b, c, d = 0,0,0,0
        x_val = self.curr_itg_pt[0]
        y_val = self.curr_itg_pt[1]

        if (i == 1):
            a = b = -1/4
            c = d = 1/4
        elif (i == 2):
            a = d = 1/4
            b = c = -1/4
        elif (i == 3):
            a = b = c = d = 1/4
        elif (i == 4):
            a = c = -1/4
            b = d = 1/4

        return (a*x_val) + (b*y_val) + (c*x_val*y_val) + d

# Derivative of basis function in integration domain coordinates
    def dN(self, i, dir):

        x_val = self.curr_itg_pt[0]
        y_val = self.curr_itg_pt[1]

        if (i == 1 and dir == 1):
            return -1/4 + y_val/4
        elif (i == 1 and dir == 2):
            return -1/4 + x_val/4
        elif (i == 2 and dir == 1):
            return 1/4 - y_val/4
        elif (i == 2 and dir == 2):
            return -1/4 - x_val/4
        elif (i == 3 and dir == 1):
            return 1/4 + y_val/4
        elif (i == 3 and dir == 2):
            return 1/4 + x_val/4
        elif (i == 4 and dir == 1):
            return -1/4 - y_val/4
        elif (i == 4 and dir == 2):
            return 1/4 - x_val/4

        print("dN : Invalid direction and number: ", dir, i)
        return 0

    # Initialize for numerical integration in domain
    def initAtItgPt(self, pt):

        a = 1.0/(3**(0.5))

        if (pt == 1):
            self.curr_itg_pt = [a, a]         
        elif (pt == 2):
            self.curr_itg_pt = [-a, a]
        elif (pt == 3):
            self.curr_itg_pt = [-a, -a]
        elif (pt == 4):
            self.curr_itg_pt = [a, -a]
    
        x_0 = self.l_l_c[0]
        y_0 = self.l_l_c[1]

        # Global coordinates transformed into integration domain
        # for use in right hand side
        self.coor_at_itg_pt = [x_0 + ((self.dx/2.0) * (self.curr_itg_pt[0] + 1.0)), 
        y_0 + ((self.dy/2.0) * (self.curr_itg_pt[1] + 1.0))]



    # Initialize for numerical integration along boundary
    # Not yet used
    def initAtItgPtBound(self, pt):

        a = 1/(3**(0.5))

        if (pt == 1):
            self.curr_itg_pt = [-a]
        elif (pt == 2):
            self.curr_itg_pt = [a]

        x_0 = self.l_l_c[0]
        y_0 = self.l_l_c[1]
     
        self.coor_at_itg_pt = [((self.dx/2) * (self.curr_itg_pt[0] +1))]


# Represents a 2d grid. 
class Grid2d:
    
    def __init__(self, prms, elm_type="4bf4gf2d"):
        self.dscprms = prms
        self.n_elms = (prms.nno_x - 1)*(prms.nno_y - 1)
        self.elem_type = elm_type
        self.boinds_with_essBC = []
        self.ess_bc_nodes = {}

    # Perhaps read from file here, to create a more complicated mesh

        dx = prms.dx
        dy = prms.dy

        self.elems = []
        for e in range(self.n_elms):
            elm = Elm4bn4gf()
            x0 = [(e % (prms.nno_x-1))*dx, (e // (prms.nno_x-1))*dy]
            elm.setPrms(x0, dx, dy, e, -1)
            self.elems.append(elm)
           
        self.A = np.zeros((prms.nno_x * prms.nno_y, prms.nno_y * prms.nno_x))
        self.b = np.zeros(prms.nno_x * prms.nno_y)
        self.phi = np.zeros(prms.nno_x * prms.nno_y)

        self.dofmap = []

        for j in range(0, prms.nno_y-1):
            for i in range(0, prms.nno_x-1):
                self.dofmap.append((j*prms.nno_y + i, j*prms.nno_y + i + 1, (j+1)*prms.nno_y + i + 1, (j+1)*prms.nno_y + i))

       
        self.boNodeMap = [[[] for i in range(prms.nno_x)] for j in range(prms.nno_y)]
        for j in range(prms.nno_y):
            self.boNodeMap[0][j] = [0, j*dy, 3]
            self.boNodeMap[prms.nno_y-1][j] = [1, j*dy, 1]

        for i in range(prms.nno_x):
            self.boNodeMap[i][0] = [i*dx, 0, 4]
            self.boNodeMap[i][prms.nno_x-1] = [i*dx, 1, 2]

    def addBC(self, i, j, val):
        globdof = j*self.dscprms.nno_x + i
        self.ess_bc_nodes.update({globdof: val})

    # local node number in element to global node number
    def loc2glob(self, loc_no, e):
        return self.dofmap[e][loc_no]

    # List of boundaries along which essbc's are imposed 
    def setBoindWithEssBC(self, list):
        self.boinds_with_essBC = list

    def isBoNode(self, i, j):
        if (self.boNodeMap[i][j]):
            if (self.boNodeMap[i][j][2] in self.boinds_with_essBC):
                return True
        return False

    # Solution at nodal points
    def interpolSolution(self):
        solu = np.zeros((self.dscprms.nno_x, self.dscprms.nno_y))
        for i in range(self.dscprms.nno_x):
            for j in range(self.dscprms.nno_y):
                solu[i,j] = self.phi[i*self.dscprms.nno_x + j]
        return solu

if __name__ == '__main__':

    parameters = DiscPrms(nnx=17, nny =17, dt=1000, t_max=2000)
    grid = Grid2d(parameters)
    grid.setBoindWithEssBC([2, 4])
    for i in range(parameters.nno_x):
        for j in range(parameters.nno_y):
            if not grid.isBoNode(i, j):
                print("Node number ", j*parameters.nno_x + i, " is not on any boundary")
            else:
                print("Node number ", j*parameters.nno_x + i, " has values: ", grid.boNodeMap[i][j])
