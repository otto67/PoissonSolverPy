import numpy as np


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



class Elm2d:

    def __init__(self, elem_type="4bf4gf2d"):
        self.elm_type = elem_type
        self.curr_itg_pt = [0, 0]
        self.coor_at_itg_pt = [0, 0]
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



class Elm4bn4gf(Elm2d):

    def __init__(self):
        super().__init__("4bf4gf2d")
        self.nbf = 4
        self.ngf = 4
        self.l_l_c = [0,0]
        self.dx = 0
        self.dy = 0
        self.lb = 0


    def setPrms(self, x0, lenx, leny, elm_no, lb):
        self.l_l_c = x0
        self.dx = lenx
        self.dy = leny
        self.detJ = lenx*leny/4
        self.detJxW = lenx*leny/4
        self.elem_no = elm_no
        self.lb = lb

    def N(self, i):
        a, b, c, d = 0,0,0,0
        x_val = self.curr_itg_pt[0]
        y_val = self.curr_itg_pt[1]

        if (i == 1):
            a = b = -1/4
            c = d = 1/4
        if (i == 2):
            a = d = 1/4
            b = c = -1/4
        if (i == 3):
            a = b = c = d = 1/4
        if (i == 4):
            a = c = -1/4
            b = d = 1/4

        return (a*x_val) + (b*y_val) + (c*x_val*y_val) + d

    def dN(self, i, dir):

        x_val = self.curr_itg_pt[0]
        y_val = self.curr_itg_pt[1]

        if (i == 1 and dir == 1):
            return -1/4 + y_val/4
        if (i == 1 and dir == 2):
            return -1/4 + x_val/4

        if (i == 2 and dir == 1):
            return 1/4 - y_val/4
        if (i == 2 and dir == 2):
            return -1/4 - x_val/4

        if (i == 3 and dir == 1):
            return 1/4 + y_val/4
        if (i == 3 and dir == 2):
            return 1/4 + x_val/4

        if (i == 4 and dir == 1):
            return -1/4 - y_val/4
        if (i == 4 and dir == 2):
            return 1/4 - x_val/4

        print("dN : Invalid direction and number: ", dir, i)
        return 0


    def initAtItgPt(self, pt):

        a = 1/(3**0.5)

        if (pt == 1):
            self.curr_itg_pt = [a, a]
        if (pt == 2):
            self.curr_itg_pt = [-a, a]
        if (pt == 3):
            self.curr_itg_pt = [-a, -a]
        if (pt == 4):
            self.curr_itg_pt = [a, -a]

        x_0 = self.l_l_c[0]
        y_0 = self.l_l_c[1]
        x_1 = x_0 + self.dx
        y_1 = y_0 + self.dy

        self.coor_at_itg_pt = [((self.dx/2) * self.curr_itg_pt[0]) + ((x_0 + x_1)/2), ((self.dy/2) * (self.curr_itg_pt[1])) + ((y_0 + y_1) / 2)]

    def initAtItgPtBound(self, pt):

        a = 1/(3**0.5)

        if (pt == 1):
            self.curr_itg_pt = [-a]
        if (pt == 2):
            self.curr_itg_pt = [a]

        x_0 = self.l_l_c[0]
        y_0 = self.l_l_c[1]
        x_1 = x_0 + self.dx
        y_1 = y_0 + self.dy

        self.coor_at_itg_pt = [((self.dx/2) * self.curr_itg_pt[0]) + (x_0 + x_1)/2]


# Only one type of elements in this base class
class Grid2d:
    n_elms = 100

    def __init__(self, prms, elm_type="4bf4gf2d"):
        self.dscprms = prms
        self.n_elms = (prms.nno_x - 1)*(prms.nno_y - 1)
        self.elem_type = elm_type
        self.boindsWithEssBC = []
        self.essbcnodes = {}
        self.nno = prms.nno_x

    # possibly read from file here, to create a more complicated mesh

        dx = prms.x_max / (self.nno - 1)
        dy = prms.y_max / (self.nno - 1)

        self.elems = []
        for i in range(self.n_elms):
            elm = Elm4bn4gf()
            x0 = [i % (self.nno-1), i // (self.nno-1)]
            x0[0] *= dx
            x0[1] *= dy

            elm.setPrms(x0, dx, dy, i, -1)
            self.elems.append(elm)

        self.A = np.zeros((self.nno ** 2, self.nno ** 2))
        self.b = np.zeros(self.nno ** 2)
        self.phi = np.zeros(self.nno ** 2)

        self.dofmap = []


        for j in range(0, self.nno-1):
            for i in range(0, self.nno-1):
                self.dofmap.append((j*self.nno + i, j*self.nno + i + 1, (j+1)*self.nno + i + 1, (j+1)*self.nno + i))


        self.boNodeMap = [[[] for i in range(self.nno)] for j in range(self.nno)]
        for j in range(self.nno):
            self.boNodeMap[0][j] = [0, j*dy, 3]
            self.boNodeMap[self.nno-1][j] = [self.nno, j*dy, 1]

        for i in range(self.nno):
            self.boNodeMap[i][0] = [i*dx, 0, 4]
            self.boNodeMap[i][self.nno-1] = [i*dx, self.nno, 2]

    def addBC(self, i, j, val):
        globdof = j*self.nno + i
        self.essbcnodes.update({globdof: val})

    # local node number in element to global node number
    def loc2glob(self, loc_no, e):
        return self.dofmap[e][loc_no]

#        return (self.dofmap[e][loc_no] % (self.dscprms.nno_x), self.dofmap[e][loc_no] // (self.dscprms.nno_x))

    def setBoindWithEssBC(self, list):
        self.boindsWithEssBC = list

    def isBoNode(self, i, j):
        if (self.boNodeMap[i][j]):
            if (self.boNodeMap[i][j][2] in self.boindsWithEssBC):
                return True
        return False

    def interpolSolution(self):
        solu = np.zeros((self.nno, self.nno))

        for i in range(self.nno):
            for j in range(self.nno):
                solu[i,j] = self.phi[i*self.nno + j]

        """
        for elmno in range(self.n_elms):
            for r in range(self.elems[elmno].ngf):
                for s in range(self.elems[elmno].nbf):
                    solu[self.loc2glob(r, elmno) // self.nno, self.loc2glob(s, elmno) % self.nno] += self.phi[self.loc2glob(s, elmno) % self.nno]
        """

        return solu


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    parameters = DiscPrms(nnx =3, nny =3, dt=1000, t_max=2000)
    print("parameters", parameters.nno_x)
    grid = Grid2d(parameters)
    grid.setBoindWithEssBC([2, 4])
    for i in range(6):
        for j in range(6):
            if not grid.isBoNode(i, j):
                print("Node number ", j*parameters.nno_x + i, " is not on any boundary")
            else:
                print("Node number ", j*parameters.nno_x + i, " has values: ", grid.boNodeMap[i][j])
