# Class for handling essential boundary conditions
# for both FDM and FEM

class BC:
    def __init__(self) -> None:
        pass

    def essBC(self, x, y, boind):
        return 0
    
    def attachBC(bcs_):
        pass

class constBC(BC):
    def __init__(self) -> None:
        self.bcs = []

    def attachBC(self, bcs_):
        self.bcs = bcs_

    def essBC(self, x, y, boind):

        if not self.bcs:
            return 0.0 

        if boind == 3:
            return self.bcs[2]
        if boind == 4:
            return self.bcs[3]
        if boind == 2:
            return self.bcs[1]
        if boind == 1:
            print(self.bcs[0])
            return self.bcs[0]
        
        return 0.0
