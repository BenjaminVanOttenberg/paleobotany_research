'''
This file contains code to reproduce evolutionary algorithms devised and created by Kark J. Niklas.
The following papers are referenced:
    1.  Niklas 1994: MORPHOLOGICAL EVOLUTION THROUGH COMPLEX DOMAINS OF FITNESS
    2.  Niklas 1997: ADAPTIVE WALKS THROUGH FITNESS LANDSCAPES FOR EARLY VASCULAR LAND PLANTS
'''

import matplotlib.pyplot as plt
import numpy as np
import math

'''
This class defines a single segment of a plant called a telome
it refers to its parent and maintains the properties of length,
diameter, level, branch (1 or 2) and its angle to the vertical
'''
class Telome:
    def __init__(self, plant, level, branch=1, parent=0, 
                 inheritedVerticalAngle=0,
                 L=10, d=1):
        
        if branch == 1:
            bifurcation, rotation = plant.b1, plant.r1
            plant.telomes[parent].child1 = plant.telomes.__len__()
        elif branch == 2:
            bifurcation, rotation = plant.b2, plant.r2 + math.pi
            plant.telomes[parent].child2 = plant.telomes.__len__()
        elif branch == 0:
            bifurcation, rotation = 0, 0

        # defining shape/size
        # self.length: float = L
        self.length: float = L * level # for tapering
        self.diameter: float = d

        # defining hierarchy
        self.level: int = level
        self.parent: int = parent
        self.branch: int = branch
        self.child1: int = 0
        self.child2: int = 0
        

        # defining the vector
        # self.absoluteAngle: float = math.cos(rotation) * bifurcation + inheritedVerticalAngle
        self.absoluteAngle: float = bifurcation + inheritedVerticalAngle
        self.vector2D = [float(math.sin(self.absoluteAngle))*self.length,
                         float(math.cos(self.absoluteAngle))*self.length]
        
        # self.vector3D = [float(math.sin(self.absoluteAngle))*float(math.cos(self.absoluteRotation))*self.length,
        #                  float(math.cos(self.absoluteAngle))*self.length,
        #                  float(math.sin(self.absoluteAngle))*float(math.sin(self.absoluteRotation))*self.length]
        self.vector3D: list = plant.transformCoordinates(level, branch, self.length)
        print(self.vector3D)
        if branch > 0: # branch
            self.base2D = plant.telomes[parent].terminus2D
            self.base3D = plant.telomes[parent].terminus3D
        else: # trunk condition
            self.base2D = [float(0), float(0)]
            self.base3D = [float(0), float(0), float(0)]
        self.terminus2D = [self.base2D[0] + self.vector2D[0],
                           self.base2D[1] + self.vector2D[1]]
        self.terminus3D = [self.base3D[0] + self.vector3D[0],
                           self.base3D[1] + self.vector3D[1],
                           self.base3D[2] + self.vector3D[2]]
        


    def getAbsoluteAngle(self):
        return self.absoluteAngle
    
    def getDiameter(self):
        return self.diameter
    
    def getLength(self):
        return self.diameter
    
    def getBranch(self):
        return self.branch
    
    def getBase2D(self):
        return self.base2D
    
    def getTerminus2D(self):
        return self.terminus2D
    
    def getBase3D(self):
        return self.base3D
    
    def getTerminus3D(self):
        return self.terminus3D
    

    def __str__(self):
        return f"Level {self.level}, branch {self.branch} telome. Child of  {self.parent}. Absolute angle: {self.absoluteAngle}. Moment: {self.maxBendingStress}"
    
    '''
    Three layer function for calculating the maximum bending stresses
    on a telome from Niklas 1994. 
        Eqn 2a: maxBendingStress()
        Eqn 2b: bendingMoment() is modified according my recalulation
    In generateMoment, I added functiion that reduces the moment according to its starting angle
    In Plant.generateTelomesR(), the moments are recusively added together descending down the plant which 
        is my own extension for Niklas 1994, inconsistent with his more thorough moment explanation in Niklas 1997
    '''
    def bendingMoment(self, rotationAngle, bifurcationAngle, d, L, rho=1, g=10):
        result = math.pi / rotationAngle
        result *= d * d * L * L * rho * g
        result += math.sin(bifurcationAngle)
        return result
    def maxBendingStress(self, rotationAngle, bifurcationAngle, d, L):
        result = 32 * self.bendingMoment(rotationAngle, bifurcationAngle, d, L)
        result /= math.pi * d * d * d
        return result
    def generateMoment(self, rotationAngle, bifurcationAngle):
        self.maxBendingStress = math.sin(self.absoluteAngle) * self.maxBendingStress(rotationAngle, bifurcationAngle, self.diameter, self.length)

'''
This class defines an entire plant in the morphology. It is defined by
its 6 parameters of branching probabilities, rotation angles, and bifurcation
angles. It generates from these a list of telomes
'''
class Plant:
    def __init__(self, p1, p2, r1, r2, b1, b2, levels: int):
        # these parameters are from Niklas 1994
        # branching probability
        self.p1: float = p1
        self.p2: float = p2
        # rotation angle
        self.r1: float = math.radians(r1)
        self.r2: float = -math.radians(r2)
        # bifurcation angle
        self.b1: float = math.radians(b1)
        self.b2: float = -math.radians(b2)

        self.levels = levels

        self.generateRotationMatrices1(levels)
        self.generateRotationMatrices2(levels)

    def prob1(self):
        return self.p1
    def prob2(self):
        return self.p2
    
    # angle getters return in degrees
    def r1(self):
        return math.degrees(self.r1)
    def r2(self):
        return math.degrees(self.r2)
    def b1(self):
        return math.degrees(self.b1)
    def b2(self):
        return math.degrees(self.b2)
    
    def transformCoordinates(self, level, branch, L=1):
        if branch == 1:
            transform = np.matrix(self.rotationMatrices1[self.levels - level] @ self.xPrime1)
            return [transform[0,0]*L,transform[1,0]*L,transform[2,0]*L]
        if branch == 2:
            transform = np.matrix(self.rotationMatrices2[self.levels - level] @ self.xPrime2)
            return [transform[0,0]*L,transform[1,0]*L,transform[2,0]*L]
        if branch == 0:
            return [0, L, 0]
    
    def appendTelomeCoordinates2D(self, coords):
        self.x2D += [coords[0]]
        self.y2D += [coords[1]]

    def appendTelomeCoordinates3D(self, coords):
        self.x3D += [coords[0]]
        self.y3D += [coords[1]]
        self.z3D += [coords[2]]

    def initializeRender2D(self, baseCoords):
        self.x2D: list = []
        self.y2D: list = []
        self.appendTelomeCoordinates2D(baseCoords)

    def initializeRender3D(self, baseCoords):
        self.x3D: list = []
        self.y3D: list = []
        self.z3D: list = []
        self.appendTelomeCoordinates3D(baseCoords)

    def constructRender2Drec(self, child1=1, child2=2):
        if child1 > 0:
            print(f"child {child1} of parent {self.telomes[child1].parent} of abs {self.telomes[child1].getAbsoluteAngle()}")
            self.appendTelomeCoordinates2D(self.telomes[child1].getTerminus2D())
            self.constructRender2Drec(self.telomes[child1].child1, self.telomes[child1].child2)
            self.appendTelomeCoordinates2D(self.telomes[child1].getBase2D())
        if child2 > 0:
            print(f"child {child2} of parent {self.telomes[child2].parent} of abs {self.telomes[child1].getAbsoluteAngle()}")
            self.appendTelomeCoordinates2D(self.telomes[child2].getTerminus2D())
            self.constructRender2Drec(self.telomes[child2].child1, self.telomes[child2].child2)
            self.appendTelomeCoordinates2D(self.telomes[child2].getBase2D())
        return
    
    def constructRender3Drec(self, child1=1, child2=2):
        if child1 > 0:
            self.appendTelomeCoordinates3D(self.telomes[child1].getTerminus3D())
            self.constructRender3Drec(self.telomes[child1].child1, self.telomes[child1].child2)
            self.appendTelomeCoordinates3D(self.telomes[child1].getBase3D())
        if child2 > 0:
            self.appendTelomeCoordinates3D(self.telomes[child2].getTerminus3D())
            self.constructRender3Drec(self.telomes[child2].child1, self.telomes[child2].child2)
            self.appendTelomeCoordinates3D(self.telomes[child2].getBase3D())
        return

    def constructRender2D(self):
        self.initializeRender2D(self.telomes[0].getBase2D())
        self.appendTelomeCoordinates2D(self.telomes[0].getTerminus2D())
        self.constructRender2Drec(self.telomes[0].child1, self.telomes[0].child2)

    def constructRender3D(self):
        self.initializeRender3D(self.telomes[0].getBase3D())
        self.appendTelomeCoordinates3D(self.telomes[0].getTerminus3D())
        self.constructRender3Drec(self.telomes[0].child1, self.telomes[0].child2)

    def renderPlant2D(self):
        fig, ax = plt.subplots()
        ax.plot(self.x2D, self.y2D, linewidth=2.0)
        plt.show()

    def renderPlant3D(self):
        ax = plt.figure().add_subplot(projection='3d')
        ax.plot(self.x3D, self.z3D, self.y3D, linewidth=1.5)
        plt.show()

    # template for a recursive crawl of the telomes
    def recursiveCrawlRec(self, child1=1, child2=2):
        if child1 > 0:
            self.recursiveCrawlRec(self.telomes[child1].child1, self.telomes[child1].child2)
        
        if child2 > 0:
            self.recursiveCrawlRec(self.telomes[child2].child1, self.telomes[child2].child2)
    def recursiveCrawl(self):
        self.recursiveCrawlRec(self.telomes[0].child1, self.telomes[0].child2)


    def generateRotationMatrix(self, bifurcation_angle, rotation_angle):
        # Rotation around the X-axis for bifurcation
        R_x = np.matrix([
            [1, 0, 0],
            [0, math.cos(bifurcation_angle), -math.sin(bifurcation_angle)],
            [0, math.sin(bifurcation_angle), math.cos(bifurcation_angle)]
        ])

        # Rotation around the Z-axis for branch rotation
        R_z = np.matrix([
            [math.cos(rotation_angle), -math.sin(rotation_angle), 0],
            [math.sin(rotation_angle), math.cos(rotation_angle), 0],
            [0, 0, 1]
        ])

        # Combine rotations: Apply bifurcation first, then rotation
        return R_z @ R_x

    # this function computes the coordinate transformation for 3D reference frame
    # changes at each branching event
    # branch #1
    '''
    def generateRotationMatrices1R(self, level):
        if level == 0:
            return
        rotM = self.rotationMatrices1[self.rotationMatrices1.__len__()-1]  @ self.rotationMatrices1[0]
        print(rotM)
        self.rotationMatrices1.append(rotM)
        self.generateRotationMatrices1R(level-1)

    def generateRotationMatrices1(self, levels):
        self.xPrime1 = np.matrix([[float(math.sin(self.b1) * math.cos(self.r1))],
                                  [float(math.cos(self.b1))], 
                                  [float(math.sin(self.b1) * math.sin(self.r1))]])
        rotM = np.matrix([[float( math.sin(self.b1) * math.cos(self.r1)),
                           float( math.cos(self.b1)),
                           float( math.sin(self.b1) * math.sin(self.r1))],
                          [float( math.cos(self.b1) * math.cos(self.r1)),
                           float(-math.sin(self.b1)),
                           float( math.cos(self.b1) * math.sin(self.r1))],
                          [float(-math.sin(self.r1)),
                           float(0),
                           float( math.cos(self.r1))]]).transpose()
        print(rotM)
        self.rotationMatrices1: list = [rotM]
        self.generateRotationMatrices1R(levels)

    # now for branch #2
    def generateRotationMatrices2R(self, level):
        if level == 0:
            return
        rotM = self.rotationMatrices2[self.rotationMatrices2.__len__()-1] @ self.rotationMatrices2[0]
        print(rotM)
        self.rotationMatrices2.append(rotM)
        self.generateRotationMatrices2R(level-1)

    def generateRotationMatrices2(self, levels):
        self.xPrime2 = np.matrix([[float(math.sin(self.b2) * math.cos(self.r2))], 
                                  [float(math.cos(self.b2))], 
                                  [float(math.sin(self.b2) * math.sin(self.r2))]])
        rotM = np.matrix([[float( math.sin(self.b2) * math.cos(self.r2)),
                           float( math.cos(self.b2)),
                           float( math.sin(self.b2) * math.sin(self.r2))],
                          [float( math.cos(self.b2) * math.cos(self.r2)), 
                           float(-math.sin(self.b2)), 
                           float( math.cos(self.b2) * math.sin(self.r2))],
                          [float(-math.sin(self.r2)), 
                           float(0), 
                           float( math.cos(self.r2))]]).transpose()
        print(rotM)
        self.rotationMatrices2: list = [rotM]
        self.generateRotationMatrices2R(levels)
    '''
    def generateRotationMatrices1(self, levels):
        self.xPrime1 = np.matrix([[float(math.sin(self.b1) * math.cos(self.r1))],
                                  [float(math.cos(self.b1))], 
                                  [float(math.sin(self.b1) * math.sin(self.r1))]])
        self.rotationMatrices1 = []
        for i in range(levels):
            rotM = self.generateRotationMatrix(self.b1, self.r1)
            if i > 0:
                # Compound previous rotation with the new rotation
                rotM = self.rotationMatrices1[-1] @ rotM
            self.rotationMatrices1.append(rotM)

    def generateRotationMatrices2(self, levels):
        self.xPrime2 = np.matrix([[float(math.sin(-self.b2) * math.cos(self.r2))], 
                                  [float(math.cos(-self.b2))], 
                                  [float(math.sin(-self.b2) * math.sin(self.r2))]])
        self.rotationMatrices2 = []
        for i in range(levels):
            rotM = self.generateRotationMatrix(-self.b2, self.r2)
            if i > 0:
                # Compound previous rotation with the new rotation
                rotM = self.rotationMatrices2[-1] @ rotM
            self.rotationMatrices2.append(rotM)

    '''
    the generateTelomesR function is a recursive call to generate lower branches.
    As arguments, it receives the current level and the parent telome that called the function
    Each call spawns up to two more. 
    This version contineus to branch until a numnber of levels specified in the wrapper is reached
    I want to add a filter based on the plant probability/branching frequency
    Each bending moment is returned to that they can compound down the plant
    The plant coordinates are added to the 2d render lists up and down such that matplotlib can graph it out
    '''
    def generateTelomesR(self, level, parent):
        if level == 0:
            return self.telomes[self.telomes.__len__()-1].maxBendingStress

        child1 = Telome(self, level, 1, parent, self.telomes[parent].getAbsoluteAngle())
        child1.generateMoment(self.r1, self.b1)
        # if child1.maxBendingStress > self.bendingStress:
        #     self.bendingStress = child1.maxBendingStress
        self.telomes.append(child1)
        # self.appendTelomeCoordinates2D(child1.getTerminus2D()) # plot terminus of new branch
        child1.maxBendingStress += self.generateTelomesR(level-1, self.telomes.__len__()-1)
        # self.appendTelomeCoordinates2D(child1.getBase2D()) # plot return to base of current branch

        child2 = Telome(self, level, 2, parent, self.telomes[parent].getAbsoluteAngle())
        child2.generateMoment(self.r2, self.b2)
        # if child2.maxBendingStress > self.bendingStress:
        #     self.bendingStress = child2.maxBendingStress
        self.telomes.append(child2)
        # self.appendTelomeCoordinates2D(child2.getTerminus2D()) # plot terminus of new branch
        child2.maxBendingStress += self.generateTelomesR(level-1, self.telomes.__len__()-1)
        # self.appendTelomeCoordinates2D(child1.getTerminus2D()) # plot return to base of current branch

        return child1.maxBendingStress + child2.maxBendingStress

    def generateTelomes(self, levels):
        # wrapper
        self.telomes = [Telome(self, levels, 0, 0)]

        # initialize render
        # self.initializeRender2D(self.telomes[0].base2D)
        # self.appendTelomeCoordinates2D(self.telomes[0].terminus2D)

        # establish starting conditions for moment calculations
        # then start the recursive call
        self.bendingStress = 0
        self.telomes[0].maxBendingStress = self.generateTelomesR(levels-1, 0)

        # maybe move this to a getter function
        self.constructRender2D()
        self.constructRender3D()
        return self.telomes
    


            
levels: int = 7
plant = Plant(0.5, 0.5, 20, 20, 30, 30, levels)
Telomes = plant.generateTelomes(levels)
plant.renderPlant2D()
plant.renderPlant3D()

# for telome in Telomes:
#     print(telome)