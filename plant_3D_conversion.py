'''
This file contains code to reproduce evolutionary algorithms devised and created by Kark J. Niklas.
The following papers are referenced:
    1.  Niklas 1994: MORPHOLOGICAL EVOLUTION THROUGH COMPLEX DOMAINS OF FITNESS
    2.  Niklas 1997: ADAPTIVE WALKS THROUGH FITNESS LANDSCAPES FOR EARLY VASCULAR LAND PLANTS
'''

import matplotlib.pyplot as plt
import matplotlib as mlp
import numpy as np
import math


'''
This class defines a single segment of a plant called a telome
it refers to its parent and maintains the properties of length,
diameter, level, branch (1 or 2) and its angle to the vertical
'''
class Telome:
    def __init__(self, plant, level, branch=1, parent=0, 
                 inheritedVerticalAngle=math.radians(90),
                 L=10, d=1):
        
        # go into parents and indicate this as its child; two way "pointers"
        if branch == 1:
            bifurcation, rotation = plant.b1, plant.r1
            plant.telomes[parent].child1 = plant.telomes.__len__()
        elif branch == 2:
            bifurcation, rotation = plant.b2, plant.r2 + math.pi
            plant.telomes[parent].child2 = plant.telomes.__len__()
        elif branch == 0:
            bifurcation, rotation = 0, 0

        # defining shape/size
        self.length: float = L * level # for tapering
        self.diameter: float = d * level # scaling with level, not visualized
        self.radius: float = self.diameter / 2

        # defining hierarchy
        self.level: int = level
        self.parent: int = parent
        self.branch: int = branch
        self.child1: int = 0
        self.child2: int = 0
        
        # defining the vector
        # self.absoluteAngle: float = math.cos(rotation) * bifurcation + inheritedVerticalAngle
        self.absoluteAngle: float = bifurcation + inheritedVerticalAngle
        self.vector2D = [float(math.cos(self.absoluteAngle))*self.length,
                         float(math.sin(self.absoluteAngle))*self.length]
        
        if branch > 0: # branch
            self.base2D = plant.telomes[parent].terminus2D
        else: # trunk condition
            self.base2D = [float(0), float(0)]
        self.terminus2D = [self.base2D[0] + self.vector2D[0],
                           self.base2D[1] + self.vector2D[1]]
        
        # for MOMENTS:
        self.xProj: float = math.cos(self.absoluteAngle) * self.length
        self.weight: float = self.length * (self.diameter ** 2) # length * diameter^2, dropping scalars
        self.baseMoment: float = 0.5 * self.xProj * self.weight / (self.diameter ** 4) # stiffness increases with r^4
        
        # for HEIGHT:
        self.yProj: float = math.sin(self.absoluteAngle) * self.length


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
        self.r1: float = -math.radians(r1)
        self.r2: float = math.radians(r2)
        # bifurcation angle
        self.b1: float = -math.radians(b1)
        self.b2: float = math.radians(b2)

        self.levels = levels
    
    # angle getters return in degrees
    def r1(self):
        return math.degrees(self.r1)
    def r2(self):
        return math.degrees(self.r2)
    def b1(self):
        return math.degrees(self.b1)
    def b2(self):
        return math.degrees(self.b2)
    
    def appendTelomeCoordinates2D(self, coords):
        self.x2D += [coords[0]]
        self.y2D += [coords[1]]
    def initializeRender2D(self, baseCoords):
        self.x2D: list = []
        self.y2D: list = []
        self.appendTelomeCoordinates2D(baseCoords)
    def constructRender2Drec(self, child1=1, child2=2):
        if child1 > 0:
            print(f"child {child1} of parent {self.telomes[child1].parent}, level {self.telomes[child1].level} of abs {math.degrees(self.telomes[child1].absoluteAngle)}")
            self.appendTelomeCoordinates2D(self.telomes[child1].terminus2D)
            self.constructRender2Drec(self.telomes[child1].child1, self.telomes[child1].child2)
            self.appendTelomeCoordinates2D(self.telomes[child1].base2D)
        if child2 > 0:
            print(f"child {child2} of parent {self.telomes[child2].parent}, level {self.telomes[child2].level}, of abs {math.degrees(self.telomes[child1].absoluteAngle)}")
            self.appendTelomeCoordinates2D(self.telomes[child2].terminus2D)
            self.constructRender2Drec(self.telomes[child2].child1, self.telomes[child2].child2)
            self.appendTelomeCoordinates2D(self.telomes[child2].base2D)
        return
    def constructRender2D(self):
        self.initializeRender2D(self.telomes[0].base2D)
        self.appendTelomeCoordinates2D(self.telomes[0].terminus2D)
        self.constructRender2Drec(self.telomes[0].child1, self.telomes[0].child2)
    def renderPlant2D(self, withMoments=False):
        fig, ax = plt.subplots()
        ax.plot(self.x2D, self.y2D, linewidth=2.0)
        if withMoments:
            self.visualizeMoments()
            ax.scatter(self.momentScatterX, self.momentScatterY, c=self.momentScatterMoms)
            for i, mom in enumerate(self.momentScatterMoms):
                ax.annotate(mom, (self.momentScatterX[i], self.momentScatterY[i]))
        plt.show()

    # template for a recursive crawl of the telomes
    def recursiveCrawlRec(self, this=0):
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        if child1 > 0: self.recursiveCrawlRec(child1)
        if child2 > 0: self.recursiveCrawlRec(child2)
    def recursiveCrawl(self):
        self.recursiveCrawlRec()

    # MOMENTS
    # this function returns the center of mass of the subBranches to calculate base moments
    # on the way out
    def massCenterX(self, X1, X2, W1, W2):
        if W1+W2 == 0: return float(0)
        return (W1*X1 + W2*X2) / (W1+W2)
    def computeMomentsR(self, this=0):
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        thisWeight = self.telomes[this].weight
        thisX = self.telomes[this].base2D[0]

        if child1 > 0: xW1 = self.computeMomentsR(child1)
        else: xW1 = [float(0), float(0)]
        if child2 > 0: xW2 = self.computeMomentsR(child2)
        else: xW2 = [float(0), float(0)]

        Xcom: float = self.massCenterX(xW1[0], xW2[0], xW1[1], xW2[1])
        W: float = xW1[1] + xW2[1]
        self.telomes[this].baseMoment += W*(self.telomes[this].base2D[0] - Xcom)
        if abs(self.telomes[this].baseMoment) > abs(self.maxMoment): # track greatest for constraint
            self.maxMoment = self.telomes[this].baseMoment
        return [self.massCenterX(Xcom, thisX, W, thisWeight), thisWeight+W]
    def computeMoments(self): # WRAPPER
        self.maxMoment = 0
        self.computeMomentsR()
    
    def visualizeMomentsR(self, this=0):
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        if child1 > 0: self.visualizeMomentsR(child1)
        if child2 > 0: self.visualizeMomentsR(child2)
        offset = 0
        # branch = self.telomes[this].branch
        # if branch == 1: offset = 3
        # if branch == 2: offset = -3
        self.momentScatterX += [self.telomes[this].terminus2D[0]]
        self.momentScatterY += [self.telomes[this].terminus2D[1]+offset]
        self.momentScatterMoms += [int(self.telomes[this].baseMoment)]
        return
    def visualizeMoments(self): # WRAPPER
        self.momentScatterX: list = []
        self.momentScatterY: list = []
        self.momentScatterMoms: list = []
        self.visualizeMomentsR(0)


    # template for a recursive crawl of the telomes
    def computeHeighRec(self, height, this=0):
        if height > self.maxHeight: self.maxHeight = height
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        if child1 > 0: self.computeHeighRec(height + self.telomes[child1].yProj, child1)
        if child2 > 0: self.computeHeighRec(height + self.telomes[child2].yProj, child2)
    def computeHeight(self):
        self.maxHeight = 0
        self.computeHeighRec(float(self.telomes[0].yProj))
        return self.maxHeight

    # OPTIMIZE SURFACE AREA / VOLUME
    # this function has a wrapper and two recursive subcrawls of the telomes
    # to sum the surface area and the volume. Helper functions for each do the calculation
    # for each telome based on its radius and length
    def calculateTelomeVolume(self, this):
        result: float = math.pi * self.telomes[this].length * (self.telomes[this].radius ** 2)
        return result
    def volumeCrawlR(self, volume, this=0):
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        volume += self.calculateTelomeVolume(this)
        if child1 > 0: volume += self.volumeCrawlR(volume, child1)
        if child2 > 0: volume += self.volumeCrawlR(volume, child2)
        return volume
    def calculateTelomeSurfaceArea(self, this, child1, child2):
        result: float = math.pi * self.telomes[this].length * self.telomes[this].diameter
        # add half of cap area for each missing child
        if child1 == 0: result += 0.5 * math.pi * (self.telomes[this].radius ** 2)
        if child2 == 0: result += 0.5 * math.pi * (self.telomes[this].radius ** 2)
        return result
    def surfaceAreaCrawlR(self, surfaceArea, this=0):
        child1 = self.telomes[this].child1
        child2 = self.telomes[this].child2
        surfaceArea += self.calculateTelomeSurfaceArea(this, child1, child2)
        if child1 > 0: surfaceArea += self.surfaceAreaCrawlR(surfaceArea, child1)
        if child2 > 0: surfaceArea += self.surfaceAreaCrawlR(surfaceArea, child2)
        return surfaceArea
    def surfaceAreaToVolume(self, R=False):
        volume: float = self.volumeCrawlR(float(0), 0)
        surfaceArea: float = self.surfaceAreaCrawlR(float(0), 0)
        self.surfaceAreaToVolume = surfaceArea / volume
        if R: return self.surfaceAreaToVolume

    '''
    The generateTelomesR function is a recursive call to generate lower branches.
    As arguments, it receives the current level and the parent telome that called the function
    Each call spawns up to two more. 
    This version contineus to branch until a numnber of levels specified in the wrapper is reached
    ADD: I want to add a filter based on the plant probability/branching frequency
    '''
    def generateTelomesR(self, level, parent):
        if level == 0:
            return
        child1 = Telome(self, level, 1, parent, self.telomes[parent].absoluteAngle)
        self.telomes.append(child1)
        self.generateTelomesR(level-1, self.telomes.__len__()-1)
        child2 = Telome(self, level, 2, parent, self.telomes[parent].absoluteAngle)
        self.telomes.append(child2)
        self.generateTelomesR(level-1, self.telomes.__len__()-1)
        return 
    def generateTelomes(self, levels): # WRAPPER
        self.telomes = [Telome(self, levels, 0, 0, math.radians(90))]
        self.generateTelomesR(levels-1, 0)
        self.constructRender2D()
    

levels: int = 10
plant = Plant(0.5, 0.5, 0, 0, 26, 42, levels)
plant.generateTelomes(levels)
plant.computeMoments()
print(f"{plant.computeHeight()}")
print(plant.surfaceAreaToVolume(True))
plant.renderPlant2D(withMoments=True)

