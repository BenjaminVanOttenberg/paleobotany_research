'''
This file contains code to reproduce evolutionary algorithms devised and created by Kark J. Niklas.
The following papers are referenced:
    1.  Niklas 1994: MORPHOLOGICAL EVOLUTION THROUGH COMPLEX DOMAINS OF FITNESS
    2.  Niklas 1997: ADAPTIVE WALKS THROUGH FITNESS LANDSCAPES FOR EARLY VASCULAR LAND PLANTS
'''

import matplotlib.pyplot as plt
#import matplotlib as mlp
from matplotlib.collections import LineCollection 
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import math
import random
import tkinter as tk
from tkinter import ttk
import time


'''
This class defines a single segment of a plant called a telome
it refers to its parent and maintains the properties of length,
diameter, level, branch (1 or 2) and its angle to the vertical
'''
class Telome:
    def __init__(self, plant, level, branch=1, parent=0, 
                 inheritedVerticalAngle=math.radians(90),
                 L=10, d=1, taperL=1, taperD=1):
        
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
        self.length: float = L * level * (taperL ** level) # for tapering
        self.diameter: float = d * level * (taperD ** level) # scaling with level, not visualized
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
        self.weight: float = self.length * (self.radius ** 2) # length * diameter^2, dropping scalars
        self.baseMoment: float = 0.5 * self.xProj * self.weight / (self.radius ** 4) # stiffness increases with r^4
        
        # for HEIGHT:
        self.yProj: float = math.sin(self.absoluteAngle) * self.length
    

'''
This class defines an entire plant in the morphology. It is defined by
its 6 parameters of branching probabilities, rotation angles, and bifurcation
angles. It generates from these a list of telomes
'''
class Plant:
    def __init__(self, p1, p2, r1, r2, b1, b2, levels: int, taperL:float=1, taperD:float=1):
        # these parameters are from Niklas 1994
        # branching probability
        self.p1: float = p1
        self.p2: float = p2
        # rotation angle
        self.r1deg = r1
        self.r2deg = r2
        self.r1: float = math.radians(r1)
        self.r2: float = -math.radians(r2)
        # bifurcation angle
        self.b1deg = b1
        self.b2deg = b2
        self.b1: float = math.radians(b1)
        self.b2: float = -math.radians(b2)

        self.levels: int = levels
    
    # angle getters return in degrees
    def r1(self):
        return math.degrees(self.r1)
    def r2(self):
        return math.degrees(self.r2)
    def b1(self):
        return math.degrees(self.b1)
    def b2(self):
        return math.degrees(self.b2)
    
    # def getProbabilitiesR(self, k):
    #     if k == 0: return 1
    #     result: float = float(8 / (self.levels + k)) * self.getProbabilitiesR(k-1)
    #     self.terminationProbs += [result]
    #     return result
    # def getProbabilities(self, R=False):
    #     self.terminationProbs: list = []
    #     self.getProbabilitiesR(self.levels)
    #     if R: return self.terminationProbs

    def appendTelomeCoordinates2D(self, coords):
        self.x2D += [coords[0]]
        self.y2D += [coords[1]]
    def initializeRender2D(self, baseCoords):
        self.x2D: list = []
        self.y2D: list = []
        self.appendTelomeCoordinates2D(baseCoords)
    def constructRender2Drec(self, child1=1, child2=2):
        if child1 > 0:
            #print(f"child {child1} of parent {self.telomes[child1].parent}, level {self.telomes[child1].level} of abs {math.degrees(self.telomes[child1].absoluteAngle)}")
            self.appendTelomeCoordinates2D(self.telomes[child1].terminus2D)
            self.constructRender2Drec(self.telomes[child1].child1, self.telomes[child1].child2)
            self.appendTelomeCoordinates2D(self.telomes[child1].base2D)
        if child2 > 0:
            #print(f"child {child2} of parent {self.telomes[child2].parent}, level {self.telomes[child2].level}, of abs {math.degrees(self.telomes[child1].absoluteAngle)}")
            self.appendTelomeCoordinates2D(self.telomes[child2].terminus2D)
            self.constructRender2Drec(self.telomes[child2].child1, self.telomes[child2].child2)
            self.appendTelomeCoordinates2D(self.telomes[child2].base2D)
        return
    def constructRender2D(self):
        self.initializeRender2D(self.telomes[0].base2D)
        self.appendTelomeCoordinates2D(self.telomes[0].terminus2D)
        self.constructRender2Drec(self.telomes[0].child1, self.telomes[0].child2)

    def renderPlant2D(self, withMoments=False, withMomentNumbers=False):
        # Step 1: Collect segment coordinates and line widths
        segments = []
        widths = []
        for telome in self.telomes:
            start = telome.base2D
            end = telome.terminus2D
            segments.append([start, end])
            widths.append(telome.diameter)  # Thickness based on diameter

        fig, ax = plt.subplot()

        # Step 2: Create LineCollection with varying line widths
        line_collection = LineCollection(segments, linewidths=widths, color='#654321', zorder=1)
        ax.add_collection(line_collection)

        # Step 3: Plot scatter if withMoments is True
        if withMoments:
            self.visualizeMoments()
            sc = ax.scatter(self.momentScatterX, self.momentScatterY, c=self.momentScatterMomsAbs, cmap='RdYlGn_r', zorder=2)
            if withMomentNumbers:
                for i, mom in enumerate(self.momentScatterMomsAbs):
                    ax.annotate(mom, (self.momentScatterX[i], self.momentScatterY[i]))
            plt.colorbar(sc, ax=ax)

        # Set title and labels
        ax.set_title("Plant Morphology with Varying Thickness", fontsize=6)
        ax.set_xlabel(f"Levels: {self.levels}   b1: {self.b1deg}   b2: {self.b2deg}", fontsize=6)
        ax.axis("equal")  # Ensures equal scaling


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
        self.telomes[this].baseMoment += W*(self.telomes[this].base2D[0] - Xcom) / (self.telomes[this].radius**4)
        self.telomes[this].baseMomentMax = abs(self.telomes[this].baseMoment)
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
        self.momentScatterX += [self.telomes[this].terminus2D[0] - self.telomes[this].vector2D[0]*0.8]
        self.momentScatterY += [self.telomes[this].terminus2D[1] - self.telomes[this].vector2D[1]*0.8]
        # self.momentScatterMoms += [int(self.telomes[this].baseMoment)]
        self.momentScatterMomsAbs += [abs(int(self.telomes[this].baseMoment))]
        return
    def visualizeMoments(self): # WRAPPER
        self.momentScatterX: list = []
        self.momentScatterY: list = []
        # self.momentScatterMoms: list = []
        self.momentScatterMomsAbs: list = []
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
    def calculateSurfaceAreaToVolume(self, R=False):
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

    # OLD CODE FOR GENERATING TELOMES WITHOUT BRANCHING PROBABILITY
    def generateTelomesR(self, level, parent):
        #self.minLevel = min(level, self.minLevel)
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
        self.minLevel = levels
        self.telomes = [Telome(self, levels, 0, 0, math.radians(90))]
        self.generateTelomesR(levels-1, 0)
        self.constructRender2D()

'''
    def calculate_branching_probability(self, P_prev, N, k):
        # Modified formula to increase branching likelihood
        return (12 * P_prev) / (N + k / 2) #* 1.2
    
    def generateTelomesR(self, level, parent, P_prev):
        # Termination check based on probability
        if level == 0:
            return
        
        # Calculate the branching probability for this level
        N = self.levels  # Initial levels can be used as N
        p1_branch_prob = self.calculate_branching_probability(P_prev, N, level)
        
        # Branch 1
        if random.random() < p1_branch_prob:
            child1 = Telome(self, level, 1, parent, self.telomes[parent].absoluteAngle)
            self.telomes.append(child1)
            self.generateTelomesR(level - 1, len(self.telomes) - 1, p1_branch_prob)
        
        # Calculate the branching probability again for branch 2 (can be different)
        p2_branch_prob = self.calculate_branching_probability(P_prev, N, level)
        
        # Branch 2
        if random.random() < p2_branch_prob:
            child2 = Telome(self, level, 2, parent, self.telomes[parent].absoluteAngle)
            self.telomes.append(child2)
            self.generateTelomesR(level - 1, len(self.telomes) - 1, p2_branch_prob)
    
    def generateTelomes(self):  # Wrapper for generating telomes
        self.telomes = [Telome(self, self.levels, 0, 0, math.radians(90))]
        initial_prob = 1  # Set initial probability of branching
        self.generateTelomesR(self.levels - 1, 0, initial_prob)
        self.constructRender2D()
'''

levels: int = 7
plant = Plant(0.5, 0.5, 0, 0, 15, 20, levels=levels)
plant.generateTelomes(levels)
# plant.computeMoments()
# print(f"{plant.computeHeight()}")
# print(plant.calculateSurfaceAreaToVolume(True))

plant.renderPlant2D(withMoments=False)

plt.show()



    # def renderPlant2D(self):
    #     fig, ax = plt.subplots()
    #     ax.plot(self.x2D, self.y2D, linewidth=2.0)
    #     plt.show()

    # def renderPlant3D(self):
    #     ax = plt.figure().add_subplot(projection='3d')
    #     ax.plot(self.x3D, self.z3D, self.y3D, linewidth=1.5)
    #     plt.show()