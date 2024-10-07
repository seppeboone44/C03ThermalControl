# =============================== THERMAL CONTROL DESIGNER PROGRAM ==================================================
# This program is designed to find thermal control settings for the JARED Venus mission
# Programmed by group C03 for the second year Systems Design course at the faculty of Aerospace Engineering, TU Delft

import matplotlib.pyplot as plt
import math
import numpy as np
import csv

Alphalist = [0.44, 0.25, 0.21, 0.24, 0.14, 0.60, 0.95, 0.37, 0.26, 0.29, 0.12, 0.88, 0.75, 0.40, 0.16, 0.08, 0.07] # Absorptivities corresponding to different materials.
Solutions = [["Daylight temperature [K]", "Night temperature [K]", "Heater power [W]", "Radiator Surface Area [m²]", "Emmissivity [-]", "Absorptivity [-]"]]
TargetSolutions = [["Daylight temperature [K]", "Night temperature [K]", "Heater power [W]", "Radiator Surface Area [m²]", "Emmissivity [-]", "Absorptivity [-]"]]

LowerLimit      = -15+273.15         # [K]
UpperLimit      = +38+273.15         # [K]
target          = +15+273.15         # [K]
tolerance       = 1                  # [K]


# Defining shape variables
    # Cylinder
Cyl_dim = [0.60, 1.80] #radius, length

    #Cube
Cube_dim = [1.26] #side

# Setting up thermal balance function
def Thermalbalance(Shape="Cube", DimSC= list, AreaRad= float, Heatint= list, alphaSC= float, alphaRad= float, epsilonSC= float, epsilonRad= float, heatpower= float):
    """
    Thermalbalance calculates the thermal balance of a spacecraft that has the form of a cube or a cylinder
    and radiator (1) that are approximated by a plate. The variables used should fallow the following:
    Shape       = "Cube", "Cylinder"            This is the shape of the main spacecraft body.
    DimSC       = [x], [r,h], or [x,y]          Referring to the dimensions of the cube, cylinder, or plate respectevly
    DimRad      = float                         The area of the radiator
    Heatint     = [HeatDay, HeatEclipse]        The internally dissipated heat during sunlight, and during eclipse
    alphaSC     = float                         The absorptance of the main spacecraft body
    alphaRad    = float                         The absorptance of the radiators
    epsilonSC   = float                         The emmitance of the main spacecraft body
    epsilonRad  = float                         The emmitance of the radiators
    heatpower   = float                         The power of the heating element (only used during eclipse)
    """

    # Important variables for this funcion for Venus only!
    sigma = 5.670374419*10**(-8)    # stefan Boltzmann constant [Wm**-2K**-4]
    T_eff = 226.6                   # effective surface temperature [K]
    h_orb = 1000*10**3              # orbital altitude [m]
    R_venus = 6051.8*10**3          # radius venus [m]
    J_s = 2601.3                    # solar constant [W/m**2]
    Alb = 0.65                      # albedo

    # For the radiators, assuming the radiators are parallel to the solar rays at all times.
    # and assuming worst (hottest) case conditions for the albedo and IR planetary radiation.
    Heat_Solar_rad  = 0
    Heat_Albedo_rad = J_s*AreaRad*alphaRad*Alb*(R_venus/(h_orb+R_venus))**2
    Heat_IR_rad     = sigma*T_eff**4*(R_venus/(h_orb+R_venus))**2*AreaRad
    
    # For the different spacecraft bodies
    if Shape == "Cube":
        Heat_Solar_SC   = J_s*DimSC[0]**2*alphaSC
        Heat_Albedo_SC  = J_s*DimSC[0]**2*alphaSC*Alb*(R_venus/(h_orb+R_venus))**2
        Heat_IR_SC      = sigma*T_eff**4*(R_venus/(h_orb+R_venus))**2*DimSC[0]**2*epsilonSC
        AreaSC          = 6*DimSC[0]**2

    elif Shape == "Cylinder":
        Heat_Solar_SC   = J_s*2*DimSC[0]*DimSC[1]*alphaSC
        Heat_Albedo_SC  = J_s*2*DimSC[0]*DimSC[1]*alphaSC*Alb*(R_venus/(h_orb+R_venus))**2
        Heat_IR_SC      = sigma*T_eff**4*(R_venus/(h_orb+R_venus))**2*2*DimSC[0]*DimSC[1]*epsilonSC
        AreaSC          = 2*DimSC[0]**2*math.pi+2*math.pi*DimSC[0]*DimSC[1]

    # Add up all effects to find temperature
    T_day    = ((Heatint[0] + Heat_Solar_rad + Heat_Albedo_rad + Heat_IR_rad + Heat_Solar_SC+Heat_Albedo_SC+Heat_IR_SC)/(sigma*(AreaRad*epsilonRad+AreaSC*epsilonSC)))**(1/4)
    T_night  = ((Heatint[1] + Heat_IR_SC)/(sigma*(AreaSC*epsilonSC)))**(1/4)

    return[T_day, T_night]
# We use the funciton to iterate over the available materials, and find out if any of them can keep the temperature
# within the set limits, while iterating over all epsilon values, since these can change depending on insulation,
# ass well as different radiator areas, and heater power values.

Radiatorareas = []
Heaterpowers = []

iteration = 1

for i in range(len(Alphalist)):              # Absorptivity iterations
    for e in np.arange(0.005,0.2,0.005):     # Emittance iterations
        for A in np.arange(0.0,10,0.2):      # Radiator surface area
            for h in range(0,40,5):          # The initial power budget available
                temp = Thermalbalance("Cylinder", Cyl_dim, A, [362,464], Alphalist[i], 0.07, e, 0.74, h)
                print("iteration", iteration,":",[round(temp[0],2), round(temp[1],2), round(h,1), round(A,4), round(e,4), Alphalist[i]])
                iteration += 1
                if temp[0] < UpperLimit and temp[0] > LowerLimit and temp[1] < UpperLimit and temp[1] > LowerLimit:
                    Solutions.append([round(temp[0],2), round(temp[1],2), round(h,1), round(A,4), round(e,4), Alphalist[i]])              # Stores all possible solutions
                    if abs(temp[0]-target)<tolerance and abs(temp[1]-target)<tolerance:
                        TargetSolutions.append([round(temp[0],2), round(temp[1],2), round(h,1), round(A,4), round(e,4), Alphalist[i]])    # Stores only solutions that are very close to the target temperature.
                        Radiatorareas.append(round(A,2))
                        Heaterpowers.append(round(h,1))


# Writing the solutions to csv files, this is what can be used to choose a design point.
with open("Solutions.csv",'w', newline = '') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(Solutions)
csvfile.close()

with open("EqualSolutions.csv", 'w', newline = '') as newcsv:
    writer2 = csv.writer(newcsv)
    writer2.writerows(TargetSolutions)
newcsv.close()

# The most important variables are the radiator surface area (since this adds mass) and the heater power since the heating element uses power from the battery.
# Therefore we want to know which design point minimizes both values, and then we take the surface properties (alpha and epsilon) for this desing point.
# We plot these values in a graph for visual reference.
plt.scatter(Radiatorareas, Heaterpowers)
plt.title("Design space - Radiator Area vs Heater Power for a thermal equilibrium of 15C")
plt.xlabel("Radiator Area [m²]")
plt.ylabel("Heater Power [W]")
plt.show()