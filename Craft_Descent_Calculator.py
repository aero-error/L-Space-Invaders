"""
L'SPACE INVADERS
Mars Descent Profile Calculator

Written by Michael Gromski

NOTES:
    CALCULATOR ASSUMES THESE FACTORS:
        Spacecraft is initally in a perfect circular orbit
        Rotation of planet is negligible
        Landing site is at "sealevel"
        No weather on Mars (wind is not accounted)
        Thermodynamic approximates the medium as 100% CO2
        
REF:
    1.https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
    2.http://www.braeunig.us/space/orbmech.htm
    3.https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19640016000.pdf
    4.https://www.faa.gov/about/office_org/headquarters_offices/avs/offices/aam/cami/library/online_libraries/aerospace_medicine/tutorial/media/iii.4.1.7_returning_from_space.pdf
    5.https://www.grc.nasa.gov/WWW/K-12/airplane/sound.html
   *6.https://apps.dtic.mil/dtic/tr/fulltext/u2/a505342.pdf
   
   * direct link to large .pdf

FEATURES TO ADD:
    Atmospheric heating (need to find data on thermal conductivity)
    More accurate Coefficient of lift/drag (need more data on properties of spacecraft)
 X  Verify EOM accuracy (figure out what is wrong with EOM)
 X  Maximum deceleration (not too bad)
    Include multiple stages (currently only ballistic entry)
    Speed of sound/Mach calculations (need gas properties at low pressures)
    Improve user interface (restructure code to be easily interpreted)
"""
import numpy as np
import matplotlib.pyplot as plt

########## INPUTS ##########
#Constants
G = 6.67408e-11 #Gravitaional constant

#Properties of the spacecraft
CraftMass = 100 #Mass in KG
CoD = 1.4       #Coefficient of drag
CoL = 0.2       #Coefficent of lift
CrossArea = 0.4 #Cross-sectional area of spacecraft M^2

altitude = 400e3 #Starting altitude in M #400e3
V0 = 2.38e3 #Inital velocity in M/S
angle = 1 #entry angle in degrees (Below horizontal)

#Properties of Planet (Mars)
MarsMass = 6.39e23 #Mass in KG
MarsRadius = 3.3895e6 #Radius in M

########## FUNCTIONS ##########
def Gravity(altitude): #returns acceleration due to gravity in meters/sec (EQNS REF:2)
    total_alt = altitude + MarsRadius
    GRAV = G*MarsMass / total_alt**2
    return GRAV

def PressAtm(altitude): #returns atmospheric pressure in K-Pascals (EQNS REF:1)
    if altitude > 7000: 
        T = -23.4 - 0.00222*altitude
        P = 0.699 * np.exp(-0.00009*altitude) 
    else:     
        T = -31 - 0.000998*altitude
        P = 0.699 * np.exp(-0.00009*altitude)
    Density = P / (0.1921 * (T + 273.15))
    if T < -273.1:
        T = -273.1
    if Density < 0:
        Density = 1e-10
    return P, T, Density

### DATA STORAGE ARRAYS ###
TIME = [] #Seconds
ALT = [] #Meters
VEL = [] #M/Sec
ANG = [] #degrees
QHEAT = [] #Joules
PRESSURE = [] #Pascals
GRAVITY = [] #M/s^2
AIR_DENSITY = [] #KG/M^3
TEMP = [] #Kelvin
Gs = [] #no units

########## MAIN CODE ###########
V = V0
dist_horiz = 0
angle = -angle * (np.pi/180) #Convert degrees to rad
time = 0 #sec

### MAIN LOOP ###
gamma = angle
v = V0
dt = 0.01
B = altitude
m = CraftMass
S = CrossArea
r = MarsRadius
h = altitude
t = 0
while h > 0:
    y = np.exp(-h/B)
    g = Gravity(h)
    p, T, rho = PressAtm(h)
    k1 = (rho * S * CoD) / (2*m)
    k2 = (rho * S * CoL) / (2*m)
    dh = v*np.sin(gamma)*dt
    dv = (-k1*(v**2)*y - g*np.sin(gamma))*dt
    d_gamma = (k2*v*y - (g/v)*np.cos(gamma*(1-((v**2)/(g*r)))))*dt
    h += dh
    v += dv
    gamma += d_gamma
    
    dv_horiz = v * np.cos(gamma) 
    dist_horiz += dv_horiz * dt

    Gs.append(-dv/(9.81*dt)) #Note these are Earth Gs
    TIME.append(t)
    VEL.append(v)
    ANG.append(gamma * (180/np.pi))
    ALT.append(h)
    
    t += dt


########### OUTPUT ##########
print("*"*50)
print("RESULTS:")
print()    
print("The descent will take %.2f seconds" %(t))
print("The craft will travel %.2f meters horizontally." %(dist_horiz))

plt.figure(1)
plt.scatter(TIME, VEL)
plt.title("Velocity of SpaceCraft")
plt.ylabel("Velocity (Meters/Sec)")
plt.xlabel("Time (Sec)")

plt.figure(2)
plt.scatter(TIME, ANG)
plt.title("Angle from Horizontal of Nominal Orbit")
plt.ylabel("Angle (Degrees)")
plt.xlabel("Time (Sec)")

plt.figure(3)
plt.scatter(TIME, ALT)
plt.title("Altitude of SpaceCraft")
plt.ylabel("Altitude (Meters)")
plt.xlabel("Time (Sec)")    

plt.figure(4)
plt.scatter(TIME, Gs)
plt.title("Deceleration of Spacecraft")
plt.ylabel("Deceleration (Earth Gs)")
plt.xlabel("Time (Sec)")

