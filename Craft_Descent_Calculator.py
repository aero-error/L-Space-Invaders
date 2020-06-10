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
        
REF:
    1.https://www.grc.nasa.gov/WWW/K-12/airplane/atmosmrm.html
    2.http://www.braeunig.us/space/orbmech.htm
    3.https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19640016000.pdf
    4.https://www.faa.gov/about/office_org/headquarters_offices/avs/offices/aam/cami/library/online_libraries/aerospace_medicine/tutorial/media/iii.4.1.7_returning_from_space.pdf
    5.https://www.grc.nasa.gov/WWW/K-12/airplane/sound.html
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

altitude = 400e3 #Starting altitude in M
V0 = 2.38e3 #Inital velocity in M/S
angle = 10 #entry angle in degrees (Below horizontal)
DForce = 500 #Deceleration force in newtons (not used yet)

#Properties of Planet (Mars)
MarsMass = 6.39e23 #Mass in KG
MarsRadius = 3.3895e6 #Radius in M


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

def VelocityCheck(): #This is not really needed but is a good diagonostic tool
    V_orbit = np.sqrt((G*MarsMass)/(altitude+MarsRadius))
    if V_orbit > V0:
        print("The space craft is too slow for the set altitude. It will enter the planet.")
    if V_orbit < V0:
        print("The space craft is too fast for the orbit. It is not a circular orbit or it will miss the planet.")
    else:
        print("The spacecraft is in the perfect circular orbit.")
    return

#DATA STORAGE ARRAYS
TIME = [] #Seconds
ALT = [] #Meters
VEL = [] #M/Sec
ANG = [] #degrees
QHEAT = [] #Joules
PRESSURE = [] #Pascals
GRAVITY = [] #M/s^2
AIR_DENSITY = [] #KG/M^3
TEMP = [] #Kelvin

'''
SOUND_VEL = [] #M/Sec #Feature yet to be added
MACH = [] #Unitless
'''

########## MAIN CODE ###########
V = V0
h = altitude
angle = -angle * (np.pi/180) #Convert degrees to rad
time = 0 #sec
q_heat = 0 # This is the ammount of heat energy applied

#VelocityCheck() #Mostly Used as a diagonostic tool for determining orbit speed

### MAIN LOOP ###
while h > 0: #stops when craft has reached the surface   
    P,T,B = PressAtm(h)
    P = P*1000 #convert kPa to Pa
    T = T + 273.1 #Convert C to K
    g = Gravity(h)
    K1 = (B*CrossArea*CoD) / (2*CraftMass) #EQNS REF:3 EQN:3a
    K2 = (B*CrossArea*CoL) / (2*CraftMass) #EQNS REF:3 EQN:3b
    h_dot = V * angle #Units(M/Sec) EQNS REF:3 EQN:1
    v_dot = -K1*(V**2)*np.exp(-h/B) - g*angle #Units(M/Sec^2) EQNS REF:3 EQN:2
    angle_dot = K2*(V**2)*np.exp(-h/B) - g*(1-(V**2)/(g*MarsRadius))
    angle_dot = angle_dot/V
    h = h + h_dot
    V = V - v_dot
    angle = angle + angle_dot
    
    ### DATA ACQUISITION ###
    TIME.append(time)
    PRESSURE.append(P)
    GRAVITY.append(g)
    AIR_DENSITY.append(B)
    TEMP.append(T)
    ALT.append(h)
    VEL.append(V)
    ANG.append(angle*180/np.pi)
    time += 1


########## DISPLAY DATA ##########

plt.figure(1)
plt.scatter(TIME, TEMP)
plt.title("Atmospheric Temperature")
plt.ylabel("Temperature (Kelvin)")
plt.xlabel("Time (Sec)")

plt.figure(2)
plt.scatter(TIME, PRESSURE)
plt.title("Atmospheric Pressure")
plt.ylabel("Pressure (Pascal)")
plt.xlabel("Time (Sec)")

plt.figure(3)
plt.scatter(TIME, AIR_DENSITY)
plt.title("Atmospheric Density")
plt.ylabel("Density (KG/M^3)")
plt.xlabel("Time (Sec)")

plt.figure(4)
plt.scatter(TIME, GRAVITY)
plt.title("Gravitational Force")
plt.ylabel("Gravity (M/S^2)")
plt.xlabel("Time (Sec)")

plt.figure(5)
plt.scatter(TIME, ALT)
plt.title("Altitude of SpaceCraft")
plt.ylabel("Altitude (Meters)")
plt.xlabel("Time (Sec)")

plt.figure(6)
plt.scatter(TIME, VEL)
plt.title("Velocity of SpaceCraft")
plt.ylabel("Velocity (Meters/Sec)")
plt.xlabel("Time (Sec)")

plt.figure(7)
plt.scatter(TIME, ANG)
plt.title("Angle from Horizontal of Nominal Orbit")
plt.ylabel("Angle (Degrees)")
plt.xlabel("Time (Sec)")





