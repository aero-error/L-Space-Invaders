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
    Verify EOM accuracy (figure out what is wrong with EOM)
    Maximum deceleration (not too bad)
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

altitude = 400e3 #Starting altitude in M
V0 = 2.38e3 #Inital velocity in M/S
angle = 10 #entry angle in degrees (Below horizontal)
DForce = 500 #Deceleration force in newtons (not used yet)

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

def VelocityCheck(): #This is not really needed but is a good diagonostic tool
    V_orbit = np.sqrt((G*MarsMass)/(altitude+MarsRadius))
    alt_optimal = ((G*MarsMass)/(V0**2)) - MarsRadius
    print("The optimal circular orbit for this altitude is %6.2f meters." %(alt_optimal))
    if V_orbit > V0:
        print("The spacecraft is too slow for the set altitude. It will enter the planet.")
    elif V_orbit < V0:
        print("The spacecraft is too fast for the orbit. It is not a circular orbit or it will miss the planet.")
    elif V_orbit == V0:
        print("The spacecraft is in a perfect circular orbit.")
    return

def Calc_Cf(): #calculates average skin coefficient (which can be approximated as twice the Stanton number)
    cf = (2*h_star)/(Cp*B*V)
    return cf

### DATA STORAGE ARRAYS ###
TIME = [] #Seconds
ALT = [] #Meters
alt0 = [] #Meters
VEL = [] #M/Sec
ANG = [] #degrees
QHEAT = [] #Joules
PRESSURE = [] #Pascals
GRAVITY = [] #M/s^2
AIR_DENSITY = [] #KG/M^3
TEMP = [] #Kelvin
DIST_HORIZ = [] #Meters

'''
SOUND_VEL = [] #M/Sec #Feature yet to be added
MACH = [] #Unitless
'''

########## MAIN CODE ###########
V = V0
dist_horiz = 0
h = altitude
alt = altitude
angle = -angle * (np.pi/180) #Convert degrees to rad
time = 0 #sec
q_heat = 0 # This is the ammount of heat energy applied

VelocityCheck() #Mostly Used as a diagonostic tool for determining orbit speed

### MAIN LOOP ###

while alt > 0: #Used to visualize properites of atmosphere at altitude
    P,T,B = PressAtm(alt)
    PRESSURE.append(P*1000)
    TEMP.append(T+273.1)
    AIR_DENSITY.append(B)
    alt0.append(alt)
    alt -= 100

while h > 0: #stops when craft has reached the surface   
    P,T,B = PressAtm(h)
    #Cf = Calc_Cf()
    P = P*1000 #convert kPa to Pa
    T = T + 273.1 #Convert C to K
    g = Gravity(h)
    K1 = (B*CrossArea*CoD) / (2*CraftMass) #EQNS REF:3 EQN:3a
    K2 = (B*CrossArea*CoL) / (2*CraftMass) #EQNS REF:3 EQN:3b
    vel_horiz = V*np.cos(angle) #find horizontal velocity    
    h_dot = V * angle #Units(M/Sec) EQNS REF:3 EQN:1
    v_dot = -K1*(V**2)*np.exp(-h/B) - g*angle #Units(M/Sec^2) EQNS REF:3 EQN:2
    angle_dot = K2*(V**2)*np.exp(-h/B) - g*(1-(V**2)/(g*MarsRadius))
    #q_dot = (1/4)*Cf*B*(V**3)
    #q_heat = q_heat + q_dot
    angle_dot = angle_dot/V
    h = h + h_dot
    V = V - v_dot
    angle = angle + angle_dot
    dist_horiz += vel_horiz*1

    
    ### DATA ACQUISITION ###
    TIME.append(time)
    #PRESSURE.append(P)
    GRAVITY.append(g)
    #AIR_DENSITY.append(B)
    #TEMP.append(T)
    ALT.append(h)
    VEL.append(V)
    ANG.append(angle*180/np.pi)
    #QHEAT.append(q_heat)
    DIST_HORIZ.append(dist_horiz)
    time += 1


########## DISPLAY DATA ##########
print("*"*50)
print("RESULTS:")
print()
print("Total horizontal distance covered: %6.2f meters" %(dist_horiz))
'''
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

plt.figure(8)
plt.scatter(TIME, DIST_HORIZ)
plt.title("Horizontal distance")
plt.ylabel("Distance (Meters)")
plt.xlabel("Time (Sec)")
'''
plt.figure(9)
plt.scatter(alt0, TEMP)
plt.title("Temperature at Altitude")
plt.ylabel("Temperature (Kelvin)")
plt.xlabel("Altitude (M)")

plt.figure(10)
plt.scatter(alt0, PRESSURE)
plt.title("Atmospeheric Pressure at Altitude")
plt.ylabel("Pressure (Pa)")
plt.xlabel("Altitude (M)")

plt.figure(11)
plt.scatter(alt0, AIR_DENSITY)
plt.title("Atmospeheric Density at Altitude")
plt.ylabel("Density of Atmosphere (Kg/M^3)")
plt.xlabel("Altitude (M)")

