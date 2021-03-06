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
        All rarefied fluids are approximated as ideal gases
        
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
 X  Include multiple stages (currently only ballistic entry)
 X  Speed of sound/Mach calculations (need gas properties at low pressures)
    Improve user interface (restructure code to be easily interpreted)
"""
import numpy as np
import matplotlib.pyplot as plt

########## INPUTS ##########
altitude = 350e3 #Starting altitude in M #400e3
V0 = 2.42e3 #Inital velocity in M/S
angle = 10 #entry angle in degrees (Below horizontal)
Deploy1 = 30000 #altitude to deploy parachute
Deploy2 = 2000  #altitude to stage parachute and drop heat shield

#Properties of the spacecraft
### STAGE 1 ###
CraftMass1 = 150 #Mass in KG
CoD1 = 1.1       #Coefficient of drag
CoL1 = 0.5       #Coefficent of lift
CrossArea1 = 0.4 #Cross-sectional area of spacecraft M^2

### STAGE 2 ###
CraftMass2 = 150 #Mass in KG
CoD2 = 1.6       #Coefficient of drag
CoL2 = 0.64       #Coefficent of lift
CrossArea2 = 1.4 #Cross-sectional area of spacecraft M^2

### STAGE 3 ###
CraftMass3 = 135 #Mass in KG
CoD3 = 1.75       #Coefficient of drag
CoL3 = 0.8       #Coefficent of lift
CrossArea3 = 3.14 #Cross-sectional area of spacecraft M^2


#Properties of Planet (Mars)
MarsMass = 6.39e23 #Mass in KG
MarsRadius = 3.3895e6 #Radius in M

#Constants
G = 6.67408e-11 #Gravitaional constant
R = 188.92 #Gas Constant of CO2 (J/(kg*K))
GAMMA = 1.28 #Specific heat of CO2 at 300k (NOTE: WILL NEED TO CHANGE THIS LATER)

#Program Specific
testing = False #Mark True if you want to see all graphs

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
    
    denominator = (0.1921 * (T + 273.15)) #note: this is a rough patch to fix density issiue 
    if denominator < 10 and altitude > 100000:
        denominator = 100
    Density = P / denominator
    
    if T < -273.1:
        T = -273.1
        
    if Density < 0:
        Density = 1e-10
    return P, T, Density

def Mach(T): #returns the speed of sound (m/s) (NEEDS WORK ON GAMMA)
    T += 273.15
    M = np.sqrt(GAMMA*R*T)
    return M

def VelocityCheck(): #This is not really needed but is a good diagonostic tool
    V_orbit = np.sqrt((G*MarsMass)/(altitude+MarsRadius))
    alt_optimal = ((G*MarsMass)/(V0**2)) - MarsRadius
    print("Testing Mode:")
    print()
    print("The optimal circular orbit for this altitude is %6.2f meters." %(alt_optimal))
    if V_orbit > V0:
        print("The spacecraft is too slow for the set altitude. It will enter the planet.")
    elif V_orbit < V0:
        print("The spacecraft is too fast for the orbit. It is not a circular orbit or it will miss the planet.")
    elif V_orbit == V0:
        print("The spacecraft is in a perfect circular orbit.")
    print()
    return

def Calc_Cf(): #calculates average skin coefficient (which can be approximated as twice the Stanton number) (NEEDS WORK)
    cf = (2*h_star)/(Cp*B*V)
    return cf

### DATA STORAGE ARRAYS ###
TIME = [] #Seconds
ALT = [] #Meters
Y_POS = [] #Meters (Testing Loop)
VEL = [] #M/Sec
ANG = [] #degrees
QHEAT = [] #Joules
PRESSURE = [] #Pascals  (Testing Loop)
GRAVITY = [] #M/s^2   (Testing Loop)
AIR_DENSITY = [] #KG/M^3  (Testing Loop)
TEMP = [] #Kelvin   (Testing Loop)
Gs = [] #no units
MACH = [] #no units

########## MAIN CODE ###########
V = V0
dist_horiz = 0
angle = -angle * (np.pi/180) #Convert degrees to rad
time = 0 #sec

### TESTING LOOP ###
if testing == True:
    VelocityCheck()
    for i in range(0,int(altitude)):
        pres, temp, density = PressAtm(i)
        grav = Gravity(i)
        PRESSURE.append(pres*1000)
        TEMP.append(temp+273.15)
        AIR_DENSITY.append(density)
        GRAVITY.append(grav)
        Y_POS.append(i)
        MACH.append(Mach(temp))
    

### MAIN LOOP ###
gamma = angle
v = V0
dt = 0.01
B = altitude
r = MarsRadius
h = altitude
t = 0
while h > 0:
    if h > Deploy1:
        S = CrossArea1
        CoL = CoL1
        CoD = CoD1
        m = CraftMass1
    elif h > Deploy2 and h <= Deploy1:
        S = CrossArea2
        CoL = CoL2
        CoD = CoD2
        m = CraftMass2
    elif h <= Deploy2:
        S = CrossArea3
        CoL = CoL3
        CoD = CoD3
        m = CraftMass3
        
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
print("The craft will impact Mars at %.2f m/s." %(VEL[len(VEL)-1]))

plt.figure(1)
plt.plot(TIME, VEL, "-k")
plt.title("Velocity of SpaceCraft")
plt.ylabel("Velocity (Meters/Sec)")
plt.xlabel("Time (Sec)")
plt.grid(True)

plt.figure(2)
plt.plot(TIME, ANG, "-k")
plt.title("Angle from Horizontal of Nominal Orbit")
plt.ylabel("Angle (Degrees)")
plt.xlabel("Time (Sec)")
plt.grid(True)

plt.figure(3)
plt.plot(TIME, ALT, "-k")
plt.title("Altitude of SpaceCraft")
plt.ylabel("Altitude (Meters)")
plt.xlabel("Time (Sec)")
plt.grid(True)

plt.figure(4)
plt.plot(TIME, Gs, "-r")
plt.title("Deceleration of Spacecraft")
plt.ylabel("Deceleration (Earth Gs)")
plt.xlabel("Time (Sec)")
plt.grid(True)

### TESTING PLOTS ###
if testing == True:
    plt.figure(5)
    plt.plot(Y_POS,GRAVITY,color="orange")
    plt.title("Mars Gravity with Respect to Altitude")
    plt.xlabel("Altitude (Meters)")
    plt.ylabel("Gravity (M/S^2)")
    plt.grid(True)
    
    plt.figure(6)
    plt.plot(Y_POS,TEMP,color="orange")
    plt.title("Temperature with Respect to Altitude")
    plt.xlabel("Altitude (Meters)")
    plt.ylabel("Temperature (Kelvin)")
    plt.grid(True)
    
    plt.figure(7)
    plt.plot(Y_POS,PRESSURE,color="orange")
    plt.title("Mars Atmospheric Pressure with Respect to Altitude")
    plt.xlabel("Altitude (Meters)")
    plt.ylabel("Pressure (Pascals)")
    plt.grid(True)
    
    plt.figure(8)
    plt.plot(Y_POS,AIR_DENSITY,color="orange")
    plt.title("Mars Atmospheric Pressure with Respect to Altitude")
    plt.xlabel("Altitude (Meters)")
    plt.ylabel("Air Density (Kg/M^3)")
    plt.grid(True)
    
    plt.figure(9)
    plt.plot(Y_POS,MACH,color="orange")
    plt.title("Mach Speed with Respect to Altitude")
    plt.xlabel("Altitude (Meters)")
    plt.ylabel("Velocity (M/S)")
    plt.grid(True)
