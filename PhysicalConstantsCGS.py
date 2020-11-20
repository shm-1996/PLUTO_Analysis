import numpy as np
#Lengths
AstronomicalUnit = 1.49597870691e+13
Parsec = 206265.0 * AstronomicalUnit
LightYear = 63241.0 * AstronomicalUnit
SolarRadius = 6.955e+10
EarthRadius = 637279760.0
JupiterRadius = 7.149e+9

#Time

Minute = 60.0
Hour = 60.0 * Minute
Day = 24.0 * Hour
Year = 365.25 * Day

#Masses

SolarMass   = 1.9891e+33
EarthMass   = 1.0 /   332946.0  * SolarMass
LunarMass   = 1.0 / 27068510.0  * SolarMass
JupiterMass = 1.0 /     1047.56 * SolarMass

#Radiation

c_SpeedOfLight = 2.99792458e+10
k_BoltzmannConstant = 1.3807e-16
h_PlanckConstant = 6.62606896e-27
R_UniversalGasConstant = 8.314472e+7
NA_AvogadroConstant = 6.02214129e+23
eV = 1.0 / 624150974451.150

#Ionization

mC_CarbonMass             = 12.0 / NA_AvogadroConstant
u_AtomicMassUnit          =  1.0 / NA_AvogadroConstant
mH_HydrogenMass           = 1.00794         * u_AtomicMassUnit
me_ElectronMass           = 5.4857990946e-4 * u_AtomicMassUnit

#Gravity

G_GravityConstant         = 6.67428e-8
g_AccelerationOfGravity   = 980.665
SolarLuminosity           = 3.839e+33

