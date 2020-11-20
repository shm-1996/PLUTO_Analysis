import PhysicalConstantsCGS as const
import numpy as np

################## FUNDAMENTAL UNITS ##############################################################################
unit_Length = const.Parsec      # Unit of length is 1 parsec

unit_Density = 1.e-20           # Unit of density is 1e-20 gcm^-3

unit_Time = const.Year*1.e6     # Unit of Time is 1 Myr

################## DERIVED UNITS ##################################################################################
unit_Mass = unit_Density*(unit_Length**3)
unit_Velocity = (unit_Length/unit_Time) #Unit of velocity in cm/s
unit_Energy = unit_Velocity*unit_Velocity*unit_Mass #Unit of unit_Energy in ergs per g 
unit_Pressure = unit_Density*(unit_Velocity**2)
unit_MagneticField = unit_Velocity*np.sqrt(4*np.pi*unit_Density)
unit_GravPotential = unit_Energy/unit_Mass
