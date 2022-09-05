from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np

g = 32.1740 # lbf of gravity
rollingMu = .45
EstimatedTime = 120 # seconds

class Performance:
    
    def __init__(self, thrustData, cruiseVel, maxVelocity, mVolt, motorAmps, cells, batAmps, batVolt, weight, CL_df, CD_df, Sref):
        self.thrustData = thrustData
        self.cruiseVel = cruiseVel
        self.maxVelocity = maxVelocity
        self.mVolt = mVolt
        self.motorAmps = motorAmps
        self.cells = cells
        self.batAmps = batAmps
        self.batVolt = batVolt
        self.weight = weight
        self.CL_df = CL_df
        self.CD_df = CD_df
        self.Sref = Sref
        
        volts = 4.2 * cells
        self.motorWatts = motorAmps * mVolt
        self.battWatts = batAmps * batVolt
        
        self.thrustArray = np.array([],[])
        self.pwrWat = np.array([],[])
        
        for i in thrustData.head().index:
            np.concatenate(self.thrustArray, [i, [thrustData[i]["Thrust"]]])
            np.concatenate(self.pwrWat, [i, [thrustData[i]["Power (W)"]]])
        
        # boo, volts = cruise(cruiseVel, EstimatedTime)
        time, Usedvolts = self.takeOff(CL_df[0], CD_df[0], V, To)
        
        if Usedvolts > volts:
            raise Exception("This plane will die before it flies")
            
        print("Time: " + str(time))
        print("Take off Power Use: " + str(Usedvolts))
        print("Useable Volts Left: " + str(volts - Usedvolts))
        
    
    # def cruise(cruiseVel, EstimatedTime):
        
    #     volts = 4.2 * cells
    #     motorWatts = motorAmps * mVolt
    #     battWatts = batAmps * batVolt
        
    #     if motorWatts < battWatts:
    #         wattUse = motorWatts
    #     else:
    #         wattUse = battWatts
    #     Actualtime = 0
    #     while volts < 3.2:
    #         volts -= wattUse
    #         Actualtime += 1
    #     #Wh = motorAmps * volts
        
    #     return EstimatedTime > Actualtime, volts
    
    def takeOff(self, Clg, Cdg, V, To):
        Vlof, thrust = self.takeOffThrustVelocity()
        assert(thrust != None)
        
        
        h = 0
        drag = 0
        lift = 0
        time = 0
        v = Vlof
                
        # LBS of thrust/drag when first start spinning 
        pwrUse = interp1d(self.pwrWat[0], self.pwrWat[1])
        tUse = interp1d(self.thrustArray[0], self.thrustArray[1])
        
        totPwrWatts = 0
    
        
        while (h < 200):
            v += thrust - drag
            h += lift - self.weight
            
            # LBS of thrust/drag at v
            
            totPwrWatts += pwrUse(v)
            thrust = tUse(v)
            drag = self.CD_df[v]
            
            # Lift at v
            lift = self.CL_df[v] #this assumes CL_df is 
            time += 1
        t2 = self.timeGroundRun(To, Clg, Cdg, Vlof, self.Sref, V, thrust)
        time += t2
        return time, totPwrWatts
    
    # def groundDistTakeOff(self, StaticThrust, Clg, Cdg, Vlof, T):
    #     """ Obsolete """
    #     # Vlof is velocity where lift > weight
    #     # thrust at Vlof/sqrt(2)
        
    #     D = self.CD_df[Vlof] # drag at Vlof/sqrt(2)
    #     L = self.CL_df[Vlof] # lift at Vlof/sqrt(2)
        
    #     Sq = (Vlof / 2**.5) / (g / self.weight * (T - D - rollingMu * (self.weight - L)))
        
    #     return Sq
        
    def timeGroundRun(self, To, Clg, Cdg, Vlof, V, T):
        """ Time and Distance of Ground Run to Takeoff """
        massDensity = .0000416666 #kbs/ft^3
        
        a = (T-To)/V^2 # this is a constant determined by thrust stuff
        
        # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://www.dept.aoe.vt.edu/~lutze/AOE3104/takeoff&landing.pdf
        A = self.g * (To / self.weight - self.rollingMu)
        B = g / self.weight * (1/2 * massDensity * self.Sref * (Cdg - rollingMu * Clg) + a)
        
        t = 1/(2 *(A * B)**.5) * np.ln((A**.5 + Vlof * B**.5) / (A**.5 - Vlof * B**.5))
        S = 1 / (2*B)* np.ln(A / (A - B*Vlof**2) )

        
        return t, S
    
    def takeOffThrustVelocity(self):
        velocity = 0
        thrust = 0
        drag = 0
        lift = 0
        
        for i in self.thrustData.head().index:
            try:
                thrust = self.thrustData[i]["Thrust (N)"]
            except:
                print("Thrust Req outside of Allowable Range")
                return None, None
            drag = self.CD_df[i]
            lift = self.CL_df[i]
            velocity = i
            
            if (thrust > drag and lift > self.weight*1.2):
                return velocity, thrust*1.3
                # chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.worldscientific.com/doi/pdf/10.1142/S2010194516601745
        return None, None