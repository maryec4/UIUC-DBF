#import numpy as np
#import pandas as pd

def Performance(lift, drag, thrustData, maxVelocity, mVolt, motorAmps, cells, batAmps, batVolt):
    
    def cruise(cruiseVel, EstimatedTime):
        
        volts = 4.2 * cells
        motorWatts = motorAmps * mVolt
        battWatts = batAmps * batVolt
        
        wattUse = max(motorWatts, battWatts)    
        Actualtime = 0
        while volts < 3.2:
            volts -= wattUse
            Actualtime += 1
        #Wh = motorAmps * volts
        
        return EstimatedTime > Actualtime
    
    def takeOff(maxVelocity):
        
        
        print("yee")
        