import numpy as np
import pandas as pd
#rom openpyxl import load_workbook
import os
import re

def PropulsionCalc(fileName, Throttle, maxVel, drag, eCurrMax, eRes, cells, batVolt, kv, batRes, nLCurr, mRes, Density = 0.0023769, alt = 0):
    
    def RPMfinder(motorAmpGuess, batAmpGuess, Output, velocity):

        # Main Calc for Iteration of pMotorOut
        def pMotorOutCalc(motorAmp, batAmp):

            vBat = cells * batVolt - cells * batRes * batAmp         # Total Bat Volts w/ Volt Loss
            vTerm = (vBat - eRes * batAmp) * Throttle / 100          # Volt Drop from esc/Throttle Val - voltage at terminal
            vDropMotor = motorAmp * mRes                             # Volt Drop from Motor
            vEMF = vTerm - vDropMotor                                # Voltage Motor after losses
            pMotorOut = vEMF * (motorAmp - nLCurr)                   # Power Motor Generates (Voltage of Motor * (Amps of Motor - No Load Motor Amps))
            pMotorIn = vTerm * batAmp                            # Power Motor Receives

            return pMotorIn, pMotorOut, vTerm, motorAmp, vBat, batAmp, vEMF

        pMotorIn, pMotorOut, vTerm, motorAmp, vBat, batAmp, vEMF = pMotorOutCalc(motorAmpGuess, batAmpGuess)
        #pMotorOut = 1210.83

        RPM = vEMF * kv                # Revolutions/Min of Motor
        RPS = RPM / 60                 # Revolutions/s of Motor

        # Find Closest Data to RPM Values, above and below
        RPMgroups = np.array(list(RPMdf.groups))

        if (RPM > RPMgroups.max()):
            return pMotorIn, pMotorOut, 0, vTerm, motorAmp, vBat, batAmp, 0, vEMF, 0, 0

        bottomRPM = RPMgroups[RPMgroups < RPM].max()
        topRPM = RPMgroups[RPMgroups > RPM].min()

        # Get Bottom/Top RPM Groups
        bottomRPMdata = RPMdf.get_group(bottomRPM).reset_index()
        topRPMdata = RPMdf.get_group(topRPM).reset_index()

        velocityMPH = velocity * 0.68181818      # Conversion to MPH

        # Finding Relevant Indices
        bottomRPMClosestVel = abs(bottomRPMdata['V'] - velocityMPH).idxmin()
        topRPMClosestVel = abs(topRPMdata['V'] - velocityMPH).idxmin()

        # Finding Ct
        bottomCt = bottomRPMdata['Ct'].iloc[bottomRPMClosestVel]  # Coefficeient of Thrust
        topCt = topRPMdata['Ct'].iloc[topRPMClosestVel]

        # Finding Cp
        bottomCp = bottomRPMdata['Cp'].iloc[bottomRPMClosestVel]
        topCp = topRPMdata['Cp'].iloc[topRPMClosestVel]

        # Finding Torque
        bottomTorque = bottomRPMdata['Torque'].iloc[:,0].iloc[bottomRPMClosestVel]
        topTorque = topRPMdata['Torque'].iloc[:,0].iloc[topRPMClosestVel]

        # Interpolate to get Coefficient of Power/Thrust for Motor RPM
        Cp = bottomCp + (RPM - bottomRPM) * (topCp - bottomCp) / (topRPM - bottomRPM)
        Ct = bottomCt + (RPM - bottomRPM) * (topCt - bottomCt) / (topRPM - bottomRPM)
        torque = bottomTorque + (RPM - bottomRPM) * (topTorque - bottomTorque) / (topRPM - bottomRPM)

        # Prop Power/Thrust Calc
        pProp = Density * RPS**3 * propDia**5 * Cp       # Power Prop required ft-lb/s
        ThrustProp = Density * RPS**2 * propDia**4 * Ct             # Thrust Prop generates in lbf
        # print('pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, torque, Cp, Ct')
        # print(pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, torque, Cp, Ct)
        return pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, torque, Cp, Ct

    # Pulling Prop Dia and Type
    firstword = "PER3_"
    divider = "x"

    #rx_to_last = r'^.*{}'.format(re.escape(divider))

    #remove PER3_ from filename
    removeBeginning = re.sub(r'^.*?'+firstword, '', fileName)

    #remove .dat from filename
    propType = removeBeginning[:len(removeBeginning) - 4]
    
    #remove paranthesis from propType, messes up regular expressions
    propType = re.sub(r'\([^)]*\)', '', propType)
    
    propPitch = re.sub(r'^.*?'+divider, '', propType)
    
    propDia = re.sub(divider+propPitch, '', propType)
    propDia = int(propDia)
    propDia *= 0.08333333333
    
    
    # Set Column Index for Importing Data
    header_list = [str(x) for x in range(13)]

    # Read file
    try:
        propData = pd.read_csv(fileName, sep = '\s+', names = header_list)
    except:
        header_list = [str(x) for x in range(15)]
        propData = pd.read_csv(fileName, sep = '\s+', names = header_list)
        
    propData = propData.to_numpy()
    df = pd.DataFrame(propData)
    df.to_csv(index = False)

    # Find All Rows w/ RPM Data
    searchString = "RPM"
    mask = np.column_stack([df[col].str.contains(r'' + searchString, na = False)
                            for col in df])
    RPMdata = pd.DataFrame.dropna(df.loc[mask.any(axis = 1)], axis = 1)
    
    
    
    """
    Create New Table to bin RPM Values into Groups
    Remove Characters from PROP RPM
    Change to Numeric Values for Analysis
    """

    # First Find Column w/ Numbers
    RPMs = RPMdata.apply(pd.to_numeric, errors = 'ignore')

    # Second Remove All non-int64 Columns
    RPMs = RPMs.select_dtypes(include = ['int64'])

    # Third Get First RPM Value Index
    a = RPMs.columns.values
    indexinRPM = RPMs[a[0]].index[0]

    # Delete Unnecessary Top Rows
    RPMdf = df.iloc[indexinRPM:]

    # Get Index of Each RPM Value to Bin Data
    RPMindex = RPMs.index

    # Clear Columns w/ Mostly NA
    RPMdf = pd.DataFrame.dropna(RPMdf, axis = 'columns', thresh = 10)

    # Add Index of Last Row of Database ot RPM DF to Bin
    lastIndex = RPMdf.iloc[[-1]].index.values[0]
    
    # Get Index of RPMs ot Group by RPM
    RPMindex = np.append(RPMindex, lastIndex)
    RPMlabel = RPMs[a[0]].to_list()

    # Label Columns - Assumes consistantcy in value-to-columns, V in first column
    Units = RPMdf[RPMdf[0].str.match('V')]

    # First Row of "V"
    unitsRow = Units[0].index[0]

    # Takes First Row, Sets as Column Headers
    RPMdf.columns = RPMdf.loc[unitsRow]

    # Removes All Non-Numeric Rows
    RPMdf = RPMdf.apply(lambda x: pd.to_numeric(x, errors = 'coerce')).dropna()


    # Group Each RPM
    RPMdf = RPMdf.groupby(pd.cut(RPMdf.index, RPMindex,
                                 right = False, labels = RPMlabel))
    
    RPMgroups = np.array(list(RPMdf.groups))

    # allAmps array w/ Expected Amp Draw at Max Static Thrust at each velocity
    allAmps = []

    for groups in RPMgroups:
        RPMgroup = RPMdf.get_group(groups)
        PWR = RPMgroup['PWR'].iloc[0] * 745.7  # Converting HP to Watts
        allAmps.append(PWR / (cells * batVolt))

    """ Calculations """

    velocityArray = [0]
    velocityStep = 5        #ft/s

    InitMotorGuess = 2    # Make First Motor Amp Draw Guess
    InitBatGuess = 2     # No Load Amp Draw for Motor


    # Dataframe for Storing Final Results
    resultDF = []

    powerTol = 0.0001
    relaxPara = .45
    
    velocity = 0
    velocityArray.append(velocity)

    ThrustProp = 1
    c2 = 0

    nLEMF = (cells*batVolt - nLCurr*(batRes*cells + eRes))*Throttle/100 - mRes*nLCurr
    nLRPM = kv*nLEMF

    bottomRPM = RPMgroups[RPMgroups < nLRPM].max()
    topRPM = RPMgroups[RPMgroups > nLRPM].min()

    # Get Bottom/Top RPM Groups
    bottomRPMdata = RPMdf.get_group(bottomRPM).reset_index()
    topRPMdata = RPMdf.get_group(topRPM).reset_index()

    velocityMPH = velocity * 2.23694      # Conversion to MPH

    # Finding Relevant Indices
    bottomRPMClosestVel = abs(bottomRPMdata['V'] - velocityMPH).idxmin()
    topRPMClosestVel = abs(topRPMdata['V'] - velocityMPH).idxmin()

    # Finding Cp
    bottomCp = bottomRPMdata['Cp'].iloc[bottomRPMClosestVel]
    topCp = topRPMdata['Cp'].iloc[topRPMClosestVel]

    # Interpolate to get Coefficient of Power/Thrust for Motor RPM
    nLCp = bottomCp + (nLRPM - bottomRPM) * (topCp - bottomCp) / (nLRPM - bottomRPM)

    nLpProp = Density * (nLRPM/60)**3 * (propDia)**5 * nLCp

    nLpProp *= 1.35581795

    while ThrustProp > 0:

        counter = 0

        pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, Torque, Cp, Ct = RPMfinder(InitMotorGuess, InitBatGuess, 0, velocity)

        InitMotorGuess = 2
        InitBatGuess = 2

        while abs(pProp * 1.35581795 - pMotorOut) > powerTol:
            
            pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, Torque, Cp, Ct = RPMfinder(InitMotorGuess, InitBatGuess, 0, velocity)

            aMotorGuess = pProp/vEMF + nLCurr
            aBatGuess = vTerm * motorAmp / (batVolt * cells - batAmp * (eRes + batRes))

            InitMotorGuess = (aMotorGuess - InitMotorGuess) * relaxPara + InitMotorGuess
            InitBatGuess = (aBatGuess - InitBatGuess) * relaxPara + InitBatGuess

            InitBatGuess = vTerm * motorAmp / (batVolt * cells - batAmp * (eRes + batRes*cells))

            convergence = pProp*1.35581795 - pMotorOut
            InitMotorGuess = ((motorAmp-nLCurr)/(convergence-nLpProp)*(0-nLpProp)+nLCurr)*relaxPara+(1-relaxPara)*motorAmp

            if InitMotorGuess < nLCurr:
                InitMotorGuess = nLCurr + 0.01
            if InitBatGuess < nLCurr:
                InitBatGuess = nLCurr + 0.01

            if counter > 100:
                break
            else:
                counter += 1

        efficiency = pMotorOut / pMotorIn

        resultDF.append([Throttle, Density, velocity, ThrustProp, Torque, RPM, batAmp, pProp, efficiency])
        print(resultDF)
        velocity += velocityStep

    resultDF = pd.DataFrame(resultDF, columns=['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
    print(resultDF.to_string())
    print("Done")
    return resultDF
    

print("Starting")
# PackageGen(33, 15) # input drag at lbs
propName = "C:/Users/marym/OneDrive/Documents/UIUC-DBF-main/Propulsion/propellors/PER3_19x12E2.dat"
path = "C:/Users/marym/OneDrive/Documents/UIUC-DBF-main/Propulsion/Prototype 0 - Stephanie"



# name = "Scorpion  SII-3026-710KV"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 3, 3.7, 710, 0.006, 1.56, 0.022)])
# resultProp.to_csv(path + "/" + name)



name = "Hacker A60 5S (2023)"
resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
resultProp = pd.concat([resultProp, PropulsionCalc(propName, 100, 100, 0, 120, 0.004, 6, 3.7, 295, 0.0028, 1.7, 0.015)])
resultProp.to_csv(path + "/" + name)



# name = "Scorpion A-4225-500kv"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 6, 3.7, 500, 0.006, 1.54, 0.014)])
# resultProp.to_csv(path + "/" + name)



# name = "Scorpion SII-4020-540KV"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 6, 3.7, 540, 0.006, 1.22, 0.020)])
# resultProp.to_csv(path + "/" + name)



# name = "Scorpion SII-4020-420KV"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 6, 3.7, 420, 0.006, 0.91, 0.032)])
# resultProp.to_csv(path + "/" + name)



# name = "Scorpion SII-4020-630KV"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 6, 3.7, 630, 0.006, 1.54, 0.015)])
# resultProp.to_csv(path + "/" + name)



# name = "Scorpion SII-4025-330KV"
# resultProp = pd.DataFrame(columns = ['Throttle','Air Density (slugs/ft^3)', 'Velocity (ft/s)', 'Thrust (lbs)', 'Torque (in-lbf)', 'RPM (rev/min)', 'Battery Current (A)', 'Power (W)', 'Efficiency'])
# for i in range(20,105,5):
#     resultProp = pd.concat([resultProp, PropulsionCalc(propName, i, 100, 0, 120, 0.004, 6, 3.7, 330, 0.006, 0.74, 0.037)])
# resultProp.to_csv(path + "/" + name)

print("Done")
# path = "C:/Users/marym/OneDrive/Documents/UIUC-DBF-main/Propulsion/Prototype 0 - Stephanie"
# propName = "C:/Users/marym/OneDrive/Documents/UIUC-DBF-main/Propulsion/propellors/PER3_8x5.dat"
# resultProp = PropulsionCalc(propName, 100, 200, 0, 200, 0.003, 3, 3.7, 1380, 0.0077, 0.4, 0.328)
# resultProp.to_csv(path+"/"+"HW3 Comparison")