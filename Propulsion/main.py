import numpy as np
import pandas as pd
#rom openpyxl import load_workbook
import os
import re

def PropulsionCalc(fileName, Throttle, maxVel, drag, eCurrMax, eRes, cells, batVolt, kv, batRes, nLCurr, mRes, Density = 1.225, alt = 0):
    
    def RPMfinder(motorAmpGuess, batAmpGuess, Output, velocity):

        # Main Calc for Iteration of pMotorOut
        def pMotorOutCalc(motorAmp, batAmp):

            vBat = cells * batVolt - cells * batRes * batAmp         # Total Bat Volts w/ Volt Loss
            vTerm = (vBat - eRes * batAmp) * Throttle / 100          # Volt Drop from esc/Throttle Val - voltage at terminal
            vDropMotor = motorAmp * mRes                             # Volt Drop from Motor
            pMotorIn = vTerm * batAmp                                # Power Motor Uses during Opt
            vEMF = vTerm - vDropMotor                                # Voltage Motor after losses
            pMotorOut = vEMF * (motorAmp - nLCurr)                   # Power Motor Generates (Voltage of Motor * (Amps of Motor - No Load Motor Amps))
            pMotorIn = vTerm * (motorAmp - nLCurr)                   # Power Motor Receives

            return pMotorIn, pMotorOut, vTerm, motorAmp, vBat, batAmp, vEMF

        pMotorIn, pMotorOut, vTerm, motorAmp, vBat, batAmp, vEMF = pMotorOutCalc(motorAmpGuess, batAmpGuess)

        RPM = vEMF * kv                # Revolutions/Min of Motor
        #print(str(vEMF))
        RPS = RPM / 60                 # Revolutions/s of Motor

        assert(batAmp < batVolt/batRes), "Batteries can't supply this much current: " + str(batAmp)
        assert(batAmpGuess < eCurrMax), "ESC will be overloaded, chose larger ESC: " + str(batAmpGuess)
        #print("RPM:" + str(RPM))
        assert(RPM > 0), "RPM is negative, something is wrong: " + str(RPM)


        # Find Closest Data to RPM Values, above and below
        RPMgroups = np.array(list(RPMdf.groups))

        #print("RPM max limit has been reached " + str(RPM) + " and Max is " + str(RPMgroups.max()))
        assert(RPM < RPMgroups.max()), "RPM max limit has been reached " + str(RPM) + " and Max is " + str(RPMgroups.max())
        assert(RPM > RPMgroups.min()), 'RPM minimum limit reached'

        bottomRPM = RPMgroups[RPMgroups < RPM].max()
        topRPM = RPMgroups[RPMgroups > RPM].min()

        # Get Bottom/Top RPM Groups
        bottomRPMdata = RPMdf.get_group(bottomRPM).reset_index()
        topRPMdata = RPMdf.get_group(topRPM).reset_index()

        velocityMPH = velocity * 2.23694      # Conversion to MPH

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
        bottomTorque = bottomRPMdata['Torque'].iloc[bottomRPMClosestVel]
        topTorque = topRPMdata['Torque'].iloc[topRPMClosestVel]

        # Interpolate to get Coefficient of Power/Thrust for Motor RPM
        Cp = bottomCp + (RPM - bottomRPM) * (topCp - bottomCp) / (topRPM - bottomRPM)
        Ct = bottomCt + (RPM - bottomRPM) * (topCt - bottomCt) / (topRPM - bottomRPM)
        torque = bottomTorque + (RPM - bottomRPM) * (topTorque - bottomTorque) / (topRPM - bottomRPM)

        # Prop Power/Thrust Calc
        pProp = Density * RPS**3 * propDia**5 * Cp         # Power Prop requires in kg * m^2 / s^3 (W)
        ThrustProp = Density * RPS**2 * propDia**4 * Ct    # Thrust Prop generates in kg*m/2  (N)

        return pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, torque

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
    propDia *= 0.0254
    
    
    # Set Column Index for Importing Data
    header_list = [str(x) for x in range(13)]

    # Read file
    try:
        propData = pd.read_csv(fileName, sep = '\s+', names = header_list)
    except:
        header_list = [str(x) for x in range(14)]
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
    
    
    """ First Calculations """
    RPMgroups = np.array(list(RPMdf.groups))

    # allAmps array w/ Expected Amp Draw at Max Static Thrust at each velocity
    allAmps = []

    for groups in RPMgroups:
        RPMgroup = RPMdf.get_group(groups)
        PWR = RPMgroup['PWR'].iloc[0] * 745.7  # Converting HP to Watts
        allAmps.append(PWR / (cells * batVolt))

    velocityArray = [0]
    velocityStep = 1.524
    
    
    """ Final Calculations """

    InitMotorGuess = nLCurr + 0.01    # Make First Motor Amp Draw Guess
    InitBatGuess = nLCurr + 0.01      # No Load Amp Draw for Motor

    # Dataframe for Storing Final Results
    resultProp = pd.DataFrame(columns = [ "Battery Current (A)", "RPM (rev/min)", "Thrust (N)", "Power (W)", "Efficiency", "Torque (in-lbf)"])
    resultProp.columns.name = 'Velocity (m/s)'

    powerTolerance = 0.01 # Tolerance of power to Motor to Prop

    velocityPosition = 0 # VelocityArray Position
    velocity = velocityArray[velocityPosition] # Current Velocity in m/s

    try:
        pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, Torque = RPMfinder(InitMotorGuess, InitBatGuess, 0, velocity)
    
    except AssertionError:
        return None
    
    
    drag *= 4.4482216
    
    boo = False
    boo2 = False

    while (ThrustProp > 0):

        if (velocity != 0):
            velocityArray.append(velocity)

        try:
            pMotorIn, pMotorOut, pProp, vTerm, motorAmp, vBat, batAmp, ThrustProp, vEMF, RPM, Torque = RPMfinder(InitMotorGuess, InitBatGuess, 0, velocity)
        
        except AssertionError:
            return None
        
        counter = 0
        while (abs(pProp - pMotorOut) > powerTolerance):
            
            try:
                pMotorIn,pMotorOut,pProp,vTerm,motorAmp,vBat,batAmp,ThrustProp,vEMF,RPM,Torque = RPMfinder(InitMotorGuess,InitBatGuess,0,velocity)
            
            except AssertionError:
                return None

            if (pMotorOut < 0):
                print("Motor Power Out is negative, something is wrong")
                return None

            pESCOut = pMotorIn
            vESCin = (cells * batVolt) - (cells * batRes * batAmp)
            vESCout = (cells * batVolt) - (cells * batRes + mRes) * batAmp
            InitMotorGuess = pProp / vEMF + nLCurr
            InitBatGuess = pESCOut / vESCout
            
            if counter > 1000:
                return None
            else:
                counter += 1

        # Runs one More Time (Could Keep previous motorAmp/BatAmp values, they are close to new)
        try:
            pMotorIn,pMotorOut,pProp,vTerm,motorAmp,vBat,batAmp,ThrustProp,vEMF,RPM,Torque = RPMfinder(InitMotorGuess,InitBatGuess,0,velocity)
        except AssertionError:
            return None
        
        if (ThrustProp*0.224809) > drag:
            boo = True
        
        if velocity > maxVel:
            boo2 = True

        efficiency = pMotorOut / pMotorIn

        resultProp.loc[velocityArray[velocityPosition]] = [motorAmp, RPM, ThrustProp, pProp, efficiency, Torque]

        velocityPosition += 1
        velocity += velocityStep
        
        if velocity > 100:
            return resultProp
    
    #if (velocity < maxVel):
        #print("Velocity: " + str(velocity))
        #return None
    
    thrustLB = resultProp.loc[:, 'Thrust (N)'] * 0.224809  # Newtons to LBS
    resultProp['Thrust (lbf)'] = thrustLB

    if boo == False or boo2 == False:
        return None
    
    # Removing Last Row b/c FN is 0 so results may be Invalid
    resultProp.drop(resultProp.tail(1).index, inplace = True)
    
    print(resultProp.to_string())
    return resultProp

def PackageGen(maxVel, drag):
    print("started main")
    propdir = 'propellers'
    
    # Pulling Dataframess
    mDF = pd.read_excel('config_specs.xlsx', sheet_name = 'Motors')
    batDF = pd.read_excel('config_specs.xlsx', sheet_name = 'Batteries')
    
    # Finding Number of Motors/Bats for Iteration
    mDFSize = len(mDF)
    batDFSize = len(batDF)    

    #esc information
    eCurrMax = 120     #max current of esc
    eRes = 0.004      #esc resistance in ohms
    
    # General Data
    batRes = 0.0028
    Throttle = 99
    
    '''
    Looping through motors, bats,
    and props for package data
    '''
    
    # Looping Motors
    for i in range(mDFSize):

        # Motor Data
        mName, kv, nLCurr, mRes = mDF.iloc[i, 0:4]
        
        parent_dir = "C:/Users/marym/OneDrive/Documents/Propulsion"
        directory = mName
        
        path = os.path.join(parent_dir, directory)
        os.mkdir(path)
        
        # Looping Batteries
        for i in range(batDFSize):

            # Battery Data
            batName, cells, batVolt = batDF.iloc[i, 0:3]
            
            parent_dir = directory
            directoryBat = batName
            
            path = os.path.join(parent_dir, directoryBat)
            os.mkdir(path)
            
            # Looping through Props
            for prop in os.listdir(propdir):
                print(str(prop))
                propName = "C:/Users/marym/OneDrive/Documents/Propulsion/propellers/" + prop  
                
                resultProp = PropulsionCalc(propName, Throttle, maxVel, drag, eCurrMax, eRes, cells, batVolt, kv, batRes, nLCurr, mRes, Density = 1.225*.84, alt = 0)
                
                resultProp.to_csv(path+"/"+prop)
        
        
print("Starting")
PackageGen(33, 15) # input drag at lbs