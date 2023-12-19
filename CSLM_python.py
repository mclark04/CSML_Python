import numpy as np
import math
import time
import copy
import os
#
# --------- USEFUL FORTRAN FUNCTIONS --------
# 
from FORTRANtoPyFunctions import MAX,MIN,SIGN

#
# -------- CSLM DISCRIPTION ---------------
#
#These are general functions call throughout the codebase
#   XIT   - an exit function to crash safely, used in FORTRAN to log buggs,
#           here it is just kept for compatabilty
#   EQNST - Equation of state called by CLASSL and FREECONV, based on 
#           Farmer and Carmack, 1981
#
from GenearlCSMLFunctions import XIT,EQNST


def RunLake():
    baseDir = os.getcwd()
    OverrideLakeTempsWithObs = False
    #Set up all the values of things that would come from CLASS and the GCM.
    NumOfGridCells=int(1) #Number of grid cells to be modelled
    NumOfMosaicTiles=int(3) #Maximum number of mosaic tiles per grid cell being modelled
    NumOfTilesInAllGrids=int(NumOfGridCells*NumOfMosaicTiles)
    MaxNumOfLakes=int(200)
    NBS=int(4)
    IGL=int(1)
    CurrentHour=int(0)
    CurrentMin=int(0)
    CurrentDay=int(0)
    CurrentYear=int(0)
    # -- GATHER-SCATTER INDEX ARRAYS
    IndGridLand=np.zeros((NumOfTilesInAllGrids),int)
    IndMosaicLand=np.zeros((NumOfTilesInAllGrids),int)
    IndGridWater=np.zeros((NumOfTilesInAllGrids),int)
    IndMosaicWater=np.zeros((NumOfTilesInAllGrids),int)
    #-- LAKE TILE PARAMETERS           
    NumOfLakeLayers_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),int)
    LakeDepth_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    LakeLength_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    ExtinctionCoeff_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagSensHeatOnLake_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagLatentHeatOnLake_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagNetShortwaveAtLakeSurface_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagNetLongwaveAtLakeSurface_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagPhaseChangeWaterAtLakeSurface_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowAlbedoVis_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowAlbedoNIR_In=np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    NumOfLakeLayers=np.zeros((NumOfTilesInAllGrids),int)
    LakeDepth=np.zeros((NumOfTilesInAllGrids),int)
    LakeLength=np.zeros((NumOfTilesInAllGrids),int)
    ExtinctionCoeff=np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowAlbedoVis=np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowAlbedoNIR=np.zeros((NumOfTilesInAllGrids),np.float32)
    AnemometerHeight=np.zeros((NumOfTilesInAllGrids),np.float32)
    VarHeight=np.zeros((NumOfTilesInAllGrids),np.float32)
    # -- LAKE PROGNOSTIC VARIABLES
    LakeTempProfile_In=np.zeros((NumOfGridCells,NumOfMosaicTiles,MaxNumOfLakes),np.float32)
    LakeTempProfileObs_In = np.zeros((NumOfGridCells,NumOfMosaicTiles,MaxNumOfLakes),np.float32)
    LakeTempProfile = np.zeros((NumOfTilesInAllGrids,MaxNumOfLakes),np.float32)
    LakeTempProfileObs = np.zeros((NumOfTilesInAllGrids,MaxNumOfLakes),np.float32)
    ExpansivityOfWater = np.zeros((NumOfTilesInAllGrids),np.float32)
    ThermoclineTempDelta= np.zeros((NumOfTilesInAllGrids),np.float32)
    DepthToThermocline= np.zeros((NumOfTilesInAllGrids),np.float32)
    MeanMixLyrMomentum= np.zeros((NumOfTilesInAllGrids),np.float32)
    ReducedGravity= np.zeros((NumOfTilesInAllGrids),np.float32)
    TotalKineticEnergy= np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeSkinTemp= np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeIceHeight= np.zeros((NumOfTilesInAllGrids),np.float32)
    MixLyrDensity= np.zeros((NumOfTilesInAllGrids),np.float32)
    SedimentTemp= np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowIceHeight= np.zeros((NumOfTilesInAllGrids),np.float32)
    RunoffFromIce = np.zeros((NumOfTilesInAllGrids),np.float32)
    # -- LAKE DIAGNOSTIC VARIABLES
    LakeSurfaceEvap = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagAirTempAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagZonalWindSpdAtAnemometerLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagMeridionalWindSpdAtAnemometerLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSpecificHumidAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagRelavtiveHumidAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPackMass = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPackDensity = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPackTemp = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPackAlbedo = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPackLiqudWaterContent  = np.zeros((NumOfTilesInAllGrids),np.float32)
    # -- Met VARIABLES
    PrecipOnLake_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    PrecipOnSnow_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    LakeSurfaceEvap_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowSublimation_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowRunoff_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagNetShortwaveAtSnowSurface_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagNetLongwaveAtSnowSurface_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagSensHeatOnSnow_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagLatentHeatOnSnow_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagPhaseChangeWaterInSnow_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagDelEnergyDueConducOrMassInSnow_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagDelEnergyDueConducOrMassInLake_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagAirTempAtScreenLevel_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagZonalWindSpdAtAnemometerLevel_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagMeridionalWindSpdAtAnemometerLevel_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagSpecificHumidAtScreenLevel_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagRelavtiveHumidAtScreenLevel_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    LongwaveUpwelling_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagTotalVisAtLakeSurface_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagTotalNIRAtLakeSurface_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    EvapEfficiencyAtLakeSurf_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagSurfBlkBodyTemp_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagSurfSpecificHumidity_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    LakeSurfaceDragAtNeutral_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    DiagPotentialEvap_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SensibleHeatFromLakeSubarea_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    ProdLakeDragWindSpdDelT_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    LatentHeatFromLakeSubarea_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SensibleHeatFromSnowSubarea_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    ProdLakeDragWindSpdDelSpHumid_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowPackMass_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowPackDensity_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowPackTemp_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowPackAlbedo_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SnowPackLiqudWaterContent_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    ChangeInInternalLakeEnergy = np.zeros((NumOfTilesInAllGrids),np.float64)
    ShortwaveDownwellingIRandVIS = np.zeros((NumOfTilesInAllGrids,NBS),np.float32)
    ShortwaveDownwellingIRandVIS_In = np.zeros((NumOfGridCells,NBS),np.float32)
    # -- ATMOSPHERIC AND GRID-CONSTANT INPUT VARIABLES
    FractionalCoverageOfMosaicTileOnArea_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    RefHeightWindSpeedForcing_In = np.zeros((NumOfGridCells),np.float32)
    RefHeightAirTempAndHumidForcing_In = np.zeros((NumOfGridCells),np.float32)
    AnemometerHeight_In = np.zeros((NumOfGridCells),np.float32)
    VarHeight_In = np.zeros((NumOfGridCells),np.float32)
    AtmoBlendingHeightForSurfaceRoughness_In = np.zeros((NumOfGridCells),np.float32)
    VisibleIncidentOnSurf_In = np.zeros((NumOfGridCells),np.float32)
    NIRIncidentOnSurf_In = np.zeros((NumOfGridCells),np.float32)
    LatitudeInRad_In = np.zeros((NumOfGridCells),np.float32)
    CosineOfSolarZenith_In = np.zeros((NumOfGridCells),np.float32)
    ShortwaveDownwelling_In = np.zeros((NumOfGridCells),np.float32)
    LongwaveDownwelling_In = np.zeros((NumOfGridCells),np.float32)
    ZonalWindSpeed_In = np.zeros((NumOfGridCells),np.float32)
    MeridionalWindSpeed_In = np.zeros((NumOfGridCells),np.float32)
    AirTempAtRef_In = np.zeros((NumOfGridCells),np.float32)
    SpecificHumidityAtRef_In = np.zeros((NumOfGridCells),np.float32)
    WindSpeed = np.zeros((NumOfGridCells),np.float32)
    SurfaceAirPressure_In = np.zeros((NumOfGridCells),np.float32)
    SurfacePrecipRate_In = np.zeros((NumOfGridCells),np.float32)
    PartialPressureOfDryAir_In = np.zeros((NumOfGridCells),np.float32)
    VapourPressureDeficit_In = np.zeros((NumOfGridCells),np.float32)
    DewPointTemp_In = np.zeros((NumOfGridCells),np.float32)
    DensityOfAir_In = np.zeros((NumOfGridCells),np.float32)
    SubareaFractionalCoverage_In = np.zeros((NumOfGridCells),np.float32)
    WetPrecip_In = np.zeros((NumOfGridCells),np.float32)
    WetPrecipTemp_In = np.zeros((NumOfGridCells),np.float32)
    FrozenPrecip_In = np.zeros((NumOfGridCells),np.float32)
    FrozenPrecipTemp_In = np.zeros((NumOfGridCells),np.float32)
    DensityOfFreshSnow_In = np.zeros((NumOfGridCells),np.float32)
    RefHeightWindSpeedForcing = np.zeros((NumOfTilesInAllGrids),np.float32)
    RefHeightAirTempAndHumidForcing = np.zeros((NumOfTilesInAllGrids),np.float32)
    VisibleIncidentOnSurf = np.zeros((NumOfTilesInAllGrids),np.float32)
    NIRIncidentOnSurf = np.zeros((NumOfTilesInAllGrids),np.float32)
    CosineOfSolarZenith = np.zeros((NumOfTilesInAllGrids),np.float32)
    LongwaveDownwelling = np.zeros((NumOfTilesInAllGrids),np.float32)
    ZonalWindSpeed = np.zeros((NumOfTilesInAllGrids),np.float32)
    MeridionalWindSpeed = np.zeros((NumOfTilesInAllGrids),np.float32)
    AirTempAtRef = np.zeros((NumOfTilesInAllGrids),np.float32)
    SpecificHumidityAtRef = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceAirPressure = np.zeros((NumOfTilesInAllGrids),np.float32)
    DensityOfAir = np.zeros((NumOfTilesInAllGrids),np.float32)
    WetPrecip = np.zeros((NumOfTilesInAllGrids),np.float32)
    WetPrecipTemp = np.zeros((NumOfTilesInAllGrids),np.float32)
    FrozenPrecip = np.zeros((NumOfTilesInAllGrids),np.float32)
    FrozenPrecipTemp = np.zeros((NumOfTilesInAllGrids),np.float32)
    DensityOfFreshSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    LatitudeInRad = np.zeros((NumOfTilesInAllGrids),np.float32)
    PartialPressureOfDryAir = np.zeros((NumOfTilesInAllGrids),np.float32)
    # -- LAND SURFACE DIAGNOSTIC VARIABLES
    SurfDragCoeffForHeat_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    SurfDragCoeffForMomentum_In = np.zeros((NumOfGridCells,NumOfMosaicTiles),np.float32)
    # -- PARAMETERS ORIGINALLY FROM CLASS COMMON BLOCK DATA     
    VonKarmanConst = np.float32(0.40)
    MinWindSpd = np.float32(0.1)
    ThermalEmissivityOfH2O = np.float32(0.97)
    ThermalConductOfH2O = np.float32(0.57)
    ThermalConductOfIce = np.float32(2.24)
    HeatCapacityOfH2O   = np.float32(4.187E6)
    HeatCapacityOfIce = np.float32(1.9257E6)  
    SpecificHeatOfH2O   = np.float32(4.186E3)
    SpecificHeatOfIce = np.float32(2.10E3)
    DensityOfH2O   = np.float32(1.0E3)
    DensityOfIce = np.float32(0.917E3)
    LatentHeatOfFreezingH2O = np.float32(0.334E6)
    LatentHeatOfVaporizationH2O = np.float32(2.501E6)
    # -- Lake TKE process efficiencies
    TotalKineticEnergyEff_CN = np.float32(1.33)
    TotalKineticEnergyEff_CF = np.float32(0.25)
    TotalKineticEnergyEff_CE = np.float32(1.15)
    TotalKineticEnergyEff_CS = np.float32(0.20)
    TotalKineticEnergyEff_CL = np.float32(0.2350) #This has a lot of uncertainty, can try adjusting.  Ryaner (1980) suggests zero but MacKay 2012 found this excessively enhanced deeping of the thermocline
    # -- Lake process parameter limits, and grid spacing
    MinimumMixingDepth = np.float32(0.5)
    MinimumKineticEnergy   = np.float32(1.0E-12)
    MaxThermoclineThickness   = np.float32(5.0)
    MinThermoclineThickness   = np.float32(0.5)
    LakeLayerThickness   = np.float32(0.5)
    ThicknessOfSurfaceSkin  = np.float32(0.050)
    MaxDeltaThermoclineDepth    = np.float32(2.0)
    MaxDeltaMixingLayer    = np.float32(0.1)
    # -- ASSIGN VALUES NORMALLY SPECIFIED WITHIN THE GCM.
    Pi     = np.float32(3.1415926535898) #Pi,numpy has built in defintions for pi() but redefined to keep consitant with original FORTRAN script
    GravConstant   = np.float32(9.80616) #gravity
    GasConstant   = np.float32(287.04)
    SpecificHeatOfAir = np.float32(1.00464E3)
    GasConstantOfH2OVapour  = np.float32(461.50)
    FreezingPointOfH2O  = np.float32(273.16) #freezing temperature
    StefanBoltzmannConst    = np.float32(5.66796E-8)
    DeltaTimeStep   = np.float32(600.0)
    #Preallocate some stuff
    SurfDragCoeffForHeat = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfDragCoeffForMomentum = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSensHeatOnLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagLatentHeatOnLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetShortwaveAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetLongwaveAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPhaseChangeWaterAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    PrecipOnLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    PrecipOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeSurfaceEvap = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowSublimation = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowRunoff = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSensHeatOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagLatentHeatOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetShortwaveAtSnowSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetLongwaveAtSnowSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPhaseChangeWaterInSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagDelEnergyDueConducOrMassInSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagDelEnergyDueConducOrMassInLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagAirTempAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagZonalWindSpdAtAnemometerLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagMeridionalWindSpdAtAnemometerLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSpecificHumidAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagRelavtiveHumidAtScreenLevel = np.zeros((NumOfTilesInAllGrids),np.float32)
    LongwaveUpwelling = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagTotalVisAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPotentialEvap = np.zeros((NumOfTilesInAllGrids),np.float32)
    EvapEfficiencyAtLakeSurf = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSurfBlkBodyTemp = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSurfSpecificHumidity = np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeSurfaceDragAtNeutral = np.zeros((NumOfTilesInAllGrids),np.float32)
    SensibleHeatFromLakeSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProdLakeDragWindSpdDelT = np.zeros((NumOfTilesInAllGrids),np.float32)
    LatentHeatFromLakeSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    SensibleHeatFromSnowSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProdLakeDragWindSpdDelSpHumid = np.zeros((NumOfTilesInAllGrids),np.float32)

    # -- CLASS SWITCHES select different runtime options
    # 
    #      * IF IDISP=0, VEGETATION DISPLACEMENT HEIGHTS ARE IGNORED,
    #      * BECAUSE THE ATMOSPHERIC MODEL CONSIDERS THESE TO BE PART
    #      * OF THE "TERRAIN".
    #      * IF IDISP=1, VEGETATION DISPLACEMENT HEIGHTS ARE CALCULATED.
    # 
    #      * IF IZREF=1, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
    #      * TO LIE AT THE GROUND SURFACE.
    #      * IF IZREF=2, THE BOTTOM OF THE ATMOSPHERIC MODEL IS TAKEN
    #      * TO LIE AT THE LOCAL ROUGHNESS HEIGHT.
    # 
    #      * IF ISLFD=0, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
    #      * AND THE ORIGINAL GCM SET OF SCREEN-LEVEL DIAGNOSTIC CALCULATIONS 
    #      * IS DONE.
    #      * IF ISLFD=1, DRCOEF IS CALLED FOR SURFACE STABILITY CORRECTIONS
    #      * AND SLDIAG IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 
    #      * IF ISLFD=2, FLXSURFZ IS CALLED FOR SURFACE STABILITY CORRECTIONS
    #      * AND DIASURF IS CALLED FOR SCREEN-LEVEL DIAGNOSTIC CALCULATIONS. 
    # 
    #      * IF IPCP=1, THE RAINFALL-SNOWFALL CUTOFF IS TAKEN TO LIE AT 0 C.
    #      * IF IPCP=2, A LINEAR PARTITIONING OF PRECIPITATION BETWEEEN 
    #      * RAINFALL AND SNOWFALL IS DONE BETWEEN 0 C AND 2 C.
    #      * IF IPCP=3, RAINFALL AND SNOWFALL ARE PARTITIONED ACCORDING TO
    #      * A POLYNOMIAL CURVE BETWEEN 0 C AND 6 C.
    #      * IF IPCP=4, THE RAINFALL, SNOWFALL AND TOTAL PRECIPITATION RATES
    #      * ARE READ IN DIRECTLY.
    # 
    #      * ITC, ITCG AND ITG ARE SWITCHES TO CHOOSE THE ITERATION SCHEME TO
    #      * BE USED IN CALCULATING THE CANOPY OR GROUND SURFACE TEMPERATURE
    #      * RESPECTIVELY.  IF THE SWITCH IS SET TO 1, A BISECTION METHOD IS
    #      * USED IF TO 2, THE NEWTON-RAPHSON METHOD IS USED.
    #      
    #      * IF IPAI, IHGT, IALC, IALS AND IALG ARE ZERO, THE VALUES OF 
    #      * PLANT AREA INDEX, VEGETATION HEIGHT, CANOPY ALBEDO, SNOW ALBEDO
    #      * AND SOIL ALBEDO RESPECTIVELY CALCULATED BY CLASS ARE USED.
    #      * IF ANY OF THESE SWITCHES IS SET TO 1, THE VALUE OF THE
    #      * CORRESPONDING PARAMETER CALCULATED BY CLASS IS OVERRIDDEN BY
    #      * A USER-SUPPLIED INPUT VALUE.
    # 
    IZREF= np.float32(1)
    ISLFD= np.float32(0)
    IPCP= np.float32(1)
    ITG= np.float32(1)           
    IALS= np.float32(0)
    ISNOALB= np.float32(0)
    # -- Debugging Vars (there should be nothing here)
    #TestVar = np.float32(0)
    #

    # -- OPEN FILES FOR READING AND WRITING.
    #fID_71=open(baseDir + '\Output\LAKE.of1','w')
    fID_71=open(baseDir + '\Output\LAKE_Temps.txt','w')
    fID_72=open(baseDir + '\Output\LAKE.of2','w')
    fID_73=open(baseDir + '\Output\LAKE.of3','w')
    fID_74=open(baseDir + '\Output\LAKE.of4','w')
    fID_75=open(baseDir + '\Output\LAKE.of5','w')
    fID_76=open(baseDir + '\Output\LAKE.of6','w')
    fID_77=open(baseDir + '\Output\LAKE.of7','w')
    #fID_82=open(baseDir + '\Output\ExtraEditing.txt','w')
    fID_82=open(baseDir + '\Output\MainOutput.txt','w')

    # -- Read in ini file
    f = open('LAKE.ini','r')
    IniTitle = f.readline()
    IniUserName  = f.readline()
    InIPlace = f.readline()
    line   = f.readline()
    arr = str.split(line)
    A = np.float32(arr)
    LatitudeInDeg_In=A[0]*1
    LongitudeInDeg_In=A[1]*1 #unused
    RefHeightWindSpeedForcing_In[0]=A[2]*1
    RefHeightAirTempAndHumidForcing_In[0]=A[3]*1
    AtmoBlendingHeightForSurfaceRoughness_In[0]=A[4]*1 #unused
    SubareaFractionalCoverage_In[0]=A[5]*1
    NumOfGridCellsInRun=int(A[6])*1
    NumOfMosaicTilesPerGridCellInRun=int(A[7])*1
    MosaicTileType = np.zeros((NumOfGridCellsInRun,NumOfMosaicTilesPerGridCellInRun),np.float32)
    for I in range(0,NumOfGridCellsInRun):          
        for M in range(0,NumOfMosaicTilesPerGridCellInRun):            
            tline = f.readline()
            TITLE1 = tline         #first line is the lake info
            tline = f.readline()    #second line is two variabiles MIDROT and FAREROT
            strFromFile=str.split(tline)
            A=np.float32(strFromFile)
            MosaicTileType[I,M]=A[0]*1
            FractionalCoverageOfMosaicTileOnArea_In[I,M]=A[1]*1
            if MosaicTileType[I,M] > 0.0: 
                raise NameError('RUNLAKE -> Issues loading Lake data, MIDROT is > 0')#if MIDROT is >0 safely crash and throw error            
            else:
                # -- READ IN WATER TILE INFORMATION                              
                tline = f.readline()#thrid line is HLAKROT, LLAKROT and BLAKROT
                strFromFile=str.split(tline)
                A=np.float32(strFromFile)
                LakeDepth_In[I,M]=A[0]*1
                LakeLength_In[I,M]=A[1]*1
                ExtinctionCoeff_In[I,M]=A[2]*1
                NumOfLakeLayers_In[I,M]=int(round(LakeDepth_In[I,M]/LakeLayerThickness))#calculate NLAKROT
                tline = f.readline()
                strFromFile=str.split(tline)
                A=np.float32(strFromFile)
                for J in range(0,NumOfLakeLayers_In[I,M]):
                    LakeTempProfile_In[I,M,J]=A[J]*1#,J=1,NLAKROT(I,M))  

    #Flip over lake input values, if in deg C convert to deg K by adding
    #freezing temp.
    for I in range(0,NumOfGridCellsInRun):  
        for M in range(0,NumOfMosaicTilesPerGridCellInRun):
            if NumOfLakeLayers_In[I,M]>=1: 
                for J in range(0,NumOfLakeLayers_In[I,M]):                                    
                    if (LakeTempProfile_In[I,M,J]<=200.0):
                        LakeTempProfile_In[I,M,J]=LakeTempProfile_In[I,M,J]+FreezingPointOfH2O
        # -- INITIALIZE LAKE SNOW CONDITIONS
        SnowPackMass_In[I,M]=0.0
        SnowPackAlbedo_In[I,M]=0.0
        SnowPackDensity_In[I,M]=0.0
        SnowPackTemp_In[I,M]=0.0
        SnowPackLiqudWaterContent_In[I,M]=0.0

    # -- THE CLASS GATHER-SCATTER MACHINERY IS RETAINED HERE FOR COMPATIBILITY
    TotalNumberOfMosaicTilesThatAreSurf=int(0)
    TotalNumberOfMosaicTilesThatAreH2O=int(0)
    for I in range(0,NumOfGridCellsInRun):
        for J in range(0,NumOfMosaicTilesPerGridCellInRun):
            if FractionalCoverageOfMosaicTileOnArea_In[I,J]>0.0:
                if MosaicTileType[I,J]>0:
                    TotalNumberOfMosaicTilesThatAreSurf=TotalNumberOfMosaicTilesThatAreSurf+1
                    IndGridLand[TotalNumberOfMosaicTilesThatAreSurf]=I
                    IndMosaicLand[TotalNumberOfMosaicTilesThatAreSurf]=J
                else:
                    TotalNumberOfMosaicTilesThatAreH2O=TotalNumberOfMosaicTilesThatAreH2O+1
                    IndGridWater[TotalNumberOfMosaicTilesThatAreH2O]=I
                    IndMosaicWater[TotalNumberOfMosaicTilesThatAreH2O]=J
                
    print('Total Number Of Mosaic Tiles That Are Surf = ' + str(TotalNumberOfMosaicTilesThatAreSurf) + ' and H2O = ' + str(TotalNumberOfMosaicTilesThatAreH2O))
    
    #---------------------------------------------------------------
    #Finished Initialization,
    #Starting main run loop.
    IterationCounter=int(0)
    DayCounter=int(1)
    NumberOfIterationsIn24Hours=int(86400/round(DeltaTimeStep))
    #print('    N, Nfrac, deltaT, totalT, IYEAR, IDAY, IHOUR, IMIN, QSUML-CTLSTP[I], QSUML')
    print('    N, Nfrac, deltaT, totalT, IYEAR, IDAY, IHOUR, IMIN, LE, H, Ts, ThermoClDpth, ThermoClDelta')   
    
    oldTime=time.time()
    firstTime=oldTime
    fID_51 = open('LAKE.met','r')
    endOfClimateDat=False
    DeltaOfObsAndModeledEnergyInLakeLayers =np.zeros((NumOfGridCells,NumOfLakeLayers_In[I,M]),np.float32)
    useLakeData=OverrideLakeTempsWithObs
    if useLakeData:
        raise Exception("Sorry, overriding temp profile functionality has not been added to Python yet")
        #This is the MATLAB way of doing it, should be easy to edit:
        #DHour=LakeProfileIn.IHOUR/24
        #DMin=LakeProfileIn.IMIN/(24*60)
        #TV=datenum(LakeProfileIn.IYEAR,0,0)+LakeProfileIn.IDAY+DHour+DMin
        #ObserverdLakeTempProfileMatrix=LakeProfileIn.TempProfile'+FreezingPointOfH2O    
        #[~,ObserverdLakeTempProfileLength]=size(ObserverdLakeTempProfileMatrix)    
        #ObserverdLakeTempProfileMatrixTimeStamp=[TV]
        #if isfield(LakeProfileIn,'excitationVal')
        #    VariableExcitationValue=LakeProfileIn.excitationVal
        #    UseVariableExication=true
        #else
        #    UseVariableExication=false
        #end
        #%clear flNm DHour DMin TV TempObserverdLakeTempProfileMatrix
        #clear DHour DMin TV TempObserverdLakeTempProfileMatrix LakeProfileIn 
        #%fID_52=fopen('D:\Dropbox\PostDoc\UofS_Course\Project\UnpackingNCAR\LakeProfile_Test.lkT')    
    else:
        UseVariableExication=False
    

    while not(endOfClimateDat):
        #
        #========================================================================
        #     * READ IN METEOROLOGICAL FORCING DATA FOR CURRENT TIME STEP
        #     * CALCULATE SOLAR ZENITH ANGLE AND COMPONENTS OF INCOMING SHORT-
        #     * WAVE RADIATION FLUX ESTIMATE FLUX PARTITIONS IF NECESSARY.
        #
        IterationCounter=IterationCounter+1
        for I in range(0,NumOfGridCellsInRun):
            tline = fID_51.readline()
            if len(tline)==0:
                endOfClimateDat=True
                break
            strFromFile=str.split(tline)
            metDat=strFromFile
            CurrentHour=int(metDat[0])*1
            CurrentMin=int(metDat[1])*1
            CurrentDay=int(metDat[2])*1
            CurrentYear=int(metDat[3])*1   # 
            ShortwaveDownwelling_In[I]=np.float32(metDat[4])*1  # downwelling shortwave
            LongwaveDownwelling_In[I]=np.float32(metDat[5])*1  # downwelling longwave
            SurfacePrecipRate_In[I]=np.float32(metDat[6])*1  # precipitation (m I think)
            AirTempAtRef_In[I]=np.float32(metDat[7])*1   # Air Temp
            SpecificHumidityAtRef_In[I]=np.float32(metDat[8])*1   # spacific humidity
            WindSpeed[I]=np.float32(metDat[9])*1  # windspeed
            SurfaceAirPressure_In[I]=np.float32(metDat[10])*1 #pressure in Pa
            VisibleIncidentOnSurf_In[I]=0.5*ShortwaveDownwelling_In[I]
            NIRIncidentOnSurf_In[I]=0.5*ShortwaveDownwelling_In[I]
            AirTempAtRef_In[I]=AirTempAtRef_In[I]+FreezingPointOfH2O
            ZonalWindSpeed_In[I]= WindSpeed[I] #all wind is assumed downwind for the first timestep.
            MeridionalWindSpeed_In[I]=np.float32(0.0)# %cross-wind speed                
            AnemometerHeight_In[I]=np.float32(10.0)
            VarHeight_In[I]=np.float32(2.0)
            ShortwaveDownwellingIRandVIS_In[I,0]=VisibleIncidentOnSurf_In[I]*1
            ShortwaveDownwellingIRandVIS_In[I,1]=NIRIncidentOnSurf_In[I]*1
            #Below is the MATLAB function for loading in lake profile data
            if useLakeData:
                raise Exception("Sorry, overriding temp profile functionality has not been added to Python yet")
            #    for K=1:TotalNumberOfMosaicTilesThatAreH2O                  
            #        %load Lake file and store in LakeTempProfile_In... check for potential deviations from previous loop as well!                
            #        for L=1:NumOfLakeLayers_In(IndGridWater(K),IndMosaicWater(K))
            #            currentTimeStamp=datenum(double([CurrentYear,0,CurrentDay,CurrentHour,CurrentMin,0]))                    
            #            ind=ObserverdLakeTempProfileMatrixTimeStamp(:,1)==currentTimeStamp
            #            if ~isempty(find(ind))&&L<=ObserverdLakeTempProfileLength
            #                LakeTempProfileObs_In(IndGridWater(K),IndMosaicWater(K),L)=ObserverdLakeTempProfileMatrix(ind,L)
            #            else
            #                LakeTempProfileObs_In(IndGridWater(K),IndMosaicWater(K),L)=nan
            #            end
            #        end
            #    end
            #    if UseVariableExication
            #        currentTimeStamp=datenum(double([CurrentYear,0,CurrentDay,CurrentHour,CurrentMin,0]))                    
            #        ind=ObserverdLakeTempProfileMatrixTimeStamp(:,1)==currentTimeStamp
            #        if ~isnan(VariableExcitationValue(ind))
            #            ExtinctionCoeff_In(I,M)=VariableExcitationValue(ind)
            #        end
            #    end
            #else
            #    for K=1:TotalNumberOfMosaicTilesThatAreH2O 
            #        currentTimeStamp=datenum(double([CurrentYear,0,CurrentDay,CurrentHour,CurrentMin,0]))   
            #        %Set obs to nan if no lake data is to be loaded
            #        for L=1:NumOfLakeLayers(K+TotalNumberOfMosaicTilesThatAreSurf)
            #            LakeTempProfileObs_In(IndGridWater(K),IndMosaicWater(K),L)=nan           
            #        end
            #    end
            #end              
        if endOfClimateDat:
            break        
        

        CurrentDecimalDay=CurrentDay+(CurrentHour+CurrentMin/60)/24
        SolarDeclanation=np.float32(math.sin(2.*Pi*(284+CurrentDecimalDay)/365)*23.45*Pi/180)
        SolarHour=(CurrentHour+CurrentMin/60)*Pi/12-Pi
        for I in range(0,NumOfGridCellsInRun):
            LatitudeInRad_In[I]=LatitudeInDeg_In*Pi/180
            CosineOfSolarZenith_Temp=np.float32(math.sin(LatitudeInRad_In[I])*math.cos(SolarDeclanation)+math.cos(LatitudeInRad_In[I])*math.cos(SolarDeclanation)*math.cos(SolarHour))                                                
            CosineOfSolarZenith_In[I]=SIGN(MAX(abs(CosineOfSolarZenith_Temp),1.0E-3),CosineOfSolarZenith_Temp)
        
        # -- CALCULATION OF ATMOSPHERIC INPUT VARIABLES.        
        [VapourPressureDeficit_In,DewPointTemp_In,PartialPressureOfDryAir_In,DensityOfAir_In,\
            DensityOfFreshSnow_In,WetPrecip_In,WetPrecipTemp_In,FrozenPrecip_In,\
            FrozenPrecipTemp_In,AirTempAtRef_In,SpecificHumidityAtRef_In,SurfacePrecipRate_In,\
            SurfaceAirPressure_In,IPCP,NumOfGridCells,XXX,NumOfGridCellsInRun] =\
            CLASSI(VapourPressureDeficit_In,DewPointTemp_In,PartialPressureOfDryAir_In,DensityOfAir_In,\
            DensityOfFreshSnow_In,WetPrecip_In,WetPrecipTemp_In,FrozenPrecip_In,\
            FrozenPrecipTemp_In,AirTempAtRef_In,SpecificHumidityAtRef_In,SurfacePrecipRate_In,\
            SurfaceAirPressure_In,IPCP,NumOfGridCells,1,\
            NumOfGridCellsInRun,FreezingPointOfH2O,GasConstant,GasConstantOfH2OVapour,DensityOfH2O) 
        #-- GATHER LAKE RELEVANT VARIABLES INTO WATER-TILE GAT ARRAYS          
        for K in range(0,TotalNumberOfMosaicTilesThatAreH2O):
            LakeDepth[K+TotalNumberOfMosaicTilesThatAreSurf]=LakeDepth_In[IndGridWater[K],IndMosaicWater[K]]*1
            LakeLength[K+TotalNumberOfMosaicTilesThatAreSurf]=LakeLength_In[IndGridWater[K],IndMosaicWater[K]]*1
            ExtinctionCoeff[K+TotalNumberOfMosaicTilesThatAreSurf]=ExtinctionCoeff_In[IndGridWater[K],IndMosaicWater[K]]            
            NumOfLakeLayers[K+TotalNumberOfMosaicTilesThatAreSurf]=NumOfLakeLayers_In[IndGridWater[K],IndMosaicWater[K]]*1  
            for L in range(0,int(NumOfLakeLayers[K+TotalNumberOfMosaicTilesThatAreH2O-1])): 
                LakeTempProfile[K+TotalNumberOfMosaicTilesThatAreSurf,L]=LakeTempProfile_In[IndGridWater[K],IndMosaicWater[K],L]
                LakeTempProfileObs[K+TotalNumberOfMosaicTilesThatAreSurf,L]=LakeTempProfileObs_In[IndGridWater[K],IndMosaicWater[K],L]                           
            SnowAlbedoVis[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowAlbedoVis_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowAlbedoNIR[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowAlbedoNIR_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowPackMass[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowPackMass_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowPackDensity[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowPackDensity_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowPackTemp[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowPackTemp_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowPackAlbedo[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowPackAlbedo_In[IndGridWater[K],IndMosaicWater[K]]*1
            SnowPackLiqudWaterContent[K+TotalNumberOfMosaicTilesThatAreSurf]=SnowPackLiqudWaterContent_In[IndGridWater[K],IndMosaicWater[K]]*1
        #-- ATMOSPHERIC FORCING VARIABLES NEEDED FOR LAKE TILES GATHERED ON TOP OF LAND TILES 
        for  K in range(0,TotalNumberOfMosaicTilesThatAreH2O):            
            VisibleIncidentOnSurf[K+TotalNumberOfMosaicTilesThatAreSurf]=VisibleIncidentOnSurf_In[IndGridWater[K]]                       
            NIRIncidentOnSurf[K+TotalNumberOfMosaicTilesThatAreSurf]=NIRIncidentOnSurf_In[IndGridWater[K]]                          
            CosineOfSolarZenith [K+TotalNumberOfMosaicTilesThatAreSurf]=CosineOfSolarZenith_In [IndGridWater[K]]                          
            LongwaveDownwelling [K+TotalNumberOfMosaicTilesThatAreSurf]=LongwaveDownwelling_In [IndGridWater[K]]                          
            ZonalWindSpeed  [K+TotalNumberOfMosaicTilesThatAreSurf]=ZonalWindSpeed_In  [IndGridWater[K]]                          
            MeridionalWindSpeed  [K+TotalNumberOfMosaicTilesThatAreSurf]=MeridionalWindSpeed_In  [IndGridWater[K]]                           
            AirTempAtRef  [K+TotalNumberOfMosaicTilesThatAreSurf]=AirTempAtRef_In  [IndGridWater[K]]                          
            SpecificHumidityAtRef  [K+TotalNumberOfMosaicTilesThatAreSurf]=SpecificHumidityAtRef_In  [IndGridWater[K]]                          
            SurfaceAirPressure[K+TotalNumberOfMosaicTilesThatAreSurf]=SurfaceAirPressure_In[IndGridWater[K]]                          
            DensityOfAir[K+TotalNumberOfMosaicTilesThatAreSurf]=DensityOfAir_In[IndGridWater[K]]                          
            RefHeightWindSpeedForcing[K+TotalNumberOfMosaicTilesThatAreSurf]=RefHeightWindSpeedForcing_In[IndGridWater[K]]                          
            RefHeightAirTempAndHumidForcing[K+TotalNumberOfMosaicTilesThatAreSurf]=RefHeightAirTempAndHumidForcing_In[IndGridWater[K]]                          
            for L in range(0,NBS):
                ShortwaveDownwellingIRandVIS[K+TotalNumberOfMosaicTilesThatAreSurf,L]=ShortwaveDownwellingIRandVIS_In[IndGridWater[K],L] 
            AnemometerHeight[K+TotalNumberOfMosaicTilesThatAreSurf]=AnemometerHeight_In[IndGridWater[K]]
            VarHeight[K+TotalNumberOfMosaicTilesThatAreSurf]=VarHeight_In[IndGridWater[K]]
            WetPrecip[K+TotalNumberOfMosaicTilesThatAreSurf]=WetPrecip_In[IndGridWater[K]]
            WetPrecipTemp[K+TotalNumberOfMosaicTilesThatAreSurf]=WetPrecipTemp_In[IndGridWater[K]]
            FrozenPrecip[K+TotalNumberOfMosaicTilesThatAreSurf]=FrozenPrecip_In[IndGridWater[K]]
            FrozenPrecipTemp[K+TotalNumberOfMosaicTilesThatAreSurf]=FrozenPrecipTemp_In[IndGridWater[K]]
            DensityOfFreshSnow[K+TotalNumberOfMosaicTilesThatAreSurf]=DensityOfFreshSnow_In[IndGridWater[K]]
            LatitudeInRad[K+TotalNumberOfMosaicTilesThatAreSurf]=LatitudeInRad_In[IndGridWater[K]]
            PartialPressureOfDryAir[K+TotalNumberOfMosaicTilesThatAreSurf]=PartialPressureOfDryAir_In[IndGridWater[K]]
        #Gather Other stuff        
        # Flux diagnoistic & vars
        for K in range(0,TotalNumberOfMosaicTilesThatAreH2O):
            SurfDragCoeffForHeat[K+TotalNumberOfMosaicTilesThatAreSurf] = SurfDragCoeffForHeat_In[IndGridWater[K],IndMosaicWater[K]]*1               
            SurfDragCoeffForMomentum[K+TotalNumberOfMosaicTilesThatAreSurf] = SurfDragCoeffForMomentum_In[IndGridWater[K],IndMosaicWater[K]]*1                 
            DiagSensHeatOnLake[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagSensHeatOnLake_In[IndGridWater[K],IndMosaicWater[K]]*1               
            DiagLatentHeatOnLake[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagLatentHeatOnLake_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagNetShortwaveAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagNetShortwaveAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagNetLongwaveAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagNetLongwaveAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagPhaseChangeWaterAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagPhaseChangeWaterAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]*1                
            PrecipOnLake[K+TotalNumberOfMosaicTilesThatAreSurf] = PrecipOnLake_In[IndGridWater[K],IndMosaicWater[K]]*1                
            PrecipOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf] = PrecipOnSnow_In[IndGridWater[K],IndMosaicWater[K]]*1               
            LakeSurfaceEvap[K+TotalNumberOfMosaicTilesThatAreSurf] = LakeSurfaceEvap_In[IndGridWater[K],IndMosaicWater[K]]*1                
            SnowSublimation[K+TotalNumberOfMosaicTilesThatAreSurf] = SnowSublimation_In[IndGridWater[K],IndMosaicWater[K]]*1                
            SnowRunoff[K+TotalNumberOfMosaicTilesThatAreSurf] = SnowRunoff_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagSensHeatOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagSensHeatOnSnow_In[IndGridWater[K],IndMosaicWater[K]]*1               
            DiagLatentHeatOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagLatentHeatOnSnow_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagNetShortwaveAtSnowSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagNetShortwaveAtSnowSurface_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagNetLongwaveAtSnowSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagNetLongwaveAtSnowSurface_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagPhaseChangeWaterInSnow[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagPhaseChangeWaterInSnow_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagDelEnergyDueConducOrMassInSnow[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagDelEnergyDueConducOrMassInSnow_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagDelEnergyDueConducOrMassInLake[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagDelEnergyDueConducOrMassInLake_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagAirTempAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagAirTempAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagZonalWindSpdAtAnemometerLevel[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagZonalWindSpdAtAnemometerLevel_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagMeridionalWindSpdAtAnemometerLevel[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagMeridionalWindSpdAtAnemometerLevel_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagSpecificHumidAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagSpecificHumidAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagRelavtiveHumidAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagRelavtiveHumidAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]*1                
            LongwaveUpwelling[K+TotalNumberOfMosaicTilesThatAreSurf] = LongwaveUpwelling_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagTotalVisAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagTotalVisAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]*1               
            DiagPotentialEvap[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagPotentialEvap_In[IndGridWater[K],IndMosaicWater[K]]*1               
            EvapEfficiencyAtLakeSurf[K+TotalNumberOfMosaicTilesThatAreSurf] = EvapEfficiencyAtLakeSurf_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagSurfBlkBodyTemp[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagSurfBlkBodyTemp_In[IndGridWater[K],IndMosaicWater[K]]*1                
            DiagSurfSpecificHumidity[K+TotalNumberOfMosaicTilesThatAreSurf] = DiagSurfSpecificHumidity_In[IndGridWater[K],IndMosaicWater[K]]*1                
            LakeSurfaceDragAtNeutral[K+TotalNumberOfMosaicTilesThatAreSurf] = LakeSurfaceDragAtNeutral_In[IndGridWater[K],IndMosaicWater[K]]*1                
            SensibleHeatFromLakeSubarea[K+TotalNumberOfMosaicTilesThatAreSurf] = SensibleHeatFromLakeSubarea_In[IndGridWater[K],IndMosaicWater[K]]*1                
            ProdLakeDragWindSpdDelT[K+TotalNumberOfMosaicTilesThatAreSurf] = ProdLakeDragWindSpdDelT_In[IndGridWater[K],IndMosaicWater[K]]*1                
            LatentHeatFromLakeSubarea[K+TotalNumberOfMosaicTilesThatAreSurf] = LatentHeatFromLakeSubarea_In[IndGridWater[K],IndMosaicWater[K]]*1
            SensibleHeatFromSnowSubarea[K+TotalNumberOfMosaicTilesThatAreSurf] = SensibleHeatFromSnowSubarea_In[IndGridWater[K],IndMosaicWater[K]]*1                
            ProdLakeDragWindSpdDelSpHumid[K+TotalNumberOfMosaicTilesThatAreSurf] = ProdLakeDragWindSpdDelSpHumid_In[IndGridWater[K],IndMosaicWater[K]]*1                            
        # CHECK ENERGY BALANCE - start of timestep
        # IceBottomDepth doesn't include weight of snow
        # 
        # -> Added mechanism to override lake temperature profile with
        # observations here since it has the machinery to handle ice conditions    
        DensityRatioOfIceToWater=DensityOfIce/DensityOfH2O
        JL1=1+TotalNumberOfMosaicTilesThatAreSurf                                            
        JL2=TotalNumberOfMosaicTilesThatAreH2O+TotalNumberOfMosaicTilesThatAreSurf
        for I in range(JL1-1,JL2):
            IceBottomDepth=DensityRatioOfIceToWater*LakeIceHeight[I]
            IceTopDepth=LakeIceHeight[I]-IceBottomDepth
            if IceBottomDepth >= ThicknessOfSurfaceSkin:
                HeatCapacityOfLakeSurface=HeatCapacityOfIce*1
            elif LakeIceHeight[I] <= 0.0: 
                HeatCapacityOfLakeSurface=HeatCapacityOfH2O*1
            else:
                HeatCapacityOfLakeSurface=(LakeIceHeight[I]*HeatCapacityOfIce + (ThicknessOfSurfaceSkin-IceBottomDepth)*HeatCapacityOfH2O)/ThicknessOfSurfaceSkin
            for J in range(0,int(NumOfLakeLayers[I])):
                CurrentLayerTop=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
                CurrentLayerBottom=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
                if IceBottomDepth >= CurrentLayerBottom:
                    HeatCapacityOfLakeSurface=HeatCapacityOfIce
                elif IceBottomDepth <= CurrentLayerTop:
                        HeatCapacityOfLakeSurface=HeatCapacityOfH2O
                else:
                    Z=IceBottomDepth-CurrentLayerTop
                    HeatCapacityOfLakeSurface=(Z*HeatCapacityOfIce + (LakeLayerThickness-Z)*HeatCapacityOfH2O)/LakeLayerThickness                
                ChangeInInternalLakeEnergy[I]=ChangeInInternalLakeEnergy[I] - HeatCapacityOfLakeSurface*LakeTempProfile[I,J]*LakeLayerThickness
                #This is the function to handle difference in Obs and Modeled
                #lake temps (from previous iteration)
                if useLakeData:
                    raise Exception("Sorry, overriding temp profile functionality has not been added to Python yet")
                # Below is the MATLAB way to overwrite, shuold be easy to adapt here to python.
                # #if useLakeData &&~isnan(LakeTempProfileObs(I,J))                
                #    TempDelta=LakeTempProfileObs(I,J)-LakeTempProfile(I,J)
                #    DeltaOfObsAndModeledEnergyInLakeLayers(I,J)=DeltaOfObsAndModeledEnergyInLakeLayers[I] + (HeatCapacityOfLakeSurface*TempDelta*LakeLayerThickness)
                #    LakeTempProfile(I,J)=LakeTempProfileObs(I,J)%overwrite modelled temp with actual temperature                
        
        #CLASSL is the actual Lake Model, so a lot of stuff is passed in and out of this function              
        [LakeDepth, LakeLength, ExtinctionCoeff, NumOfLakeLayers,\
            LakeTempProfile, LakeSkinTemp, DepthToThermocline, LakeIceHeight,\
            SnowIceHeight, RunoffFromIce, SnowPackMass, SnowPackDensity,\
            SnowPackTemp, SnowPackAlbedo, SnowPackLiqudWaterContent, SurfDragCoeffForHeat,\
            SurfDragCoeffForMomentum, SensibleHeatFromLakeSubarea,ProdLakeDragWindSpdDelT,LatentHeatFromLakeSubarea,\
            SensibleHeatFromSnowSubarea,ProdLakeDragWindSpdDelSpHumid,DiagPotentialEvap, EvapEfficiencyAtLakeSurf,\
            DiagSurfBlkBodyTemp, DiagSurfSpecificHumidity, LakeSurfaceDragAtNeutral,DiagAirTempAtScreenLevel,\
            DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagSpecificHumidAtScreenLevel,DiagRelavtiveHumidAtScreenLevel,\
            LongwaveUpwelling,DiagTotalVisAtLakeSurface, DiagNetShortwaveAtLakeSurface, DiagNetLongwaveAtLakeSurface,\
            DiagSensHeatOnLake, DiagLatentHeatOnLake, DiagPhaseChangeWaterAtLakeSurface, DiagDelEnergyDueConducOrMassInLake,\
            DiagNetShortwaveAtSnowSurface, DiagNetLongwaveAtSnowSurface, DiagSensHeatOnSnow, DiagLatentHeatOnSnow,\
            DiagPhaseChangeWaterInSnow, DiagDelEnergyDueConducOrMassInSnow, PrecipOnLake,PrecipOnSnow,\
            LakeSurfaceEvap,SnowSublimation,SnowRunoff,ExpansivityOfWater,\
            ThermoclineTempDelta, TotalKineticEnergy, MeanMixLyrMomentum, ReducedGravity,\
            MixLyrDensity, LongwaveDownwelling, ZonalWindSpeed, MeridionalWindSpeed,\
            AirTempAtRef, SpecificHumidityAtRef, DensityOfAir, PartialPressureOfDryAir,\
            SurfaceAirPressure, CosineOfSolarZenith, RefHeightWindSpeedForcing, RefHeightAirTempAndHumidForcing,\
            AnemometerHeight,VarHeight,WetPrecipTemp,FrozenPrecipTemp,\
            DensityOfFreshSnow, LatitudeInRad, SnowAlbedoVis,SnowAlbedoNIR,\
            ShortwaveDownwellingIRandVIS, NumOfGridCells, MaxNumOfLakes, ISLFD,\
            IZREF, ITG, IALS, NBS, ISNOALB, IGL, SedimentTemp] =\
        CLASSL(LakeDepth, LakeLength, ExtinctionCoeff, NumOfLakeLayers,\
            LakeTempProfile,LakeSkinTemp, DepthToThermocline, LakeIceHeight,\
            SnowIceHeight, RunoffFromIce,SnowPackMass,SnowPackDensity,\
            SnowPackTemp,SnowPackAlbedo,SnowPackLiqudWaterContent,DiagAirTempAtScreenLevel,\
            DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagSpecificHumidAtScreenLevel,DiagRelavtiveHumidAtScreenLevel,\
            LakeSurfaceEvap,ExpansivityOfWater, ThermoclineTempDelta, TotalKineticEnergy,\
            MeanMixLyrMomentum, ReducedGravity, MixLyrDensity,VisibleIncidentOnSurf,\
            NIRIncidentOnSurf, LongwaveDownwelling, ZonalWindSpeed, MeridionalWindSpeed,\
            AirTempAtRef, SpecificHumidityAtRef,DensityOfAir, PartialPressureOfDryAir,\
            SurfaceAirPressure, CosineOfSolarZenith, RefHeightWindSpeedForcing, RefHeightAirTempAndHumidForcing,\
            AnemometerHeight,VarHeight,WetPrecip,WetPrecipTemp,\
            FrozenPrecip,FrozenPrecipTemp,DensityOfFreshSnow, LatitudeInRad,\
            SnowAlbedoVis,SnowAlbedoNIR, ShortwaveDownwellingIRandVIS,NumOfTilesInAllGrids,\
            JL1, JL2, NumOfGridCells, MaxNumOfLakes,\
            ISLFD, IZREF,ITG, IALS,\
            NBS, ISNOALB, IGL, IterationCounter,\
            CurrentYear, CurrentDay, CurrentHour, CurrentMin,\
            SedimentTemp,FreezingPointOfH2O,DensityOfIce,DensityOfH2O,\
            HeatCapacityOfH2O,SpecificHeatOfIce,LatentHeatOfFreezingH2O,MinWindSpd,\
            MinimumKineticEnergy,LakeLayerThickness,GravConstant,ThicknessOfSurfaceSkin,\
            VonKarmanConst,LatentHeatOfVaporizationH2O,SpecificHeatOfAir,ThermalEmissivityOfH2O,\
            StefanBoltzmannConst,ThermalConductOfH2O,DeltaTimeStep,TotalKineticEnergyEff_CN,\
            TotalKineticEnergyEff_CF,TotalKineticEnergyEff_CE,TotalKineticEnergyEff_CS,TotalKineticEnergyEff_CL,\
            SpecificHeatOfH2O, MaxDeltaThermoclineDepth,MinimumMixingDepth,MaxDeltaMixingLayer,\
            MaxThermoclineThickness, MinThermoclineThickness,ThermalConductOfIce,HeatCapacityOfIce,\
            fID_71,fID_72,fID_73,fID_74,\
            fID_75,fID_76,fID_77,fID_82)
        
        # -- CHECK ENERGY BALANCE - again here at the end of timestep                
        for I in range(JL1-1,JL2):
            IceBottomDepth=DensityRatioOfIceToWater*LakeIceHeight[I]
            ICETOP=LakeIceHeight[I]-IceBottomDepth
            if IceBottomDepth >= ThicknessOfSurfaceSkin:
                HCAP=HeatCapacityOfIce*1
            elif LakeIceHeight[I] <= 0.0:
                HCAP=HeatCapacityOfH2O*1
            else:
                HCAP=(LakeIceHeight[I]*HeatCapacityOfIce + (ThicknessOfSurfaceSkin-IceBottomDepth)*HeatCapacityOfH2O)/ThicknessOfSurfaceSkin
            
            ChangeInInternalLakeEnergy[I]= ChangeInInternalLakeEnergy[I] + HCAP*LakeSkinTemp[I]*ThicknessOfSurfaceSkin
            
            for J in range(0,int(NumOfLakeLayers[I])):
                ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
                ZBOT=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
                if IceBottomDepth >= ZBOT:
                    HCAP=HeatCapacityOfIce*1
                elif IceBottomDepth <= ZTOP:
                    HCAP=HeatCapacityOfH2O*1
                else:
                    Z=IceBottomDepth-ZTOP
                    HCAP=(Z*HeatCapacityOfIce + (LakeLayerThickness-Z)*HeatCapacityOfH2O)/LakeLayerThickness            
                ChangeInInternalLakeEnergy[I]=ChangeInInternalLakeEnergy[I] + HCAP*LakeTempProfile[I,J]*LakeLayerThickness        
            ChangeInInternalLakeEnergy[I]=ChangeInInternalLakeEnergy[I]/DeltaTimeStep        
            QSUML=DiagNetShortwaveAtLakeSurface[I]+DiagNetLongwaveAtLakeSurface[I]-DiagSensHeatOnLake[I]-DiagLatentHeatOnLake[I]-DiagPhaseChangeWaterAtLakeSurface[I]
            
            if IterationCounter%1000==0:            
                newTime=time.time()
                deltaT=(newTime-oldTime)
                oldTime=newTime
                elapseTime=(newTime-firstTime)#*60*24
                LngthOfMetDat=46000 #This is just holdover from matlab            
                diff=QSUML-ChangeInInternalLakeEnergy[I]   
                print(f'{IterationCounter:5d},{IterationCounter/LngthOfMetDat:6.2f},{deltaT:7.3f},{elapseTime:7.1f},{CurrentYear:6d},{CurrentDay:5d},{CurrentHour:6d},{CurrentMin:5d},{LatentHeatFromLakeSubarea[I]:5.2e},{SensibleHeatFromLakeSubarea[I]:5.2e},{LakeSkinTemp[I]:3.2e},{DepthToThermocline[I]:3.2f},{ThermoclineTempDelta[I]:3.2f}')
        # -- Scatter code (maintained for CLASS compatabilty)
        # Lake diagnostic vars
        for K in range(0,TotalNumberOfMosaicTilesThatAreH2O):
            LakeDepth_In[IndGridWater[K],IndMosaicWater[K]]=LakeDepth[K+TotalNumberOfMosaicTilesThatAreSurf]*1                 
            LakeLength_In[IndGridWater[K],IndMosaicWater[K]]=LakeLength[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            ExtinctionCoeff_In[IndGridWater[K],IndMosaicWater[K]]=ExtinctionCoeff[K+TotalNumberOfMosaicTilesThatAreSurf]*1                 
            NumOfLakeLayers_In[IndGridWater[K],IndMosaicWater[K]]=NumOfLakeLayers[K+TotalNumberOfMosaicTilesThatAreSurf]*1                 
            for L in range(0,int(NumOfLakeLayers[K+TotalNumberOfMosaicTilesThatAreH2O-1])):
                LakeTempProfile_In[IndGridWater[K],IndMosaicWater[K],L]=LakeTempProfile[K+TotalNumberOfMosaicTilesThatAreSurf,L]*1
            SnowAlbedoVis_In[IndGridWater[K],IndMosaicWater[K]]=SnowAlbedoVis[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowAlbedoNIR_In[IndGridWater[K],IndMosaicWater[K]]=SnowAlbedoNIR[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowPackMass_In[IndGridWater[K],IndMosaicWater[K]]=SnowPackMass[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowPackDensity_In[IndGridWater[K],IndMosaicWater[K]]=SnowPackDensity[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowPackTemp_In[IndGridWater[K],IndMosaicWater[K]]=SnowPackTemp[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowPackAlbedo_In[IndGridWater[K],IndMosaicWater[K]]=SnowPackAlbedo[K+TotalNumberOfMosaicTilesThatAreSurf]*1
            SnowPackLiqudWaterContent_In[IndGridWater[K],IndMosaicWater[K]]=SnowPackLiqudWaterContent[K+TotalNumberOfMosaicTilesThatAreSurf]*1
        # Flux diagnoistic vars
        for K in range(0,TotalNumberOfMosaicTilesThatAreH2O):
            SurfDragCoeffForHeat_In[IndGridWater[K],IndMosaicWater[K]]=SurfDragCoeffForHeat[K+TotalNumberOfMosaicTilesThatAreSurf]*1                 
            SurfDragCoeffForMomentum_In[IndGridWater[K],IndMosaicWater[K]]=SurfDragCoeffForMomentum[K+TotalNumberOfMosaicTilesThatAreSurf]*1                 
            DiagSensHeatOnLake_In[IndGridWater[K],IndMosaicWater[K]]=DiagSensHeatOnLake[K+TotalNumberOfMosaicTilesThatAreSurf]*1               
            DiagLatentHeatOnLake_In[IndGridWater[K],IndMosaicWater[K]]=DiagLatentHeatOnLake[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagNetShortwaveAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagNetShortwaveAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagNetLongwaveAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagNetLongwaveAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagPhaseChangeWaterAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagPhaseChangeWaterAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            PrecipOnLake_In[IndGridWater[K],IndMosaicWater[K]]=PrecipOnLake[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            PrecipOnSnow_In[IndGridWater[K],IndMosaicWater[K]]=PrecipOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf]*1               
            LakeSurfaceEvap_In[IndGridWater[K],IndMosaicWater[K]]=LakeSurfaceEvap[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            SnowSublimation_In[IndGridWater[K],IndMosaicWater[K]]=SnowSublimation[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            SnowRunoff_In[IndGridWater[K],IndMosaicWater[K]]=SnowRunoff[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagSensHeatOnSnow_In[IndGridWater[K],IndMosaicWater[K]]=DiagSensHeatOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf]*1               
            DiagLatentHeatOnSnow_In[IndGridWater[K],IndMosaicWater[K]]=DiagLatentHeatOnSnow[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagNetShortwaveAtSnowSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagNetShortwaveAtSnowSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagNetLongwaveAtSnowSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagNetLongwaveAtSnowSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagPhaseChangeWaterInSnow_In[IndGridWater[K],IndMosaicWater[K]]=DiagPhaseChangeWaterInSnow[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagDelEnergyDueConducOrMassInSnow_In[IndGridWater[K],IndMosaicWater[K]]=DiagDelEnergyDueConducOrMassInSnow[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagDelEnergyDueConducOrMassInLake_In[IndGridWater[K],IndMosaicWater[K]]=DiagDelEnergyDueConducOrMassInLake[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagAirTempAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]=DiagAirTempAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagZonalWindSpdAtAnemometerLevel_In[IndGridWater[K],IndMosaicWater[K]]=DiagZonalWindSpdAtAnemometerLevel[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagMeridionalWindSpdAtAnemometerLevel_In[IndGridWater[K],IndMosaicWater[K]]=DiagMeridionalWindSpdAtAnemometerLevel[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagSpecificHumidAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]=DiagSpecificHumidAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagRelavtiveHumidAtScreenLevel_In[IndGridWater[K],IndMosaicWater[K]]=DiagRelavtiveHumidAtScreenLevel[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            LongwaveUpwelling_In[IndGridWater[K],IndMosaicWater[K]]=LongwaveUpwelling[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagTotalVisAtLakeSurface_In[IndGridWater[K],IndMosaicWater[K]]=DiagTotalVisAtLakeSurface[K+TotalNumberOfMosaicTilesThatAreSurf]*1               
            DiagPotentialEvap_In[IndGridWater[K],IndMosaicWater[K]]=DiagPotentialEvap[K+TotalNumberOfMosaicTilesThatAreSurf]*1               
            EvapEfficiencyAtLakeSurf_In[IndGridWater[K],IndMosaicWater[K]]=EvapEfficiencyAtLakeSurf[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagSurfBlkBodyTemp_In[IndGridWater[K],IndMosaicWater[K]]=DiagSurfBlkBodyTemp[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            DiagSurfSpecificHumidity_In[IndGridWater[K],IndMosaicWater[K]]=DiagSurfSpecificHumidity[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            LakeSurfaceDragAtNeutral_In[IndGridWater[K],IndMosaicWater[K]]=LakeSurfaceDragAtNeutral[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            SensibleHeatFromLakeSubarea_In[IndGridWater[K],IndMosaicWater[K]]=SensibleHeatFromLakeSubarea[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            ProdLakeDragWindSpdDelT_In[IndGridWater[K],IndMosaicWater[K]]=ProdLakeDragWindSpdDelT[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            LatentHeatFromLakeSubarea_In[IndGridWater[K],IndMosaicWater[K]]=LatentHeatFromLakeSubarea[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            SensibleHeatFromSnowSubarea_In[IndGridWater[K],IndMosaicWater[K]]=SensibleHeatFromSnowSubarea[K+TotalNumberOfMosaicTilesThatAreSurf]*1                
            ProdLakeDragWindSpdDelSpHumid_In[IndGridWater[K],IndMosaicWater[K]]=ProdLakeDragWindSpdDelSpHumid[K+TotalNumberOfMosaicTilesThatAreSurf]*1                            
        #increase daily counter    
        DayCounter=DayCounter+1
        if DayCounter > NumberOfIterationsIn24Hours:
            DayCounter=1 
        

        


def CLASSL(LakeDepth_Ref, LakeLength_Ref, ExtinctionCoeff_Ref, NumOfLakeLayers_Ref,\
            LakeTempProfile_Ref, LakeSkinTemp_Ref, DepthToThermocline_Ref, LakeIceHeight_Ref,\
            SnowIceHeight_Ref, RunoffFromIce_Ref, SnowPackMass_Ref, SnowPackDensity_Ref,\
            SnowPackTemp_Ref, SnowPackAlbedo_Ref, SnowPackLiqudWaterContent_Ref, DiagAirTempAtScreenLevel_Ref,\
            DiagZonalWindSpdAtAnemometerLevel_Ref,DiagMeridionalWindSpdAtAnemometerLevel_Ref, DiagSpecificHumidAtScreenLevel_Ref, DiagRelavtiveHumidAtScreenLevel_Ref,\
            LakeSurfaceEvap_Ref,ExpansivityOfWater_Ref, ThermoclineTempDelta_Ref, TotalKineticEnergy_Ref,\
            MeanMixLyrMomentum_Ref, ReducedGravity_Ref, MixLyrDensity_Ref, VisibleIncidentOnSurf_Ref,\
            NIRIncidentOnSurf_Ref,LongwaveDownwelling_Ref, ZonalWindSpeed_Ref, MeridionalWindSpeed_Ref,\
            AirTempAtRef_Ref, SpecificHumidityAtRef_Ref,DensityOfAir_Ref, PartialPressureOfDryAir_Ref,\
            SurfaceAirPressure_Ref, CosineOfSolarZenith_Ref, RefHeightWindSpeedForcing_Ref, RefHeightAirTempAndHumidForcing_Ref,\
            AnemometerHeight_Ref, VarHeight_Ref, WetPrecip_Ref, WetPrecipTemp_Ref,\
            FrozenPrecip_Ref, FrozenPrecipTemp_Ref, DensityOfFreshSnow_Ref, LatitudeInRad_Ref,\
            SnowAlbedoVis_Ref, SnowAlbedoNIR_Ref, ShortwaveDownwellingIRandVIS_Ref,NumOfTilesInAllGrids,\
            IL1, IL2, NumOfGridCells, MaxNumOfLakes,\
            ISLFD, IZREF,ITG, IALS,\
            NBS, ISNOALB, IGL, IterationCounter,\
            CurrentYear, CurrentDay, CurrentHour, CurrentMin,\
            SedimentTemp,FreezingPointOfH2O,DensityOfIce,DensityOfH2O,\
            HeatCapacityOfH2O,SpecificHeatOfIce,LatentHeatOfFreezingH2O,MinWindSpd,\
            MinimumKineticEnergy,LakeLayerThickness,GravConstant,ThicknessOfSurfaceSkin,\
            VonKarmanConst,LatentHeatOfVaporizationH2O,SpecificHeatOfAir,ThermalEmissivityOfH2O,\
            StefanBoltzmannConst,ThermalConductOfH2O,DeltaTimeStep,TotalKineticEnergyEff_CN,\
            TotalKineticEnergyEff_CF,TotalKineticEnergyEff_CE,TotalKineticEnergyEff_CS,TotalKineticEnergyEff_CL,\
            SpecificHeatOfH2O,MaxDeltaThermoclineDepth,MinimumMixingDepth,MaxDeltaMixingLayer,\
            MaxThermoclineThickness,MinThermoclineThickness,ThermalConductOfIce,HeatCapacityOfIce,\
            fID_71,fID_72,fID_73,fID_74,\
            fID_75,fID_76,fID_77,fID_82):
    #Break the refrences to the numpy arrays and just pass along values    
    LakeDepth = EnsurePassByValue(LakeDepth_Ref)
    LakeLength = EnsurePassByValue(LakeLength_Ref)
    ExtinctionCoeff = EnsurePassByValue(ExtinctionCoeff_Ref)
    NumOfLakeLayers = EnsurePassByValue(NumOfLakeLayers_Ref)    
    LakeTempProfile = EnsurePassByValue(LakeTempProfile_Ref)
    LakeSkinTemp = EnsurePassByValue(LakeSkinTemp_Ref)
    DepthToThermocline = EnsurePassByValue(DepthToThermocline_Ref)
    LakeIceHeight = EnsurePassByValue(LakeIceHeight_Ref)
    SnowIceHeight = EnsurePassByValue(SnowIceHeight_Ref)
    RunoffFromIce = EnsurePassByValue(RunoffFromIce_Ref)
    SnowPackMass = EnsurePassByValue(SnowPackMass_Ref)
    SnowPackDensity = EnsurePassByValue(SnowPackDensity_Ref)
    SnowPackTemp = EnsurePassByValue(SnowPackTemp_Ref)
    SnowPackAlbedo = EnsurePassByValue(SnowPackAlbedo_Ref)
    SnowPackLiqudWaterContent = EnsurePassByValue(SnowPackLiqudWaterContent_Ref)
    DiagAirTempAtScreenLevel = EnsurePassByValue(DiagAirTempAtScreenLevel_Ref)
    DiagZonalWindSpdAtAnemometerLevel = EnsurePassByValue(DiagZonalWindSpdAtAnemometerLevel_Ref)
    DiagMeridionalWindSpdAtAnemometerLevel = EnsurePassByValue(DiagMeridionalWindSpdAtAnemometerLevel_Ref)
    DiagSpecificHumidAtScreenLevel = EnsurePassByValue(DiagSpecificHumidAtScreenLevel_Ref)
    DiagRelavtiveHumidAtScreenLevel = EnsurePassByValue(DiagRelavtiveHumidAtScreenLevel_Ref)
    LakeSurfaceEvap = EnsurePassByValue(LakeSurfaceEvap_Ref)
    ExpansivityOfWater = EnsurePassByValue(ExpansivityOfWater_Ref)
    ThermoclineTempDelta = EnsurePassByValue(ThermoclineTempDelta_Ref)
    TotalKineticEnergy = EnsurePassByValue(TotalKineticEnergy_Ref)
    MeanMixLyrMomentum = EnsurePassByValue(MeanMixLyrMomentum_Ref)
    ReducedGravity = EnsurePassByValue(ReducedGravity_Ref)
    MixLyrDensity = EnsurePassByValue(MixLyrDensity_Ref)
    VisibleIncidentOnSurf = EnsurePassByValue(VisibleIncidentOnSurf_Ref)
    NIRIncidentOnSurf = EnsurePassByValue(NIRIncidentOnSurf_Ref)
    LongwaveDownwelling = EnsurePassByValue(LongwaveDownwelling_Ref)
    ZonalWindSpeed = EnsurePassByValue(ZonalWindSpeed_Ref)
    MeridionalWindSpeed = EnsurePassByValue(MeridionalWindSpeed_Ref)
    AirTempAtRef = EnsurePassByValue(AirTempAtRef_Ref)
    SpecificHumidityAtRef = EnsurePassByValue(SpecificHumidityAtRef_Ref)
    DensityOfAir = EnsurePassByValue(DensityOfAir_Ref)
    PartialPressureOfDryAir = EnsurePassByValue(PartialPressureOfDryAir_Ref)
    SurfaceAirPressure = EnsurePassByValue(SurfaceAirPressure_Ref)
    CosineOfSolarZenith = EnsurePassByValue(CosineOfSolarZenith_Ref)
    RefHeightWindSpeedForcing = EnsurePassByValue(RefHeightWindSpeedForcing_Ref)
    RefHeightAirTempAndHumidForcing = EnsurePassByValue(RefHeightAirTempAndHumidForcing_Ref)
    AnemometerHeight = EnsurePassByValue(AnemometerHeight_Ref)
    VarHeight = EnsurePassByValue(VarHeight_Ref)
    WetPrecip = EnsurePassByValue(WetPrecip_Ref)
    WetPrecipTemp = EnsurePassByValue(WetPrecipTemp_Ref)
    FrozenPrecip = EnsurePassByValue(FrozenPrecip_Ref)
    FrozenPrecipTemp = EnsurePassByValue(FrozenPrecipTemp_Ref)
    DensityOfFreshSnow = EnsurePassByValue(DensityOfFreshSnow_Ref)
    LatitudeInRad = EnsurePassByValue(LatitudeInRad_Ref)
    SnowAlbedoVis = EnsurePassByValue(SnowAlbedoVis_Ref)
    SnowAlbedoNIR = EnsurePassByValue(SnowAlbedoNIR_Ref)
    ShortwaveDownwellingIRandVIS = EnsurePassByValue(ShortwaveDownwellingIRandVIS_Ref)
    
    #Define Internal Arrays
    DelEnergyDueToWetPrecipIntoLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DelEnergyDueToSnowIntoLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DelEnergyDueToSnowIntoLakeFromMelting = np.zeros((NumOfTilesInAllGrids),np.float32)
    PondedWater = np.zeros((NumOfTilesInAllGrids),np.float32)
    FractionalIceCover = np.zeros((NumOfTilesInAllGrids),np.float32)
    TemperatureOfMixingLayerInCelcius = np.zeros((NumOfTilesInAllGrids),np.float32)
    WindSpeed2D = np.zeros((NumOfTilesInAllGrids),np.float32)
    DownwellingShortwaveTotal = np.zeros((NumOfTilesInAllGrids),np.float32)
    FrictionVelocity = np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeIceHeightAtStartOfTimestep = np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeTempProfileAtStartOfTimestep = np.zeros((MaxNumOfLakes),np.float32)    
    HeatCapacityOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)   
    SnowHeight = np.zeros((NumOfTilesInAllGrids),np.float32)  
    AlbedoWater = np.zeros((NumOfTilesInAllGrids),np.float32)  
    AlbedoIce = np.zeros((NumOfTilesInAllGrids),np.float32)  
    FractionOfSurfaceWithSnow = np.zeros((NumOfTilesInAllGrids),np.float32)  
    #SnowAlbedoVis = np.zeros((NumOfTilesInAllGrids,1),np.float32)  
    SnowAlbedoIR = np.zeros((NumOfTilesInAllGrids),np.float32)  
    GroundAlbedoVis = np.zeros((NumOfTilesInAllGrids),np.float32)  
    LakeAlbedoVis = np.zeros((NumOfTilesInAllGrids),np.float32)  
    EquationInterceptForSnowHeatFlux = np.zeros((NumOfTilesInAllGrids),np.float32)
    EquationMultiplierForSnowHeatFlux = np.zeros((NumOfTilesInAllGrids),np.float32)
    LatentHeatOfSublimationH2O = np.zeros((NumOfTilesInAllGrids),np.float32)
    ThermalConductOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    IntervalBtwAtmoBtmAndWindRef = np.zeros((NumOfTilesInAllGrids),np.float32)
    IntervalBtwAtmoBtmAndTempRef = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeightOfWindRefAboveZeroWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeightOfTempRefAboveZeroWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeightAboveAtmoBtmAndWindDiagnostic = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeightAboveAtmoBtmAndTempDiagnostic = np.zeros((NumOfTilesInAllGrids),np.float32)
    RatioOfRoughnessLengthToRefHightWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    RatioOfRoughnessLengthToRefTempWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    LogOfRoughnessLengthForMomentumOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    LogOfRoughnessLengthForHeatOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceRoughnessLengthForMomentum = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceRoughnessLengthForHeat = np.zeros((NumOfTilesInAllGrids),np.float32)
    VirtualPotentialTempOfAirAtRefHeight = np.zeros((NumOfTilesInAllGrids),np.float32)
    PotentialTempOfAir = np.zeros((NumOfTilesInAllGrids),np.float32)
    RichardsonNumberCoefficient = np.zeros((NumOfTilesInAllGrids),np.float32)
    StartPointForSurfTempIteration = np.zeros((NumOfTilesInAllGrids),np.float32)
    CoriolisForce = np.zeros((NumOfTilesInAllGrids),np.float32)
    FractionOfLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    MeanSurfaceTemperature = np.zeros((NumOfTilesInAllGrids),np.float32)
    TempAtBottomOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    FractionOfFrozenSurfaceAndSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    WetPrecipOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    WetPrecipOnSnowTemperature = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPrecipOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowPrecipOnSnowTemperature = np.zeros((NumOfTilesInAllGrids),np.float32)
    MaxSnowDensity = np.zeros((NumOfTilesInAllGrids),np.float32)
    RainfallRateForSnowAlbedoCalculations = np.zeros((NumOfTilesInAllGrids),np.float32)
    ResidualWaterUnavailableOnSurfaceForHousekeeping = np.zeros((NumOfTilesInAllGrids),np.float32)
    SN = np.zeros((NumOfTilesInAllGrids),np.float32)
    SL = np.zeros((NumOfTilesInAllGrids),np.float32)
    RN = np.zeros((NumOfTilesInAllGrids),np.float32)
    RL = np.zeros((NumOfTilesInAllGrids),np.float32)
    TempOfRunoff = np.zeros((NumOfTilesInAllGrids),np.float32)
    Runoff = np.zeros((NumOfTilesInAllGrids),np.float32)
    TempOfOverlandFlow = np.zeros((NumOfTilesInAllGrids),np.float32)
    OverlandFlow = np.zeros((NumOfTilesInAllGrids),np.float32)
    VolumetricLiquidWaterContentOfSoil = np.zeros((NumOfTilesInAllGrids,IGL),np.float32)
    ResidualSoilWaterPostFreezeOrEvap = np.zeros((NumOfTilesInAllGrids,IGL),np.float32)
    PermeableThicknessOfSoilLayer = np.zeros((NumOfTilesInAllGrids,IGL),np.float32)
    SurfaceConditionFlag = np.zeros((NumOfTilesInAllGrids),np.float32)
    # -- TSOLVE Output Arrays
    UpwellingLongwaveFromSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    SensibleHeatFromSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    LatentHeatFromSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    EvaporationRateAtSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    TemperatureAtSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    SpacificHumidityAtSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeatConductionIntoSnowPack = np.zeros((NumOfTilesInAllGrids),np.float32)
    AvailableEnergyForMeltingOfSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceDragCoefficentForHeat = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceDragCoefficentForMomentum = np.zeros((NumOfTilesInAllGrids),np.float32)
    BulkRichardsonNumber = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProductOfSurfaceDragAndWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProductOfSurfaceAirTempDragAndWind = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProductOfSurfaceAirHumidityDragAndWindSpeed = np.zeros((NumOfTilesInAllGrids),np.float32)
    InverseOfMoninObukhovRoughness = np.zeros((NumOfTilesInAllGrids),np.float32)
    FrictionVelocityOfAir = np.zeros((NumOfTilesInAllGrids),np.float32)
    HeightOfAtmoBoundaryLayer = np.zeros((NumOfTilesInAllGrids),np.float32)       
    CounterOfIterationsForEnergyBalance = np.zeros((NumOfTilesInAllGrids,6,50),np.float32)
    # -- Radiation Band-Dependent Arrays
    TransmissivityOfSnowForShortwave = np.zeros((NumOfTilesInAllGrids,NBS),np.float32)   
    AllWaveAlbedoOfSnowpack = np.zeros((NumOfTilesInAllGrids,NBS),np.float32)   
    ShortwaveTransmissivityOfSnow = np.zeros((NumOfTilesInAllGrids,NBS),np.float32)   
    # -- Internal Work Arrays
    TemperatureStepDeltaForTSOLVE = np.zeros((NumOfTilesInAllGrids),np.float32)
    VirtualPotentialTemperatureOfAirAtSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfaceEvaporationEfficiency = np.zeros((NumOfTilesInAllGrids),np.float32)
    InitialSpecificHumidityAtSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    TheResidualOfTheEnergyBalance = np.zeros((NumOfTilesInAllGrids),np.float32)
    DCFLXM = np.zeros((NumOfTilesInAllGrids),np.float32)
    CFLUXM = np.zeros((NumOfTilesInAllGrids),np.float32)
    MixingRatioAtSaturation = np.zeros((NumOfTilesInAllGrids),np.float32)
    FlagForLoopExit = np.zeros((NumOfTilesInAllGrids),np.float32)
    IterationCounterForTSOLVE = np.zeros((NumOfTilesInAllGrids),int)
    IndexOfCounterInTSOLVE = np.zeros((NumOfTilesInAllGrids),int)
    #Diagnostic Output Fields
    DRAGS = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagDelEnergyDueConducOrMassInLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    PrecipOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    PrecipOnLake = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagDelEnergyDueConducOrMassInSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPhaseChangeWaterInSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowSublimation = np.zeros((NumOfTilesInAllGrids),np.float32)
    SnowRunoff = np.zeros((NumOfTilesInAllGrids),np.float32)
    GZEROSL = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPhaseChangeWaterAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    OSW = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfDragCoeffForMomentum = np.zeros((NumOfTilesInAllGrids),np.float32)
    SurfDragCoeffForHeat = np.zeros((NumOfTilesInAllGrids),np.float32)
    LakeSurfaceDragAtNeutral = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSurfSpecificHumidity = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagTotalVisAtLakeSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetShortwaveAtSnowSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagNetLongwaveAtSnowSurface = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSensHeatOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagLatentHeatOnSnow = np.zeros((NumOfTilesInAllGrids),np.float32)
    LongwaveUpwelling = np.zeros((NumOfTilesInAllGrids),np.float32)
    SensibleHeatFromLakeSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProdLakeDragWindSpdDelT = np.zeros((NumOfTilesInAllGrids),np.float32)
    LatentHeatFromLakeSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    SensibleHeatFromSnowSubarea = np.zeros((NumOfTilesInAllGrids),np.float32)
    ProdLakeDragWindSpdDelSpHumid = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagPotentialEvap = np.zeros((NumOfTilesInAllGrids),np.float32)
    EvapEfficiencyAtLakeSurf = np.zeros((NumOfTilesInAllGrids),np.float32)
    DiagSurfBlkBodyTemp = np.zeros((NumOfTilesInAllGrids),np.float32)
    #Other
    SolarFluxInLayer = np.zeros((NumOfTilesInAllGrids,MaxNumOfLakes),np.float32)
    ConductiveHeatFluxInLayer = np.zeros((NumOfTilesInAllGrids,MaxNumOfLakes),np.float32)
    #Defintions
    ThicknessOfSediment= np.float32(10.0)                           # !THICKNESS OF SEDIMENT LAYER (m)  (STD)
    TemperatureConstantAtBottomOfSediment=FreezingPointOfH2O+6.0    # !CONSTANT TEMP AT BOTTOM OF SEDIMENT (K)  (STD)
    VolumetricHeatOfSediment= np.float32(2.13E6)                    # !VOLUMETRIC HEAT CAPACITY OF SEDIMENT     (STD)
    ThermalConductivityOfSediment= np.float32(0.0)                           # !THERMAL CONDUCTIVITY OF SEDIMENT (0 FOR ADIABATIC LOWER BC)
    FractionOfShortwaveAtTopOfSediment= np.float32(0.0)                      # !FRACTION OF NET SW THAT REACHES SEDIMENTS (0 to turn off)
    WindShelteringFactor= np.float32(1.0)                                    # !Wind sheltering factor reduces surface drag (1.0 turns off)
    IsSnow= np.float32(1)
    MaxSnowHeightBeforeFullCover= np.float32(0.10)
    DensityRatioOfIceToWater=DensityOfIce/DensityOfH2O 
    
    # -- INITIAL CONDITIONS (IterationCounter=1)
    #    * EQUATION OF STATE FROM FARMER AND CARMACK (1981,JPO), BUT NEGLECTS SALINITY AND PRESSURE EFFECTS   
    if IterationCounter == 1:        
        for I in range(IL1-1,IL2):
            LakeIceHeight[I]=0.0 #                  !INITIAL ICE COVER SET TO ZERO
            SnowIceHeight[I]=0.0 #                  !INITIAL SNOW ICE  SET TO ZERO
            RunoffFromIce[I]=0.0 #                  !INITIAL runoff ICE  SET TO ZERO
            LakeSkinTemp[I]=LakeTempProfile[I,0]#   !INITIAL SKIN TEMP SET TO FIRST LAYER TEMP
            NLEV=NumOfLakeLayers[I]
            SedimentTemp[I]=LakeTempProfile[I,NLEV-1]#!INITIAL SEDIMENT TEMP SET TO LOWEST LAYER TEMP
            TotalKineticEnergy[I]=MinimumKineticEnergy
            MeanMixLyrMomentum[I]=0.0
            #---* INITIAL MIXED LAYER DEPTH ESTIMATED BASED ON INITIAL T PROFILE
            for J in range(0,int(NumOfLakeLayers[I])-1):
                JMIX=J
                TTOP=LakeTempProfile[I,J]
                TBOT=LakeTempProfile[I,J+1]
                if TTOP-TBOT > 1.0:
                    #Exit loop
                    break
            DepthToThermocline[I]=LakeLayerThickness*np.float32(JMIX+1)#zero indexing requires plus one here.
            ThermoclineTempDelta[I]=LakeTempProfile[I,JMIX]-LakeTempProfile[I,JMIX+1]# !see Spigel et al 1986
        # ---* INITIAL MIXED LAYER TEMP (CELSIUS), DENSITY, EXPANSIVITY, AND REDUCED GRAVITY        
            TemperatureOfMixingLayerInCelcius[I]=( (LakeTempProfile[I,0]+LakeTempProfile[I,JMIX])/2.0) - FreezingPointOfH2O
            [ExpansivityOfWater[I],MixLyrDensity[I]]=EQNST(TemperatureOfMixingLayerInCelcius[I],0.0)
            ReducedGravity[I] = GravConstant*abs(DensityOfH2O-MixLyrDensity[I])/DensityOfH2O
            if (ReducedGravity[I] < 0): #this is impossible because of the abs, unless G is pointed upsidedown.
                raise Exception("RUNLAKE -- In CLASSL ReducedGravity < 0")
    # -- ICE COVER AND SNOW PARAMETERS AT BEGINNING OF TIMESTEP    
    for I in range(IL1-1,IL2):
        PondedWater[I]=0.0
        LakeIceHeightAtStartOfTimestep[I]=LakeIceHeight[I]
        ICELIM=LakeLength[I]*1.0E-5 #ICELIM is "limiting ice thickness" in Mackay et al 2017       
        FractionalIceCover[I]=MIN(ICELIM,LakeIceHeight[I])/ICELIM   
        FractionOfLakeSurface[I]=1.0
        if SnowPackMass[I]>0.0:
            SnowHeight[I]=SnowPackMass[I]/SnowPackDensity[I]
            if SnowHeight[I]>=(MaxSnowHeightBeforeFullCover-0.00001):
                FractionOfSurfaceWithSnow[I]=1.0            
            else:
                FractionOfSurfaceWithSnow[I]=SnowHeight[I]/MaxSnowHeightBeforeFullCover
                SnowHeight[I]=MaxSnowHeightBeforeFullCover
                SnowPackLiqudWaterContent[I]=SnowPackLiqudWaterContent[I]/FractionOfSurfaceWithSnow[I]
        else:
            SnowHeight[I]=0.0
            FractionOfSurfaceWithSnow[I]=0.0
        if FrozenPrecip[I]>0.0:
            if LakeIceHeight[I]>=ThicknessOfSurfaceSkin:
                if FractionOfSurfaceWithSnow[I]>0.0:
                    SN[I]=FrozenPrecip[I]*FractionalIceCover[I]/FractionOfSurfaceWithSnow[I]
                else:
                    SN[I]=FrozenPrecip[I]
                SL[I]=(1.0-FractionalIceCover[I])*FrozenPrecip[I]
            else:
                SN[I]=0.0
                SL[I]=FrozenPrecip[I]
        else:
            SN[I]=0.0
            SL[I]=0.0        
        if WetPrecip[I]>0.0:
            if FractionOfSurfaceWithSnow[I]>0.0:
                RN[I]=WetPrecip[I]
                RL[I]=(1.0-FractionOfSurfaceWithSnow[I])*WetPrecip[I]
            else:
                RN[I]=0.0
                RL[I]=WetPrecip[I]
        else:
            RN[I]=0.0
            RL[I]=0.0
        if FractionOfSurfaceWithSnow[I]>0:
            FractionOfFrozenSurfaceAndSnow[I]=FractionOfSurfaceWithSnow[I]
            StartPointForSurfTempIteration[I]=SnowPackTemp[I]
        else:
            FractionOfFrozenSurfaceAndSnow[I]=FractionalIceCover[I]
            StartPointForSurfTempIteration[I]=FreezingPointOfH2O
        PrecipOnSnow[I] =FractionOfSurfaceWithSnow[I]*DensityOfH2O*RN[I]+FractionOfFrozenSurfaceAndSnow[I]*DensityOfFreshSnow[I]*SN[I]
        PrecipOnLake[I] =DensityOfH2O*RL[I]+DensityOfFreshSnow[I]*SL[I]
        DiagDelEnergyDueConducOrMassInSnow[I]=0.0
        DiagPhaseChangeWaterInSnow[I]=0.0
        SnowSublimation[I]=0.0
        SnowRunoff[I]=0.0
        GZEROSL[I]=0.0
        AvailableEnergyForMeltingOfSnow[I]=0.0
        DelEnergyDueToWetPrecipIntoLake[I]=(WetPrecipTemp[I]+FreezingPointOfH2O-LakeSkinTemp[I])*RL[I]*HeatCapacityOfH2O
        DelEnergyDueToSnowIntoLake[I]=(FrozenPrecipTemp[I]+FreezingPointOfH2O-LakeSkinTemp[I])*SL[I]*SpecificHeatOfIce*DensityOfFreshSnow[I]
        DelEnergyDueToSnowIntoLakeFromMelting[I]=LatentHeatOfFreezingH2O*SL[I]*DensityOfFreshSnow[I]
        DiagDelEnergyDueConducOrMassInLake[I]=DelEnergyDueToWetPrecipIntoLake[I]+DelEnergyDueToSnowIntoLake[I]-DelEnergyDueToSnowIntoLakeFromMelting[I]
        # -- INITIALIZE OTHER SNOW-RELATED VARIABLES.
        SnowPrecipOnSnow  [I]=0
        SnowPrecipOnSnowTemperature [I]=0
        WetPrecipOnSnow  [I]= np.float32(0)
        WetPrecipOnSnowTemperature [I]= 0
        SurfaceDragCoefficentForMomentum  [I]=0
        SurfaceDragCoefficentForHeat  [I]=0
        DRAGS [I]=0
        TemperatureAtSurface[I]=0
        SpacificHumidityAtSurface[I]=0
        SnowAlbedoVis[I]=np.float32(0)
        SnowAlbedoIR[I]=0
        UpwellingLongwaveFromSurface [I]=0
        SensibleHeatFromSurface[I]=0
        LatentHeatFromSurface[I]=0
        EvaporationRateAtSurface [I]=0
        HeatCapacityOfSnow[I]=0
        ThermalConductOfSnow[I]=0
        HeatConductionIntoSnowPack [I]=0         
    # --  WIND SPEED AND PRECIP
    for I in range(IL1-1,IL2):
        WindSpeed2D[I]=MAX(MinWindSpd,math.sqrt(ZonalWindSpeed[I]*ZonalWindSpeed[I]+MeridionalWindSpeed[I]*MeridionalWindSpeed[I]))
        CoriolisForce[I]=2.0*7.29E-5*math.sin(LatitudeInRad[I])
    # -- DEFINE ICE AND WATER ALBEDOS, USED IN TSOLVL.
    #     * STARTING "GROUND ALBEDOS" FOR ROUTINE SNOALBA.
    for I in range(IL1-1,IL2):
        AlbedoWater[I]=0.045/MAX(CosineOfSolarZenith[I],0.1)# !std value for water        
        AlbedoIce[I]=0.08+0.44*(LakeIceHeight[I])**0.28   #    !thin ice albedo (Vavrus et al 1996)
        AlbedoIce[I]=MIN(AlbedoIce[I],0.44)        
        if LakeIceHeight[I]>0.0:
            GroundAlbedoVis[I]=AlbedoIce[I]
        else:
            GroundAlbedoVis[I]=0.0
    # -- SNOW PACK ENERGY AND WATER BALANCE
    [SnowAlbedoVis,SnowAlbedoIR,AllWaveAlbedoOfSnowpack, TransmissivityOfSnowForShortwave] = \
        SNOALBA(SnowAlbedoVis,SnowAlbedoIR,SnowPackAlbedo, AllWaveAlbedoOfSnowpack, TransmissivityOfSnowForShortwave,\
        SnowHeight,FractionOfSurfaceWithSnow,SnowAlbedoVis,SnowAlbedoNIR,IL1,IL2,NumOfGridCells,IALS,NBS,ISNOALB)
    [EquationMultiplierForSnowHeatFlux,EquationInterceptForSnowHeatFlux,LatentHeatOfSublimationH2O,ThermalConductOfSnow,HeatCapacityOfSnow,\
        SurfaceConditionFlag,IntervalBtwAtmoBtmAndWindRef,IntervalBtwAtmoBtmAndTempRef,HeightOfWindRefAboveZeroWind,HeightOfTempRefAboveZeroWind,\
        HeightAboveAtmoBtmAndWindDiagnostic,HeightAboveAtmoBtmAndTempDiagnostic,RatioOfRoughnessLengthToRefHightWind,RatioOfRoughnessLengthToRefTempWind,\
        LogOfRoughnessLengthForMomentumOfSnow,LogOfRoughnessLengthForHeatOfSnow,SurfaceRoughnessLengthForMomentum,SurfaceRoughnessLengthForHeat,\
        VirtualPotentialTempOfAirAtRefHeight,PotentialTempOfAir,RichardsonNumberCoefficient,DRAGS,CEVAP,EvaporationFlag,SandContentFlag] = \
        TLSPREP(EquationMultiplierForSnowHeatFlux,EquationInterceptForSnowHeatFlux,LatentHeatOfSublimationH2O,ThermalConductOfSnow,HeatCapacityOfSnow,\
        SurfaceConditionFlag,IntervalBtwAtmoBtmAndWindRef,IntervalBtwAtmoBtmAndTempRef,HeightOfWindRefAboveZeroWind,
        HeightOfTempRefAboveZeroWind,HeightAboveAtmoBtmAndWindDiagnostic,HeightAboveAtmoBtmAndTempDiagnostic,RatioOfRoughnessLengthToRefHightWind,\
        RatioOfRoughnessLengthToRefTempWind,LogOfRoughnessLengthForMomentumOfSnow,LogOfRoughnessLengthForHeatOfSnow,SurfaceRoughnessLengthForMomentum,\
        SurfaceRoughnessLengthForHeat,VirtualPotentialTempOfAirAtRefHeight,PotentialTempOfAir,RichardsonNumberCoefficient,DRAGS,FractionOfSurfaceWithSnow,\
        SnowHeight,SnowPackTemp,SnowPackDensity,SnowPackLiqudWaterContent,RefHeightWindSpeedForcing,RefHeightAirTempAndHumidForcing,AnemometerHeight,\
        VarHeight,AirTempAtRef,SpecificHumidityAtRef,WindSpeed2D,IZREF,IL1,IL2,IGL,GravConstant,SpecificHeatOfAir,VonKarmanConst,HeatCapacityOfIce,\
        DensityOfIce,HeatCapacityOfH2O,DensityOfH2O,LatentHeatOfVaporizationH2O,LatentHeatOfFreezingH2O)
    [QSWNS,UpwellingLongwaveFromSurface,QTRANSL,SensibleHeatFromSurface,LatentHeatFromSurface,EvaporationRateAtSurface,TemperatureAtSurface,\
        SpacificHumidityAtSurface,HeatConductionIntoSnowPack,AvailableEnergyForMeltingOfSnow,SurfaceDragCoefficentForHeat,SurfaceDragCoefficentForMomentum,\
        BulkRichardsonNumber,ProductOfSurfaceDragAndWind,ProductOfSurfaceAirTempDragAndWind,ProductOfSurfaceAirHumidityDragAndWindSpeed,\
        InverseOfMoninObukhovRoughness,FrictionVelocityOfAir,HeightOfAtmoBoundaryLayer,HeightOfTempRefAboveZeroWind,HeightOfWindRefAboveZeroWind,\
        SurfaceRoughnessLengthForHeat,SurfaceRoughnessLengthForMomentum,CoriolisForce,NumOfTilesInAllGrids,CounterOfIterationsForEnergyBalance,\
        TemperatureStepDeltaForTSOLVE,VirtualPotentialTemperatureOfAirAtSurface,SurfaceEvaporationEfficiency,InitialSpecificHumidityAtSurface,\
        TheResidualOfTheEnergyBalance,DCFLXM,CFLUXM,MixingRatioAtSaturation,ShortwaveTransmissivityOfSnow,FlagForLoopExit,IterationCounterForTSOLVE,\
        JEVAP,IndexOfCounterInTSOLVE] = \
        TSOLVE(IsSnow,FractionOfSurfaceWithSnow,UpwellingLongwaveFromSurface,SensibleHeatFromSurface,LatentHeatFromSurface,EvaporationRateAtSurface,\
        TemperatureAtSurface,SpacificHumidityAtSurface,HeatConductionIntoSnowPack,AvailableEnergyForMeltingOfSnow,SurfaceDragCoefficentForHeat,\
        SurfaceDragCoefficentForMomentum,BulkRichardsonNumber,ProductOfSurfaceDragAndWind,ProductOfSurfaceAirTempDragAndWind,\
        ProductOfSurfaceAirHumidityDragAndWindSpeed,InverseOfMoninObukhovRoughness,FrictionVelocityOfAir,HeightOfAtmoBoundaryLayer,\
        HeightOfTempRefAboveZeroWind,HeightOfWindRefAboveZeroWind,SurfaceRoughnessLengthForHeat,SurfaceRoughnessLengthForMomentum,CoriolisForce,\
        LongwaveDownwelling,PotentialTempOfAir,SpecificHumidityAtRef,WindSpeed2D,PartialPressureOfDryAir,DensityOfAir,SnowAlbedoVis,SnowAlbedoIR,\
        RichardsonNumberCoefficient,LatentHeatOfSublimationH2O,CEVAP,VirtualPotentialTempOfAirAtRefHeight,RatioOfRoughnessLengthToRefTempWind,\
        RatioOfRoughnessLengthToRefHightWind,EquationInterceptForSnowHeatFlux,EquationMultiplierForSnowHeatFlux,StartPointForSurfTempIteration,\
        TransmissivityOfSnowForShortwave,ShortwaveDownwellingIRandVIS,AllWaveAlbedoOfSnowpack,VolumetricLiquidWaterContentOfSoil,\
        ResidualSoilWaterPostFreezeOrEvap,PermeableThicknessOfSoilLayer,SnowPackDensity,SnowHeight,PondedWater,SurfaceConditionFlag,\
        EvaporationFlag,CounterOfIterationsForEnergyBalance,SandContentFlag,ISLFD,ITG,NumOfTilesInAllGrids,IL1,IL2,NBS,ISNOALB,\
        TemperatureStepDeltaForTSOLVE,VirtualPotentialTemperatureOfAirAtSurface,SurfaceEvaporationEfficiency,InitialSpecificHumidityAtSurface,\
        TheResidualOfTheEnergyBalance,DCFLXM,CFLUXM,MixingRatioAtSaturation,ShortwaveTransmissivityOfSnow,FlagForLoopExit,IterationCounterForTSOLVE,\
        IndexOfCounterInTSOLVE,DeltaTimeStep,FreezingPointOfH2O,DensityOfH2O,StefanBoltzmannConst,SpecificHeatOfAir,GravConstant,VonKarmanConst)
    [HeatConductionIntoSnowPack,SnowPackTemp,SnowPackLiqudWaterContent,SnowPackDensity,AvailableEnergyForMeltingOfSnow,GZEROSL,TempAtBottomOfSnow,\
        DiagDelEnergyDueConducOrMassInSnow,DiagPhaseChangeWaterInSnow,SnowSublimation,EvaporationRateAtSurface,WetPrecipOnSnow,WetPrecipOnSnowTemperature,\
        SnowPrecipOnSnow,SnowPrecipOnSnowTemperature,HeatCapacityOfSnow] = \
        TLSPOST(HeatConductionIntoSnowPack,SnowPackTemp,SnowPackLiqudWaterContent,SnowPackDensity,AvailableEnergyForMeltingOfSnow,GZEROSL,TempAtBottomOfSnow,\
        DiagDelEnergyDueConducOrMassInSnow,DiagPhaseChangeWaterInSnow,EvaporationRateAtSurface,LakeSkinTemp,SnowHeight,ThermalConductOfSnow,HeatCapacityOfSnow,\
        RN,WetPrecipTemp,SN,FrozenPrecipTemp,TemperatureAtSurface,DensityOfFreshSnow,FractionOfSurfaceWithSnow,ThicknessOfSurfaceSkin,IL1,IL2,DensityOfH2O,\
        DeltaTimeStep,FreezingPointOfH2O,LatentHeatOfFreezingH2O,HeatCapacityOfIce,DensityOfIce,HeatCapacityOfH2O)
    [SnowPackDensity,SnowHeight,HeatCapacityOfSnow,SnowPackTemp,EvaporationRateAtSurface,SnowSublimation,LakeSurfaceEvap,DiagDelEnergyDueConducOrMassInSnow,\
        ResidualWaterUnavailableOnSurfaceForHousekeeping,TempOfRunoff,Runoff,TempOfOverlandFlow,OverlandFlow, SnowPackLiqudWaterContent] = \
        SNOVAP(SnowPackDensity,SnowHeight,HeatCapacityOfSnow,SnowPackTemp,EvaporationRateAtSurface,SnowSublimation,LakeSurfaceEvap,\
        DiagDelEnergyDueConducOrMassInSnow,ResidualWaterUnavailableOnSurfaceForHousekeeping,TempOfRunoff,Runoff,TempOfOverlandFlow,OverlandFlow,\
        FractionOfSurfaceWithSnow,WetPrecipOnSnow,SnowPrecipOnSnow,DensityOfFreshSnow,SnowPackLiqudWaterContent,IL1,IL2,\
        DeltaTimeStep,FreezingPointOfH2O,DensityOfH2O,HeatCapacityOfIce,DensityOfIce,HeatCapacityOfH2O)
    [SnowHeight,SnowPackTemp,AvailableEnergyForMeltingOfSnow,WetPrecipOnSnow,WetPrecipOnSnowTemperature,GZEROSL,RainfallRateForSnowAlbedoCalculations,\
        DiagPhaseChangeWaterInSnow,DiagDelEnergyDueConducOrMassInSnow,HeatCapacityOfSnow,SnowPackLiqudWaterContent] = \
        TMELT(SnowHeight,SnowPackTemp,AvailableEnergyForMeltingOfSnow,WetPrecipOnSnow,WetPrecipOnSnowTemperature,GZEROSL,RainfallRateForSnowAlbedoCalculations,\
        DiagPhaseChangeWaterInSnow,DiagDelEnergyDueConducOrMassInSnow,FractionOfSurfaceWithSnow,HeatCapacityOfSnow,SnowPackDensity,SnowPackLiqudWaterContent,\
        SandContentFlag,IL1,IL2,DeltaTimeStep,FreezingPointOfH2O,DensityOfH2O,HeatCapacityOfIce,DensityOfIce,HeatCapacityOfH2O,LatentHeatOfFreezingH2O)
    [WetPrecipOnSnow,WetPrecipOnSnowTemperature,SnowHeight,SnowPackTemp,SnowPackDensity,HeatCapacityOfSnow,SnowPackLiqudWaterContent,\
        DiagDelEnergyDueConducOrMassInSnow,DiagPhaseChangeWaterInSnow,PrecipOnLake,SnowRunoff] = \
        SNINFL(WetPrecipOnSnow,WetPrecipOnSnowTemperature,SnowHeight,SnowPackTemp,SnowPackDensity,HeatCapacityOfSnow,SnowPackLiqudWaterContent,\
        DiagDelEnergyDueConducOrMassInSnow,DiagPhaseChangeWaterInSnow,PrecipOnLake,SnowRunoff,FractionOfSurfaceWithSnow,NumOfTilesInAllGrids,IL1,IL2,\
        NumOfGridCells,FreezingPointOfH2O,DeltaTimeStep,HeatCapacityOfIce,DensityOfIce,HeatCapacityOfH2O,DensityOfH2O,LatentHeatOfFreezingH2O)
    [SnowPackAlbedo,SnowPackDensity,SnowHeight,HeatCapacityOfSnow] = \
        SNOALBW(SnowPackAlbedo,SnowPackDensity,SnowHeight,HeatCapacityOfSnow,SnowPackTemp,FractionOfSurfaceWithSnow,SnowPrecipOnSnow,\
        RainfallRateForSnowAlbedoCalculations,SnowPackLiqudWaterContent,MaxSnowDensity,IL1,IL2,DeltaTimeStep,DensityOfH2O,HeatCapacityOfIce,DensityOfIce,\
        HeatCapacityOfH2O)            
    [SnowPackAlbedo,SnowPackTemp,SnowPackDensity,SnowHeight,HeatCapacityOfSnow,DiagDelEnergyDueConducOrMassInSnow] = \
        SNOADD(SnowPackAlbedo,SnowPackTemp,SnowPackDensity,SnowHeight,HeatCapacityOfSnow,DiagDelEnergyDueConducOrMassInSnow,FractionOfFrozenSurfaceAndSnow,\
        SnowPrecipOnSnow,SnowPrecipOnSnowTemperature,DensityOfFreshSnow,SnowPackLiqudWaterContent,IL1,IL2,FreezingPointOfH2O,DeltaTimeStep,HeatCapacityOfIce,\
        DensityOfIce,HeatCapacityOfH2O,DensityOfH2O)
    #Deal with snow
    for I in range(IL1-1,IL2):
        if SnowHeight[I]>0.0:
            SnowPackTemp[I]=SnowPackTemp[I]+FreezingPointOfH2O
            SnowPackMass[I]=SnowHeight[I]*SnowPackDensity[I]*FractionOfFrozenSurfaceAndSnow[I]
            SnowPackLiqudWaterContent[I]=SnowPackLiqudWaterContent[I]*FractionOfSurfaceWithSnow[I]
        else:
            SnowPackTemp[I]=0.0
            SnowPackDensity[I]=0.0
            SnowPackMass[I]=0.0
            SnowPackLiqudWaterContent[I]=0.0
        # -- Melt last bit of snow if below threshold or no ice exists        
        if SnowPackMass[I]>0.0 and (SnowPackMass[I]<1.0E-6 or LakeIceHeight[I] <=0.01):
            SnowRunoff[I]=SnowRunoff[I]+(SnowPackMass[I]+SnowPackLiqudWaterContent[I])/DeltaTimeStep
            PrecipOnLake[I]=PrecipOnLake[I]+(SnowPackMass[I]+SnowPackLiqudWaterContent[I])/DeltaTimeStep
            SL[I]=SL[I]+SnowPackMass[I]/(SnowPackDensity[I]*DeltaTimeStep)
            RL[I]=RL[I]+SnowPackLiqudWaterContent[I]/(DensityOfH2O*DeltaTimeStep)
            DiagDelEnergyDueConducOrMassInSnow[I]=DiagDelEnergyDueConducOrMassInSnow[I]-SnowPackTemp[I]*(SpecificHeatOfIce*SnowPackMass[I]+SpecificHeatOfH2O*SnowPackLiqudWaterContent[I])/DeltaTimeStep
            SnowPackTemp[I]=0.0
            SnowPackDensity[I]=0.0
            SnowPackMass[I]=0.0
            SnowPackLiqudWaterContent[I]=0.0
        # Heat capacity calculation for use if necessary by: (1) heat 
        # added from snowmelt/rain runoff and (2) snow ice production
        # Any heat added goes into layer 1, not the skin layer.
        ICEBOT=DensityRatioOfIceToWater*LakeIceHeightAtStartOfTimestep[I]# !bottom of floating ice
        ZTOP=ThicknessOfSurfaceSkin                                      # !top of first layer
        ZBOT=ThicknessOfSurfaceSkin + LakeLayerThickness                 # !bottom of first layer
        if ICEBOT >= ZBOT: 
            HCAP=HeatCapacityOfIce
        elif ICEBOT <= ZTOP:
            HCAP=HeatCapacityOfH2O
        else: 
            Z=ICEBOT-ZTOP
            HCAP=(Z*HeatCapacityOfIce + (LakeLayerThickness-Z)*HeatCapacityOfH2O)/LakeLayerThickness
        
        # Add heat of snow melt/rain reaching bottom of snowcover to ice
        # Heat added to first layer below skin
        # SnowRunoff cools then runs off - ie does not add to LakeIceHeight
        #
        #-----------------------------------------------------------
        # Runoff ICE PRODUCTION/HEATING TURNED OFF: JUNE 2016, M.MACKAY
        #
        if (SnowRunoff[I]>0.0 and LakeIceHeightAtStartOfTimestep[I]>0.0):
            ROFICE=0.0 #  !runoff ice production turned off
            # Cool SnowRunoff to FreezingPointOfH2O, then add heat to ice layer below skin
            HWARM=HeatCapacityOfH2O*(WetPrecipOnSnowTemperature[I]-0.0)*DeltaTimeStep*SnowRunoff[I]/DensityOfH2O
            HFREZ=DensityOfIce*LatentHeatOfFreezingH2O*ROFICE
            LakeTempProfile[I,0]=LakeTempProfile[I,0] + (HFREZ+HWARM)/(HCAP*LakeLayerThickness)
            #mdm         LakeIceHeight[I]=LakeIceHeight[I]+ROFICE
        #======================================================================
        # SNOW ICE
        #---------
        # Produce snow-ice if weight of snow cover sufficient
        # Pore volume in snow fills with water then freezes
        #        
        if(SnowPackMass[I]>0.0): 
            SnowRunoff[I]=SnowRunoff[I]+SnowPackMass[I]/DeltaTimeStep
            DiagDelEnergyDueConducOrMassInSnow[I]=DiagDelEnergyDueConducOrMassInSnow[I]-HeatCapacityOfSnow[I]*SnowPackTemp[I]*SnowPackMass[I]/(SnowPackDensity[I]*DeltaTimeStep)        
        BBLFAC=np.float32(1.0)  # !non-bubble fraction in snow ice
        ICECAP =LakeIceHeightAtStartOfTimestep[I]*(DensityOfH2O-DensityOfIce) # !snow holding capacity of ice
        ZBOT = MIN(10.0, (NumOfLakeLayers[I]-1.)*LakeLayerThickness)# !limit snow ice production to 10m or depth of lake (less 0.5m)  
            
        if (SnowPackMass[I]>ICECAP and LakeIceHeightAtStartOfTimestep[I]>0.0 and LakeIceHeightAtStartOfTimestep[I]<=ZBOT):
            ICE0=LakeIceHeightAtStartOfTimestep[I]  #   !initial ice thickness
            SNO0=SnowPackMass[I]                    #   !initial snow mass
            ZSNOW0=SnowPackMass[I]/SnowPackDensity[I]#  !initial mean snow depth 
            PORE=(DensityOfIce-SnowPackDensity[I])/DensityOfIce#  !pore volume in snow
            ALPHA=(DensityOfH2O*PORE*BBLFAC + DensityOfIce*(1.0-PORE))/DensityOfIce
            ETA=(DensityOfH2O-DensityOfIce)/SnowPackDensity[I]
            # Final snow cover mass equals new ice holding capacity
            SnowPackMass[I]=SnowPackDensity[I]*ETA*(ICE0+ALPHA*ZSNOW0)/(1.0+ALPHA*ETA)
            LakeIceHeight[I]=SnowPackMass[I]/(SnowPackDensity[I]*ETA)
            # Mass of liquid water that freezes in slush layer
            MASSL=DensityOfH2O*PORE*BBLFAC*(SNO0-SnowPackMass[I])/SnowPackDensity[I]
            # Cumulative snow-ice diagnostic
            SnowIceHeight[I]=SnowIceHeight[I]+LakeIceHeight[I]-ICE0
            # Add heat of fusion to snow pack
            # First warm snow to FreezingPointOfH2O, then freeze lake water
            # Lake water is assumed to be at FreezingPointOfH2O already (ie comes from layer that
            # contains ICEBOT, which is at FreezingPointOfH2O)
            HWARM=HeatCapacityOfSnow[I]*(FreezingPointOfH2O-SnowPackTemp[I])*(SNO0-SnowPackMass[I])/SnowPackDensity[I]
            HFREZ=LatentHeatOfFreezingH2O*MASSL
            HFLX=HFREZ-HWARM
            # Latent heat goes goes into snowpack unless snow too thin        
            if (SnowPackMass[I] >= 5.0E-3):
            # Partition latent heat between ice and snow based on heat capacities
                HFLXS=HFLX  #          !all latent heat goes into snow cover
                HFLXI=HFLX-HFLXS                
                LakeTempProfile[I,0]=LakeTempProfile[I,0] + HFLXI/(HCAP*LakeLayerThickness)#!mdm test
                SnowPackTemp[I]=SnowPackTemp[I]+HFLXS/(HeatCapacityOfSnow[I]*SnowPackMass[I]/SnowPackDensity[I])
                if SnowPackTemp[I]>FreezingPointOfH2O: 
                    QMLTS=(SnowPackTemp[I]-FreezingPointOfH2O)*HeatCapacityOfSnow[I]*SnowPackMass[I]/SnowPackDensity[I]# !J/m2
                    if QMLTS > SnowPackMass[I]*LatentHeatOfFreezingH2O:
                        # All snow melts. Excess heat put in lake
                        print('\nall snow melts during flooding')
                        HFLX=QMLTS-SnowPackMass[I]*LatentHeatOfFreezingH2O                        
                        LakeTempProfile[I,0]=LakeTempProfile[I,0] + HFLX/(HCAP*LakeLayerThickness)
                        SnowPackMass[I]=0.0
                        SnowPackTemp[I]=0.0
                        SnowPackDensity[I]=0.0
                        SnowPackLiqudWaterContent [I]=0.0
                        SnowHeight [I]=0.0
                    else:
                        # Some snow melts
                        SnowPackMass[I]=SnowPackMass[I]-QMLTS/LatentHeatOfFreezingH2O
                        SnowPackTemp[I]=FreezingPointOfH2O           
            else:
                print('\nno snow: all heat into lake')                
                LakeTempProfile[I,0]=LakeTempProfile[I,0] + (HFREZ-HWARM)/(HCAP*LakeLayerThickness)
        if SnowPackMass[I]>0.0:
            SnowRunoff[I]=SnowRunoff[I]-SnowPackMass[I]/DeltaTimeStep
            DiagDelEnergyDueConducOrMassInSnow[I]=DiagDelEnergyDueConducOrMassInSnow[I]+HeatCapacityOfSnow[I]*SnowPackTemp[I]*SnowPackMass[I]/(SnowPackDensity[I]*DeltaTimeStep)
           
    # -- LAKE LakeSurfaceDragAtNeutral COEFFICIENTS 
    [CDML,CDHL,DRAGL] =\
        DRCOEFL(WindSpeed2D,LakeSkinTemp,AirTempAtRef,SpecificHumidityAtRef,SurfaceAirPressure,RefHeightWindSpeedForcing,\
        RefHeightAirTempAndHumidForcing,FractionalIceCover,FractionOfSurfaceWithSnow,IL1,IL2,GravConstant,VonKarmanConst,FreezingPointOfH2O)
    # -- TOTAL SOLAR RADIATION, AND WATER FRICTION VELOCITY 
    #       Reduce momentum transfer due to sheltering
    for I in range(IL1-1,IL2):
        DownwellingShortwaveTotal [I] = VisibleIncidentOnSurf[I] + NIRIncidentOnSurf[I]
        CDML[I]=CDML[I]*WindShelteringFactor
        if LakeIceHeight[I] <= 0: 
            FrictionVelocity[I]=WindSpeed2D[I]*math.sqrt(CDML[I]*DensityOfAir[I]/DensityOfH2O)
        else:
            FrictionVelocity[I]=0.0#            !no wind stress through ice
        
    [LakeSkinTemp,LakeTempProfile]=\
        FREECONV(LakeIceHeight,LakeSkinTemp,LakeTempProfile,DensityRatioOfIceToWater,NumOfLakeLayers,IL1,IL2,\
        FreezingPointOfH2O,ThicknessOfSurfaceSkin,LakeLayerThickness)
      
    # -- COMPUTE WATER EXTINCTION COEFFICIENTS
    [CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI] = LKTRANS(ExtinctionCoeff,IL1,IL2)
    # -- COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP AND STEP
    # -- FORWARD SKIN TEMP LakeSkinTemp.  COMPUTE ICE COVER IN SKIN LAYER.
    [LakeSkinTemp,LakeIceHeight,Q0,F0,QSURFL,DiagNetShortwaveAtLakeSurface,DiagNetLongwaveAtLakeSurface,DiagSensHeatOnLake,\
        DiagLatentHeatOnLake,LakeSurfaceEvap,LakeAlbedoVis,G0] = \
        TSOLVL(LakeTempProfile,LakeSkinTemp,LakeIceHeight,LakeSurfaceEvap,LakeAlbedoVis,\
        DownwellingShortwaveTotal ,LongwaveDownwelling,AirTempAtRef,SpecificHumidityAtRef,WindSpeed2D,SurfaceAirPressure,DensityOfAir,CDHL,\
        GZEROSL,QTRANSL,DiagDelEnergyDueConducOrMassInLake,FractionOfSurfaceWithSnow,FractionalIceCover,AlbedoWater,AlbedoIce,\
        CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI,IL1,IL2,\
        LakeLayerThickness,ThicknessOfSurfaceSkin,DensityOfIce,DensityOfH2O,FreezingPointOfH2O,\
        LatentHeatOfVaporizationH2O,SpecificHeatOfAir,ThermalEmissivityOfH2O,StefanBoltzmannConst,ThermalConductOfH2O,DeltaTimeStep,HeatCapacityOfH2O,\
        LatentHeatOfFreezingH2O,ThermalConductOfIce,HeatCapacityOfIce)
    # -- COMPUTE SOLAR FLUX (SolarFluxInLayer) FOR CURRENT TIMESTEP
    #    include attenuation through ice (including leads) if present
    for I in range(IL1-1,IL2):
        ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]
        ICETOP=LakeIceHeight[I]-ICEBOT
        for J in range(0,int(NumOfLakeLayers[I])):
            Z=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
            if LakeIceHeight[I] <= 0.0: #       !NO ICE 
                ATTEN1=CQ1B[I]*Z
                ATTEN2=CQ2B[I]*Z
                ATTEN3=CQ3B[I]*Z
            else: 
                if ICEBOT > Z: #        !Z inside ice
                    ATTEN1=FractionalIceCover[I]*(CQ1BI[I]*(Z+ICETOP)) + (1.-FractionalIceCover[I])*CQ1B[I]*Z
                    ATTEN2=FractionalIceCover[I]*(CQ2BI[I]*(Z+ICETOP)) + (1.-FractionalIceCover[I])*CQ2B[I]*Z
                    ATTEN3=FractionalIceCover[I]*(CQ3BI[I]*(Z+ICETOP)) + (1.-FractionalIceCover[I])*CQ3B[I]*Z
                else:
                    ATTEN1=FractionalIceCover[I]*(CQ1BI[I]*LakeIceHeight[I]+CQ1B[I]*(Z-ICEBOT)) + (1.-FractionalIceCover[I])*CQ1B[I]*Z
                    ATTEN2=FractionalIceCover[I]*(CQ2BI[I]*LakeIceHeight[I]+CQ2B[I]*(Z-ICEBOT)) + (1.-FractionalIceCover[I])*CQ2B[I]*Z
                    ATTEN3=FractionalIceCover[I]*(CQ3BI[I]*LakeIceHeight[I]+CQ3B[I]*(Z-ICEBOT)) + (1.-FractionalIceCover[I])*CQ3B[I]*Z
                
            #----------- Remove fraction of shortwave flux from water column to be used ---
            #            to heat sediment if sediment heating option activated.
            #            Only for ice-free conditions            
            if LakeIceHeight[I] <= 0.0:#        !NO ICE 
                SolarFluxInLayer[I,J]=(1.0-FractionOfShortwaveAtTopOfSediment)*DiagNetShortwaveAtLakeSurface[I]*(CQ1A[I]*math.exp(-ATTEN1) + CQ2A[I]*math.exp(-ATTEN2) + CQ3A[I]*math.exp(-ATTEN3) )
            else:
                SolarFluxInLayer[I,J]=DiagNetShortwaveAtLakeSurface[I]*(CQ1A[I]*math.exp(-ATTEN1) + CQ2A[I]*math.exp(-ATTEN2) + CQ3A[I]*math.exp(-ATTEN3) )
    # --- * COMPUTE CONDUCTIVE HEAT FLUX (ConductiveHeatFluxInLayer) FOR CURRENT TIMESTEP
    # --- * THERMAL CONDUCTIVITY IS WEIGHTED AVERAGE OF ICE AND WATER
    # --- * if ICE IS PRESENT IN LAYER
    for I in range(IL1-1,IL2):
        ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]
        ICETOP=LakeIceHeight[I]-ICEBOT
        for J in range(0,int(NumOfLakeLayers[I])):
            ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
            ZBOT=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
            if ICEBOT >= ZBOT: 
                ConductiveHeatFluxInLayer[I,J]=(-ThermalConductOfIce/LakeLayerThickness)*(LakeTempProfile[I,J+1]-LakeTempProfile[I,J])
            elif ICEBOT < ZBOT and ICEBOT > ZTOP:
                TC=((ICEBOT-ZTOP)*ThermalConductOfIce + (ZBOT-ICEBOT)*ThermalConductOfH2O)/LakeLayerThickness
                ConductiveHeatFluxInLayer[I,J]=(-TC/LakeLayerThickness)*(LakeTempProfile[I,J+1]-LakeTempProfile[I,J])
            else:
                ConductiveHeatFluxInLayer[I,J]=(-ThermalConductOfH2O/LakeLayerThickness)*(LakeTempProfile[I,J+1]-LakeTempProfile[I,J])            
        NLEV=NumOfLakeLayers[I]
        # -- COMPUTE THERMAL FLUX AT LAKE BOTTOM INTO SEDIMENT
        #    ASSUMES LAKE DOES NOT FREEZE TO BOTTOM
        #    (ADIABATIC BC RECOVERED WHEN ThermalConductivityOfSediment=0)
        TSEDT=( (ThermalConductivityOfSediment*SedimentTemp[I]/ThicknessOfSediment)+(ThermalConductOfH2O*LakeTempProfile[I,NLEV-1]/LakeLayerThickness) )/( (ThermalConductivityOfSediment/ThicknessOfSediment)+(ThermalConductOfH2O/LakeLayerThickness) )
        ConductiveHeatFluxInLayer[I,NLEV-1]= (-2.0*ThermalConductOfH2O/LakeLayerThickness)*(TSEDT-LakeTempProfile[I,NLEV-1])
    # -- COMPUTE TEMPERATURE PROFILE (LakeTempProfile) BEFORE TURBULENT MIXING 
    #    FOR NEXT TIMESTEP.  
    #    HEAT CAPACITY OF LAYER IS WEIGHTED AVERAGE OF WATER AND ICE if
    #    NECESSARY
    for I in range(IL1-1,IL2):
        ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]
        ICETOP=LakeIceHeight[I]-ICEBOT
        ICEBOT0=DensityRatioOfIceToWater*LakeIceHeightAtStartOfTimestep[I]
        ICETOP0=LakeIceHeightAtStartOfTimestep[I]-ICEBOT
        # --COLUMN  --- TEMP BEFORE MELT/FREEZE 
              
        for J in range(0,int(NumOfLakeLayers[I])):
            LakeTempProfileAtStartOfTimestep[J]=LakeTempProfile[I,J]
            ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
            ZBOT=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
            if ICEBOT >= ZBOT: 
                HCAP=HeatCapacityOfIce
            elif ICEBOT <= ZTOP: 
                HCAP=HeatCapacityOfH2O
            else: 
                Z=ICEBOT-ZTOP
                HCAP=(Z*HeatCapacityOfIce + (LakeLayerThickness-Z)*HeatCapacityOfH2O)/LakeLayerThickness
            if J == 0:
                EFLUX=F0[I]-ConductiveHeatFluxInLayer[I,0]+Q0[I]-SolarFluxInLayer[I,0]
            else:
                EFLUX=ConductiveHeatFluxInLayer[I,J-1]-ConductiveHeatFluxInLayer[I,J]+SolarFluxInLayer[I,J-1]-SolarFluxInLayer[I,J]
            LakeTempProfile[I,J]=LakeTempProfile[I,J] + (DeltaTimeStep/(LakeLayerThickness*HCAP))*EFLUX
        # -- UPDATE SEDIMENT TEMPERATURE
        NLEV=NumOfLakeLayers[I]
        FFLXSED=(-2.0*ThermalConductivityOfSediment/ThicknessOfSediment)*(TemperatureConstantAtBottomOfSediment-SedimentTemp[I])
        SedimentTemp[I]=SedimentTemp[I] + (DeltaTimeStep/(ThicknessOfSediment*VolumetricHeatOfSediment))*(ConductiveHeatFluxInLayer[I,NLEV-1]-FFLXSED+FractionOfShortwaveAtTopOfSediment*DiagNetShortwaveAtLakeSurface[I])
        #   ICE GROWTH OR DECAY    
        for J in range(0,int(NumOfLakeLayers[I])):
            ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)#In FORTRAN this is -1 but Py is zero indexed.
            ZBOT=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
            if J == 0: 
                EFLUX=F0[I]-ConductiveHeatFluxInLayer[I,0]+Q0[I]-SolarFluxInLayer[I,0]
            else:
                EFLUX=ConductiveHeatFluxInLayer[I,J-1]-ConductiveHeatFluxInLayer[I,J]+SolarFluxInLayer[I,J-1]-SolarFluxInLayer[I,J]
            # -- FREEZING --            
            if EFLUX < 0.0 and ICEBOT0 < ZBOT and ICEBOT0 >= ZTOP:
                # -- Net energy flux used to lower T to FreezingPointOfH2O
                if (LakeTempProfileAtStartOfTimestep[J] > FreezingPointOfH2O):
                    ECOOL=(ZBOT-ICEBOT0)*HeatCapacityOfH2O*(LakeTempProfileAtStartOfTimestep[J]-FreezingPointOfH2O)/DeltaTimeStep
                else:
                    ECOOL=0.0
                
                # -- Remaining energy flux (if any) used to freeze ice
                EAVAIL=EFLUX+ECOOL
                if (EAVAIL < 0.0):
                    NEWICE=-(DeltaTimeStep/(DensityOfIce*LatentHeatOfFreezingH2O))*EAVAIL
                    LakeIceHeight[I]=LakeIceHeight[I]+NEWICE
                    ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]
                    LakeTempProfile[I,J]=FreezingPointOfH2O
                    # -- LIMIT ICE GROWTH TO THE CURRENT LAYER
                    if (ICEBOT > ZBOT): 
                        EHEAT=(DensityOfIce*LatentHeatOfFreezingH2O*(ICEBOT-ZBOT))/DeltaTimeStep
                        LakeTempProfile[I,J]=LakeTempProfile[I,J] - (EHEAT*DeltaTimeStep)/(LakeLayerThickness*HeatCapacityOfIce)
                        LakeIceHeight[I]=ZBOT/DensityRatioOfIceToWater                    
            
            # -- MELTING
            if EFLUX > 0.0 and ICEBOT0 > ZTOP:  
            # -- Net energy flux used to raise T to FreezingPointOfH2O
                ZICE=MIN(LakeLayerThickness, ICEBOT0-ZTOP)
                EHEAT=ZICE*HeatCapacityOfIce*(FreezingPointOfH2O-LakeTempProfileAtStartOfTimestep[J])/DeltaTimeStep
                # -- Remaining energy flux (if any) used to melt ice
                EAVAIL=EFLUX-EHEAT
                if EAVAIL > 0.0:
                    NEWICE=-(DeltaTimeStep/(DensityOfIce*LatentHeatOfFreezingH2O))*EAVAIL
                    LakeIceHeight[I]=LakeIceHeight[I]+NEWICE
                    ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]
                    LakeTempProfile[I,J]=FreezingPointOfH2O
                    # -- LIMIT ICE MELT TO THE CURRENT LAYER
                    if ICEBOT < ZTOP:
                        EHEAT=DensityOfIce*LatentHeatOfFreezingH2O*(ZTOP-ICEBOT)/DeltaTimeStep
                        LakeTempProfile[I,J]=FreezingPointOfH2O + (EHEAT*DeltaTimeStep)/(LakeLayerThickness*HeatCapacityOfH2O)
                        LakeIceHeight[I]=ZTOP/DensityRatioOfIceToWater
    # -- COMPUTE MIXED LAYER DEPTH (DepthToThermocline), TotalKineticEnergy, SHEAR FLOW (MeanMixLyrMomentum)
    #    FOR NEXT TIMESTEP  
        
    [DepthToThermocline,TotalKineticEnergy,MeanMixLyrMomentum,FQU,BFLX,DISS,FSHEAR,FENTRA,TRAN] = \
        MIXLYR(ThermoclineTempDelta,NumOfLakeLayers,FrictionVelocity,IL1,IL2,DepthToThermocline,TotalKineticEnergy,MeanMixLyrMomentum,\
        ExpansivityOfWater,DiagNetShortwaveAtLakeSurface,LakeDepth,LakeLength,ReducedGravity,CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,MixLyrDensity,\
        DiagNetLongwaveAtLakeSurface,DiagSensHeatOnLake,DiagLatentHeatOnLake,LakeIceHeight,TotalKineticEnergyEff_CN,TotalKineticEnergyEff_CF,\
        TotalKineticEnergyEff_CE,TotalKineticEnergyEff_CS,TotalKineticEnergyEff_CL,GravConstant,ThicknessOfSurfaceSkin,SpecificHeatOfH2O,\
        DeltaTimeStep,MaxDeltaThermoclineDepth,LakeLayerThickness,MinimumMixingDepth,MinimumKineticEnergy,MaxDeltaMixingLayer) 
           
    # -- MIX TEMP OVER MIXED LAYER (now mass weighted)    
    for I in range(IL1-1,IL2):
        ICEBOT=DensityRatioOfIceToWater*LakeIceHeight[I]    
        Z2=LakeLayerThickness+ThicknessOfSurfaceSkin
        ZBOT=ThicknessOfSurfaceSkin+LakeLayerThickness*np.float32(NumOfLakeLayers[I])
        if DepthToThermocline[I] < Z2: #!fully stratified: no mixing
            JMIX=0
            ZMIX=np.round(LakeLayerThickness+ThicknessOfSurfaceSkin)
        elif DepthToThermocline[I] >= ZBOT: #	!fully mixed
            JMIX=NumOfLakeLayers[I]
            ZMIX=ZBOT
        else:            
            for J in range(1,int(NumOfLakeLayers[I])):
                Z=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J+1)
                if DepthToThermocline[I] < Z:
                    break
            JMIX=J-1
            ZMIX=np.round(Z-LakeLayerThickness,2)#floating point precision actually messed witht this one.
        # -- MIXING UNDER ICE ALLOWED (ie temp mixed between bottom of mixed layer and bottom of ice cover)               
        if LakeIceHeight[I] <= 0.0:
            TC1=LakeSkinTemp[I]-FreezingPointOfH2O
            [XXX,RHO1a]=EQNST(TC1,0.05)                        
            TBAR=ThicknessOfSurfaceSkin*RHO1a*LakeSkinTemp[I]            
            MASS=ThicknessOfSurfaceSkin*RHO1a
            #TBAR=MASS*LakeSkinTemp[I]
            RHO1_1=RHO1a
        else:
            TBAR=np.float32(0.0)
            MASS=np.float32(0.0)
        for J in range(0,JMIX+1):
            ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
            ZBOT=ZTOP+LakeLayerThickness
            if ICEBOT <= ZTOP: 
                TC1=LakeTempProfile[I,J]-FreezingPointOfH2O
                [XXX,RHO1b]=EQNST(TC1,ZBOT)
                TBAR=TBAR + RHO1b*LakeTempProfile[I,J]*LakeLayerThickness
                MASS=MASS + RHO1b*LakeLayerThickness
                RHO1_2=RHO1b                
        if ((JMIX >= 1) and ((ZMIX-ICEBOT)>=LakeLayerThickness)):
            TMIX=TBAR/MASS
        else:
            TMIX=LakeTempProfile[I,0]
        # MIX TEMPERATURES: include skin layer if no ice
        # -- Mix skin with first layer temperature                
        if ( (LakeIceHeight[I]<=0.0) and (JMIX==0) ):
            TC1=LakeSkinTemp[I]-FreezingPointOfH2O
            TC2=LakeTempProfile[I,0]-FreezingPointOfH2O
            [XXX,RHO1]=EQNST(TC1,0.05)
            [XXX,RHO2]=EQNST(TC2,0.5)
            LakeSkinTemp[I] = (RHO1*ThicknessOfSurfaceSkin*LakeSkinTemp[I] + RHO2*LakeLayerThickness*LakeTempProfile[I,0])/(RHO1*ThicknessOfSurfaceSkin + RHO2*LakeLayerThickness)
            LakeTempProfile[I,0]=LakeSkinTemp[I]
            RHO1_3=RHO1
            RHO2_1=RHO2
        elif (LakeIceHeight[I]<=0.0) and (JMIX>=1):
            LakeSkinTemp[I]=TMIX
            for J in range(0,JMIX+1):
                LakeTempProfile[I,J]=TMIX
        # -- MIXING UNDER ICE
        elif (JMIX >= 1) and ((ZMIX-ICEBOT)>=LakeLayerThickness):            
            for J in range(0,JMIX+1):
                ZTOP=ThicknessOfSurfaceSkin + LakeLayerThickness*np.float32(J)
                if ICEBOT <= ZTOP:
                    LakeTempProfile[I,J]=TMIX                
        # -- COMPUTE TEMPERATURE, DENSITY, EXPANSIVITY OF MIXED LAYER WATER
        #    EQN OF STATE FROM FARMER AND CARMACK, 1982
        TemperatureOfMixingLayerInCelcius[I]=TMIX-FreezingPointOfH2O
        JMIXP1=MIN(int(NumOfLakeLayers[I]),JMIX+1)
        [ExpansivityOfWater[I],MixLyrDensity[I]] = EQNST(TemperatureOfMixingLayerInCelcius[I],DepthToThermocline[I])
        [XXX,RHO1]            = EQNST(LakeTempProfile[I,JMIXP1]-FreezingPointOfH2O,DepthToThermocline[I])
        ReducedGravity[I]             = GravConstant*abs(RHO1-MixLyrDensity[I])/DensityOfH2O
        RHO1_4=RHO1
        if (ReducedGravity[I]<0.):
            XIT('CLASSL',-2)
        # ---* MIX TEMPERATURE IN THERMOCLINE LAYER
        #      DISABLED FOR NOW: JAN 2013 (M.MACKAY)
        # -- Compute thermocline layer thickness (DELTHRM)
        # -- Limit allowed values
        if LakeIceHeight[I] <= 0.0: 
            if ReducedGravity[I] > 0.0:
                DELTHRM=np.float32(0.3)*MeanMixLyrMomentum[I]*MeanMixLyrMomentum[I]/ReducedGravity[I]
            else:
                DELTHRM=np.float32(0.0)            
            if DELTHRM>MaxThermoclineThickness:
                DELTHRM=MaxThermoclineThickness
            if DELTHRM<MinThermoclineThickness:
                DELTHRM=MinThermoclineThickness
            HTOP=DepthToThermocline[I]-np.float32(DELTHRM/2.)
            HBOT=DepthToThermocline[I]+np.float32(DELTHRM/2.)
            for J in range(0,int(NumOfLakeLayers[I])):
                Z=LakeLayerThickness*np.float32(J+1)
                if HTOP < Z:
                    break
            JTOP=J-1
            if JTOP<0:
                JTOP=0            
            ZTOP=LakeLayerThickness*np.float32(JTOP+1)
            TTOP=LakeTempProfile[I,JTOP]

            for J in range(0,int(NumOfLakeLayers[I])):
                Z=LakeLayerThickness*np.float32(J+1)
                if HBOT < Z:
                    break
            JBOT=J-1
            if JBOT<0:
                JBOT=0            
            ZBOT=LakeLayerThickness*np.float32(JBOT+1)
            TBOT=LakeTempProfile[I,JBOT]
            if JBOT<JTOP:
                XIT('CLASSL',-3)            
            if JMIX < NumOfLakeLayers[I]:
                ThermoclineTempDelta[I]=TMIX-LakeTempProfile[I,JBOT]# 	!see Spigel et al 1986
            else:
                ThermoclineTempDelta[I]=0.
        else:
            ThermoclineTempDelta[I]=0.
            DELTHRM=0.
    # -- Compute Wedderburn number diagnostic
    for I in range(IL1-1,IL2):
        if LakeIceHeight[I] <= 0.0 and FrictionVelocity[I] > 0.0: 
            WEDB=ReducedGravity[I]*DepthToThermocline[I]*DepthToThermocline[I]/(np.float32(LakeLength[I])*FrictionVelocity[I]*FrictionVelocity[I])
        else:
            WEDB=np.float32(-999) #	!because FrictionVelocity=0
        # -- Compute heat of melting/freezing for energy balance diagnostic        
        DiagPhaseChangeWaterAtLakeSurface[I] = (LakeIceHeight[I]-LakeIceHeightAtStartOfTimestep[I])*DensityOfIce*LatentHeatOfFreezingH2O/DeltaTimeStep
        OSW[I] = DownwellingShortwaveTotal [I] - FractionOfSurfaceWithSnow[I]*(QSWNS[I]-QTRANSL[I]) - DiagNetShortwaveAtLakeSurface[I]
        # -- RESET SnowIceHeight if ICE HAS COMPLETELY MELTED; ENSURE CONSISTENCY BETWEEN 
        #  FractionalIceCover AND LakeIceHeight
        if LakeIceHeight[I] <= 0.0:
            SnowIceHeight[I]=0.0
        ICELIM=np.float32(LakeLength[I]*1.0E-5)
        FractionalIceCover[I]=MIN(ICELIM,LakeIceHeight[I])/ICELIM
        # -- WRITE OUTPUT FOR THIS TIMESTEP        
        '''
        fID_71.write(f"{CurrentYear:5d} {CurrentDay:4d} {CurrentHour:3d} {CurrentMin:3d} {G0[I,0]:8.2e}"\
            +f" {F0[I,0]:8.2f} {DiagNetShortwaveAtLakeSurface[I,0]:8.2f} {Q0[I,0]:8.2f} {DiagNetShortwaveAtLakeSurface[I,0]:8.2f}"\
            +f" {HFSL[I,0]:8.2f} {HEVL[I,0]:8.2f} {T0[I,0]:8.2f} {LKICEH[I,0]:7.3f}"\
            +f" {SNO[I,0]:7.3f} {RHOSNI[I,0]:7.3f} {RHOW:10.3e} {OSW[I,0]:10.3e}"\
            +f" {SNO[I,0]:10.3e} {SNICEH[I,0]:10.3e}\n")
        
        lakeTempStr=np.array2string(TLAK[I,0:NLAK[I,0]]-TFREZ, formatter={'float_kind':lambda x: "%7.2f" % x},threshold=np.inf, max_line_width=np.inf)
        lakeTempStr=lakeTempStr[1:-1] #I'm not really familar with python formatting so this is what I have to work with.
        fID_72.write(f"{CurrentYear:5d} {CurrentYear:4d} {CurrentHour:3d} {CurrentMin:3d} {lakeTempStr} \n")
        fID_73.write(f"{CurrentYear:5d} {CurrentYear:4d} {CurrentHour:3d} {CurrentMin:3d} {FQU[I,0]:10.2e}"\
            +f"{BFLX[I,0]:10.2e} {DISS[I,0]:10.2e} {TKE[I,0]:10.2e} {DELU[I,0]:5.2f}"\
            +f"{HDPTH[I,0]:6.2f} {JMIX:3d} {ZMIX:7.2f} {TMIX:7.2f}"\
            +f"{DTEMP[I,0]:7.2f} {FSHEAR[I,0]:10.2e} {FENTRA[I,0]:10.2e} \n")
        '''
        fID_71.write(str(AirTempAtRef[I]) + ' ')
        fID_71.write(str(LakeSkinTemp[I]) + ' ')
        for J in range(0,NumOfLakeLayers[I]):
            T_Lake_String=str(LakeTempProfile[I,J])
            fID_71.write(T_Lake_String + ' ')
        fID_71.write(str(SedimentTemp[I]) + ' ')
        fID_71.write('\n')
        
        fID_82.write(str(CurrentYear) + ' ')
        fID_82.write(str(CurrentDay) + ' ')
        fID_82.write(str(CurrentHour) + ' ')
        fID_82.write(str(CurrentMin) + ' ')
        fID_82.write(str(LatentHeatFromLakeSubarea[I]) + ' ')
        fID_82.write(str(SensibleHeatFromLakeSubarea[I]) + ' ')
        fID_82.write(str(LakeSkinTemp[I]) + ' ')
        fID_82.write(str(DepthToThermocline[I]) + ' ')
        fID_82.write(str(ThermoclineTempDelta[I]) + ' ')
        fID_82.write('\n')
        '''fID_82.write(f"{CurrentYear:5d} {CurrentDay:4d} {CurrentHour:3d} {CurrentMin:3d}"\
            +f"{LatentHeatFromLakeSubarea[I]:5.2e} {SensibleHeatFromLakeSubarea[I]:5.2e} {LakeSkinTemp[I]:3.2e} "\
            +f"{DepthToThermocline[I]:3.2f} {ThermoclineTempDelta[I]:3.2f} \n")'''
    # -- SCREEN LEVEL DIAGNOSTICS.        
    for I in range(IL1-1,IL2):
        SurfDragCoeffForMomentum[I]=FractionOfSurfaceWithSnow[I]*SurfaceDragCoefficentForMomentum[I]+(1.0-FractionOfSurfaceWithSnow[I])*CDML[I]
        SurfDragCoeffForHeat[I]=FractionOfSurfaceWithSnow[I]*SurfaceDragCoefficentForHeat[I]+(1.0-FractionOfSurfaceWithSnow[I])*CDHL[I]
        LakeSurfaceDragAtNeutral[I]=FractionOfSurfaceWithSnow[I]*DRAGS[I]+(1.0-FractionOfSurfaceWithSnow[I])*DRAGL[I]
        SurfaceRoughnessLengthForMomentum[I]=RefHeightWindSpeedForcing[I]/math.exp(VonKarmanConst/math.sqrt(LakeSurfaceDragAtNeutral[I]))
        SurfaceRoughnessLengthForHeat[I]=SurfaceRoughnessLengthForMomentum[I]/3.0
        MeanSurfaceTemperature[I]=FractionOfSurfaceWithSnow[I]*TemperatureAtSurface[I]+(1.0-FractionOfSurfaceWithSnow[I])*LakeSkinTemp[I]
        DiagSurfSpecificHumidity[I]=FractionOfSurfaceWithSnow[I]*SpacificHumidityAtSurface[I]+(1.0-FractionOfSurfaceWithSnow[I])*QSURFL[I]
        DiagTotalVisAtLakeSurface[I]=FractionOfSurfaceWithSnow[I]*SnowAlbedoVis[I]+(1.0-FractionOfSurfaceWithSnow[I])*LakeAlbedoVis[I]
    if ISLFD==0:
        for I in range(IL1-1,IL2):
            FACTM=np.float32(AnemometerHeight[I])+SurfaceRoughnessLengthForMomentum[I]  
            FACTH=np.float32(VarHeight[I])+SurfaceRoughnessLengthForMomentum[I]  
            RATIOM=math.sqrt(SurfDragCoeffForMomentum[I])*math.log(FACTM/SurfaceRoughnessLengthForMomentum[I])/VonKarmanConst            
            RATIOM=MIN(RATIOM,1.)                                  
            RATIOH=math.sqrt(SurfDragCoeffForMomentum[I])*math.log(FACTH/SurfaceRoughnessLengthForHeat[I])/VonKarmanConst           
            RATIOH=MIN(RATIOH,1.)                                   
            if MeanSurfaceTemperature[I]>AirTempAtRef[I]:
                RATIOH=RATIOH*SurfDragCoeffForHeat[I]/SurfDragCoeffForMomentum[I]
                RATIOH=MIN(RATIOH,(FACTH/RefHeightAirTempAndHumidForcing[I])**(1./3.))
            DiagAirTempAtScreenLevel[I]=MeanSurfaceTemperature[I]-(MIN(RATIOH,1.))*(MeanSurfaceTemperature[I]-AirTempAtRef[I])
            DiagSpecificHumidAtScreenLevel[I]=DiagSurfSpecificHumidity[I]-(MIN(RATIOH,1.))*(DiagSurfSpecificHumidity[I]-SpecificHumidityAtRef[I])
            DiagZonalWindSpdAtAnemometerLevel[I]=RATIOM*ZonalWindSpeed[I]
            DiagMeridionalWindSpdAtAnemometerLevel[I]=RATIOM*MeridionalWindSpeed[I]
                                                                            
        #After fiddling around with this function, Murray MacKay told me it is
        #not needed.  It only changes DiagRelavtiveHumidAtScreenLevel, which is passed back out of CLASSL
        #and assigned as a diagnostic variabile. 
        #[DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2]=SCREENRH(DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2);           
        #[DiagRelavtiveHumidAtScreenLevel]=SCREENRH(DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2);           
        #[DiagRelavtiveHumidAtScreenLevel]=SCREENRH(IL1,IL2);        
    elif ISLFD==1:
        [DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,\
            SurfDragCoeffForMomentum,SurfDragCoeffForHeat,ZonalWindSpeed,MeridionalWindSpeed,AirTempAtRef,SpecificHumidityAtRef,\
            MeanSurfaceTemperature,DiagSurfSpecificHumidity,SurfaceRoughnessLengthForMomentum,SurfaceRoughnessLengthForHeat,\
            FractionOfLakeSurface,RefHeightWindSpeedForcing,AnemometerHeight,VarHeight,NumOfTilesInAllGrids,IL1,IL2,NumOfGridCells]=\
            SLDIAG(DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,\
            SurfDragCoeffForMomentum,SurfDragCoeffForHeat,ZonalWindSpeed,MeridionalWindSpeed,AirTempAtRef,SpecificHumidityAtRef,\
            MeanSurfaceTemperature,DiagSurfSpecificHumidity,SurfaceRoughnessLengthForMomentum,SurfaceRoughnessLengthForHeat,\
            FractionOfLakeSurface,RefHeightWindSpeedForcing,AnemometerHeight,VarHeight,NumOfTilesInAllGrids,IL1,IL2,NumOfGridCells)
        #See note above where this is called about this function not being needed.
        #[DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2]=SCREENRH(DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2);           
        #[DiagRelavtiveHumidAtScreenLevel]=SCREENRH(DiagRelavtiveHumidAtScreenLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,SurfaceAirPressure,FractionOfLakeSurface,NumOfTilesInAllGrids,IL1,IL2);           
        #[DiagRelavtiveHumidAtScreenLevel]=SCREENRH(IL1,IL2);        
    elif ISLFD==2:
        raise NameError('DIASURFZ is not enabled for python version of the CSLM') 
        '''
            [DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel]=\
                DIASURFZ(DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagAirTempAtScreenLevel,DiagSpecificHumidAtScreenLevel,\
                ZonalWindSpeed,MeridionalWindSpeed,MeanSurfaceTemperature,DiagSurfSpecificHumidity,SurfaceRoughnessLengthForMomentum,\
                SurfaceRoughnessLengthForHeat,ILMO,RefHeightWindSpeedForcing,HBL,UE,FTEMP,FVAP,AnemometerHeight,VarHeight,LatitudeInRad,\
                FractionOfLakeSurface,IL1,IL2)
        '''
    for I in range(IL1-1,IL2):
        DiagNetShortwaveAtSnowSurface[I] =FractionOfSurfaceWithSnow[I]*(QSWNS[I]-QTRANSL[I])
        DiagNetLongwaveAtSnowSurface[I] =FractionOfSurfaceWithSnow[I]*(LongwaveDownwelling[I]-UpwellingLongwaveFromSurface[I])
        DiagSensHeatOnSnow[I] =FractionOfSurfaceWithSnow[I]*SensibleHeatFromSurface[I]                     
        DiagLatentHeatOnSnow[I] =FractionOfSurfaceWithSnow[I]*LatentHeatFromSurface[I]
        DiagDelEnergyDueConducOrMassInSnow[I] =DiagDelEnergyDueConducOrMassInSnow[I]-(FractionOfSurfaceWithSnow[I]*GZEROSL[I])
        DiagDelEnergyDueConducOrMassInLake[I] =DiagDelEnergyDueConducOrMassInLake[I]+(FractionOfSurfaceWithSnow[I]*GZEROSL[I])
        LongwaveUpwelling[I]=FractionOfSurfaceWithSnow[I]*UpwellingLongwaveFromSurface[I]+(1.0-FractionOfSurfaceWithSnow[I])*LongwaveDownwelling[I]-DiagNetLongwaveAtLakeSurface[I]
        SensibleHeatFromLakeSubarea[I]=DiagSensHeatOnLake[I]+DiagSensHeatOnSnow[I]
        ProdLakeDragWindSpdDelT[I]=-SensibleHeatFromLakeSubarea[I]/(DensityOfAir[I]*SpecificHeatOfAir)
        LatentHeatFromLakeSubarea[I]=DiagLatentHeatOnLake[I]+DiagLatentHeatOnSnow[I]
        SensibleHeatFromSnowSubarea[I]=LakeSurfaceEvap[I]+SnowSublimation[I]
        ProdLakeDragWindSpdDelSpHumid[I]=-SensibleHeatFromSnowSubarea[I]/DensityOfAir[I]
        DiagPotentialEvap[I]=SensibleHeatFromSnowSubarea[I]
        EvapEfficiencyAtLakeSurf[I]=1.0    
        DiagSurfBlkBodyTemp[I]=FractionOfSurfaceWithSnow[I]*TemperatureAtSurface[I]+(1.0-FractionOfSurfaceWithSnow[I])*LakeSkinTemp[I]  
                   
    return [LakeDepth, LakeLength, ExtinctionCoeff, NumOfLakeLayers,\
            LakeTempProfile, LakeSkinTemp, DepthToThermocline, LakeIceHeight,\
            SnowIceHeight, RunoffFromIce, SnowPackMass, SnowPackDensity,\
            SnowPackTemp, SnowPackAlbedo, SnowPackLiqudWaterContent, SurfDragCoeffForHeat,\
            SurfDragCoeffForMomentum, SensibleHeatFromLakeSubarea,ProdLakeDragWindSpdDelT,LatentHeatFromLakeSubarea,\
            SensibleHeatFromSnowSubarea,ProdLakeDragWindSpdDelSpHumid,DiagPotentialEvap, EvapEfficiencyAtLakeSurf,\
            DiagSurfBlkBodyTemp, DiagSurfSpecificHumidity, LakeSurfaceDragAtNeutral,DiagAirTempAtScreenLevel,\
            DiagZonalWindSpdAtAnemometerLevel,DiagMeridionalWindSpdAtAnemometerLevel,DiagSpecificHumidAtScreenLevel,DiagRelavtiveHumidAtScreenLevel,\
            LongwaveUpwelling,DiagTotalVisAtLakeSurface, DiagNetShortwaveAtLakeSurface, DiagNetLongwaveAtLakeSurface,\
            DiagSensHeatOnLake, DiagLatentHeatOnLake, DiagPhaseChangeWaterAtLakeSurface, DiagDelEnergyDueConducOrMassInLake,\
            DiagNetShortwaveAtSnowSurface, DiagNetLongwaveAtSnowSurface, DiagSensHeatOnSnow, DiagLatentHeatOnSnow,\
            DiagPhaseChangeWaterInSnow, DiagDelEnergyDueConducOrMassInSnow, PrecipOnLake,PrecipOnSnow,\
            LakeSurfaceEvap,SnowSublimation,SnowRunoff,ExpansivityOfWater,\
            ThermoclineTempDelta, TotalKineticEnergy, MeanMixLyrMomentum, ReducedGravity,\
            MixLyrDensity, LongwaveDownwelling, ZonalWindSpeed, MeridionalWindSpeed,\
            AirTempAtRef, SpecificHumidityAtRef, DensityOfAir, PartialPressureOfDryAir,\
            SurfaceAirPressure, CosineOfSolarZenith, RefHeightWindSpeedForcing, RefHeightAirTempAndHumidForcing,\
            AnemometerHeight,VarHeight,WetPrecipTemp,FrozenPrecipTemp,\
            DensityOfFreshSnow, LatitudeInRad, SnowAlbedoVis,SnowAlbedoNIR,\
            ShortwaveDownwellingIRandVIS, NumOfGridCells, MaxNumOfLakes, ISLFD,\
            IZREF, ITG, IALS, NBS, ISNOALB, IGL, SedimentTemp]

#
# -------- GENERAL FUNCTIONS ---------------
#
def XIT(NAME,N):
    #Exit function, kept for compatablity from FORTRAN
    if N<0:        
        raise NameError(f'ERROR: {NAME} N={N:3d}')
    else:
        print(f'WARNING! Have a proper shut down here: {NAME} N={N:3d}')
        
def EQNST(TCEL,H):
    #Equation of state, returns expansion coeffecent and density of water
    Z = np.float32(0.0) 
    S = np.float32(0.0) 
    # Cmdm  S=0.3
    # Cmdm  S=0.014
    # Cmdm  S=0.03
    # Cmdm  Z=MIN(10.0,H)
    GRAV = np.float32(9.80616)
    C0 = np.float32(4.9388E-5)
    ALPHA = np.float32(3.3039E-7)
    BETA = np.float32(8.2545E-6)
    XI = np.float32(8.0477E-1)
    KAPPA = np.float32(0.2151)
    T0 = np.float32(3.98275)
    RHO0=np.float32(999.975 + XI*S)
    P=np.float32(RHO0*GRAV*Z*1.0E-5)
    T=np.float32(TCEL-T0-KAPPA*S)
    RHO= np.float32(RHO0*(1. + P*(C0-ALPHA*T) - BETA*T*T ))
    EXPW = np.float32((2.*RHO0*BETA*T + RHO0*P*ALPHA )/RHO)
    return [EXPW,RHO]
#
# -------- FUNCTIONS IN CLASSL ---------------
#
#------- SNOWALBA -----------------------------------------------------------------------------<
def SNOALBA(ALVSSN,ALIRSN,ALBSNO,ALSNO, TRSNOWG,\
        ZSNOW,FSNOW,ASVDAT,ASIDAT,IL1,IL2,JL,IALS,NBS,ISNOALB):
    #Empty Vars
    ALVSSC=np.zeros((len(range(IL1-1,IL2))),np.float32)
    ALIRSC=np.zeros((len(range(IL1-1,IL2))),np.float32)
    TRSNOWC=np.zeros((len(range(IL1-1,IL2))),np.float32)
    #ALSNO=np.zeros((len(range(IL1-1,IL2)),1),np.float32)
    
    if ISNOALB != 0: #  !4-BAND SW DISABLED 
        print('IN SNOALBA ISNOALB!=0:',ISNOALB)        
    IPTBAD=np.float32(0)
    for I in range(IL1-1,IL2):
        if ALBSNO[I]<0.50 and ALBSNO[I]>0.499:
            ALBSNO[I]=0.50        
        if FSNOW[I]>0.0 and IALS==0:
            if ALBSNO[I]>0.70:
                ALVSSN[I]=0.79*(ALBSNO[I]-0.70)+0.84
                ALIRSN[I]=1.21*(ALBSNO[I]-0.70)+0.56                   
            else:                                                       
                ALVSSN[I]=0.97*(ALBSNO[I]-0.50)+0.62
                ALIRSN[I]=1.03*(ALBSNO[I]-0.50)+0.38           
            if ALVSSN[I]>0.999 or ALVSSN[I]<0.001:
                IPTBAD=I+1
            if ALIRSN[I]>0.999 or ALIRSN[I]<0.001: 
                IPTBAD=I+1        
        elif FSNOW[I]>0.0 and IALS==1:
            ALVSSN[I]=ASVDAT[I]*1
            ALIRSN[I]=ASIDAT[I]*1            
        ALVSSC[I]=ALVSSN[I]*1
        ALIRSC[I]=ALIRSN[I]*1
        TRSNOWC[I]=math.exp(-25.0*ZSNOW[I])      
    if IPTBAD!=0:
        #WRITE(6,6100) IPTBAD,JL,ALVSSN(IPTBAD),ALIRSN(IPTBAD)          
        #6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALVSSN,ALIRSN = ',2F10.5)  
        XIT('SNOALBA IPTBAD!=0',-1) #error('SNOALBA IPTBAD~=0')                        
    if ISNOALB == 0:
        for I in range(IL1-1,IL2):
            ALSNO[I,0] = ALVSSN[I]*1
            ALSNO[I,1] = ALIRSN[I]*1
            ALSNO[I,2] = ALIRSN[I]*1
            ALSNO[I,3] = ALIRSN[I]*1
            TRSNOWG[I,0:NBS-1] = TRSNOWC[I]*1
    return [ALVSSN, ALIRSN, ALSNO, TRSNOWG]

#------- TLSPREP -------------------------------------------------------------------------<
def TLSPREP(GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,\
        ZRSLDM,ZRSLDH,ZRSLFM,ZRSLFH,ZDSLM,ZDSLH,\
        ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,\
        TVIRTA,TPOTA,CRIB,DRAGS,\
        FLS,ZSNOW,TSNOW,RHOSNO,WSNOW,ZREFM,ZREFH,\
        ZDIAGM,ZDIAGH,TA,QA,VA,IZREF,IL1,IL2,IG,\
        CGRAV,CPD,CKARM,HCPICE,RHOICE,HCPW,RHOW,CLHVAP,CLHMLT):
    
    #Empty Vars
    CEVAP=np.zeros((len(range(IL1-1,IL2))),np.float32)
    IEVAP=np.zeros((len(range(IL1-1,IL2))),np.float32)
    ISAND=np.zeros((len(range(IL1-1,IL2)),1),np.float32)

    for I in range(IL1-1,IL2):
        if FLS[I]>0.:
            ZOM[I]=0.001
            ZOMLNS[I]=math.log(ZOM[I])
            ZOH[I]=0.0003
            ZOELNS[I]=math.log(ZOH[I])
            if IZREF==1:
                ZRSLDM[I]=ZREFM[I]*1
                ZRSLDH[I]=ZREFH[I]*1
                ZRSLFM[I]=ZREFM[I]-ZOM[I]
                ZRSLFH[I]=ZREFH[I]-ZOM[I]
                ZDSLM[I]=ZDIAGM[I]-ZOM[I]
                ZDSLH[I]=ZDIAGH[I]-ZOM[I]
                TPOTA[I]=TA[I]+ZRSLFH[I]*CGRAV/CPD
            else:
                ZRSLDM[I]=ZREFM[I]+ZOM[I]
                ZRSLDH[I]=ZREFH[I]+ZOM[I]
                ZRSLFM[I]=ZREFM[I]*1
                ZRSLFH[I]=ZREFH[I]*1
                ZDSLM[I]=ZDIAGM[I]*1
                ZDSLH[I]=ZDIAGH[I]*1
                TPOTA[I]=TA[I]*1
            
            ZOSCLM[I]=ZOM[I]/ZRSLDM[I]
            ZOSCLH[I]=ZOH[I]/ZRSLDH[I]
            TVIRTA[I]=TPOTA[I]*(1.0+0.61*QA[I])
            CRIB[I]=-CGRAV*ZRSLDM[I]/(TVIRTA[I]*VA[I]**2)
            DRAGS[I]=(CKARM/(math.log(ZRSLDM[I])-ZOMLNS[I]))**2
        
    #
    #     * THERMAL PROPERTIES OF SNOW.
    #
    for I in range(IL1-1,IL2):
        if ZSNOW[I]>0:
            HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
            if RHOSNO[I]<156.0:
                TCSNOW[I]=0.234E-3*RHOSNO[I]+0.023
            else:
                TCSNOW[I]=3.233E-6*RHOSNO[I]*RHOSNO[I]-1.01E-3*RHOSNO[I]+0.138
            
                    #mdm          TCSNOW[I]=2.576E-6*RHOSNO[I]*RHOSNO[I]+0.074       !test: Mellor
                    #mdm          TCSNOW[I]=0.021+4.2E-4*RHOSNO[I]+
                    #mdm 1            2.2E-9*RHOSNO[I]*RHOSNO[I]*RHOSNO[I]      !test: Rogers et al 1995
                    #mdm          TCSNOW[I]=2.22362*(RHOSNO[I]/1000.0)**1.885   !test: Yen 1981
                    #mdm          TCSNOW[I]=0.0688*exp(0.0088*
                    #mcm 1                      (TSNOW[I]-273.16)+4.6682*RHOSNO[I]/1000.0) !Pitman&Zuckerman 1967        
    #
    #     * CALCULATE COEFFICIENTS.
    #
    for I in range(IL1-1,IL2):
        if FLS[I]>0:
            GCOEFFS[I]=3.0*TCSNOW[I]/ZSNOW[I]
            GCONSTS[I]=-3.0*TCSNOW[I]*TSNOW[I]/ZSNOW[I]
            CPHCHS[I]=CLHVAP+CLHMLT
            IWATER[I]=np.float32(2)
        else:
            IWATER[I]=np.float32(1)
        
        CEVAP[I]=np.float32(1.0)
        IEVAP[I]=np.float32(1)
        for J in range(0,IG):
            ISAND[I,J]=np.float32(-4)
        
    #Retrun Vars
    return [GCOEFFS,GCONSTS,CPHCHS,TCSNOW,HCPSNO,IWATER,ZRSLDM,ZRSLDH,ZRSLFM,\
        ZRSLFH,ZDSLM,ZDSLH,ZOSCLM,ZOSCLH,ZOMLNS,ZOELNS,ZOM,ZOH,TVIRTA,TPOTA,\
            CRIB,DRAGS,CEVAP,IEVAP,ISAND]

#------- TSOLVE --------------------------------------------------------------------------<
def TSOLVE(ISNOW,FI,QLWOUT,QSENS,QEVAP,EVAP,\
        TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,\
        FTEMP,FVAP,ILMO,UE,H,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,\
        QLWIN,TPOTA,QA,VA,PADRY,RHOAIR,\
        ALVISG,ALNIRG,CRIB,CPHCH,CEVAP,TVIRTA,\
        ZOSCLH,ZOSCLM,\
        GCONST,GCOEFF,TSTART,TRSNOWG,FSSB,ALSNO,\
        THLIQ,THLMIN,DELZW,RHOSNO,ZSNOW,ZPOND,\
        IWATER,IEVAP,ITERCT,ISAND,\
        ISLFD,ITG,ILG,IL1,IL2,NBS,ISNOALB,\
        TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,\
        DCFLXM,CFLUXM,WZERO,TRTOP,ITER,NITER,KF,\
        DELT,TFREZ,RHOW,SBC,SPHAIR,GRAV,VKC):
    #Empty Vars
    QSWNET=np.zeros((len(range(IL1-1,IL2))),np.float32)
    QTRANS=np.zeros((len(range(IL1-1,IL2))),np.float32)
    A=np.zeros((len(range(IL1-1,IL2))),np.float32)    
    B=np.zeros((len(range(IL1-1,IL2))),np.float32)
    #IG=np.zeros((len(range(IL1-1,IL2)),1),np.float32)
    #JL=np.zeros((len(range(IL1-1,IL2)),1),np.float32)
    JEVAP=np.zeros((len(range(IL1-1,IL2))),np.float32)    
    EVPMAX=np.zeros((len(range(IL1-1,IL2))),np.float32)
    EZERO=np.float32(0.0)
    #[JEVAP,EVPMAX]=makeZeros(length(IL1:IL2),1);    

    if ITG<2:
        ITERMX=np.float32(50)
    else:
        ITERMX=np.float32(5)

    for I in range(IL1-1,IL2):
        QSWNET[I]=np.float32(0.0)
        QTRANS[I]=np.float32(0.0)                                                 
                                                                           
    if ISNOW == 0: # Use usual snow-free bare soil formulation
        for I in range(IL1-1,IL2):
            if FI[I]>0.:
                TRTOP[I,0]=0.
                QSWNV=FSSB(I,0)*(1.0-ALVISG[I])
                if (ISNOALB == 0):
                    QSWNI=FSSB(I,1)*(1.0-ALNIRG[I])
                elif ISNOALB == 1:
                    QSWNI=0.0
                    for IB in range(1,NBS):
                        QSWNI=QSWNI+FSSB[I,IB]*(1.0-ALNIRG[I])                                          
                QSWNET[I]=QSWNV+QSWNI
                QTRANS[I]=QSWNET[I]*TRTOP[I,0]
                QSWNET[I]=QSWNET[I]-QTRANS[I]                                                                
    else:                                                              
        if ISNOALB == 0: # Use the existing snow albedo and transmission 
            for I in range(IL1-1,IL2):
                if FI[I]>0.: 
                    TRTOP[I,0]=TRSNOWG[I,0]*1                               
                    QSWNV=FSSB[I,0]*(1.0-ALSNO[I,0])                      
                    QSWNI=FSSB[I,1]*(1.0-ALSNO[I,1])                      
                    QSWNET[I]=QSWNV+QSWNI
                    QTRANS[I]=QSWNET[I]*TRTOP[I,0]
                    QSWNET[I]=QSWNET[I]-QTRANS[I]                                           
        elif ISNOALB == 1: # Use the band-by-band snow albedo and transmission
            for I in range(IL1-1,IL2):
                QTRANS[I] = 0.0
                QSWNET[I] = 0.0                                                         
            for IB in range(0,NBS):
                for I in range(IL1-1,IL2):
                    if FI[I]>0.:
                        TRTOP[I,IB]=TRSNOWG[I,IB]*1
                        QSWNV=FSSB[I,IB]*(1.0-ALSNO[I,IB])
                        QSWNET[I]=QSWNET[I]+FSSB[I,IB]*(1.0-ALSNO[I,IB])
                        QTRANS[I]=QTRANS[I]+QSWNV*TRTOP[I,IB]                           
            for I in range(IL1-1,IL2):
                if FI[I]>0.:
                    QSWNET[I]=QSWNET[I]-QTRANS[I]                                         
    #                                                                       
    for I in range(IL1-1,IL2):
        if FI[I]>0.:
            TZERO[I]=TSTART[I]*1
            TSTEP[I]=1.0
            ITER[I]=1
            NITER[I]=0
            #                                                                       
            QMELT[I]=0.0
            RESID[I]=999999.
            DCFLXM[I]=0.0
            CFLUX[I]=0.0
            if ISNOW==1:
                KF[I]=3                                               
                EVPMAX[I]=RHOSNO[I]*ZSNOW[I]/DELT
            else:
                KF[I]=6
                EVPMAX[I]=RHOW*(ZPOND[I]+(THLIQ[I,0]-THLMIN[I,0])*DELZW[I,0])/DELT            
    #                                                                       
    #     * ITERATION SECTION.                                              
    #     * LOOP IS REPEATED UNTIL SOLUTIONS HAVE BEEN FOUND FOR ALL POINTS 
    #     * ON THE CURRENT LATITUDE CIRCLE(S).                              
    #  
    continueLoop=True
    while continueLoop:    
        NUMIT=0                                                           
        NIT=0                                                             
        for I in range(IL1-1,IL2):
            if FI[I]>0. and ITER[I]==1:
                NIT=NIT+1
                CFLUXM[I]=CFLUX[I]*1
                if TZERO[I]>=TFREZ:
                    A[I]=17.269
                    B[I]=35.86
                else:                                                      
                    A[I]=21.874
                    B[I]=7.66                                                             
                WZERO[I]=0.622*611.0*math.exp(A[I]*(TZERO[I]-TFREZ)/(TZERO[I]-B[I]))/PADRY[I]                           
                Q0SAT[I]=WZERO[I]/(1.0+WZERO[I])
                if(IWATER[I]>0):
                    EVBETA[I]=1.0
                    QZERO[I]=Q0SAT[I]*1
                else:
                    EVBETA[I]=CEVAP[I]*1                                    
                    QZERO[I]=EVBETA[I]*Q0SAT[I]+(1.0-EVBETA[I])*QA[I]
                    if QZERO[I]>QA[I] and IEVAP[I]==0:
                        EVBETA[I]=0.0
                        QZERO[I]=QA[I]*1                                                                        
            TVIRTS[I]=TZERO[I]*(1.0+0.61*QZERO[I])                                                              
        if NIT>0:
            #                                                                       
            #     * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND   
            #     * OTHER RELATED QUANTITIES.                                       
            #                                                                       
            if ISLFD<2:                
                [CDM,CDH,RIB,CFLUX] = \
                    DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,\
                        TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,\
                        GRAV,VKC)
            else:
                raise NameError('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead')
                #CALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
                # 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
                # 2                    TZERO,QZERO,H,ZOM,ZOH,                        
                # 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )                                                                         
            #                                                                       
            #     * REMAINING CALCULATIONS.                                         
            #                                                                       
            for I in range(IL1-1,IL2):
                if FI[I]>0. and ITER[I]==1:
                    QLWOUT[I]=SBC*TZERO[I]*TZERO[I]*TZERO[I]*TZERO[I]
                    if TZERO[I]<TPOTA[I]:
                        QSENS[I]=(RHOAIR[I]*SPHAIR*CFLUX[I]+EZERO)*(TZERO[I]-TPOTA[I])
                    else:
                        QSENS[I]=RHOAIR[I]*SPHAIR*CFLUX[I]*(TZERO[I]-TPOTA[I])                    
                    EVAP[I]=RHOAIR[I]*CFLUX[I]*(QZERO[I]-QA[I])
                    if EVAP[I]>EVPMAX[I]:
                        EVAP[I]=EVPMAX[I]*1                    
                    QEVAP[I]=CPHCH[I]*EVAP[I]
                    GZERO[I]=GCOEFF[I]*TZERO[I]+GCONST[I]
                    RESID[I]=QSWNET[I]+QLWIN[I]-QLWOUT[I]-QSENS[I]-QEVAP[I]-GZERO[I]
                    if abs(RESID[I])<5.0:
                        ITER[I]=0                    
                    if abs(TSTEP[I])< 1.0E-2:
                        ITER[I]=0                    
                    if NITER[I]==ITERMX-1 and ITER[I]==1:
                        ITER[I]=-1                                          
        if ITG<2:
            #                                                                       
            #     * OPTION #1: BISECTION ITERATION METHOD.                          
            #                                                                       
            if NIT>0:
                for I in range(IL1-1,IL2):
                    if FI[I]>0. and ITER[I]==1:
                        if NITER[I]==0:
                            if RESID[I]>0.0:
                                TZERO[I]=TZERO[I]+1.0
                            else:
                                TZERO[I]=TZERO[I]-1.0                                                       
                        else:
                            if (RESID[I]>0. and TSTEP[I]<0.) or (RESID[I]<0. and TSTEP[I]>0.):
                                TSTEP[I]=-TSTEP[I]/2.0                                                      
                            TZERO[I]=TZERO[I]+TSTEP[I]
                    if FI[I]>0. and ITER[I]==1:
                        NITER[I]=NITER[I]+1
                        NUMIT=NUMIT+1                                                           
    #                                                                       
    #     DO 185 I=IL1,IL2                                                  
    #         if(FI[I]>0. && ITER[I]==-1)                      %THEN 
    #             WRITE(6,6250) I,JL,RESID[I],TZERO[I],RIB[I]               
    #6250         FORMAT('0GROUND ITERATION LIMIT',3X,2I3,3(F8.2,E12.4))    
    #         end                                                         
    # 185 CONTINUE                                                          
    #                                                                       
            if NUMIT>0:
                #%GO TO 100       
                continue                                                    
        else: #Hard to tell if this else belongs to line 396 (ITG<2) or line 431 (NUMIT>0) but I think the logic makes sense here                                                             
    #                                                                       
    #     * OPTION #2: NEWTON-RAPHSON ITERATION METHOD.                     
    #                                                                       
            if NIT>0:
                for I in range(IL1-1,IL2):
                    if FI[I]>0. and ITER[I]==1:
                        if NITER[I]>0:
                            DCFLUX=(CFLUX[I]-CFLUXM[I])/SIGN(MAX(.001,abs(TSTEP[I])),TSTEP[I])
                            if(abs(TVIRTA[I]-TVIRTS[I])<0.4):
                                DCFLUX=MAX(DCFLUX,0.8*DCFLXM[I])                            
                            DCFLXM[I]=DCFLUX*1
                        else:
                            DCFLUX=0.                                                                    
                        DRDT0= -4.0*SBC*TZERO[I]**3-RHOAIR[I]*SPHAIR\
                            *(CFLUX[I]+MAX(0.,TZERO[I]-TPOTA[I])*DCFLUX)\
                            -GCOEFF[I]+ CPHCH[I]*RHOAIR[I]*(CFLUX[I]\
                            *WZERO[I]*A[I]*EVBETA[I]*(B[I]-TFREZ)/(\
                            (TZERO[I]-B[I])*(1.0+WZERO[I]))**2\
                            -(QZERO[I]-QA[I])*DCFLUX)            
                        TSTEP[I]=-RESID[I]/DRDT0
                        TSTEP[I]=MAX(-10.,MIN(5.,TSTEP[I]))
                        TZERO[I]=TZERO[I]+TSTEP[I]
                        NITER[I]=NITER[I]+1
                        NUMIT=NUMIT+1                                                                            
            if NUMIT>0:
                #GO TO 100       
                continue            
            #                                                                       
            #     * if CONVERGENCE HAS NOT BEEN REACHED, CALCULATE TEMPERATURE AND  
            #     * FLUXES ASSUMING NEUTRAL STABILITY.                              
            #                                                                       
            for I in range(IL1-1,IL2):
                NUMIT=0
                JEVAP[I]=0
                if FI[I]>0. and ITER[I]==-1:
                    TZEROT=TVIRTA[I]/(1.0+0.61*QZERO[I])
                    if abs(RESID[I])>50.:
                        TZERO[I]=TZEROT*1
                        if TZERO[I]>=TFREZ:
                            A[I]=17.269                                       
                            B[I]=35.86                                        
                        else:
                            A[I]=21.874                                       
                            B[I]=7.66                        
                        WZERO[I]=0.622*611.0*math.exp(A[I]*(TZERO[I]-TFREZ)/\
                            (TZERO[I]-B[I]))/PADRY[I]
                        Q0SAT[I]=WZERO[I]/(1.0+WZERO[I])
                        QZERO[I]=EVBETA[I]*Q0SAT[I]+(1.0-EVBETA[I])*QA[I]
                        QLWOUT[I]=SBC*TZERO[I]*TZERO[I]*TZERO[I]*TZERO[I]
                        GZERO[I]=GCOEFF[I]*TZERO[I]+GCONST[I]
                        RESID[I]=QSWNET[I]+QLWIN[I]-QLWOUT[I]-GZERO[I]
                        if RESID[I]>0.:
                            QEVAP[I]=RESID[I]*1
                        else:
                            QEVAP[I]=RESID[I]*0.5                        
                        if IEVAP[I]==0:
                            QEVAP[I]=0.0                        
                        QSENS[I]=RESID[I]-QEVAP[I]
                        RESID[I]=0.
                        EVAP[I]=QEVAP[I]/CPHCH[I]
                        TVIRTS[I]=TZERO[I]*(1.0+0.61*QZERO[I])
                        JEVAP[I]=1
                        NUMIT=NUMIT+1                                       
            if NUMIT>0:
                if ISLFD<2:
                    #ALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      
                    # 1                   CRIB,TVIRTS,TVIRTA,VA,FI,JEVAP,                
                    # 2                   ILG,IL1,IL2)
                    [CDM,CDH,RIB,CFLUX] = \
                        DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,\
                            TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,\
                            GRAV,VKC)
                else:
                    raise NameError('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead')
                    #ALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
                    # 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
                    # 2                    TZERO,QZERO,H,ZOM,ZOH,                        
                    # 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,JEVAP,JL )    
        continueLoop=False #Break out of the while loop. Could use break.
    #                                                                       
    #     * CHECK FOR BAD ITERATION TEMPERATURES.                           
    #                                                                       
    IBAD=0
    for I in range(IL1-1,IL2):
        if FI[I]>0. and (TZERO[I]<123.16 or TZERO[I]>373.16):
            IBAD=I+1                                                          
    if IBAD!=0:
        # WRITE(6,6275) IBAD,JL,TZERO(IBAD),NITER(IBAD),ISNOW           
        # 6275     FORMAT('0BAD ITERATION TEMPERATURE',3X,2I3,F16.2,2I4)         
        # WRITE(6,6280) QSWNET(IBAD),QLWIN(IBAD),QSENS(IBAD),           
        # 1        QEVAP(IBAD),GZERO(IBAD),CFLUX(IBAD),RIB(IBAD)             
        # 6280     FORMAT(2X,7F12.4)                                             
        #ALL XIT('TSOLVE',-1)                                         
        XIT('TSOLVE',-1)                                                
    #                                                                       
    #     * POST-ITERATION CLEAN-UP.                                        
    #                                                                       
    NIT=0
    for I in range(IL1-1,IL2):
        if FI[I]>0.:
            if ((IWATER[I]==1 and TZERO[I]<TFREZ) or\
                (IWATER[I]==2 and TZERO[I]>TFREZ)) or\
                (ISAND[I,0]==-4 and TZERO[I]>TFREZ):
                TZERO[I]=TFREZ*1
                WZERO[I]=0.622*611.0/PADRY[I]
                QZERO[I]=WZERO[I]/(1.0+WZERO[I])
                TVIRTS[I]=TZERO[I]*(1.0+0.61*QZERO[I])
                ITER[I]=1
                NIT=NIT+1
            else:
                ITER[I]=0                                                   
    if NIT>0:
        #                                                                       
        #       * CALCULATE SURFACE DRAG COEFFICIENTS (STABILITY-DEPENDENT) AND 
        #       * OTHER RELATED QUANTITIES.                                     
        #                                                                       
        if ISLFD<2:
            #ALL DRCOEF (CDM,CDH,RIB,CFLUX,QZERO,QA,ZOSCLM,ZOSCLH,      
            # 1                   CRIB,TVIRTS,TVIRTA,VA,FI,ITER,                 
            # 2                   ILG,IL1,IL2)   
            [CDM,CDH,RIB,CFLUX] = \
                DRCOEF (CDM,CDH,RIB,CFLUX,ZOSCLM,ZOSCLH,CRIB,\
                    TVIRTS,TVIRTA,VA,FI,ITER,IL1,IL2,\
                    GRAV,VKC)
        else:
            raise NameError('FLXSURFZ is not converted from FORTRAN, work in progress.  Set ISLFD switch to 1 to use DRCOEF instead')
        #ALL FLXSURFZ(CDM,CDH,CFLUX,RIB,FTEMP,FVAP,ILMO,            
        # 1                    UE,FCOR,TPOTA,QA,ZRSLFM,ZRSLFH,VA,            
        # 2                    TZERO,QZERO,H,ZOM,ZOH,                        
        # 3                    LZZ0,LZZ0T,FM,FH,ILG,IL1,IL2,FI,ITER,JL )   
    #                                                                       
    #     * REMAINING CALCULATIONS.                                         
    #                                                                       
    for I in range(IL1-1,IL2):
        if FI[I]>0. and ITER[I]==1:
            QLWOUT[I]=SBC*TZERO[I]*TZERO[I]*TZERO[I]*TZERO[I]
            if TZERO[I]<TPOTA[I]:
                QSENS[I]=(RHOAIR[I]*SPHAIR*CFLUX[I]+EZERO)*(TZERO[I]-TPOTA[I])                                         
            else:
                QSENS[I]=RHOAIR[I]*SPHAIR*CFLUX[I]*(TZERO[I]-TPOTA[I])            
            EVAP[I]=RHOAIR[I]*CFLUX[I]*(QZERO[I]-QA[I])
            if EVAP[I]>EVPMAX[I]:
                EVAP[I]=EVPMAX[I]*1            
            QEVAP[I]=CPHCH[I]*EVAP[I]
            GZERO[I]=GCOEFF[I]*TZERO[I]+GCONST[I]
            QMELT[I]=QSWNET[I]+QLWIN[I]-QLWOUT[I]-QSENS[I]-QEVAP[I]-GZERO[I]
            RESID[I]=0.0                                              
            if QMELT[I]<0.0:
                QMELT[I]=QMELT[I]+QEVAP[I]
                QEVAP[I]=0.0
                EVAP[I] =0.0                                                    
        if FI[I]>0.:
            if abs(EVAP[I])<1.0E-8:
                RESID[I]=RESID[I]+QEVAP[I]
                EVAP[I]=0.0
                QEVAP[I]=0.0                                                    
            if (ISNOW==1 and QMELT[I]<0.0) or (ISNOW==0 and QMELT[I]>0.0):
                GZERO[I]=GZERO[I]+QMELT[I]                            
                QMELT[I]=0.0                        
        #              QSENS[I]=QSENS[I]+0.5*RESID[I]                           
        #              GZERO[I]=GZERO[I]+0.5*RESID[I]                           
            QSENS[I]=QSENS[I]+RESID[I]
            QSWNET[I]=QSWNET[I]+QTRANS[I]
            EVAP[I]=EVAP[I]/RHOW
            ITERCT[I,KF[I]-1,NITER[I]]=ITERCT[I,KF[I]-1,NITER[I]]+1                 
    #Retrun Vars
    return [QSWNET,QLWOUT,QTRANS,QSENS,QEVAP,EVAP,\
        TZERO,QZERO,GZERO,QMELT,CDH,CDM,RIB,CFLUX,\
        FTEMP,FVAP,ILMO,UE,H,ZRSLFH,ZRSLFM,ZOH,ZOM,FCOR,ILG,\
        ITERCT,TSTEP,TVIRTS,EVBETA,Q0SAT,RESID,\
        DCFLXM,CFLUXM,WZERO,TRTOP,ITER,NITER,JEVAP,KF]

#------- DRCOEF ----------------------------------------------------------------------
def  DRCOEF(CDM,CDH,RIB,CFLUX,ZOMIN,ZOHIN,\
          CRIB,TVIRTG,TVIRTA,VA,FI,ITER,\
          IL1,IL2,GRAV,VKC):
    #Empty Vars
    ZOM=np.zeros((len(range(IL1-1,IL2))),np.float32) 
    ZOH=np.zeros((len(range(IL1-1,IL2))),np.float32)     
    AA=np.float32(9.5285714)
    AA1=np.float32(14.285714)
    BETA=np.float32(1.2)
    PR = np.float32(1.)
    for I in range(IL1-1,IL2):
        if FI[I]>0. and ITER[I]==1:
            RIB[I]=CRIB[I]*(TVIRTG[I]-TVIRTA[I])
            if RIB[I]>=0.0:
                ZLEV=-CRIB[I]*TVIRTA[I]*(VA[I]**2)/GRAV  
                ZS=MAX(10.,5.*MAX(ZOMIN[I]*ZLEV, ZOHIN[I]*ZLEV))
                ZS=ZLEV*(1.+RIB[I])/(1.+(ZLEV/ZS)*RIB[I])
                ZOM[I]=ZOMIN[I]*ZLEV/ZS
                ZOH[I]=ZOHIN[I]*ZLEV/ZS
                RIB[I]=RIB[I]*ZS/ZLEV
            else:
                ZOM[I]=ZOMIN[I]*1
                ZOH[I]=ZOHIN[I]*1            
            ZOLN=math.log(ZOH[I])
            ZMLN=math.log(ZOM[I])
            if RIB[I]<0.0:
                CPR=MAX(ZOLN/ZMLN,0.74)
                CPR=MIN(CPR,1.0)
                ZI=1000.0
                OLSF=BETA**3*ZI*VKC**2/ZMLN**3
                OLFACT=1.7*(math.log(1.+ZOM[I]/ZOH[I]))**0.5+0.9
                OLSF=OLSF*OLFACT
                ZL = -CRIB[I]*TVIRTA[I]*(VA[I]**2)/GRAV
                ZMOL=ZOM[I]*ZL/OLSF
                ZHOL=ZOH[I]*ZL/OLSF
                XM=(1.00-15.0*ZMOL)**(0.250)
                XH=(1.00-9.0*ZHOL)**0.25
                BH1=-math.log(-2.41*ZMOL)+math.log(((1.+XM)/2.)**2*(1.+XM**2)/2.)
                BH1=BH1-2.*math.atan(XM)+math.atan(1.)*2.
                BH1=BH1**1.5
                BH2=-math.log(-0.25*ZHOL)+2.*math.log(((1.00+XH**2)/2.00))
                BH=VKC**3.*BETA**1.5/(BH1*(BH2)**1.5)
                WB=math.sqrt(GRAV*(TVIRTG[I]-TVIRTA[I])*ZI/TVIRTG[I])                
                WSTAR=BH**(0.333333)*WB
                RIB0=RIB[I]*1

                WSPEED=math.sqrt(VA[I]**2+(BETA*WSTAR)**2)
                RIB[I]=RIB0*VA[I]**2/WSPEED**2
                AU1=1.+5.0*(ZOLN-ZMLN)*RIB[I]*(ZOH[I]/ZOM[I])**0.25
                OLS=-RIB[I]*ZMLN**2/(CPR*ZOLN)*(1.0+AU1/(1.0-RIB[I]/(ZOM[I]*ZOH[I])**0.25))
                PSIM1=math.log(((1.00+(1.00-15.0*OLS)**0.250)/2.00)**2*\
                    (1.0+(1.00-15.0*OLS)**0.5)/2.0)-2.0*math.atan((1.00-15.0*OLS)**0.250)+math.atan(1.00)*2.00
                PSIM0=math.log(((1.00+(1.00-15.0*OLS*ZOM[I])**0.250)/2.00)**2*\
                    (1.0+(1.00-15.0*OLS*ZOM[I])**0.5)/2.0)-2.0*math.atan((1.00-15.0*OLS*ZOM[I])**0.250)+math.atan(1.00)*2.0
                PSIH1=math.log(((1.00+(1.00-9.0*OLS)**0.50)/2.00)**2)
                PSIH0=math.log(((1.00+(1.00-9.0*OLS*ZOH[I])**0.50)/2.00)**2)

                USTAR=VKC/(-ZMLN-PSIM1+PSIM0)
                TSTAR=VKC/(-ZOLN-PSIH1+PSIH0)
                CDH[I]=USTAR*TSTAR/PR
                WTS=CDH[I]*WSPEED*(TVIRTG[I]-TVIRTA[I])
                WSTAR=(GRAV*ZI/TVIRTG[I]*WTS)**(0.333333)

                WSPEED=math.sqrt(VA[I]**2+(BETA*WSTAR)**2)
                RIB[I]=RIB0*VA[I]**2/WSPEED**2
                AU1=1.+5.0*(ZOLN-ZMLN)*RIB[I]*(ZOH[I]/ZOM[I])**0.25
                OLS=-RIB[I]*ZMLN**2/(CPR*ZOLN)*(1.0+AU1/(1.0-RIB[I]/(ZOM[I]*ZOH[I])**0.25))
                PSIM1=math.log(((1.00+(1.00-15.0*OLS)**0.250)/2.00)**2*\
                    (1.0+(1.00-15.0*OLS)**0.5)/2.0)-2.0*math.atan((1.00-15.0*OLS)**0.250)+math.atan(1.00)*2.00
                PSIM0=math.log(((1.00+(1.00-15.0*OLS*ZOM[I])**0.250)/2.00)**2*\
                    (1.0+(1.00-15.0*OLS*ZOM[I])**0.5)/2.0)-2.0*math.atan((1.00-15.0*OLS*ZOM[I])**0.250)+math.atan(1.00)*2.0
                PSIH1=math.log(((1.00+(1.00-9.0*OLS)**0.50)/2.00)**2)
                PSIH0=math.log(((1.00+(1.00-9.0*OLS*ZOH[I])**0.50)/2.00)**2)

            else:
                WSPEED=VA[I]*1
                AS1=10.0*ZMLN*(ZOM[I]-1.0)
                AS2=5.00/(2.0-8.53*RIB[I]*math.exp(-3.35*RIB[I])+0.05*RIB[I]**2)
                #<<<
                AS2=AS2*PR*math.sqrt(-ZMLN)/2.
                AS3=27./(8.*PR*PR)
                #>>>
                OLS=RIB[I]*(ZMLN**2+AS3*AS1*(RIB[I]**2+AS2*RIB[I]))/(AS1*RIB[I]-PR*ZOLN)
                PSIM1=-0.667*(OLS-AA1)*math.exp(-0.35*OLS)-AA-OLS
                PSIM0=-0.667*(OLS*ZOM[I]-AA1)*math.exp(-0.35*OLS*ZOM[I])-AA-OLS*ZOM[I]
                PSIH1=-(1.0+2.0*OLS/3.0)**1.5-0.667*(OLS-AA1)*math.exp(-0.35*OLS)-AA+1.0
                PSIH0=-(1.0+2.0*OLS*ZOH[I]/3.0)**1.5-0.667*(OLS*ZOH[I]-AA1)*math.exp(-0.35*OLS*ZOH[I])-AA+1.0

            USTAR=VKC/(-ZMLN-PSIM1+PSIM0)
            TSTAR=VKC/(-ZOLN-PSIH1+PSIH0)

            CDM[I]=USTAR**2.0
            CDH[I]=USTAR*TSTAR/PR
            #
            #         * CALCULATE CD*MOD(V) UNDER FREE-CONVECTIVE LIMIT.
            #
            if TVIRTG[I]>TVIRTA[I]:
                CLIMIT=1.9E-3*(TVIRTG[I]-TVIRTA[I])**0.333333
            else:
                CLIMIT=0.
               
            CFLUX[I]=MAX(CDH[I]*WSPEED,CLIMIT)        

    #Retrun Vars
    return [CDM,CDH,RIB,CFLUX]

#------- FLXSURFZ ----------------------------------------------------------------------
def FLXSURFZ():
    raise NameError('ERROR FLXSURFZ has not been translated from Python')

#------- TLSPREP ----------------------------------------------------------------------
def TLSPOST(GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZERO,TSNBOT,\
        HTCS,HMFN,EVAPS,\
        T0,ZSNOW,TCSNOW,HCPSNO,\
        RPCP,TRPCP,SPCP,TSPCP,TZEROS,RHOSNI,\
        FLS,DELSKIN,IL1,IL2,\
        RHOW,DELT,TFREZ,CLHMLT,HCPICE,RHOICE,HCPW):
    #Empty Vars
    QFN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    SPCN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    RPCN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    TRPCN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    TSPCN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    RADD=np.zeros((len(range(IL1-1,IL2))),np.float32)
    SADD=np.zeros((len(range(IL1-1,IL2))),np.float32)
    for I in range(IL1-1,IL2):
        if FLS[I]>0:
            TSNBOT[I]=(ZSNOW[I]*TSNOW[I]+DELSKIN*T0[I])/(ZSNOW[I]+DELSKIN)
            GZERO[I]=-2.0*TCSNOW[I]*(TSNBOT[I]-TSNOW[I])/ZSNOW[I]
    #     1                 +TCICE*(T0[I]-TSNBOT[I])/DELSKIN)
            if QMELTS[I]<0.:
                GSNOW[I]=GSNOW[I]+QMELTS[I]
                QMELTS[I]=0
            
            TSNOW[I]=TSNOW[I]+(GSNOW[I]-GZERO[I])*DELT/(HCPSNO[I]*ZSNOW[I])-TFREZ
            if TSNOW[I]>0:
                QMELTS[I]=QMELTS[I]+TSNOW[I]*HCPSNO[I]*ZSNOW[I]/DELT
                GSNOW[I]=GSNOW[I]-TSNOW[I]*HCPSNO[I]*ZSNOW[I]/DELT
                TSNOW[I]=0
            
    #         C              GZERO[I]=GZERO[I]+QMELTS[I]
    #         C              QMELTS[I]=0.0
      
    for I in range(IL1-1,IL2):
        if FLS[I]>0 and TSNOW[I]<0 and WSNOW[I]>0:
            HTCS[I]=HTCS[I]-FLS[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
            HADD=-TSNOW[I]*HCPSNO[I]*ZSNOW[I]
            HCONV=CLHMLT*WSNOW[I]
            if HADD<=HCONV:
                WFREZ=HADD/CLHMLT
                HADD=np.float32(0.0)
                WSNOW[I]=MAX(0.0,WSNOW[I]-WFREZ)
                TSNOW[I]=0.0
                RHOSNO[I]=RHOSNO[I]+WFREZ/ZSNOW[I]
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
            else:
                HADD=HADD-HCONV
                WFREZ=WSNOW[I]*1
                WSNOW[I]=0.0
                RHOSNO[I]=RHOSNO[I]+WFREZ/ZSNOW[I]
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE
                TSNOW[I]=-HADD/(HCPSNO[I]*ZSNOW[I])
            
            HMFN[I]=HMFN[I]-FLS[I]*CLHMLT*WFREZ/DELT
            HTCS[I]=HTCS[I]-FLS[I]*CLHMLT*WFREZ/DELT
            HTCS[I]=HTCS[I]+FLS[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
        
    #
    for I in range(IL1-1,IL2):
        QFN[I]=FLS[I]*EVAPS[I]*RHOW
        if SPCP[I]>0 or EVAPS[I]<0:
            SADD[I]=SPCP[I]-EVAPS[I]*RHOW/RHOSNI[I]
            if abs(SADD[I]) < 1.0E-12: 
                SADD[I]=np.float32(0.0)
            
            if SADD[I] > 0.0:
                SPCN[I]=SADD[I]*1
                if SPCP[I] > 0.0:
                    TSPCN[I]=TSPCP[I]*1
                else:
                    TSPCN[I]=MIN((TZEROS[I]-TFREZ),0.0)
                
                EVAPS[I]=np.float32(0.0)
            else:                                           
                EVAPS[I]=-SADD[I]*RHOSNI[I]/RHOW
                SPCN [I]=np.float32(0.0)
                TSPCN[I]=np.float32(0.0)
            
        else:
            SPCN [I]=np.float32(0.0)
            TSPCN[I]=np.float32(0.0)
        
        if RPCP[I] > 0:
            RADD[I]=RPCP[I]-EVAPS[I]
            if abs(RADD[I]) < 1.0E-12: 
                RADD[I]=np.float32(0.0)
            
            if RADD[I] > 0.:
                RPCN [I]=RADD[I]*1
                TRPCN[I]=TRPCP[I]*1
                EVAPS[I]=np.float32(0.0)
            else:
                EVAPS[I]=-RADD[I]*1
                RPCN [I]=np.float32(0.0)
                TRPCN[I]=np.float32(0.0)
                          
        else:
            RPCN [I]=np.float32(0.0)
            TRPCN[I]=np.float32(0.0)
        
    #Retrun Vars
    return [GSNOW,TSNOW,WSNOW,RHOSNO,QMELTS,GZERO,TSNBOT,\
        HTCS,HMFN,QFN,EVAPS,RPCN,TRPCN,SPCN,TSPCN,\
        HCPSNO]   

#------- SNOVAP ----------------------------------------------------------------------
def SNOVAP(RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS, \
    WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW, \
    FI,R,S,RHOSNI,WSNOW,IL1,IL2, \
    DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW):
    CLHMLT = 0.334E6 #Latent heat of freezing of water ($J kg^{-1}$)
    CLHVAP = 2.501E6 #Latent heat of vaporization of water ($J kg^{-1}$)
    #Empty Vars (none)
    for I in range(IL1-1,IL2):
        if FI[I]>0. and (S[I]<1.0E-11 or R[I]<1.0E-11) and ZSNOW[I]>0.:
            HTCS[I]=HTCS[I]-FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
            if EVAP[I]<0.:
                ZADD=-EVAP[I]*DELT*RHOW/RHOSNI[I]
                RHOSNO[I]=(ZSNOW[I]*RHOSNO[I]+ZADD*RHOSNI[I])/(ZSNOW[I]+ZADD);                         
                ZSNOW [I]=ZSNOW[I]+ZADD
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                EVAP  [I]=0.0
            else:
                ZLOST=EVAP[I]*DELT*RHOW/RHOSNO[I]
                if ZLOST<=ZSNOW[I]:
                    ZSNOW[I]=ZSNOW[I]-ZLOST
                    EVAP [I]=0.0
                    HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                else:
                    ZREM=(ZLOST-ZSNOW[I])*RHOSNO[I]/RHOW
                    ZSNOW[I]=0.0
                    HCPSNO[I]=0.0
                    EVAP[I]=ZREM*(CLHMLT+CLHVAP)/(CLHVAP*DELT)
                    WLOST[I]=WLOST[I]-ZREM*RHOW*CLHMLT/CLHVAP
                    if RUNOFF[I]>0. or WSNOW[I]>0.:
                        TRUNOF[I]=(TRUNOF[I]*RUNOFF[I]+(TSNOW[I]+TFREZ)*WSNOW[I]/RHOW)/(RUNOFF[I]+WSNOW[I]/RHOW)
                    
                    RUNOFF[I]=RUNOFF[I]+WSNOW[I]/RHOW
                    if OVRFLW[I]>0. or WSNOW[I]>0.:
                        TOVRFL[I]=(TOVRFL[I]*OVRFLW[I]+(TSNOW[I]+TFREZ)*FI[I]*WSNOW[I]/RHOW)/(OVRFLW[I]+FI[I]*WSNOW[I]/RHOW)
                    
                    OVRFLW[I]=OVRFLW[I]+FI[I]*WSNOW[I]/RHOW
                    TSNOW[I]=0.0
                    WSNOW[I]=0.0
                    QFN[I]=QFN[I]-FI[I]*ZREM*RHOW/DELT
                    QFG[I]=QFG[I]+FI[I]*EVAP[I]*RHOW
                
            HTCS[I]=HTCS[I]+FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
    #Return Vars
    return [RHOSNO,ZSNOW,HCPSNO,TSNOW,EVAP,QFN,QFG,HTCS,\
        WLOST,TRUNOF,RUNOFF,TOVRFL,OVRFLW,WSNOW]

#------- TMELT ----------------------------------------------------------------------
def TMELT (ZSNOW,TSNOW,QMELT,R,TR,GZERO,RALB,\
        HMFN,HTCS,FI,HCPSNO,RHOSNO,WSNOW,\
        ISAND,IL1,IL2,\
        DELT,TFREZ,RHOW,HCPICE,RHOICE,HCPW,CLHMLT):
    #HTC
    #Empty Vars (none)
    for I in range(IL1-1,IL2):
        if FI[I]>0.:
            if QMELT[I]>0. and ZSNOW[I]>0.:
                HTCS[I]=HTCS[I]-FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
                HADD =QMELT[I]*DELT                                                            
                HCONV=(0.0-TSNOW[I])*HCPSNO[I]*ZSNOW[I] + CLHMLT*RHOSNO[I]*ZSNOW[I]                         
                if HADD<=HCONV:
                    ZMELT=HADD/((0.0-TSNOW[I])*HCPSNO[I]+CLHMLT*RHOSNO[I])
                    RMELTS=ZMELT*RHOSNO[I]/(RHOW*DELT)
                    RMELT=RMELTS*1
                    TRMELT=0.0
                    ZSNOW[I]=ZSNOW[I]-ZMELT
                    HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                    HTCS [I]=HTCS[I]-FI[I]*(QMELT[I]-CLHMLT*RMELT*RHOW)
                else:                                                                        
                    RMELTS=ZSNOW[I]*RHOSNO[I]/RHOW
                    RMELT=RMELTS+WSNOW[I]/RHOW
                    HADD=HADD-HCONV
                    TRMELT=HADD/(HCPW*RMELT)
                    RMELT=RMELT/DELT
                    RMELTS=RMELTS/DELT
                    ZSNOW [I]=0.0
                    HCPSNO[I]=0.0
                    TSNOW [I]=0.0
                    WSNOW [I]=0.0
                    HTCS [I]=HTCS[I]-FI[I]*(QMELT[I]-CLHMLT*RMELTS*RHOW-HADD/DELT)
                         
                HMFN [I]=HMFN[I]+FI[I]*CLHMLT*RMELTS*RHOW
                TR   [I]=(R[I]*TR[I]+RMELT*TRMELT)/(R[I]+RMELT)
                R    [I]=R[I]+RMELT
                QMELT[I]=0.0
                HTCS[I]=HTCS[I]+FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)* ZSNOW[I]/DELT
            
            RALB[I]=R[I]*1
    
    for I in range(IL1-1,IL2):
        if FI[I]>0.:
            if QMELT[I]>0 and ISAND[I,0]>-4:
                GZERO[I]=GZERO[I]+QMELT[I]
                HTCS [I]=HTCS[I]-FI[I]*QMELT[I]
                #HTC[I,0]=HTC[I,0]+FI[I]*QMELT[I]
            
            RALB[I]=R[I]*1
        
    #Return Vars
    return [ZSNOW,TSNOW,QMELT,R,TR,GZERO,RALB,\
            HMFN,HTCS,HCPSNO,WSNOW]
    #HTC,

#------- SNINFL ----------------------------------------------------------------------
def SNINFL(R,TR,ZSNOW,TSNOW,RHOSNO,HCPSNO,WSNOW,HTCS,HMFN,PCPG,\
        ROFN,FI,ILG,IL1,IL2,JL,TFREZ,DELT,HCPICE,RHOICE,HCPW,RHOW,CLHMLT):
    #Empty Vars    
    #PCPG=np.zeros((len(range(IL1-1,IL2)),1),np.float32)
    #ROFN=np.zeros((len(range(IL1-1,IL2)),1),np.float32)

    WSNCAP=np.float32(0.04)

    for I in range(IL1-1,IL2):
        if FI[I]>0. and R[I]>0. and ZSNOW[I]>0.:
            HTCS[I]=HTCS[I]-FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
            RAIN=R[I]*DELT
            HRCOOL=TR[I]*HCPW*RAIN
            HRFREZ=CLHMLT*RHOW*RAIN
            HSNWRM=(0.0-TSNOW[I])*HCPSNO[I]*ZSNOW[I]
            HSNMLT=CLHMLT*RHOSNO[I]*ZSNOW[I]
            if HRCOOL>=(HSNWRM+HSNMLT): 
                HRCOOL=HRCOOL-(HSNWRM+HSNMLT)
                ZMELT=ZSNOW[I]*RHOSNO[I]/RHOW
                HMFN[I]=HMFN[I]+FI[I]*CLHMLT*ZMELT*RHOW/DELT
                HTCS[I]=HTCS[I]+FI[I]*CLHMLT*ZMELT*RHOW/DELT
                TR[I]=HRCOOL/(HCPW*(ZMELT+RAIN+WSNOW[I]/RHOW))
                R[I]=R[I]+(ZMELT+WSNOW[I]/RHOW)/DELT
                ZSNOW[I]=0.0
                TSNOW[I]=0.0
                RHOSNO[I]=0.0
                HCPSNO[I]=0.0
                WSNOW[I]=0.0
            elif HRCOOL>=HSNWRM and HRCOOL<(HSNWRM+HSNMLT):
                HSNMLT=HRCOOL-HSNWRM
                ZMELT=HSNMLT/(CLHMLT*RHOSNO[I])
                HMFN[I]=HMFN[I]+FI[I]*CLHMLT*ZMELT*RHOSNO[I]/DELT
                HTCS[I]=HTCS[I]+FI[I]*CLHMLT*ZMELT*RHOSNO[I]/DELT
                ZSNOW[I]=ZSNOW[I]-ZMELT
                WAVAIL=ZMELT*RHOSNO[I]+WSNOW[I]
                if WAVAIL>(WSNCAP*ZSNOW[I]*RHOSNO[I]):
                    WSNOW[I]=WSNCAP*ZSNOW[I]*RHOSNO[I]
                    ZMELT=(WAVAIL-WSNOW[I])/RHOW
                else:
                    WSNOW[I]=WAVAIL*1
                    ZMELT=0.0
                
                TSNOW[I]=0.0
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                TR[I]=0.0
                R[I]=R[I]+ZMELT/DELT
            elif HSNWRM>=(HRCOOL+HRFREZ):
                HSNWRM=(HRCOOL+HRFREZ)-HSNWRM
                HMFN[I]=HMFN[I]-FI[I]*HRFREZ/DELT
                HTCS[I]=HTCS[I]-FI[I]*HRFREZ/DELT
                RHOSNO[I]=(RHOSNO[I]*ZSNOW[I]+RHOW*RAIN)/ZSNOW[I]
                if RHOSNO[I]>RHOICE:
                    ZSNOW[I]=RHOSNO[I]*ZSNOW[I]/RHOICE
                    RHOSNO[I]=RHOICE*1
                
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                TSNOW[I]=HSNWRM/(HCPSNO[I]*ZSNOW[I])
                TR[I]=0.0
                R[I]=0.0
            elif HSNWRM>=HRCOOL and HSNWRM<(HRCOOL+HRFREZ):
                HRFREZ=HSNWRM-HRCOOL
                ZFREZ=HRFREZ/(CLHMLT*RHOW)
                HMFN[I]=HMFN[I]-FI[I]*CLHMLT*ZFREZ*RHOW/DELT
                HTCS[I]=HTCS[I]-FI[I]*CLHMLT*ZFREZ*RHOW/DELT
                RHOSNO[I]=(RHOSNO[I]*ZSNOW[I]+RHOW*ZFREZ)/ZSNOW[I]
                if RHOSNO[I]>RHOICE:
                    ZSNOW[I]=RHOSNO[I]*ZSNOW[I]/RHOICE
                    RHOSNO[I]=RHOICE*1
                
                WAVAIL=(RAIN-ZFREZ)*RHOW+WSNOW[I]
                if WAVAIL>(WSNCAP*ZSNOW[I]*RHOSNO[I]):
                    WSNOW[I]=WSNCAP*ZSNOW[I]*RHOSNO[I]
                    WAVAIL=WAVAIL-WSNOW[I]
                else:
                    WSNOW[I]=WAVAIL*1
                    WAVAIL=0.0
                
                HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                R[I]=WAVAIL/(RHOW*DELT)
                TR[I]=0.0
                TSNOW[I]=0.0
            
            HTCS[I]=HTCS[I]+FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
            PCPG[I]=PCPG[I]+FI[I]*R[I]*RHOW
            ROFN[I]=ROFN[I]+FI[I]*R[I]*RHOW        
    
    #Return Vars
    return [R,TR,ZSNOW,TSNOW,RHOSNO,HCPSNO,WSNOW,HTCS,\
            HMFN,PCPG,ROFN]

#------- SNOALBW ----------------------------------------------------------------------
def SNOALBW(ALBSNO,RHOSNO,ZSNOW,HCPSNO,TSNOW,FI,S,RMELT,WSNOW,RHOMAX,\
        IL1,IL2,DELT,RHOW,HCPICE,RHOICE,HCPW):
    #Empty Vars (none)
    IPTBAD=np.float32(0)                                                          
    for I in range(IL1-1,IL2):
        if ZSNOW[I]>0. and FI  [I]>0. and S[I]*DELT<1.0E-4:
            if ALBSNO[I]>0.5001 and (RMELT[I]>1.0E-7 or TSNOW[I]>=-0.01):
                ALBSNO[I]=(ALBSNO[I]-0.50)*math.exp(-0.01*DELT/3600.0)+ 0.50
            elif ALBSNO[I]>0.7001 and RMELT[I]<=1.0E-7:
                ALBSNO[I]=(ALBSNO[I]-0.70)*math.exp(-0.01*DELT/3600.0)+0.70                                                                    
    #                                                                       
        if FI[I]>0. and ZSNOW[I]>0.0001:
            if TSNOW[I]<-0.01:
                RHOMAX[I]=450.0-(204.7/ZSNOW[I])*(1.0-math.exp(-ZSNOW[I]/0.673))
            else:
                RHOMAX[I]=700.0-(204.7/ZSNOW[I])*(1.0-math.exp(-ZSNOW[I]/0.673))                        
                                                            
    #                                                                       
        if FI[I]>0. and ZSNOW[I]>0.0001 and RHOSNO[I]<(RHOMAX[I]-0.01):
            RHOOLD=RHOSNO[I]
            RHOSNO[I]=(RHOSNO[I]-RHOMAX[I])*math.exp(-0.01*DELT/3600.0)+RHOMAX[I]
            ZSNOW[I]=ZSNOW[I]*RHOOLD/RHOSNO[I]
            HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
                                                                
        if (ALBSNO[I]<0.49 or ALBSNO[I]>1.0) and ZSNOW [I]>0. and FI[I]>0.:
            IPTBAD=I+1        
    #                                                                       
    if IPTBAD!=0:
        #WRITE(6,6100) IPTBAD,JL,ALBSNO(IPTBAD)                         
        #6100    FORMAT('0AT (I,J)= (',I3,',',I3,'), ALBSNO = ',F10.5)          
        #write error code, but for now just exit
        XIT('SNOALBW',-1)                                         
    
    #Return Vars
    return [ALBSNO,RHOSNO,ZSNOW,HCPSNO]

#------- SNOADD ----------------------------------------------------------------------
def SNOADD(ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS,FI,S,TS,RHOSNI,WSNOW,IL1,IL2,\
        TFREZ,DELT,HCPICE,RHOICE,HCPW,RHOW):
    #Empty Vars (none)
    for I in range(IL1-1,IL2):
        if FI[I]>0. and S[I]>0.:
            HTCS  [I]=HTCS[I]-FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT
            SNOFAL=S[I]*DELT
            if SNOFAL>=1.E-4:
                #   ipy test
                #              if(SNOFAL>=0.005)                               THEN
                ALBSNO[I]=0.84
            elif not ZSNOW[I]>0.:
                ALBSNO[I]=0.50
            
            HCPSNP=HCPICE*RHOSNI[I]/RHOICE
            TSNOW [I]=((TSNOW[I]+TFREZ)*ZSNOW[I]*HCPSNO[I] + (TS   [I]+TFREZ)*SNOFAL  *HCPSNP)\
                /(ZSNOW[I]*HCPSNO[I] + SNOFAL*HCPSNP) - TFREZ
            RHOSNO[I]=(ZSNOW[I]*RHOSNO[I] + SNOFAL*RHOSNI[I])/(ZSNOW[I]+SNOFAL)
            ZSNOW [I]=ZSNOW[I]+SNOFAL
            HCPSNO[I]=HCPICE*RHOSNO[I]/RHOICE+HCPW*WSNOW[I]/(RHOW*ZSNOW[I])
            HTCS  [I]=HTCS[I]+FI[I]*HCPSNO[I]*(TSNOW[I]+TFREZ)*ZSNOW[I]/DELT

    #Return Vars
    return [ALBSNO,TSNOW,RHOSNO,ZSNOW,HCPSNO,HTCS]

#------- DRCOEFL ----------------------------------------------------------------------
def DRCOEFL(VA,T0,TA,QA,PRES,ZREFM,ZREFH,FICE,FLS,IL1,IL2,\
        GRAV,VKC,TFREZ):
    #Empty Vars (none)
    G=GRAV
    TOL=np.float32(1.0E-8)
    TOL2=np.float32(1.0E-2)
    ITMAX=np.float32(100)
    CHARN=np.float32(0.0175)
    PI=np.float32(3.14159)
    CDM=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CDH=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CDMN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    
    #
    #======================================================================
    for I in range(IL1-1,IL2):
        #-----------------------------------------------------------------------
        # NEUTRAL DRAG COEFFICIENTS  
        #  Iterative routine to solve for CDN based on Charnock relation for 
        #  roughness
        #
        CDHNW=np.float32(1.35E-3)
        CDMNW=np.float32(1.0E-3)    #initial trial value
        VSQ=VA[I]*VA[I]
        RESID=999.0
        while abs(RESID) > TOL:
            #if (abs(RESID) <= TOL ) 
            #    break;%EXIT
            #end
            Cold=CDMNW*1
            CDMNW=VKC*VKC/(math.log(ZREFM[I]*G/(CHARN*Cold*VSQ))*math.log(ZREFM[I]*G/(CHARN*Cold*VSQ)))
            RESID=CDMNW-Cold
        
        #
        CDMNI=(VKC/(math.log(ZREFM[I]/0.002)))**2
        CDHNI=(VKC/(math.log(ZREFH[I]/0.00067)))**2
        if FICE[I]>(FLS[I]+0.001):
            XICE=MAX(FICE[I]-FLS[I],0.0)
            CDMN[I]=(XICE*CDMNI+(1.0-FICE[I])*CDMNW)/(1.0-FLS[I])
            CDHN=(XICE*CDHNI+(1.0-FICE[I])*CDHNW)/(1.0-FLS[I])
        else:
            CDMN[I]=CDMNW*1
            CDHN=CDHNW*1
        
        #-----------------------------------------------------------------------
        # INITIAL TRIAL VALUES FOR TRANSFER COEFFICIENTS: SET TO NEUTRAL VALUES
        #
        CDH[I]=CDHN*1
        CDM[I]=CDMN[I]*1
        #-----------------------------------------------------------------------
        # ITERATIVELY COMPUTE DRAG COEFFICIENTS UNTIL M.O. LENGTH CONVERGES
        #
        RESID=999.0
        MOL=9999.0
        ITER=0
        loopFlag=True
        while loopFlag:    
            if abs(RESID) <= TOL2 or ITER < ITMAX: 
                break            
            #-----------------------------------------------------------------------
            # HEAT FLUXES
            #----------------------------------------------------
            #     * CALCULATION OF EASAT CONSISTENT WITH CLASSI
            #     * BUT CONSTANTS DifFER FROM ROGERS&YAU
            #     * Rogers and Yau values
            #         CA=17.67
            #         CB=29.65
            #----------------------------------------------------
            if T0[I]>=TFREZ:
                CA=17.269                                                     
                CB=35.86                                                        
            else:
                CA=21.874
                CB=7.66
                                                                   
            SHF=CDH[I]*VA[I]*(T0[I]-TA[I])
            EASAT=611.0*math.exp(CA*(T0[I]-TFREZ)/(T0[I]-CB))
            QASAT=0.622*EASAT/(PRES[I]-0.378*EASAT)
            LHF=CDH[I]*VA[I]*(QASAT-QA[I])

            #-----------------------------------------------------------------------
            # VIRTUAL TEMPERATURE AND M.-O. LENGTH
            #
            TVIRT=TA[I]*(1.0+0.61*QA[I])
            MOLold=MOL*1
            MOL = -VA[I]*VA[I]*VA[I]*CDM[I]*math.sqrt(CDM[I])*TVIRT/( VKC*G*(SHF + 0.61*LHF*TA[I]) )
            Z_L = ZREFM[I]/MOL
            RESID=MOL-MOLold

            #-----------------------------------------------------------------------
            # STABILITY CORRECTIONS
            #
            #
            #- UNSTABLE CASE
            #---------------
            if Z_L < 0.0:
                X = (1.0 - (16.0*Z_L))**0.25
                PSIM = 2.0*math.log((1.0+X)/2.0) + math.log((1.0+X*X)/2.0) - 2.0*math.atan(X) + PI/2.0
                PSIH = 2.0*math.log((1.0+X*X)/2.0)
                #
                #- STABLE CASE
                #-------------
            elif Z_L >= 0 and Z_L < 0.5:
                PSIM = -5.0*Z_L
                PSIH=PSIM*1
            elif Z_L >= 0.5 and Z_L < 10.0:
                PSIM = (0.5/(Z_L*Z_L)) - (4.25/Z_L) - 7.0*math.log(Z_L) - 0.852
                PSIH=PSIM*1
            else: 
                PSIM = math.log(Z_L) - 0.76*Z_L - 12.093
                PSIH=PSIM*1           

            #-----------------------------------------------------------------------
            # RECOMPUTE DRAG COEFFICIENTS WITH STABILTY CORRECTIONS
            #
            DENOM = (1.0 + (CDMN[I]/(VKC*VKC))*(PSIM*PSIM - (2.0*VKC*PSIM/math.sqrt(CDMN[I]))) )
            if DENOM<1.0E-6:
                DENOM=1.0E-6
            
            #        if(abs(DENOM)<1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
            CDM[I] = CDMN[I]/DENOM
            DENOM = (1.0 + (CDHN/(VKC*VKC))*(PSIM*PSIH -(VKC*PSIH/math.sqrt(CDMN[I]))\
                -(VKC*PSIM*math.sqrt(CDMN[I])/CDHN)))
            if DENOM<1.0E-6:
                DENOM=1.0E-6
            
            #        if(abs(DENOM)<1.0E-6) DENOM=SIGN(1.0E-6,DENOM)
            CDH[I] = CDHN/DENOM
            ITER=ITER+1
        

        if ITER>=ITMAX:
            CDM[I]=CDMN[I]*1
            CDH[I]=CDHN*1
        

        #     if (ITER >= ITMAX) print*, "** max iters reached"
        if CDH[I] < 0.0: 
            CDH[I]=CDHN*1
        
        if CDM[I] < 0.0:
            CDM[I]=CDMN[I]*1
        

        CDHRAT=CDH[I]/CDHN
        if CDHRAT >= 8.0:
            CDH[I]=CDHN*1        

    #Return Vars
    return [CDM,CDH,CDMN]

#------- FREECONV ----------------------------------------------------------------------
def FREECONV(LKICEH,T0,TLAK,RHOIW,NLAK,IL1,IL2,\
        TFREZ,DELSKIN,DELZLK):
    #This function mixes based on density
    #If a layer is more dense than the one below it mixes them
    LKICEH = EnsurePassByValue(LKICEH)
    T0 = EnsurePassByValue(T0)
    TLAK = EnsurePassByValue(TLAK)    
    #Empty Vars (none)
    for I in range(IL1-1,IL2):
        if LKICEH[I] <= 0.0:
            TC1=T0[I]-TFREZ
            TC2=TLAK[I,0]-TFREZ
            [XXX,RHO1]=EQNST(TC1,0.05)
            [XXX,RHO2]=EQNST(TC2,0.5)
            if RHO1 > RHO2:
                TBAR=((DELSKIN*RHO1*T0[I])+(DELZLK*RHO2*TLAK[I,0]))/((DELSKIN*RHO1)+(DELZLK*RHO2))
                T0[I]=TBAR*1
                TLAK[I,0]=TBAR*1
            
        ICEBOT=RHOIW*LKICEH[I]

        NMIX=1
        for J in range(0,int(NLAK[I])-1):
            ZTOP=DELSKIN + np.float32(J)*DELZLK
            ZBOT=ZTOP+DELZLK
            if ICEBOT <= ZTOP:
                TC1=TLAK[I,J]-TFREZ
                TC2=TLAK[I,J+1]-TFREZ
                #mdm       TTEST=(TC1-3.9816)*(TC2-3.9816)
                TTEST=(TC1-np.float32(3.98275))*(TC2-np.float32(3.98275))
                [XXX,RHO1]=EQNST(TC1,ZBOT)
                [XXX,RHO2]=EQNST(TC2,ZBOT+DELZLK)
                #--------- MIX LAYERS if RHO1>RHO2 OR TEMPERATURES SPAN 
                #--------- T_MAXDENSITY=3.9816 C.
                if (RHO1 > RHO2) or (TTEST < 0.0):
                    TBAR=((NMIX*RHO1*TLAK[I,J])+(RHO2*TLAK[I,J+1]))/(NMIX*RHO1+RHO2)
                    #need to mix all layers above with > densities
                    for K in range(J-NMIX+1,J+2):#DO 430, K=J-NMIX+1,J+1                        
                        TLAK[I,K]=TBAR
                    
                    NMIX=NMIX+1
                    #           WRITE(6,6666) "static instability removed under ice:",
                    #    >                IYEAR,IDAY,IHOUR,IMIN,j*DELZLK
                else:
                    NMIX=1
                
    #Return Vars
    return [T0,TLAK]

#------- LKTRANS ----------------------------------------------------------------------
def LKTRANS (BLAK,IL1,IL2):
    CQ1A=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ2A=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ2B=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ3A=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ3B=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ1BI=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ2BI=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ3BI=np.zeros((len(range(IL1-1,IL2))),np.float32)
    CQ1B=np.zeros((len(range(IL1-1,IL2))),np.float32)    
    for I in range(IL1-1,IL2):
        # FIXED WATER VALUES
        CQ1A[I]=np.float32(0.5817)
        CQ2A[I]=np.float32(0.4183)
        CQ2B[I]=np.float32(6.89)
        CQ3A[I]=np.float32(0.0)
        CQ3B[I]=np.float32(69.0)
        # FIXED ICE VALUES (from Patterson and Hamblin, 1988, L&O)
        CQ1BI[I]=np.float32(1.5)
        # Cmdm  CQ1BI[I]=3.75    !test - high extinction for snow ice
        CQ2BI[I]=np.float32(20.0)
        CQ3BI[I]=np.float32(69.0)

        #======================================================================
        # CQ1B NOW READ IN .INI FILE
        #----------------------------------------------------------------------
        CQ1B[I]=BLAK[I]
    return [CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,CQ1BI,CQ2BI,CQ3BI]

#------- TSOLVL ----------------------------------------------------------------------
def TSOLVL(TLAK,T0,LKICEH,EVAP,ALVS,\
           QSWIN,QLWIN,TA,QA,VA,PRES,RHOAIR,CDH,\
           GZEROL,QTRANSL,HTCL,FLS,FICE,ALBW,ALBI,\
           CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,\
           CQ1BI,CQ2BI,CQ3BI,IL1,IL2,\
           DELZLK,DELSKIN,RHOICE,RHOW,TFREZ,\
           CLHVAP,SPHAIR,EMSW,SBC,TCW,DELT,HCPW,\
           CLHMLT,TCICE,HCPICE):
    #Empty Vars
    KLAK=np.zeros((len(range(IL1-1,IL2))),np.float32)
    GLAK=np.zeros((len(range(IL1-1,IL2))),np.float32)
    Q0SAT=np.zeros((len(range(IL1-1,IL2))),np.float32)
    KSTAR=np.zeros((len(range(IL1-1,IL2))),np.float32)
    LSTAR=np.zeros((len(range(IL1-1,IL2))),np.float32)
    QSENS=np.zeros((len(range(IL1-1,IL2))),np.float32)
    QEVAP=np.zeros((len(range(IL1-1,IL2))),np.float32)
    G0=np.zeros((len(range(IL1-1,IL2))),np.float32) 
    #
    # ----* LOCAL PARAMETER DEFINITIONS *-------------------------------
    #
    DZ=DELZLK
    DS=DELSKIN
    RHOIW=RHOICE/RHOW
    XICE=np.zeros((len(range(IL1-1,IL2))),np.float32)     
    #----------------------------------------------------------------------
    for I in range(IL1-1,IL2):
        #======================================================================
        # COMPUTE SURFACE ENERGY BALANCE FOR CURRENT TIMESTEP
        #======================================================================
        # COMPUTE SW FLUXES --------------------------------------------------
        #   KSTAR is net SW at surface of skin
        #   KLAK is penetrating SW at base of skin
        #
        XICE[I]=MAX(FICE[I]-FLS[I],0.0)
        if FICE[I]>(FLS[I]+0.001):
            ALBTOT=(XICE[I]*ALBI[I]+(np.float32(1)-FICE[I])*ALBW[I])/(np.float32(1.0)-FLS[I])
        else:
            ALBTOT=ALBW[I]*np.float32(1)   
        KSTAR[I]=(1.-FLS[I])*(np.float32(1)-ALBTOT)*QSWIN[I]+FLS[I]*QTRANSL[I]
        if KSTAR[I] < 0.0:
            KSTAR[I]=np.float32(0.0)        
        ALVS[I]=copy.deepcopy(ALBTOT)
        #ALIR[I]=copy.deepcopy(ALBTOT)

        #--- Attenuation through ice (now includes attenuation through leads)
        #---
        ICEBOT=RHOIW*LKICEH[I]      #   !ICE DRAFT
        ICETOP=LKICEH[I]-ICEBOT     #   !ICE FREEBOARD
        if LKICEH[I] <= 0.0:        #   !NO ICE 
            ATTEN1=CQ1B[I]*DS
            ATTEN2=CQ2B[I]*DS
            ATTEN3=CQ3B[I]*DS
        else:
            if ICEBOT > DS:         #   !Z inside ice
                ATTEN1=FICE[I]*(CQ1BI[I]*(DS+ICETOP)) + (np.float32(1)-FICE[I])*CQ1B[I]*DS
                ATTEN2=FICE[I]*(CQ2BI[I]*(DS+ICETOP)) + (np.float32(1)-FICE[I])*CQ2B[I]*DS
                ATTEN3=FICE[I]*(CQ3BI[I]*(DS+ICETOP)) + (np.float32(1)-FICE[I])*CQ3B[I]*DS
            else:
                ATTEN1=FICE[I]*(CQ1BI[I]*LKICEH[I] + CQ1B[I]*(DS-ICEBOT)) + (np.float32(1)-FICE[I])*CQ1B[I]*DS
                ATTEN2=FICE[I]*(CQ2BI[I]*LKICEH[I] + CQ2B[I]*(DS-ICEBOT)) + (np.float32(1)-FICE[I])*CQ2B[I]*DS
                ATTEN3=FICE[I]*(CQ3BI[I]*LKICEH[I] + CQ3B[I]*(DS-ICEBOT)) + (np.float32(1)-FICE[I])*CQ3B[I]*DS
            
        KLAK[I]=KSTAR[I]*(CQ1A[I]*math.exp(-ATTEN1) + CQ2A[I]*math.exp(-ATTEN2) + CQ3A[I]*math.exp(-ATTEN3) )
        #
        # COMPUTE TURBULENT FLUXES -------------------------------------------
        #     * CALCULATION OF E0SAT CONSISTENT WITH CLASSI
        #     * BUT CONSTANTS DifFER FROM ROGERS&YAU
        #     * Rogers and Yau values
        #         CA=17.67
        #         CB=29.65
        #
        if T0[I]>=TFREZ:
            CA=17.269                                       
            CB=35.86                                       
        else:
            CA=21.874
            CB=7.66                         

        if FICE[I]>(FLS[I]+0.001):
            CPHCH=(XICE[I]*(CLHMLT+CLHVAP)+(1.-FICE[I])*CLHVAP)/(np.float32(1.0)-FLS[I])
        else:
            CPHCH=CLHVAP*np.float32(1)

        QSENS[I]=(np.float32(1)-FLS[I])*RHOAIR[I]*SPHAIR*CDH[I]*VA[I]*(T0[I]-TA[I])
        E0SAT=np.float32(611.0)*math.exp(CA*(T0[I]-TFREZ)/(T0[I]-CB))
        Q0SAT[I]=np.float32(0.622)*E0SAT/(PRES[I]-np.float32(0.378)*E0SAT)
        EVAP[I]=(np.float32(1)-FLS[I])*RHOAIR[I]*CDH[I]*VA[I]*(Q0SAT[I]-QA[I])
        #mdm    QEVAP[I]=CPHCH*EVAP[I]
        QEVAP[I]=(np.float32(1)-FLS[I])*RHOAIR[I]*CPHCH*CDH[I]*VA[I]*(Q0SAT[I]-QA[I])
        #
        # COMPUTE NET LW AND NET SFC ENERGY -------------------------------
        # GLAK IS THERMAL FLUX AT BASE OF SKIN.  THERMAL CONDUCTIVITY BASED
        # ON WEIGHTED AVERAGE FOR WATER AND ICE if ICE PRESENT IN LAYER
        #
        LSTAR[I]=(np.float32(1)-FLS[I])*EMSW*(QLWIN[I]-SBC*T0[I]*T0[I]*T0[I]*T0[I])
        G0[I]=KSTAR[I]-KLAK[I]+LSTAR[I]-QSENS[I]-QEVAP[I]+HTCL[I]+FLS[I]*GZEROL[I]
        if ICEBOT >= DS: 
            GLAK[I]=(np.float32(-2.0)*TCICE/(DZ+DS))*(TLAK[I,0]-T0[I])
        elif ICEBOT < DS and LKICEH[I] > 0.0:
            TC=(ICEBOT*TCICE + (DS-ICEBOT)*TCW)/DS
            #mdm      TC=20.0          !mdm test
            GLAK[I]=(np.float32(-2.0)*TC/(DZ+DS))*(TLAK[I,0]-T0[I])
        else:
            GLAK[I]=(np.float32(-2.0)*TCW/(DZ+DS))*(TLAK[I,0]-T0[I])        
        #-----NET ENERGY FLUX INTO SKIN (W/M2)
        ESKIN= G0[I] - GLAK[I]
        #
        # STEP FORWARD SKIN TEMP T0
        #
        T0old=copy.deepcopy(T0[I])
        if LKICEH[I] <= 0.0:
            T0[I] = T0[I] + (DELT/(DS*HCPW))*ESKIN
        elif LKICEH[I] > 0.0 and ICEBOT <= DS: 
            T0[I] = T0[I] + (DELT/((LKICEH[I]*HCPICE) + (DS-ICEBOT)*HCPW))*ESKIN
        else:
            T0[I] = T0[I] + (DELT/((DS+ICETOP)*HCPICE))*ESKIN        
        #
        # ICE GROWTH OR DECAY
        #
        if ESKIN < 0.0 and ICEBOT < DS:
            #-----NET ENERGY FLUX USED TO LOWER T0 TO TFREZ 
            if T0old > TFREZ:
                ECOOL=(DS-ICEBOT)*HCPW*(T0old-TFREZ)/DELT
            else:
                ECOOL=np.float32(0.0)           
        #-----REMAINING ENERGY FLUX (if ANY) USED TO FREEZE ICE
            EAVAIL=ESKIN+ECOOL
            if EAVAIL < 0.0:
                NEWICE=-(DELT/(RHOICE*CLHMLT))*EAVAIL
                LKICEH[I]=LKICEH[I]+NEWICE
                ICEBOT=RHOIW*LKICEH[I]
                T0[I]=TFREZ*1
                #-----LIMIT ICE GROWTH TO THE CURRENT LAYER
                if ICEBOT > DS:
                    EHEAT=(RHOICE*CLHMLT*(ICEBOT-DS))/DELT
                    T0[I]=TFREZ - (EHEAT*DELT)/(DS*HCPICE)
                    LKICEH[I]=DS/RHOIW                

        if ESKIN > 0.0 and LKICEH[I] > 0.0:
            #-----NET ENERGY FLUX USED FIRST TO RAISE T TO ZERO
            if ICEBOT <= DS:
                EHEAT=LKICEH[I]*HCPICE*(TFREZ-T0old)/DELT
            else:
                EHEAT=(DS+ICETOP)*HCPICE*(TFREZ-T0old)/DELT

            #-----NET ENERGY FLUX USED TO MELT ICE
            if ESKIN > EHEAT:
                NEWICE=-(DELT/(RHOICE*CLHMLT))*(ESKIN-EHEAT)
                LKICEH[I]=LKICEH[I]+NEWICE
                T0[I]=TFREZ
                #-----LIMIT ICE MELT TO THE CURRENT LAYER
                if LKICEH[I] < 0.0:
                    EHEAT=-(RHOICE*CLHMLT*LKICEH[I])/DELT
                    T0[I]=TFREZ + (EHEAT*DELT)/(DS*HCPW)
                    LKICEH[I]=np.float32(0.0)                
        # -----------------
    #Return Vars
    return [T0,LKICEH,KLAK,GLAK,Q0SAT,KSTAR,LSTAR,QSENS,\
           QEVAP,EVAP,ALVS,G0]

#------- MIXLYR ----------------------------------------------------------------------
def MIXLYR(DTEMP,NLAK,USTAR,IL1,IL2,\
        HDPTH,TKE,DELU,EXPW,QSTAR,\
        HLAK,LLAK,GRED,\
        CQ1A,CQ1B,CQ2A,CQ2B,CQ3A,CQ3B,RHOMIX,\
        LSTAR,QSENS,QEVAP,LKICEH,\
        TKECN,TKECF,TKECE,TKECS,TKECL,GRAV,\
        DELSKIN,SPHW,DELT,DHMAX,DELZLK,HDPTHMIN,\
        TKEMIN,DUMAX):
    #Ensure no overwrite
    DTEMP = EnsurePassByValue(DTEMP)
    NLAK = EnsurePassByValue(NLAK)
    USTAR = EnsurePassByValue(USTAR)
    HDPTH = EnsurePassByValue(HDPTH)
    TKE = EnsurePassByValue(TKE)
    DELU = EnsurePassByValue(DELU)
    EXPW = EnsurePassByValue(EXPW)
    QSTAR = EnsurePassByValue(QSTAR)
    HLAK = EnsurePassByValue(HLAK)
    LLAK = EnsurePassByValue(LLAK)    
    GRED = EnsurePassByValue(GRED)
    CQ1A = EnsurePassByValue(CQ1A)    
    CQ1B = EnsurePassByValue(CQ1B)
    CQ2A = EnsurePassByValue(CQ2A)
    CQ2B = EnsurePassByValue(CQ2B)
    CQ3A = EnsurePassByValue(CQ3A)
    CQ3B = EnsurePassByValue(CQ3B)
    RHOMIX = EnsurePassByValue(RHOMIX)
    LSTAR = EnsurePassByValue(LSTAR)
    QSENS = EnsurePassByValue(QSENS)
    QEVAP = EnsurePassByValue(QEVAP)
    LKICEH = EnsurePassByValue(LKICEH)
    #DELSKIN = EnsurePassByValue(DELSKIN)
    #SPHW = EnsurePassByValue(SPHW)
    #DELT = EnsurePassByValue(DELT)
    #DHMAX = EnsurePassByValue(DHMAX)
    #DELZLK = EnsurePassByValue(DELZLK)
    #HDPTHMIN = EnsurePassByValue(HDPTHMIN)
    #TKEMIN = EnsurePassByValue(TKEMIN)
    #DUMAX = EnsurePassByValue(DUMAX)


    #Empty Vars (none)
    FQU=np.zeros((len(range(IL1-1,IL2))),np.float32)
    BFLX=np.zeros((len(range(IL1-1,IL2))),np.float32)
    DISS=np.zeros((len(range(IL1-1,IL2))),np.float32)
    FSHEAR=np.zeros((len(range(IL1-1,IL2))),np.float32)
    FENTRA=np.zeros((len(range(IL1-1,IL2))),np.float32)
    TRAN=np.zeros((len(range(IL1-1,IL2))),np.float32)
    DHDT=np.zeros((len(range(IL1-1,IL2))),np.float32)    
    CN = TKECN
    CF = TKECF
    CE = TKECE
    CS = TKECS
    CL = TKECL
    G  = GRAV
    #
    #=======================================================================
    for I in range(IL1-1,IL2):
    #-----------------------------------------------------------------------
        HDPTH_OLD=HDPTH[I]
        TKE_OLD=TKE[I]
        DELU_OLD=DELU[I]        
        #
        #=======================================================================
        # MEAN MIXING LAYER TKE 
        #-----------------------------------------------------------------------
        # (1) Buoyancy Flux 	
        #
        LTERM=LSTAR[I]-QSENS[I]-QEVAP[I]
        DEPTH = HDPTH_OLD + DELSKIN
        QH = np.float32(QSTAR[I]*( CQ1A[I]*math.exp(-CQ1B[I]*DEPTH) + CQ2A[I]*math.exp(-CQ2B[I]*DEPTH) + CQ3A[I]*math.exp(-CQ3B[I]*DEPTH) ))        
        QINT=np.float32((2.*QSTAR[I]/DEPTH)*( (CQ1A[I]/CQ1B[I])*(math.exp(-CQ1B[I]*DEPTH)-1.) + (CQ2A[I]/CQ2B[I])*(math.exp(-CQ2B[I]*DEPTH)-1.) + (CQ3A[I]/CQ3B[I])*(math.exp(-CQ3B[I]*DEPTH)-1.) ))
        QTERM=QSTAR[I]+QH+QINT
        if QTERM < 0.:
            XIT('MIXLYR',-1)
        
        BFLX[I]=np.float32(0.5)*DEPTH*(G*EXPW[I]/(SPHW*RHOMIX[I]))*(-LTERM - QTERM)# !m3/s3
        #----------
        # Suppress buoyancy production under ice
        #
        if (LKICEH[I] > 0.0):
            BFLX[I]=np.float32(0.0)

        #-----------------------------------------------------------------------
        # (2) Mechanical Forcing 
        #
        FQU[I]= np.float32(0.5)*CN*CN*CN*USTAR[I]*USTAR[I]*USTAR[I]#	!m3/s3

        #-----------------------------------------------------------------------
        # (3) Dissipation and transport of TKE to thermocline
        #       (tendency in TKE due to dissipation and transport)
        #
        DISS[I]=np.float32(0.5)*CE*math.sqrt(TKE_OLD*TKE_OLD*TKE_OLD)#		!m3/s3
        TRAN[I]=np.float32(0.5)*CF*math.sqrt(TKE_OLD*TKE_OLD*TKE_OLD)#		!m3/s3

        #-----------------------------------------------------------------------
        # (4) TKE (m2/s2)
        #
        #-- Forced tendency
        TKE1= (np.float32(2.0)*DELT/HDPTH_OLD)*( FQU[I]  + BFLX[I] )

        #-- Dissipation and transport tendency
        #-- Solved analytically using HDPTH_OLD for MLD  
        C1=(CE+CF)/HDPTH_OLD
        TKE2= (np.float32(1.0)/( (np.float32(0.5)*C1*DELT)+(np.float32(1.0/math.sqrt(TKE_OLD))) )) * (np.float32(1.0)/( (np.float32(0.5)*C1*DELT)+(np.float32(1.0/math.sqrt(TKE_OLD))) )) - TKE_OLD
        #TKE[I]=4.54
        TKE[I]= TKE_OLD + TKE1 + TKE2
        #TKE[I]= TKE1
        #
        #=======================================================================
        # MIXING LAYER DEPTH (HDPTH) AND TENDENCY (DHDT)
        #
        DENOM=TKE_OLD-CS*DELU_OLD*DELU_OLD+EXPW[I]*G*HDPTH_OLD*DTEMP[I]
        if abs(DENOM) <= 1.E-10:
            if abs(DTEMP[I]) > 1.E-10:
                # *** set H to shear penetration depth (Pollard et al 1973) ****
                # --need to add correction from Spigel at al 1986
                HDPTH[I]=CS*DELU_OLD*DELU_OLD/(EXPW[I]*G*DTEMP[I])
                #         print*, "reset to shear penetration depth=========: ",HDPTH[I]
            else:
                HDPTH[I]=HDPTH_OLD            
        else:
            HDPTH[I]=HDPTH_OLD + DELT*((CF-CL)*math.sqrt(TKE_OLD*TKE_OLD*TKE_OLD))/DENOM
        
        #
        # *** limit deepening to a maximum of DHMAX in 1 timestep
        #
        if (HDPTH[I]-HDPTH_OLD) > DHMAX:
            HDPTH[I] = HDPTH_OLD + DHMAX
        
        DPTHMAX=NLAK[I]*DELZLK
        HDPTH[I]=MIN(MAX(HDPTH[I],HDPTHMIN),DPTHMAX)
        DHDT[I]=(HDPTH[I]-HDPTH_OLD)/DELT

        #-----------------------------------------------------------------------
        # MIXED LAYER RETREAT FOLLOWING RAYNER (1980)
        #
        if TKE[I] <= TKEMIN:
            TKE[I]=TKEMIN*1
            if -BFLX[I] >= 1.0E-15:
                MOL=-HDPTH_OLD*FQU[I]/BFLX[I]#	!Monin-Obukhov length
                #mdm     HDPTH[I]=MAX(HDPTHMIN,MOL)
                HDPTH[I]=MIN(MAX(MOL,HDPTHMIN),DPTHMAX)
            else:
                HDPTH[I]=HDPTHMIN           
            if HDPTH[I] <= 0.:
                XIT('MIXLYR',-2)
            
            DHDT[I]=np.float32(0.0)
            DELU[I]=np.float32(0.0)

        #=======================================================================
        # SHEAR PRODUCTION AND ENTRAINMENT FLUXES DIAGNOSTIC
        #
        FSHEAR[I]=np.float32(0.5)*CS*DELU_OLD*DELU_OLD*DHDT[I]
        FENTRA[I]=np.float32(0.5)*EXPW[I]*G*HDPTH_OLD*DTEMP[I]*DHDT[I]
        #
        #=======================================================================
        # MEAN MIXING LAYER MOMENTUM (DELU)
        #
        DELU[I]=DELU_OLD + DELT*(  (USTAR[I]*USTAR[I]/HDPTH_OLD) - (DHDT[I]*DELU_OLD/HDPTH_OLD) )
        DELU[I]=MAX(DELU[I],np.float32(0.0))
        DELTADU=DELU[I]-DELU_OLD
        if DELTADU > DUMAX:
            DELU[I]=DELU_OLD + DUMAX
        
        #
        #-----------------------------------------------------------------------
        # RESET DELU if MAXIMUM VALUE REACHED
        # (EG Spigel and Imberger, 1980; Spigel et al 1986) 
        #
        H2=np.float32(HLAK[I])-HDPTH[I]
        if H2 >= 1.0E-12 and GRED[I] > 0.0:
            TI=np.float32(2.0*LLAK[I]/(math.sqrt(GRED[I]*HDPTH[I]*H2/HLAK[I])))
        #        TI=0.0		!test mdm to turn off shear term
        else:
            TI=np.float32(0.0)
        UMAX=USTAR[I]*USTAR[I]*TI/(np.float32(4.0)*HDPTH[I])
        if DELU[I] >= UMAX:
            DELU[I]=np.float32(0.0)       
        #-----------------------------------------------------------------------    
    #Return Vars
    return [HDPTH,TKE,DELU,FQU,BFLX,DISS,FSHEAR,FENTRA,TRAN]

#------- SCREENRH ----------------------------------------------------------------------
def ESW(TTT):
    RW1 = np.float32(53.67957)
    RW2 = np.float32(6743.769)
    RW3 = np.float32(4.8451)
    Y = math.exp(RW1+RW2/TTT)*TTT**RW3
    return Y

def ESI(TTT):
    RI1 = np.float32(23.33086)
    RI2 = np.float32(-6111.72784)
    RI3 = np.float32(0.15215)
    Y = math.exp(RI1+RI2/TTT)*TTT**RI3
    return Y

def ESTEFF(TTT,UUU):
    Y = UUU*ESW(TTT) + (np.float32(1.0)-UUU)*ESI(TTT)
    return Y

def SCREENRH(SRH,ST,SQ,PRESSG,FMASK,ILG,IL1,IL2):    
    EPS1=np.float32(0)
    EPS2=np.float32(0)
    T1S=np.float32(0)
    EPSLIM=np.float32(0.001)
    FACTE=np.float32(1.0)/EPS1-np.float32(1.0)

    for IL in range(IL1-1,IL2):
        if FMASK[IL]>0.:
        #
        #       * COMPUTE THE FRACTIONAL PROBABILITY OF WATER PHASE      
        #       * EXISTING AS A FUNCTION OF TEMPERATURE (FROM ROCKEL,     
        #       * RASCHKE AND WEYRES, BEITR. PHYS. ATMOSPH., 1991.)       
        #
            if ST[IL] >= T1S:
                FRACW=np.float32(1.0)
            else:
                FRACW=np.float32(0.0059)+np.float32(0.9941)*np.float32(math.exp(-0.003102*(T1S-ST[IL])**2))
            ETMP=ESTEFF(ST[IL],FRACW)
            ESTREF=np.float32(0.01)*PRESSG[IL]*(np.float32(1.0)-EPSLIM)/(np.float32(1.0)-EPSLIM*EPS2)
            if ETMP < ESTREF:
                ESAT=ETMP
            else:
                ESAT=ESTREF
            QSW=EPS1*ESAT/(np.float32(0.01)*PRESSG[IL]-EPS2*ESAT)
            SRH[IL]=MIN(MAX((SQ[IL]*(np.float32(1.0)+QSW*FACTE))/(QSW*(np.float32(1.0)+SQ[IL]*FACTE)),0.),np.float32(1.0))    
    return [SRH]

#------- SLDIAG ----------------------------------------------------------------------
#    * STABILITY FUNCTIONS FOR THE STABLE CASE 
def PSM(X):
    Y = np.float32(-X -.667*(X-5/.35)*math.exp(-.35*X))
    return Y

def PSE(X):
    Y = np.float32(-(1+.667*X)**1.5 -.667*(X-5/.35)*math.exp(-.35*X))
    return Y

#     * STABILITY defS FOR THE UNSTABLE CASE                       
def Y_fun(X):
    Y=np.float32((1-16*X)**.25)
    return Y

def PIM(X):
    Y= np.float32(math.log((1+X)**2*(1+X**2)) -2*math.atan(X))
    return Y

def PIE(X):
    Y= np.float32(2*math.log(1+X**2))
    return Y



def SLDIAG (SUT,SVT,STT,SQT,CDM,CDH,UA,VA,TA,QA,T0,Q0,Z0M,Z0E,F,ZA,ZU,ZT,IL1,IL2):
    raise NameError('Issues on this script, is not called though')
    PR=np.float32(1.0)                                                           
    for I in range(IL1-1,IL2):
        if F[I]>0.:
            #     * CALCULATION OF SURFACE FLUXES AND MONIN-OBUKHOV LENGTH          
            WSPD=MAX(VMIN,math.sqrt(UA[I]**2+VA[I]**2))
            CM=math.sqrt(CDM[I])            
            US=CM*WSPD
            if abs(TA[I]-T0[I])<0.01:
                TS=-0.01*CDH[I]/CM                                        
            else:                                                           
                TS=CDH[I]*(TA[I]-T0[I])/CM                                             
            if abs(QA[I]-Q0[I])<1.0E-7:
                QS=-1.0E-7*CDH[I]/CM
            else:
                QS=CDH[I]*(QA[I]-Q0[I])/CM            
            L=TA[I]*US**2/(VKC*GRAV*(TS*(1+.61*QA[I])+.61*TA[I]*QS))
            #     * CALCULATE CORRECTION FACTORS TO TAKE INTO ACCOUNT THE APPROXIMATIONS
            #     * IN DRCOEF                                                          
            if L>0.: 
            #     * STABLE CASE                                                     
                UVA=US/VKC*(math.log(ZA[I]/Z0M[I])-PSM(ZA[I]/L)+PSM(Z0M[I]/L))
                RATIO=WSPD/UVA
                UVU=US/VKC*(math.log((ZU[I]+Z0M[I])/Z0M[I])-PSM((ZU[I]+Z0M[I])/L)\
                    +PSM(Z0M[I]/L))*RATIO
                TTA=T0[I]+TS/VKC*PR*(math.log(ZA[I]/Z0E[I])-PSE(ZA[I]/L)\
                    +PSE(Z0E[I]/L))                                          
                RATIO=(TA[I]-T0[I])/SIGN(MAX(abs(TTA-T0[I]),1.E-4),TTA-T0[I])
                CE=(math.log((ZT[I]+Z0M[I])/Z0E[I])-PSE((ZT[I]+Z0M[I])/L)\
                    +PSE(Z0E[I]/L))*RATIO*PR/VKC
            else:
            #     * UNSTABLE CASE                                                   
                UVA=US/VKC*(math.log(ZA[I]/Z0M[I])-PIM(Y_fun(ZA[I]/L))+PIM(Y_fun(Z0M[I]/L)))
                RATIO=WSPD/UVA
                UVU=US/VKC*(math.log((ZU[I]+Z0M[I])/Z0M[I])-PIM(Y_fun((ZU[I]+Z0M[I])/L))\
                    +PIM(Y_fun(Z0M[I]/L)))*RATIO
                TTA=T0[I]+TS/VKC*PR*(math.log(ZA[I]/Z0E[I])-PIE(Y_fun(ZA[I]/L))\
                    +PIE(Y_fun(Z0E[I]/L)))
                RATIO=(TA[I]-T0[I])/SIGN(MAX(abs(TTA-T0[I]),1.E-4),TTA-T0[I])
                CE=(math.log((ZT[I]+Z0M[I])/Z0E[I])-PIE(Y_fun((ZT[I]+Z0M[I])/L))\
                    +PIE(Y_fun(Z0E[I]/L)))*RATIO*PR/VKC                               
            #                                                                       
            SUT[I]=UVU*UA[I]/WSPD
            SVT[I]=UVU*VA[I]/WSPD
            STT[I]=T0[I]+TS*CE
            SQT[I]=Q0[I]+QS*CE        
    return [SUT,SVT,STT,SQT]

#------- DIASURFZ ----------------------------------------------------------------------
def FMI(Z2,Z02,LZZ02,ILMO2):#,X,X0
    # implicit none
    # C
    # REAL, INTENT(IN ) :: Z2,Z02,LZZ02,ILMO2
    # REAL, INTENT(OUT) :: X,X0
    # c
    CI = 40
    BETA = 1.0
    raise NameError('This should never be reached')   
    '''X =(1-CI*Z2 *BETA*ILMO2)**(0.16666666)
    X0=(1-CI*Z02*BETA*ILMO2)**(0.16666666)
    FMIout=LZZ02+math.log((X0+1)**2*math.sqrt(X0**2-X0+1)*(X0**2+X0+1)**1.5\
        /((X+1)**2*math.sqrt(X**2-X+1)*(X**2+X+1)**1.5))\
        +RAC3*math.atan(RAC3*((X**2-1)*X0-(X0**2-1)*X)\
        /((X0**2-1)*(X**2-1)+3*X*X0))'''
    return FMIout


# c
# C   Internal def FHI
# C   Stability def for heat and moisture in the unstable regime (ILMO<0)
# c   Reference: Delage Y. and Girard C. BLM 58 (19-31) Eq. 17
# c
# REAL def FHI(Z2,Z0T2,LZZ0T2,ILMO2,Y,Y0)
def FHI(Z2,Z0T2,LZZ0T2,ILMO2):#,Y,Y0
    # implicit none
    # C
    # REAL, INTENT(IN ) :: Z2,Z0T2,LZZ0T2,ILMO2
    # REAL, INTENT(OUT) :: Y,Y0
    # c
    raise NameError('This should never be reached')   
    '''CI = 40
    BETA = 1.0
    Y =(1-CI*Z2  *BETA*ILMO2)**(0.33333333)
    Y0=(1-CI*Z0T2*BETA*ILMO2)**(0.33333333)
    FHIout=BETA*(LZZ0T2+1.5*math.log((Y0**2+Y0+1)/(Y**2+Y+1))+RAC3\
        *math.atan(RAC3*2*(Y-Y0)/((2*Y0+1)*(2*Y+1)+3)))'''
    return FHIout

# C
# C   Internal def PSI
# C   Stability def for momentum in the stable regime (unsl>0)
# c   Reference :  Y. Delage, BLM, 82 (p23-48) (Eqs.33-37)
# c
# REAL def PSI(Z2,HI2,ILMO2)
def PSI(Z2,HI2,ILMO2):
    # # implicit none
    # C
    # REAL a,b,c,d
    # REAL, INTENT(IN ) :: ILMO2,Z2,HI2
    # c
    AS = np.float32(12.0)
    BETA = np.float32(1.0)
    D = np.float32(4)*AS*BETA*ILMO2
    C = D*HI2 - HI2**np.float32(2)
    B = D - np.float32(2)*HI2
    A = np.float32(math.sqrt(1 + B*Z2 - C*Z2**2))
    PSIout = np.float32(0.5 * (A-Z2*HI2-math.log(1+B*Z2*0.5+A)\
        -B/(2*math.sqrt(C))*math.asin((B-2*C*Z2)/D)))
    return PSIout
 
def DIASURFZ (UZ,VZ,TZ,QZ,U,V,TG,QG,Z0,Z0T,ILMO,ZA,\
        H,UE,FTEMP,FVAP,ZU,ZT,LAT,F,IL1,IL2):
    raise NameError('This should never be reached')    
    RAC3=np.float32(math.sqrt(3.))

    warning('\nDIASURFZ has not been debugged post traslation from MATLAB');
    for J in range(IL1-1,IL2):
        if F[J]>0.0:
            LZZ0T=math.log((ZT[J]+Z0[J])/Z0T[J])
            LZZ0=math.log(ZU[J]/Z0[J]+1)            
            if ILMO[J]<=0.:
                # *---------------------------------------------------------------------
                # *                      UNSTABLE CASE
                # *
                HI=0.
                # *CDIR IEXPAND
                FH=FHI(ZT[J]+Z0[J],Z0T[J],LZZ0T,ILMO[J])
                # *CDIR IEXPAND
                FM=FMI(ZU[J]+Z0[J],Z0 [J],LZZ0 ,ILMO[J])
            else:
            # *---------------------------------------------------------------------
            # *                        STABLE CASE
                HI=1/MAX(HMIN,H[J],(ZA[J]+10*Z0[J])*FACTN,FACTN/(4*AS*BETA*ILMO[J]))
                # *CDIR IEXPAND
                FH=BETA*(LZZ0T+MIN( PSI(ZT[J]+Z0[J],HI,ILMO[J])-PSI(Z0T[J],HI,ILMO[J]),\
                    ASX*ILMO[J]*(ZT[J]+Z0[J]-Z0T[J])))
                # *CDIR IEXPAND
                FM=LZZ0+MIN(PSI(zu[J]+Z0[J],HI,ILMO[J])-PSI(Z0[J],HI,ILMO[J]),\
                    ASX*ILMO[J]*ZU[J])            
            # *---------------------------------------------------------------------
            CT=KARMAN/FH
            CM=KARMAN/FM
            TZ[J]=TZ[J]+F[J]*(TG[J]-FTEMP[J]/(CT*UE[J])-GRAV/CPD*ZT[J])
            QZ[J]=QZ[J]+F[J]*(QG[J]-FVAP[J]/(CT*UE[J]))
            VITS=UE[J]/CM

            # * CALCULATE WIND DIRECTION CHANGE FROM TOP OF SURFACE LAYER
            DANG= (ZA[J]-ZU[J])*HI*ANGMAX*math.sin(LAT[J])
            ANGI=math.atan2(V[J],SIGN(abs(U[J])+1.e-05,U[J]))#SIGN should have 2 inputs!        
            if ILMO[J]>0.:#    THEN
                ANG=ANGI+DANG
            else:
                ANG=ANGI           

            UZ[J]=UZ[J]+F[J]*VITS*math.cos(ANG)
            VZ[J]=VZ[J]+F[J]*VITS*math.sin(ANG)
    return [UZ,VZ,TZ,QZ]

#
# -------- CLASSI FUNCTION ---------------
# Set up the surface conditions from the atmospheric conditions
# using modified machinery from CLASS 
'''
VapourPressureDeficit_In,DewPointTemp_In,PartialPressureOfDryAir_In,DensityOfAir_In,\
VPD                     ,TADP           ,PADRY                     ,RHOAIR
DensityOfFreshSnow_In,WetPrecip_In,WetPrecipTemp_In,FrozenPrecip_In,\
RHOSNI               ,RPCP        ,TRPCP           ,SPCP
FrozenPrecipTemp_In,AirTempAtRef_In,SpecificHumidityAtRef_In,SurfacePrecipRate_In,\
TSPCP              ,TA             ,QA                      ,PCPR
SurfaceAirPressure_In,IPCP,NumOfGridCells,1,\
PRESSG               ,IPCP,NL            ,IL1
NumOfGridCellsInRun,FreezingPointOfH2O,GasConstant,GasConstantOfH2OVapour,DensityOfH2O)
IL2                ,TFREZ             ,RGAS       ,RGASV                 ,RHOW

RRATE                ,SRATE 
'''
def CLASSI(VPD,TADP,PADRY,RHOAIR,RHOSNI,RPCP,TRPCP,\
        SPCP,TSPCP,TA,QA,PCPR,PRESSG,IPCP,\
        NL,IL1,IL2,TFREZ,RGAS,RGASV,RHOW):
    #RRATE,SRATE
    PHASE=np.zeros((len(range(IL1-1,IL2)),1),np.float32)
    for I in range(IL1-1,IL2):
        EA=QA[I]*PRESSG[I]/np.float32(0.622+0.378*QA[I])
        if TA[I]>=TFREZ:
            CA=np.float32(17.269)
            CB=np.float32(35.86)
        else:
            CA=np.float32(21.874)
            CB=np.float32(7.66)
        EASAT=np.float32(611.0*math.exp(CA*(TA[I]-TFREZ)/(TA[I]-CB)))
        VPD[I]=MAX(np.float32(0.0),(EASAT-EA)/np.float32(100.0))
        PADRY[I]=PRESSG[I]-EA
        RHOAIR[I]=PADRY[I]/(RGAS*TA[I])+EA/(RGASV*TA[I])
        CONST=np.float32(math.log(EA/611.0))
        TADP[I]=(CB*CONST-CA*TFREZ)/(CONST-CA)
        #
        #     * DENSITY OF FRESH SNOW.
        #
        if TA[I]<=TFREZ:
            RHOSNI[I]=np.float32(67.92+51.25*math.exp((TA[I]-TFREZ)/2.59))
        else:
            RHOSNI[I]=np.float32(MIN((119.17+20.0*(TA[I]-TFREZ)),200.0))
        
        #
        #     * PRECIPITATION PARTITIONING BETWEEN RAIN AND SNOW.
        #
        RPCP [I]=np.float32(0.0)
        TRPCP[I]=np.float32(0.0)
        SPCP [I]=np.float32(0.0)
        TSPCP[I]=np.float32(0.0)
        if PCPR>1.0E-8:
            if IPCP==1:
                if TA[I]>TFREZ:
                    RPCP [I]=PCPR[I]/RHOW
                    TRPCP[I]=MAX((TA[I]-TFREZ),np.float32(0.0))
                else:
                    SPCP [I]=PCPR[I]/RHOSNI[I]
                    TSPCP[I]=MIN((TA[I]-TFREZ),np.float32(0.0))            
            '''elif IPCP==2:
                if TA[I]<=TFREZ:
                    PHASE[I]=1.0
                elif TA[I]>=(TFREZ+2.0):
                    PHASE[I]=0.0
                else:
                    PHASE[I]=1.0-0.5*(TA[I]-TFREZ)                
                RPCP[I]=(1.0-PHASE[I])*PCPR[I]/RHOW
                if(RPCP[I]>0.0):
                    TRPCP[I]=MAX((TA[I]-TFREZ),0.0)                
                SPCP[I]=PHASE[I]*PCPR[I]/RHOSNI[I]
                if(SPCP[I]>0.0):
                    TSPCP[I]=MIN((TA[I]-TFREZ),0.0)                
            elif(IPCP==3):
                if TA[I]<=TFREZ:
                    PHASE[I]=1.0
                elif TA[I]>=(TFREZ+6.0):
                    PHASE[I]=0.0
                else:
                    PHASE[I]=(0.0202*(TA[I]-TFREZ)**6-0.3660*(TA[I]-TFREZ)**5+2.0399*(TA[I]-TFREZ)**4-1.5089*(TA[I]-TFREZ)**3-15.038*(TA[I]-TFREZ)**2+4.6664*(TA[I]-TFREZ)+100.0)/100.0
                    PHASE[I]=MAX(0.0,MIN(1.0,PHASE[I]))                
                RPCP[I]=(1.0-PHASE[I])*PCPR[I]/RHOW
                if RPCP[I]>0.0:
                    TRPCP[I]=MAX((TA[I]-TFREZ),0.0)                
                SPCP[I]=PHASE[I]*PCPR[I]/RHOSNI[I]
                if SPCP[I]>0.0:
                    TSPCP[I]=MIN((TA[I]-TFREZ),0.0)                
            elif IPCP==4:
                RPCP[I]=RRATE[I]/RHOW
                if RPCP[I]>0.0:
                    TRPCP[I]=MAX((TA[I]-TFREZ),0.0)
                    SPCP[I]=SRATE[I]/RHOSNI[I]                
                if SPCP[I]>0.0:
                    TSPCP[I]=MIN((TA[I]-TFREZ),0.0)'''
    return [VPD,TADP,PADRY,RHOAIR,RHOSNI,RPCP,TRPCP,SPCP,TSPCP,TA,QA,PCPR,\
        PRESSG,IPCP,NL,IL1,IL2]
    #RRATE,SRATE,

def EnsurePassByValue(a_in):
    '''Python passes numpy arrays and lists by refrence, which causes some issues 
    since this code was refactored from the MATLAB version of the CSLM (not the 
    original FORTRAN).  In the future this will be addressed as it is not 
    efficent to reallocate variabiles each time a fuction is called, but
    for now it ensures that variables in the parent are not being overwritten by 
    the children unless spacifically being returned.  
    -MG Clark
    '''
    try:
        [row,col]=a_in.shape   
        a_out = np.zeros((row,col),a_in.dtype)
        for I in range(0,row):
            for J in range(0,col):
                a_out[I,J] = a_in[I,J] 
    except:
        [row]=a_in.shape    
        col=0
        a_out = np.zeros((row),a_in.dtype)    
        for I in range(0,row):        
            a_out[I] = a_in[I]
    return a_out

RunLake()