ReadJuly2020Covid = False
ReadAugust2020Covid = False
ReadJan2021Covid = False
ReadApril2021Covid = False
ReadNov2021Covid = False
ReadMay2022Covid = False
Read7dayCovid = False
ScaleProperties = False
ConvertDynamicPredictedQuantity = False
ConvertDynamicProperties = True
GenerateFutures = False
GenerateSequences = False
PredictionsfromInputs = False
RereadMay2020 = False
UseOLDCovariates = False
Dropearlydata = 0
NIHCovariates = False
UseFutures = True
Usedaystart = False
PopulationNorm = False
SymbolicWindows = False
Hydrology = False
Earthquake = False

CDSpecial = False
RootCasesDeaths = True
NumpredbasicperTime = 2
NumpredFuturedperTime = 2
NumTimeSeriesCalculated = 0
Dailyunit = 1
TimeIntervalUnitName = 'Day'
InitialDate = datetime(2000, 1, 1)
NumberofTimeunits = 0
Num_Time = 0
FinalDate = datetime(2000, 1, 1)
GlobalTrainingLoss = 0.0
GlobalValidationLoss = 0.0

# Type of Testing
LocationBasedValidation = False
LocationValidationFraction = 0.0
LocationTrainingfraction = 1.0
RestartLocationBasedValidation = False

# Plotting
SeparateValandTrainingPlots = True
Plotsplitsize = -1  # if > 1 split time in plots
Plotrealnumbers = True
ListofTestFIPS = []
PlotsOnlyinTestFIPS = True
EarthquakeImagePlots = False
AddSpecialstoSummedplots = False
UseRealDatesonplots = False
Dumpoutkeyplotsaspics = False
OutputNetworkPictures = False
JournalSimplePrint = False
PlotinDL2F = False
FONTSIZE = 20

GarbageCollect = True
GarbageCollectionLimit = 5000000

PrintTitle('Start Dataset')

SubName = RunName[0:6]
if SubName == 'BEST14' or SubName == 'BEST15' or SubName == 'BEST16':
    UseOLDCovariates = False
    ReadAugust2020Covid = True
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    NIHCovariates = True
    ConvertDynamicProperties = True
    Dropearlydata = 37
    CDSpecial = True

if SubName == 'CovidA' or SubName == 'CovidN' or SubName == 'CovidM' or SubName == 'Covid7':
    UseOLDCovariates = False
    ReadApril2021Covid = True
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    UseFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    NIHCovariates = True
    ConvertDynamicProperties = True
    CDSpecial = True
    if SubName == 'CovidN':
        ReadNov2021Covid = True
    if SubName == 'CovidM':
        ReadMay2022Covid = True
    if SubName == 'Covid7':
        ReadMay2022Covid = True
        Read7dayCovid = True

if SubName == 'C2021A' or SubName == 'C2021B':
    UseOLDCovariates = False
    ReadJan2021Covid = True
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    NIHCovariates = True
    ConvertDynamicProperties = True
    Dropearlydata = 0
    CDSpecial = True

if SubName == 'Hydrol':
    Hydrology = True

if SubName == 'EARTHQ':
    Earthquake = True

if RunName == 'BEST10' or RunName == 'BEST13-10D' or RunName == 'BEST12-10' or RunName == 'BEST12-Test' or RunName == 'BEST13' or RunName == 'BEST13-10' or RunName == 'BEST13-10A' or RunName == 'BEST13-10C':
    UseOLDCovariates = False
    ReadAugust2020Covid = True
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    CDSpecial = True

if RunName == 'BEST11' or RunName == 'BEST11A':
    UseOLDCovariates = True
    ReadAugust2020Covid = True
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    CDSpecial = True

if RunName == 'BEST12':
    UseOLDCovariates = True
    RereadMay2020 = True
    ReadAugust2020Covid = False
    ScaleProperties = True
    ConvertDynamicPredictedQuantity = True
    GenerateFutures = True
    GenerateSequences = True
    PredictionsfromInputs = True
    CDSpecial = True

if RunName == 'BEST8' or RunName == 'BEST8A' or RunName == 'BEST12-LSTM-8':
    ReadJuly2020Covid = True