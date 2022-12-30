import numpy as np
import math
from hec.util import printexit
from hec.util import NaN
from datetime import datetime
from datetime import timedelta
from textwrap import wraptotext
from textwrap import wrap
from csv import reader
from cloudmesh.common.Shell import Shell
from datetime import timenow
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
import gc

from hec.util import startbold
from hec.util import startred
from hec.util import resetfonts

global ReadApril2021Covid
global CDSpecial
global ReadNov2021Covid
global ReadMay2022Covid
global Read7dayCovid
global ReadJan2021Covid
global ReadAugust2020Covid
global RereadMay2020

global ConvertDynamicPredictedQuantity

global APPLDIR
global RunName
global RunComment
global GarbageCollect

if ReadApril2021Covid:
    Dropearlydata = 40  # 3 more than needed by covariates so as to get "round number of days"
    if ReadNov2021Covid:
        Dropearlydata = 0
    if ReadMay2022Covid:
        Dropearlydata = 40  # XXX
    if Read7dayCovid:
        Dropearlydata = 0
    NIHCovariates = True
    UseOLDCovariates = False
    LengthFutures = 0

    if ReadNov2021Covid:
        # Set Dropearlydata to any number of days >= 0
        InitialDate = datetime(2020, 2, 29) + timedelta(days=Dropearlydata)
        FinalDate = datetime(2021, 11, 29)

        # Or set InitialData >= actual Date and FinalDate
        Dropearlydata = (InitialDate - datetime(2020, 2, 29)).days
        if Dropearlydata < 0:
            printexit('Illegal start date ' + str(Dropearlydata))
        NumberofTimeunits = (FinalDate - InitialDate).days + 1
        print("Total number of Days November 2021 Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(
            Dropearlydata))

        DATASETDIR = APPLDIR + '/CovidNovember2021'
        CasesFile = DATASETDIR + '/' + 'Cases.csv'
        DeathsFile = DATASETDIR + '/' + 'Deaths.csv'

    elif ReadMay2022Covid:
        if Read7dayCovid:
            InitialDate = datetime(2020, 2, 29) + timedelta(days=Dropearlydata)
            FinalDate = datetime(2022, 2, 28)
            FinalDate = datetime(2021, 11, 29)
            RootCasesDeaths = False
            NumberofTimeunits = (FinalDate - InitialDate).days + 1
            print(
                "Total number of Days 7 day Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(Dropearlydata))

            DATASETDIR = APPLDIR + '/Covid7Day/2021-11-29'
            CasesFile = DATASETDIR + '/' + 'Cases.csv'
            DeathsFile = DATASETDIR + '/' + 'Deaths.csv'

        else:  # May 2022 Data
            # Set Dropearlydata to any number of days >= 0
            InitialDate = datetime(2020, 1, 22) + timedelta(days=Dropearlydata)
            FinalDate = datetime(2022, 5, 15)  # In current code this is minmum final date between cases and deaths
            # Or set InitialData >= actual Initial Date and FinalDate
            # Initial Date must be >= datetime(2020,1,22) and final date must be <= datetime(2022,5,15)
            Dropearlydata = (InitialDate - datetime(2020, 1, 22)).days
            if Dropearlydata < 0:
                printexit('Illegal start date ' + str(Dropearlydata))
            NumberofTimeunits = (FinalDate - InitialDate).days + 1
            print("Total number of Days May 2022 Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(
                Dropearlydata))

            DATASETDIR = APPLDIR + '/CovidMay17-2022'
            CasesFile = DATASETDIR + '/' + 'Cases.csv'
            DeathsFile = DATASETDIR + '/' + 'Deaths.csv'

    else:  # April2021 Covid Data
        InitialDate = datetime(2020, 1, 22) + timedelta(days=Dropearlydata)
        FinalDate = datetime(2021, 4, 14)
        NumberofTimeunits = (FinalDate - InitialDate).days + 1
        print("Total number of Days April 2021 Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(
            Dropearlydata))

        DATASETDIR = APPLDIR + '/CovidApril14-2021'

        CasesFile = DATASETDIR + '/' + 'US_daily_cumulative_cases_April14.csv'
        DeathsFile = DATASETDIR + '/' + 'US_daily_cumulative_deaths_April14.csv'

    LocationdataFile = DATASETDIR + '/Population.csv'
    VotingdataFile = DATASETDIR + '/2020votes.csv'
    AlaskaVotingdataFile = DATASETDIR + '/Alaskavoting2016.csv'

    Nloc = 3142
    NFIPS = 3142

    # Set up location information
    Num_Time = NumberofTimeunits
    Locationfips = np.empty(NFIPS, dtype=int)  # integer version of FIPs
    Locationcolumns = []  # String version of FIPS
    FIPSintegerlookup = {}
    FIPSstringlookup = {}
    BasicInputTimeSeries = np.empty([Num_Time, Nloc, 2], dtype=np.float32)

    # Read in  cases Data into BasicInputTimeSeries
    with open(CasesFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS' and Ftype != 'casrn':
            printexit('EXIT: Wrong file type Cases ' + Ftype)
        hformat = '%Y-%m-%d'
        TargetDate = datetime.strptime(header[1], hformat)
        if (InitialDate - TargetDate).days != Dropearlydata:
            printexit('Incorrect cases initial date ' + str(Dropearlydata) + ' ' + str(TargetDate) + ' ' + str(InitialDate))

        iloc = 0
        for nextrow in csv_reader:
            if (len(nextrow) < NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Cases ' + str(iloc) + ' ' + str(len(nextrow)))
            # skip first entry
            localfips = nextrow[0]
            Locationcolumns.append(localfips)
            Locationfips[iloc] = int(localfips)
            FIPSintegerlookup[int(localfips)] = iloc
            FIPSstringlookup[localfips] = iloc
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 0] = nextrow[itime + 1 + Dropearlydata]
                if Dropearlydata > 0:
                    floatlast = np.float(nextrow[Dropearlydata])
                    BasicInputTimeSeries[itime, iloc, 0] = BasicInputTimeSeries[itime, iloc, 0] - floatlast
            iloc += 1
    # End Reading in cases data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths Cases ' + str(iloc) + ' ' + str(Nloc))
    print('Read Cases data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    # Read in deaths Data into BasicInputTimeSeries
    with open(DeathsFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS' and Ftype != 'casrn':
            printexit('EXIT: Wrong file type Deaths ' + Ftype)
        hformat = '%Y-%m-%d'
        TargetDate = datetime.strptime(header[1], hformat)
        if (InitialDate - TargetDate).days != Dropearlydata:
            printexit(
                'Incorrect deaths initial date ' + str(Dropearlydata) + ' ' + str(TargetDate) + ' ' + str(InitialDate))

        iloc = 0
        for nextrow in csv_reader:
            if (len(nextrow) < NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Deaths ' + str(iloc) + ' ' + str(len(nextrow)))
            localfips = nextrow[0]
            if (Locationfips[iloc] != int(localfips)):
                printexit('EXIT: Unexpected FIPS Deaths ' + localfips + ' ' + str(Locationfips[iloc]))
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 1] = np.float(nextrow[itime + 1 + Dropearlydata])
                if Dropearlydata > 0:
                    floatlast = np.float(nextrow[Dropearlydata])
                    BasicInputTimeSeries[itime, iloc, 1] = BasicInputTimeSeries[itime, iloc, 1] - floatlast
            iloc += 1
    # End Reading in deaths data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths ' + str(iloc) + ' ' + str(Nloc))
    print('Read Deaths data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    Locationname = ['Empty'] * NFIPS
    Locationstate = ['Empty'] * NFIPS
    Locationpopulation = np.empty(NFIPS, dtype=int)
    with open(LocationdataFile, 'r', encoding='latin1') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type Prop Data ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            localfips = int(nextrow[0])
            if localfips in FIPSintegerlookup.keys():
                jloc = FIPSintegerlookup[localfips]
                Locationname[jloc] = nextrow[4]
                Locationstate[jloc] = nextrow[3]
                Locationpopulation[jloc] = int(nextrow[2])
                iloc += 1  # just counting lines
            else:
                printexit('EXIT Inconsistent FIPS ' + str(iloc) + ' ' + str(localfips))
            # END setting NFIPS location properties

    DemVoting = np.full(NFIPS, -1.0, dtype=np.float32)
    RepVoting = np.full(NFIPS, -1.0, dtype=np.float32)
    with open(VotingdataFile, 'r', encoding='latin1') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'state_name':
            printexit('EXIT: Wrong file type Voting Data ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            localfips = int(nextrow[1])
            if localfips > 2900 and localfips < 2941:  # Alaska not useful
                continue
            if localfips in FIPSintegerlookup.keys():
                jloc = FIPSintegerlookup[localfips]
                if DemVoting[jloc] >= 0.0:
                    printexit('EXIT Double Setting of FIPS ' + str(iloc) + ' ' + str(localfips))
                DemVoting[jloc] = nextrow[8]
                RepVoting[jloc] = nextrow[7]
                iloc += 1  # just counting lines
            else:
                printexit('EXIT Inconsistent FIPS ' + str(iloc) + ' ' + str(localfips))

    with open(AlaskaVotingdataFile, 'r', encoding='utf-8-sig') as read_obj:  # remove ufeff
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'SpecialAlaska':
            printexit('EXIT: Wrong file type Alaska Voting Data ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            localfips = int(nextrow[1])
            if localfips in FIPSintegerlookup.keys():
                jloc = FIPSintegerlookup[localfips]
                if DemVoting[jloc] >= 0.0:
                    printexit('EXIT Double Setting of FIPS ' + str(iloc) + ' ' + str(localfips))
                DemVoting[jloc] = float(nextrow[2]) * 42.77 / 36.5
                RepVoting[jloc] = float(nextrow[3]) * 52.83 / 51.3
                iloc += 1  # just counting lines
            else:
                printexit('EXIT Inconsistent FIPS ' + str(iloc) + ' ' + str(localfips))

    for iloc in range(0, NFIPS):
        if DemVoting[iloc] >= 0.0:
            continue
        print(str(iloc) + ' Missing Votes ' + str(Locationfips[iloc]) + ' ' + Locationname[iloc] + ' ' + Locationstate[
            iloc] + ' pop ' + str(Locationpopulation[iloc]))
        DemVoting[iloc] = 0.5
        RepVoting[iloc] = 0.5

    # Set Static Properties of the Nloc studied locations
    # Order is Static, Dynamic, Cases, Deaths
    # Voting added as 13th covariate
    # Add fully vaccinated in November 2021
    NpropperTimeDynamic = 13
    if ReadNov2021Covid:
        NpropperTimeDynamic = 14
    if ReadMay2022Covid:
        NpropperTimeDynamic = 15
        if Read7dayCovid:
            NpropperTimeDynamic = 7
    NpropperTimeStatic = 0

    NpropperTime = NpropperTimeStatic + NpropperTimeDynamic + 2
    InputPropertyNames = [] * NpropperTime
    Property_is_Intensive = np.full(NpropperTime, True, dtype=np.bool)
    print('Initial Date ' + str(InitialDate) + ' Final Date ' + str(FinalDate) + ' NpropperTimeStatic ' + str(
        NpropperTimeStatic) + ' NpropperTimeDynamic ' + str(NpropperTimeDynamic))

"""### Read January 2021 Covid Data"""

if ReadJan2021Covid:
    Dropearlydata = 37
    NIHCovariates = True
    UseOLDCovariates = False

    InitialDate = datetime(2020, 1, 22) + timedelta(days=Dropearlydata)
    FinalDate = datetime(2021, 1, 26)
    NumberofTimeunits = (FinalDate - InitialDate).days + 1
    print(
        "Total number of Days January 2021 Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(Dropearlydata))

    DATASETDIR = APPLDIR + '/January2021'

    CasesFile = DATASETDIR + '/' + 'US_daily_cumulative_cases.csv'
    DeathsFile = DATASETDIR + '/' + 'US_daily_cumulative_deaths.csv'
    LocationdataFile = DATASETDIR + '/Population.csv'

    Nloc = 3142
    NFIPS = 3142

    # Set up location information
    Num_Time = NumberofTimeunits
    Locationfips = np.empty(NFIPS, dtype=int)  # integer version of FIPs
    Locationcolumns = []  # String version of FIPS
    FIPSintegerlookup = {}
    FIPSstringlookup = {}
    BasicInputTimeSeries = np.empty([Num_Time, Nloc, 2], dtype=np.float32)

    # Read in  cases Data into BasicInputTimeSeries
    with open(CasesFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type Cases ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            if (len(nextrow) < NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Cases ' + str(iloc) + ' ' + str(len(nextrow)))
            # skip first entry
            localfips = nextrow[0]
            Locationcolumns.append(localfips)
            Locationfips[iloc] = int(localfips)
            FIPSintegerlookup[int(localfips)] = iloc
            FIPSstringlookup[localfips] = iloc
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 0] = nextrow[itime + 1 + Dropearlydata]
                if Dropearlydata > 0:
                    floatlast = np.float(nextrow[Dropearlydata])
                    BasicInputTimeSeries[itime, iloc, 0] = BasicInputTimeSeries[itime, iloc, 0] - floatlast
            iloc += 1
    # End Reading in cases data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths Cases ' + str(iloc) + ' ' + str(Nloc))
    print('Read Cases data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    # Read in deaths Data into BasicInputTimeSeries
    with open(DeathsFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type Deaths ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            if (len(nextrow) < NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Deaths ' + str(iloc) + ' ' + str(len(nextrow)))
            localfips = nextrow[0]
            if (Locationfips[iloc] != int(localfips)):
                printexit('EXIT: Unexpected FIPS Deaths ' + localfips + ' ' + str(Locationfips[iloc]))
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 1] = nextrow[itime + 1 + Dropearlydata]
                if Dropearlydata > 0:
                    floatlast = np.float(nextrow[Dropearlydata])
                    BasicInputTimeSeries[itime, iloc, 1] = BasicInputTimeSeries[itime, iloc, 1] - floatlast
            iloc += 1
    # End Reading in deaths data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths ' + str(iloc) + ' ' + str(Nloc))
    print('Read Deaths data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    Locationname = ['Empty'] * NFIPS
    Locationstate = ['Empty'] * NFIPS
    Locationpopulation = np.empty(NFIPS, dtype=int)
    with open(LocationdataFile, 'r', encoding='latin1') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type Prop Data ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            localfips = int(nextrow[0])
            if localfips in FIPSintegerlookup.keys():
                jloc = FIPSintegerlookup[localfips]
                Locationname[jloc] = nextrow[4]
                Locationstate[jloc] = nextrow[3]
                Locationpopulation[jloc] = int(nextrow[2])
                iloc += 1  # just counting lines
            else:
                printexit('EXIT Inconsistent FIPS ' + str(iloc) + ' ' + str(localfips))
            # END setting NFIPS location properties

    # Set Static Properties of the Nloc studied locations
    # Order is Static, Dynamic, Cases, Deaths
    NpropperTimeDynamic = 12
    NpropperTimeStatic = 0

    NpropperTime = NpropperTimeStatic + NpropperTimeDynamic + 2
    InputPropertyNames = [' '] * NpropperTime
    Property_is_Intensive = np.full(NpropperTime, True, dtype=np.bool)

# Finish this after NIH Covariate

"""### Read Data defining COVID problem August 2020 Dataset"""

if ReadAugust2020Covid:
    InitialDate = datetime(2020, 1, 22) + timedelta(days=Dropearlydata)
    FinalDate = datetime(2020, 8, 13)
    NumberofTimeunits = (FinalDate - InitialDate).days + 1
    print("Total number of Days August Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(Dropearlydata))

    DATASETDIR = APPLDIR + '/MidAugust2020Data'

    CasesFile = DATASETDIR + '/' + 'covid-cases.csv'
    DeathsFile = DATASETDIR + '/' + 'covid-deaths.csv'
    CovariatesFile = DATASETDIR + '/' + 'PVI-31July2020.csv'
    if RereadMay2020 or UseOLDCovariates:
        CovariatesFile = DATASETDIR + '/' + 'Static_316USCities_Pop.csv'
    LocationdataFile = DATASETDIR + '/' + 'Static_316USCities_Pop.csv'

    Nloc = 314
    NFIPS = 316

if RereadMay2020:
    InitialDate = datetime(2020, 1, 22) + timedelta(days=Dropearlydata)
    FinalDate = datetime(2020, 5, 25)
    NumberofTimeunits = (FinalDate - InitialDate).days + 1
    print("Total number of Days May Dataset " + str(NumberofTimeunits) + ' Dropping at start ' + str(Dropearlydata))

    DATASETDIR = APPLDIR + '/EndMay2020fromfiles'

    CasesFile = DATASETDIR + '/' + 'Covid19-cases-110USCities.csv'
    DeathsFile = DATASETDIR + '/' + 'Covid19-deaths-110USCities.csv'
    CovariatesFile = DATASETDIR + '/' + 'PVI-31July2020.csv'
    if UseOLDCovariates:
        CovariatesFile = DATASETDIR + '/' + 'Static_316USCities_Pop.csv'
    LocationdataFile = DATASETDIR + '/' + 'Static_316USCities_Pop.csv'

    Nloc = 110
    NFIPS = 112

if ReadAugust2020Covid or RereadMay2020:

    # Set up location information
    Num_Time = NumberofTimeunits
    Locationfips = np.empty(NFIPS, dtype=int)  # integer version of FIPs
    Locationcolumns = []  # String version of FIPS
    FIPSintegerlookup = {}
    FIPSstringlookup = {}
    BasicInputTimeSeries = np.empty([Num_Time, Nloc, 2], dtype=np.float32)

    # Read in  cases Data into BasicInputTimeSeries
    with open(CasesFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type ' + Ftype)

        iloc = 0

        for nextrow in csv_reader:
            if (len(nextrow) != NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Cases ' + str(iloc) + ' ' + str(len(nextrow)))
            localfips = nextrow[0]
            Locationcolumns.append(localfips)
            Locationfips[iloc] = int(localfips)
            FIPSintegerlookup[int(localfips)] = iloc
            FIPSstringlookup[localfips] = iloc
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 0] = nextrow[itime + 1 + Dropearlydata]
            iloc += 1
    # End Reading in cases data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths Cases ' + str(iloc) + ' ' + str(Nloc))
    print('Read Cases data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    # Read in deaths Data into BasicInputTimeSeries
    with open(DeathsFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            if (len(nextrow) != NumberofTimeunits + 1 + Dropearlydata):
                printexit('EXIT: Incorrect row length Deaths ' + str(iloc) + ' ' + str(len(nextrow)))
            localfips = nextrow[0]
            if (Locationfips[iloc] != int(localfips)):
                printexit('EXIT: Unexpected FIPS Deaths ' + localfips + ' ' + str(Locationfips[iloc]))
            for itime in range(0, NumberofTimeunits):
                BasicInputTimeSeries[itime, iloc, 1] = nextrow[itime + 1 + Dropearlydata]
            iloc += 1
    # End Reading in deaths data

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths ' + str(iloc) + ' ' + str(Nloc))
    print('Read Deaths data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))

    # START setting location properties -- there are NFIPS of these
    # NFIPS can be larger than Nloc. Any fips in studied group must be in fips property group
    # Add missing FIPS in 315 and not 314 set are 49057  and 49053
    # while 48203 is in 314 but not 315; 316 adds 48023
    Locationfips[Nloc] = 49057
    Locationfips[Nloc + 1] = 49053
    Locationcolumns.append(str(Locationfips[Nloc]))
    FIPSintegerlookup[Locationfips[Nloc]] = Nloc
    FIPSstringlookup[str(Locationfips[Nloc])] = Nloc
    Locationcolumns.append(str(Locationfips[Nloc + 1]))
    FIPSintegerlookup[Locationfips[Nloc + 1]] = Nloc + 1
    FIPSstringlookup[str(Locationfips[Nloc + 1])] = Nloc + 1

    Locationname = ['Empty'] * NFIPS
    Locationstate = ['Empty'] * NFIPS
    Locationpopulation = np.empty(NFIPS, dtype=int)
    with open(LocationdataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'FIPS':
            printexit('EXIT: Wrong file type ' + Ftype)

        iloc = 0
        for nextrow in csv_reader:
            localfips = int(nextrow[0])
            if localfips in FIPSintegerlookup.keys():
                jloc = FIPSintegerlookup[localfips]
                Locationname[jloc] = nextrow[2]
                Locationstate[jloc] = nextrow[1]
                Locationpopulation[jloc] = int(nextrow[5])
                iloc += 1  # just counting lines

    if iloc != Nloc + 2:
        printexit('EXIT Inconsistent old static data lengths ' + str(iloc) + ' ' + str(Nloc + 2))
    if 48203 in FIPSintegerlookup.keys():
        iloc = FIPSintegerlookup[48203]
        Locationname[iloc] = 'Harrison'
        Locationstate[iloc] = 'Texas'
        Locationpopulation[iloc] = 66553
    # END setting NFIPS location properties

    # Set Static Properties of the Nloc studied locations
    # Order is Static, Dynamic, Cases, Deaths
    if NIHCovariates:
        NpropperTimeDynamic = 11
        NpropperTimeStatic = 0
    else:
        NpropperTimeDynamic = 0
        NpropperTimeStatic = 12
        if UseOLDCovariates:
            NpropperTimeStatic = 26
    NpropperTime = NpropperTimeStatic + NpropperTimeDynamic + 2
    InputPropertyNames = [] * NpropperTime
    Property_is_Intensive = np.full(NpropperTime, True, dtype=np.bool)

    if not NIHCovariates:
        BasicInputStaticProps = np.empty([Nloc, NpropperTimeStatic], dtype=np.float32)

        with open(CovariatesFile, 'r') as read_obj:
            csv_reader = reader(read_obj)
            header = next(csv_reader)
            Ftype = header[0]
            if Ftype != 'FIPS':
                printexit('EXIT: Wrong file type ' + Ftype)
            throwaway = 2
            if UseOLDCovariates:
                throwaway = 6
            if (len(header) != (throwaway + NpropperTimeStatic)):
                printexit('EXIT: Incorrect property header length ' + str(len(header)) + ' ' + str(2 + NpropperTimeStatic))
            InputPropertyNames[:] = header[throwaway:]

            iloc = 0
            for nextrow in csv_reader:
                if (len(nextrow) != (throwaway + NpropperTimeStatic)):
                    printexit('EXIT: Incorrect row length ' + str(iloc) + ' ' + str(2 + NpropperTimeStatic) + ' ' + str(
                        len(nextrow)))
                localfips = int(nextrow[0])
                if not localfips in FIPSintegerlookup.keys():
                    continue
                #           printexit('EXIT: Missing FIPS ' + str(localfips))
                jloc = FIPSintegerlookup[localfips]
                if jloc >= Nloc:
                    print('FIPS ' + str(localfips) + ' skipped in property read')
                    continue  # skip this FIPS
                BasicInputStaticProps[jloc, :] = np.asarray(nextrow[throwaway:], dtype=np.float32)
                iloc += 1
        # End Reading in Static Properties data

        if iloc != Nloc:
            printexit('EXIT Inconsistent location lengths ' + str(iloc) + ' ' + str(Nloc))
        print('Read Static Properties for locations ' + str(Nloc) + ' Properties ' + str(NpropperTimeStatic))

        # August Covariates all intensive and no missing data
        # May Coviates have intensive properties missing for Harrison TX
        if UseOLDCovariates:
            Property_is_Intensive[20] = False
            Property_is_Intensive[21] = False
            Property_is_Intensive[22] = False

# Finish this after NIH Covariate

"""###  Clean up Extensive and Undefined Properties

### Read and setup NIH Covariates August 2020 and January, April 2021 Data

new collection of time dependent covariates (even if constant).

cases and deaths and location property from previous data
"""

import re

if NIHCovariates:
    if ReadJan2021Covid:
        Propfilenames = ["Age Distribution.csv", "Air Pollution.csv", "Comorbidities.csv", "Demographics.csv",
                         "Disease Spread.csv",
                         "Health Disparities.csv", "Hospital Beds.csv", "Intervention Testing.csv", "Mobility.csv",
                         "Residential Density.csv", "Social Distancing.csv", "Transmissible Cases.csv"]
        Propnames = ["Age Distribution", "Air Pollution", "Co-morbidities", "Demographics", "Disease Spread",
                     "Health Disparities", "Hospital Beds", "Intervention Testing", "Mobility", "Residential Density",
                     "Social Distancing", "Transmissible Cases"]

    elif ReadApril2021Covid:
        if ReadMay2022Covid:
            if Read7dayCovid:
                Propfilenames = ["Age Distribution.csv", "Disease Spread.csv",
                                 "Health Disparities.csv",
                                 "Social Distancing.csv", "Transmissible Cases.csv", "Vaccination.csv", "NOFILE"]
                Propnames = ["Age Distribution", "Disease Spread",
                             "Health Disparities",
                             "Social Distancing", "Transmissible Cases", "Vaccination", "voting"]
            else:
                Propfilenames = ["Age Distribution.csv", "Air Pollution.csv", "Co-morbidities.csv", "Pop Demographics.csv",
                                 "Disease Spread.csv",
                                 "Health Disparities.csv", "Hospital Beds.csv", "Pop Mobility.csv",
                                 "Residential Density.csv", "Social Distancing.csv", "Testing.csv",
                                 "Transmissible Cases.csv", "Vaccination.csv", "VaccinationOneDose.csv", "NOFILE"]
                Propnames = ["Age Distribution", "Air Pollution", "Co-morbidities", "Pop Demographics", "Disease Spread",
                             "Health Disparities", "Hospital Beds", "Pop Mobility", "Residential Density",
                             "Social Distancing", "Testing", "Transmissible Cases", "Vaccination", "VaccinationOneDose",
                             "voting"]
        else:
            Propfilenames = ["Age Distribution.csv", "Air Pollution.csv", "Comorbidities.csv", "Demographics.csv",
                             "Disease Spread.csv",
                             "Health Disparities.csv", "Hospital Beds.csv", "Mobility.csv",
                             "Residential Density.csv", "Social Distancing.csv", "Testing.csv", "Transmissible Cases.csv",
                             "NOFILE"]
            Propnames = ["Age Distribution", "Air Pollution", "Co-morbidities", "Demographics", "Disease Spread",
                         "Health Disparities", "Hospital Beds", "Mobility", "Residential Density",
                         "Social Distancing", "Testing", "Transmissible Cases", "voting"]
            if ReadNov2021Covid:
                Propfilenames.append("FullyVaccinated.csv")
                Propnames.append("Fully Vaccinated")
    else:
        Propfilenames = ["Age Distribution.csv", "Air Pollution.csv", "Co-morbidities.csv", "Health Disparities.csv",
                         "Hospital Beds.csv", "Pop Demographics.csv", "Pop Mobility.csv", "Residential Density.csv",
                         "Social Distancing.csv", "Testing.csv", "Transmissible Cases.csv"]
        Propnames = ["Age Distribution", "Air Pollution", "Co-morbidities", "Health Disparities", "Hospital Beds",
                     "Pop Demographics", "Pop Mobility", "Residential Density", "Social Distancing", "Testing",
                     "Transmissible Cases"]

    NIHDATADIR = DATASETDIR + '/'
    numberfiles = len(Propnames)
    NpropperTimeStatic = 0
    if NpropperTimeDynamic != numberfiles:
        printexit('EXIT: Dynamic Properties set wrong ' + str(numberfiles) + ' ' + str(NpropperTimeDynamic))
    DynamicPropertyTimeSeries = np.zeros([Num_Time, Nloc, numberfiles], dtype=np.float32)
    enddifference = NaN

    for ifiles in range(0, numberfiles):
        InputPropertyNames.append(Propnames[ifiles])
        if Propfilenames[ifiles] == 'NOFILE':  # Special case of Voting Data
            for iloc in range(0, Nloc):
                Demsize = DemVoting[iloc]
                RepSize = RepVoting[iloc]
                Votingcovariate = Demsize / (RepSize + Demsize)
                DynamicPropertyTimeSeries[:, iloc, ifiles] = Votingcovariate
            continue  # over ifile loop

        DynamicPropFile = NIHDATADIR + Propfilenames[ifiles]
        if not (ReadJan2021Covid or ReadApril2021Covid):
            DynamicPropFile = DATASETDIR + '/ThirdCovariates/' + Propfilenames[ifiles]

        # Read in  Covariate Data into DynamicPropertyTimeSeries
        with open(DynamicPropFile, 'r') as read_obj:
            csv_reader = reader(read_obj)
            header = next(csv_reader)
            skip = 1
            if ReadNov2021Covid:
                skip = 2
                Ftype = header[1]
                if Ftype != 'Name':
                    printexit('EXIT: Wrong file type ' + Ftype)
                Ftype = header[0]
                if Ftype != 'FIPS':
                    printexit('EXIT: Wrong file type ' + Ftype)
            elif ReadMay2022Covid:
                if Read7dayCovid:
                    skip = 2
                    Ftype = header[0]
                    if Ftype != 'Name':
                        printexit('EXIT: Wrong file type ' + Ftype)
                    Ftype = header[1]
                    if Ftype != 'FIPS':
                        printexit('EXIT: Wrong file type ' + Ftype)
                else:
                    if (Propfilenames[ifiles] == "Vaccination.csv") or (Propfilenames[ifiles] == "VaccinationOneDose.csv"):
                        skip = 1
                        Ftype = header[0]
                        if Ftype != 'FIPS':
                            printexit('EXIT: Wrong file type ' + Ftype)
                    else:
                        skip = 2
                        Ftype = header[1]
                        if Ftype != 'Name':
                            printexit('EXIT: Wrong file type ' + Ftype)
                        Ftype = header[0]
                        if Ftype != 'FIPS':
                            printexit('EXIT: Wrong file type ' + Ftype)
            else:
                if (ReadJan2021Covid or ReadApril2021Covid):
                    skip = 2
                    Ftype = header[0]
                    if Ftype != 'Name':
                        printexit('EXIT: Wrong file type ' + Ftype)
                Ftype = header[skip - 1]
                if Ftype != 'FIPS':
                    printexit('EXIT: Wrong file type ' + Ftype)

            # Check Date
            hformat = '%m-%d-%Y'
            if ReadJan2021Covid or ReadApril2021Covid:
                hformat = '%Y-%m-%d'
            if Propfilenames[ifiles] == "FullyVaccinated.csv" and ReadNov2021Covid:
                hformat = '%m/%d/%Y'
            print(str(ifiles) + " " + str(skip) + " " + header[0] + " " + header[1] + ' ' + Propnames[ifiles])
            stringdate = header[skip]
            stringdatelist = re.findall('(.*) .*', stringdate)
            if stringdatelist:
                stringdate = stringdatelist[0]
            firstdate = datetime.strptime(stringdate, hformat)
            tdelta = (firstdate - InitialDate).days
            if tdelta > 0:
                print(Propnames[ifiles] + ' Missing Covariate Data at start ' + str(tdelta))
            stringdate = header[len(header) - 1]
            stringdatelist = re.findall('(.*) .*', stringdate)
            if stringdatelist:
                stringdate = stringdatelist[0]
            lastdate = datetime.strptime(stringdate, hformat)
            enddifference1 = (FinalDate - lastdate).days
            if math.isnan(enddifference):
                enddifference = enddifference1
                print(Propnames[ifiles] + ' Missing days at the end ' + str(enddifference))
            else:
                if enddifference != enddifference1:
                    print('Change in time length at end ' + Propnames[ifiles] + ' expected ' + str(
                        enddifference) + ' actual ' + str(enddifference1))
            iloc = 0

            for nextrow in csv_reader:
                if (len(nextrow) != NumberofTimeunits + skip - enddifference1 - tdelta):
                    printexit('EXIT: Incorrect row length ' + Propnames[ifiles] + ' Location ' + str(iloc) + ' ' + str(
                        len(nextrow)))
                if ReadNov2021Covid or ReadMay2022Covid:
                    localfips = nextrow[0]
                    if Read7dayCovid:
                        localfips = nextrow[1]
                else:
                    localfips = nextrow[skip - 1]
                intversion = int(localfips)
                if intversion > 56045:
                    continue
                jloc = FIPSstringlookup[localfips]
                FinalTimeIndex = min(NumberofTimeunits - enddifference1, NumberofTimeunits)
                FirstTimeIndex = max(tdelta, 0)
                for itime in range(FirstTimeIndex, FinalTimeIndex):
                    DynamicPropertyTimeSeries[itime, jloc, ifiles] = nextrow[itime + skip - tdelta]
                # Use previous week value for missing data at the end
                for itime in range(FinalTimeIndex, NumberofTimeunits):
                    DynamicPropertyTimeSeries[itime, jloc, ifiles] = DynamicPropertyTimeSeries[itime - 7, jloc, ifiles]
                iloc += 1
        # End Reading in dynamic property data

        if iloc != Nloc:
            printexit('EXIT Inconsistent location lengths ' + Propnames[ifiles] + str(iloc) + ' ' + str(Nloc))
        if tdelta <= 0:
            print('Read ' + Propnames[ifiles] + ' data for locations ' + str(Nloc) + ' Time Steps ' + str(
                Num_Time) + ' Days dropped at start ' + str(-tdelta))
        else:
            print('Read ' + Propnames[ifiles] + ' data for locations ' + str(Nloc) + ' Time Steps ' + str(
                Num_Time) + ' zero value Days added at start ' + str(tdelta))

    if ReadApril2021Covid:
        CovidPopulationCut = 0  # Use this if NumberCut = 0
        NumberCut = 2642
        if Read7dayCovid:
            NumberCut = 0
        uselocation = np.full(Nloc, True, dtype=np.bool)
        if (CovidPopulationCut > 0) or (NumberCut > 0):
            if NumberCut > 0:
                smalllocations = np.argsort(Locationpopulation)
                for jloc in range(0, NumberCut):
                    uselocation[smalllocations[jloc]] = False
                CovidPopulationCut = Locationpopulation[smalllocations[NumberCut]]
            else:
                NumberCut = 0
                for iloc in range(0, Nloc):
                    if Locationpopulation[iloc] < CovidPopulationCut:
                        uselocation[iloc] = False
                        NumberCut += 1
            print(' Population Cut ' + str(CovidPopulationCut) + ' removes ' + str(NumberCut) + ' of ' + str(Nloc))
        if (NumberCut > 0):
            NewNloc = Nloc - NumberCut
            NewNFIPS = NewNloc
            NewLocationfips = np.empty(NewNFIPS, dtype=int)  # integer version of FIPs
            NewLocationcolumns = []  # String version of FIPS
            NewFIPSintegerlookup = {}
            NewFIPSstringlookup = {}
            NewBasicInputTimeSeries = np.empty([Num_Time, NewNloc, 2], dtype=np.float32)
            NewLocationname = ['Empty'] * NewNFIPS
            NewLocationstate = ['Empty'] * NewNFIPS
            NewLocationpopulation = np.empty(NewNFIPS, dtype=int)
            NewDynamicPropertyTimeSeries = np.empty([Num_Time, NewNloc, numberfiles], dtype=np.float32)

            Newiloc = 0
            for iloc in range(0, Nloc):
                if not uselocation[iloc]:
                    continue
                NewBasicInputTimeSeries[:, Newiloc, :] = BasicInputTimeSeries[:, iloc, :]
                NewDynamicPropertyTimeSeries[:, Newiloc, :] = DynamicPropertyTimeSeries[:, iloc, :]
                localfips = Locationcolumns[iloc]
                NewLocationcolumns.append(localfips)
                NewLocationfips[Newiloc] = int(localfips)
                NewFIPSintegerlookup[int(localfips)] = Newiloc
                NewFIPSstringlookup[localfips] = Newiloc
                NewLocationpopulation[Newiloc] = Locationpopulation[iloc]
                NewLocationstate[Newiloc] = Locationstate[iloc]
                NewLocationname[Newiloc] = Locationname[iloc]
                Newiloc += 1

            BasicInputTimeSeries = NewBasicInputTimeSeries
            DynamicPropertyTimeSeries = NewDynamicPropertyTimeSeries
            Locationname = NewLocationname
            Locationstate = NewLocationstate
            Locationpopulation = NewLocationpopulation
            FIPSstringlookup = NewFIPSstringlookup
            FIPSintegerlookup = NewFIPSintegerlookup
            Locationcolumns = NewLocationcolumns
            Locationfips = NewLocationfips
            NFIPS = NewNFIPS
            Nloc = NewNloc

"""## Process Input Data  in various ways

### Convert Cumulative to Daily
"""

# Convert  cumulative to Daily.
# Replace negative daily values by zero
# remove daily to sqrt(daily)  and Then normalize maximum to 1
if ConvertDynamicPredictedQuantity:
    NewBasicInputTimeSeries = np.empty_like(BasicInputTimeSeries, dtype=np.float32)
    Zeroversion = np.zeros_like(BasicInputTimeSeries, dtype=np.float32)
    Rolleddata = np.roll(BasicInputTimeSeries, 1, axis=0)
    Rolleddata[0, :, :] = Zeroversion[0, :, :]
    NewBasicInputTimeSeries = np.maximum(np.subtract(BasicInputTimeSeries, Rolleddata), Zeroversion)
    originalnumber = np.sum(BasicInputTimeSeries[NumberofTimeunits - 1, :, :], axis=0)
    newnumber = np.sum(NewBasicInputTimeSeries, axis=(0, 1))
    print('Original summed counts ' + str(originalnumber) + ' become ' + str(newnumber) + ' Cases, Deaths')

    BasicInputTimeSeries = NewBasicInputTimeSeries

"""### Static and Dynamic specials for COVID

except case where Romeo data read
"""

# Remove special status of Cases and Deaths
if CDSpecial:

    NewNpropperTimeDynamic = NpropperTimeDynamic + 2
    NewNpropperTime = NpropperTimeStatic + NewNpropperTimeDynamic

    NewProperty_is_Intensive = np.full(NewNpropperTime, True, dtype=np.bool)
    NewInputPropertyNames = []
    NewDynamicPropertyTimeSeries = np.empty([Num_Time, Nloc, NewNpropperTimeDynamic], dtype=np.float32)

    for casesdeaths in range(0, 2):
        NewDynamicPropertyTimeSeries[:, :, casesdeaths] = BasicInputTimeSeries[:, :, casesdeaths]
    BasicInputTimeSeries = None

    for iprop in range(0, NpropperTimeStatic):
        NewInputPropertyNames.append(InputPropertyNames[iprop])
        NewProperty_is_Intensive[iprop] = Property_is_Intensive[iprop]
    NewProperty_is_Intensive[NpropperTimeStatic] = False
    NewProperty_is_Intensive[NpropperTimeStatic + 1] = False
    NewInputPropertyNames.append('Cases')
    NewInputPropertyNames.append('Deaths')
    for ipropdynamic in range(0, NpropperTimeDynamic):
        Newiprop = NpropperTimeStatic + 2 + ipropdynamic
        iprop = NpropperTimeStatic + ipropdynamic
        NewDynamicPropertyTimeSeries[:, :, Newiprop] = DynamicPropertyTimeSeries[:, :, iprop]
        NewInputPropertyNames.append(InputPropertyNames[iprop])
        NewProperty_is_Intensive[Newiprop] = Property_is_Intensive[iprop]

    NpropperTimeDynamic = NewNpropperTimeDynamic
    NpropperTime = NewNpropperTime
    DynamicPropertyTimeSeries = NewDynamicPropertyTimeSeries
    InputPropertyNames = NewInputPropertyNames
    Property_is_Intensive = NewProperty_is_Intensive

"""### Static Property Manipulations for Covid Case"""

# Execute under all COVID circumstances properties generated here
if CDSpecial:
    if NpropperTimeStatic > 0:
        Num_Extensive = 0
        for iprop in range(0, NpropperTimeStatic):
            if not Property_is_Intensive[iprop]:
                Num_Extensive += 1
        print(startbold + startred + ' Number of Extensive parameters ' + str(Num_Extensive) + resetfonts)
        for iprop in range(0, NpropperTimeStatic):
            if not Property_is_Intensive[iprop]:
                print(InputPropertyNames[iprop])

        # Convert Extensive covariates to SQRT(Population normed)
        # Replace negatives by mean of positives and zeroes
        positivemean = np.zeros(NpropperTimeStatic, dtype=np.float32)
        countvalidentries = np.zeros(NpropperTimeStatic, dtype=np.float32)
        for iloc in range(0, Nloc):
            for iprop in range(0, NpropperTimeStatic):
                if not Property_is_Intensive[iprop]:
                    BasicInputStaticProps[iloc, iprop] = np.sqrt(
                        BasicInputStaticProps[iloc, iprop] / Locationpopulation[iloc])
                else:
                    if BasicInputStaticProps[iloc, iprop] >= 0:
                        positivemean[iprop] += BasicInputStaticProps[iloc, iprop]
                        countvalidentries[iprop] += 1.0

        for iprop in range(0, NpropperTimeStatic):
            if Property_is_Intensive[iprop]:
                positivemean[iprop] /= countvalidentries[iprop]

        for iloc in range(0, Nloc):
            for iprop in range(0, NpropperTimeStatic):
                if Property_is_Intensive[iprop]:
                    if BasicInputStaticProps[iloc, iprop] < 0:
                        BasicInputStaticProps[iloc, iprop] = positivemean[iprop]