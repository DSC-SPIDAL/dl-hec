from cloudmesh.common.Shell import Shell

def ReadEarthquakeData():
    global BasicInputTimeSeries, BasicInputStaticProps, CalculatedTimeSeries, InputPropertyNames, NamespredCalculated, DynamicNames

    global Locationname, Locationstate, Locationpopulation, Locationfips, Locationcolumns, FIPSintegerlookup, FIPSstringlookup

    global Nloc, NFIPS, NpropperTime, NpropperTimeDynamicCalculated, NpropperTimeDynamic, NpropperTimeDynamicInput, NpropperTimeStatic, NumTimeSeriesCalculatedBasic

    global read1950, UseEarthquakeEigenSystems, RundleEigenvectors, MagnitudeMethod, \
        numberspecialeqs, Specialuse, Specialmags, Specialdepth, Speciallong, Speciallat, Specialdate, Specialxpos, Specialypos, Specialeqname, \
        addRundleEMA, RundleLambda, RundleSteps, InputIndextogenerateEMA, FirstEMAIndex, Plottingdelay, OriginalNloc, MapLocation

    global ReadJuly2020Covid, ReadAugust2020Covid, ReadJan2021Covid, ReadApril2021Covid, ReadNov2021Covid, ReadMay2022Covid, Read7dayCovid, \
        ScaleProperties, ConvertDynamicPredictedQuantity, ConvertDynamicProperties, GenerateFutures, GenerateSequences, PredictionsfromInputs, \
        RereadMay2020, UseOLDCovariates, Dropearlydata, NIHCovariates, UseFutures, Usedaystart, PopulationNorm, SymbolicWindows, Hydrology, Earthquake, \
        CDSpecial, RootCasesDeaths, NumpredbasicperTime, NumpredFuturedperTime, NumTimeSeriesCalculated, Dailyunit, TimeIntervalUnitName, InitialDate, \
        NumberofTimeunits, Num_Time, FinalDate, GlobalTrainingLoss, GlobalValidationLoss, LocationBasedValidation, LocationValidationFraction, \
        LocationTrainingfraction, RestartLocationBasedValidation, SeparateValandTrainingPlots, Plotsplitsize, Plotrealnumbers, ListofTestFIPS, \
        PlotsOnlyinTestFIPS, EarthquakeImagePlots, AddSpecialstoSummedplots, UseRealDatesonplots, Dumpoutkeyplotsaspics, OutputNetworkPictures, \
        JournalSimplePrint, PlotinDL2F, FONTSIZE, GarbageCollect, GarbageCollectionLimit
    global NaN, PLOTNUMBER, COLABROOTDIR, APPLDIR, Directoryaddon, CHECKPOINTDIR, CovidofSomeType

    read1950 = True
    RundleEigenvectors = 2
    UseEarthquakeEigenSystems = False
    Dailyunit = 14

    content = Shell.run("ls '/content/gdrive/My Drive/Colab Datasets/EarthquakeDec2020'")
    print(content)

    if read1950:
        MagnitudeDataFile = APPLDIR + '/1950start/SC_1950-2019.freq-D-25567x2400-log_eng.multi.csv'
        DepthDataFile = APPLDIR + '/1950start/SC_1950-2019.freq-D-25567x2400-w_depth.multi.csv'
        MultiplicityDataFile = APPLDIR + '/1950start/SC_1950-2019.freq-D-25567x2400-n_shock.multi.csv'
        RundleMultiplicityDataFile = APPLDIR + '/1950start/SC_1950-2019.freq-D-25567x2400-n_shock-mag-3.29.multi.csv'
        NumberofTimeunits = 25567
        InitialDate = datetime(1950, 1, 1)
    else:
        MagnitudeDataFile = APPLDIR + '/SC_1990-2019.freq-D-10759x2400.csv'
        DepthDataFile = APPLDIR + '/SC_1990-2019.freq-D-w_depth-10759x2400.multi.csv'
        MultiplicityDataFile = APPLDIR + '/SC_1990-2019.freq-D-num_evts-10759x2400.csv'
        RundleMultiplicityDataFile = APPLDIR + '/SC_1990-2019.freq-D-10755x2400-n_shock-mag-3.29.multi.csv'
        NumberofTimeunits = 10759
        InitialDate = datetime(1990, 1, 1)

    Topearthquakesfile = APPLDIR + '/topearthquakes_20.csv'
    FaultLabelDataFile = APPLDIR + '/pix_faults_SmallJan21.csv'
    MagnitudeMethod = 0
    ReadFaultMethod = 2  # one set of x values for each input row
    Numberxpixels = 60
    Numberypixels = 40
    Numberpixels = Numberxpixels * Numberypixels
    Nloc = Numberpixels
    Nlocdimension = 2
    Nlocaxislengths = np.array((Numberxpixels, Numberypixels), ndmin=1, dtype=int)  # First row is top (north)
    vertices = cal_gilbert2d(Numberxpixels, Numberypixels)
    #    print(vertices[0], vertices[1],vertices[2399], vertices[1198], vertices[1199],vertices[1200], vertices[1201])
    sfcurvelist = vertices
    plot_gilbert2d_space_filling(sfcurvelist, Numberxpixels, Numberypixels)

    Dropearlydata = 0

    FinalDate = InitialDate + timedelta(days=NumberofTimeunits - 1)
    print(startbold + startred + InitialDate.strftime("%d/%m/%Y") + ' To ' + FinalDate.strftime("%d/%m/%Y")
          + ' days ' + str(NumberofTimeunits) + resetfonts)
    print(
        ' Pixels ' + str(Nloc) + ' x dimension ' + str(Nlocaxislengths[0]) + ' y dimension ' + str(Nlocaxislengths[1]))

    # Set up location information
    Num_Time = NumberofTimeunits
    NFIPS = Numberpixels
    Locationname = [''] * NFIPS
    Locationstate = [' '] * NFIPS
    Locationpopulation = np.ones(NFIPS, dtype=int)
    Locationfips = np.empty(NFIPS, dtype=int)  # integer version of FIPs
    Locationcolumns = []  # String version of FIPS
    FIPSintegerlookup = {}
    FIPSstringlookup = {}
    for iloc in range(0, Numberpixels):
        localfips = iloc
        xvalue = localfips % Nlocaxislengths[0]
        yvalue = np.floor(localfips / Nlocaxislengths[0])
        Stringfips = str(xvalue) + ',' + str(yvalue)
        Locationcolumns.append(Stringfips)
        Locationname[iloc] = Stringfips
        Locationfips[iloc] = localfips
        FIPSintegerlookup[localfips] = localfips
        FIPSstringlookup[Stringfips] = localfips

    # TimeSeries 0 magnitude 1 depth 2 Multiplicity 3 Rundle Multiplicity
    NpropperTimeDynamicInput = 4
    BasicInputTimeSeries = np.empty([Num_Time, Nloc, NpropperTimeDynamicInput], dtype=np.float32)
    # StaticProps 0...NumFaultLabels-1 Fault Labels
    NumFaultLabels = 4
    BasicInputStaticProps = np.empty([Nloc, NumFaultLabels], dtype=np.float32)
    RawFaultData = np.empty(Nloc, dtype=np.int)

    # Read in Magnitude Data into BasicInputTimeSeries
    with open(MagnitudeDataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != '':
            printexit('EXIT: Wrong header on line 1 ' + Ftype + ' of ' + MagnitudeDataFile)

        itime = 0
        for nextrow in csv_reader:
            if (len(nextrow) != Numberpixels + 1):
                printexit('EXIT: Incorrect row length Magnitude ' + str(itime) + ' ' + str(len(nextrow)))
            localtime = nextrow[0]
            if (itime != int(localtime)):
                printexit('EXIT: Unexpected Time in Magnitude ' + localtime + ' ' + str(itime))
            for iloc in range(0, Numberpixels):
                BasicInputTimeSeries[itime, iloc, 0] = TransformMagnitude(float(nextrow[iloc + 1]))
            itime += 1

    if itime != Num_Time:
        printexit('EXIT Inconsistent time lengths in Magnitude Data ' + str(itime) + ' ' + str(Num_Time))
    print('Read Magnitude data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))
    # End Reading in Magnitude data

    # Read in Depth Data into BasicInputTimeSeries
    with open(DepthDataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != '':
            printexit('EXIT: Wrong header on line 1 ' + Ftype + ' of ' + DepthDataFile)

        itime = 0
        for nextrow in csv_reader:
            if (len(nextrow) != Numberpixels + 1):
                printexit('EXIT: Incorrect row length Depth ' + str(itime) + ' ' + str(len(nextrow)))
            localtime = nextrow[0]
            if (itime != int(localtime)):
                printexit('EXIT: Unexpected Time in Depth ' + localtime + ' ' + str(itime))
            for iloc in range(0, Numberpixels):
                BasicInputTimeSeries[itime, iloc, 1] = nextrow[iloc + 1]
            itime += 1

    if itime != Num_Time:
        printexit('EXIT Inconsistent time lengths in Depth Data ' + str(itime) + ' ' + str(Num_Time))
    print('Read Depth data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))
    # End Reading in Depth data

    # Read in Multiplicity Data into BasicInputTimeSeries
    with open(MultiplicityDataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != '':
            printexit('EXIT: Wrong header on line 1 ' + Ftype + ' of ' + MultiplicityDataFile)

        itime = 0
        for nextrow in csv_reader:
            if (len(nextrow) != Numberpixels + 1):
                printexit('EXIT: Incorrect row length Multiplicity ' + str(itime) + ' ' + str(len(nextrow)))
            localtime = nextrow[0]
            if (itime != int(localtime)):
                printexit('EXIT: Unexpected Time in Multiplicity ' + localtime + ' ' + str(itime))
            for iloc in range(0, Numberpixels):
                BasicInputTimeSeries[itime, iloc, 2] = nextrow[iloc + 1]
            itime += 1

    if itime != Num_Time:
        printexit('EXIT Inconsistent time lengths in Multiplicity Data ' + str(itime) + ' ' + str(Num_Time))
    print('Read Multiplicity data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))
    # End Reading in Multiplicity data

    # Read in Rundle Multiplicity Data into BasicInputTimeSeries
    with open(RundleMultiplicityDataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != '':
            printexit('EXIT: Wrong header on line 1 ' + Ftype + ' of ' + RundleMultiplicityDataFile)

        itime = 0
        for nextrow in csv_reader:
            if (len(nextrow) != Numberpixels + 1):
                printexit('EXIT: Incorrect row length Rundle Multiplicity ' + str(itime) + ' ' + str(len(nextrow)))
            localtime = nextrow[0]
            if (itime != int(localtime)):
                printexit('EXIT: Unexpected Time in Rundle Multiplicity ' + localtime + ' ' + str(itime))
            for iloc in range(0, Numberpixels):
                BasicInputTimeSeries[itime, iloc, 3] = nextrow[iloc + 1]
            itime += 1

    if itime != Num_Time:
        printexit('EXIT Inconsistent time lengths in Rundle Multiplicity Data ' + str(itime) + ' ' + str(Num_Time))
    print('Read Rundle Multiplicity data locations ' + str(Nloc) + ' Time Steps ' + str(Num_Time))
    # End Reading in Rundle Multiplicity data

    # Read in Top Earthquake Data
    numberspecialeqs = 20
    Specialuse = np.full(numberspecialeqs, True, dtype=bool)
    Specialuse[14] = False
    Specialuse[15] = False
    Specialuse[18] = False
    Specialuse[19] = False
    Specialmags = np.empty(numberspecialeqs, dtype=np.float32)
    Specialdepth = np.empty(numberspecialeqs, dtype=np.float32)
    Speciallong = np.empty(numberspecialeqs, dtype=np.float32)
    Speciallat = np.empty(numberspecialeqs, dtype=np.float32)
    Specialdate = np.empty(numberspecialeqs, dtype='datetime64[D]')
    Specialxpos = np.empty(numberspecialeqs, dtype=np.int32)
    Specialypos = np.empty(numberspecialeqs, dtype=np.int32)
    Specialeqname = []

    with open(Topearthquakesfile, 'r') as read_obj:
        csv_reader = reader(read_obj)
        header = next(csv_reader)
        Ftype = header[0]
        if Ftype != 'date':
            printexit('EXIT: Wrong header on line 1 ' + Ftype + ' of ' + Topearthquakesfile)

        iquake = 0
        for nextrow in csv_reader:
            if (len(nextrow) != 6):
                printexit('EXIT: Incorrect row length Special Earthquakes ' + str(iquake) + ' ' + str(len(nextrow)))
            Specialdate[iquake] = nextrow[0]
            Speciallong[iquake] = nextrow[1]
            Speciallat[iquake] = nextrow[2]
            Specialmags[iquake] = nextrow[3]
            Specialdepth[iquake] = nextrow[4]
            Specialeqname.append(nextrow[5])
            ixpos = math.floor((Speciallong[iquake] + 120.0) * 10.0)
            ixpos = max(0, ixpos)
            ixpos = min(59, ixpos)
            iypos = math.floor((36.0 - Speciallat[iquake]) * 10.0)
            iypos = max(0, iypos)
            iypos = min(39, iypos)
            Specialxpos[iquake] = ixpos
            Specialypos[iquake] = iypos
            iquake += 1

    for iquake in range(0, numberspecialeqs):
        line = str(iquake) + ' mag ' + str(round(Specialmags[iquake], 1)) + ' Lat/Long '
        line += str(round(Speciallong[iquake], 2)) + ' ' + str(
            round(Speciallong[iquake], 2)) + ' ' + np.datetime_as_string(Specialdate[iquake])
        line += Specialeqname[iquake]
        print(line)
    printeq()

    # Possibly change Unit
    current_time = timenow()
    print(startbold + startred + current_time + ' Data read in ' + RunName + ' ' + RunComment + resetfonts)
    if Dailyunit != 1:
        if Dailyunit == 14:
            TimeIntervalUnitName = 'Fortnight'
        if Dailyunit == 28:
            TimeIntervalUnitName = 'LunarMonth'
        BasicInputTimeSeries, NumberofTimeunits, Num_Time, InitialDate, FinalDate = GatherUpData(BasicInputTimeSeries)
        current_time = timenow()
        print(startbold + startred + current_time + ' Data unit changed ' + RunName + ' ' + RunComment + resetfonts)
        Dateaxis = np.empty(Num_Time, dtype='datetime64[D]')
        Dateaxis[0] = np.datetime64(InitialDate).astype('datetime64[D]')
        for idate in range(1, Num_Time):
            Dateaxis[idate] = Dateaxis[idate - 1] + np.timedelta64(Dailyunit, 'D')
        for idate in range(0, Num_Time):
            Dateaxis[idate] = Dateaxis[idate] + np.timedelta64(int(Dailyunit / 2), 'D')
        print('Mid unit start time ' + np.datetime_as_string(Dateaxis[0]))

        Totalmag = np.zeros(Num_Time, dtype=np.float32)
        Totalefourthroot = np.zeros(Num_Time, dtype=np.float32)
        Totalesquareroot = np.zeros(Num_Time, dtype=np.float32)
        Totaleavgedmag = np.zeros(Num_Time, dtype=np.float32)
        Totalmult = np.zeros(Num_Time, dtype=np.float32)

        Totalmag[:] = BasicInputTimeSeries[:, :, 0].sum(axis=1)
        Totaleavgedmag = log_energy(BasicInputTimeSeries[:, :, 0], sumaxis=1)
        Totalmult[:] = BasicInputTimeSeries[:, :, 3].sum(axis=1)
        MagnitudeMethod = 1
        Tempseries = TransformMagnitude(BasicInputTimeSeries[:, :, 0])
        Totalefourthroot = Tempseries.sum(axis=1)
        MagnitudeMethod = 2
        Tempseries = TransformMagnitude(BasicInputTimeSeries[:, :, 0])
        Totalesquareroot = Tempseries.sum(axis=1)
        MagnitudeMethod = 0

        basenorm = Totalmult.max(axis=0)
        magnorm = Totalmag.max(axis=0)
        eavgedmagnorm = Totaleavgedmag.max(axis=0)
        efourthrootnorm = Totalefourthroot.max(axis=0)
        esquarerootnorm = Totalesquareroot.max(axis=0)
        print('Maximum Mult ' + str(round(basenorm, 2)) + ' Mag 0.15 ' + str(round(magnorm, 2))
              + ' E-avg 0.5 ' + str(round(eavgedmagnorm, 2)) + ' E^0.25 1.0 ' + str(round(efourthrootnorm, 2))
              + ' E^0.5 1.0 ' + str(round(esquarerootnorm, 2)))
        Totalmag = np.multiply(Totalmag, 0.15 * basenorm / magnorm)
        Totaleavgedmag = np.multiply(Totaleavgedmag, 0.5 * basenorm / eavgedmagnorm)
        Totalefourthroot = np.multiply(Totalefourthroot, basenorm / efourthrootnorm)
        Totalesquareroot = np.multiply(Totalesquareroot, basenorm / esquarerootnorm)

        plt.rcParams["figure.figsize"] = [16, 8]
        figure, ax = plt.subplots()
        datemin, datemax = makeadateplot(figure, ax, Dateaxis)
        ax.plot(Dateaxis, Totalmult, label='Multiplicity')
        ax.plot(Dateaxis, Totalmag, label='Summed Magnitude')
        ax.plot(Dateaxis, Totaleavgedmag, label='E-averaged Magnitude')
        ax.plot(Dateaxis, Totalefourthroot, label='Summed E^0.25')
        ax.plot(Dateaxis, Totalesquareroot, label='Summed E^0.5')
        ax.set_title('Observables summed over space')
        ax.set_xlabel("Years")
        ax.set_ylabel("Mult/Mag/Energy")
        ax.grid(True)
        ax.legend(loc='upper right')
        Addfixedearthquakes(ax, datemin, datemax)
        ax.tick_params('x', direction='in', length=15, width=2, which='major')
        ax.xaxis.set_minor_locator(mdates.YearLocator(1))
        ax.tick_params('x', direction='in', length=10, width=1, which='minor')
        figure.tight_layout()
        plt.show()

    else:
        print(' Data unit is the day and input this way')
        Dateaxis = np.empty(Num_Time, dtype='datetime64[D]')
        Dateaxis[0] = np.datetime64(InitialDate).astype('datetime64[D]')
        for idate in range(1, Num_Time):
            Dateaxis[idate] = Dateaxis[idate - 1] + np.timedelta64(Dailyunit, 'D')
        for idate in range(0, Num_Time):
            Dateaxis[idate] = Dateaxis[idate] + np.timedelta64(int(Dailyunit / 2), 'D')
        print('Mid unit start time ' + np.datetime_as_string(Dateaxis[0]))

    # Read in Fault Label Data into BasicInputStaticProps
    # No header for data
    with open(FaultLabelDataFile, 'r') as read_obj:
        csv_reader = reader(read_obj)

        iloc = 0
        if ReadFaultMethod == 1:
            for nextrow in csv_reader:
                if (len(nextrow) != 1):
                    printexit('EXIT: Incorrect row length Fault Label Data ' + str(iloc) + ' ' + str(len(nextrow)))
                RawFaultData[iloc] = nextrow[0]
                iloc += 1
        else:
            for nextrow in csv_reader:
                if (len(nextrow) != Numberxpixels):
                    printexit('EXIT: Incorrect row length Fault Label Data ' + str(iloc) + ' ' + str(
                        len(nextrow)) + ' ' + str(Numberxpixels))
                for jloc in range(0, len(nextrow)):
                    RawFaultData[iloc] = nextrow[jloc]
                    iloc += 1

    if iloc != Nloc:
        printexit('EXIT Inconsistent location lengths in Fault Label Data ' + str(iloc) + ' ' + str(Nloc))
    print('Read Fault Label data locations ' + str(Nloc))
    # End Reading in Fault Label data

    if NumFaultLabels == 1:
        BasicInputStaticProps[:, 0] = RawFaultData.astype(np.float32)
    else:  # remap fault label more reasonably
        unique, counts = np.unique(RawFaultData, return_counts=True)
        num = len(unique)
        print('Number Fault Collections ' + str(num))
        #    for i in range(0,num):
        #      print(str(unique[i]) + ' ' + str(counts[i]))

        BasicInputStaticProps[:, 0] = remapfaults(RawFaultData, Numberxpixels, Numberypixels, sfcurvelist).astype(
            np.float32)
        pix_faults = np.reshape(BasicInputStaticProps[:, 0], (40, 60)).astype(np.int)
        annotate_faults_ndarray(pix_faults, figsize=(24, 16))
        sfcurvelist2 = []
        for yloc in range(0, Numberypixels):
            for xloc in range(0, Numberxpixels):
                pixellocation = yloc * Numberxpixels + xloc
                [x, y] = sfcurvelist[pixellocation]
                sfcurvelist2.append([x, 39 - y])
        BasicInputStaticProps[:, 1] = remapfaults(RawFaultData, Numberxpixels, Numberypixels, sfcurvelist2).astype(
            np.float32)
        sfcurvelist3 = []
        for yloc in range(0, Numberypixels):
            for xloc in range(0, Numberxpixels):
                pixellocation = yloc * Numberxpixels + xloc
                [x, y] = sfcurvelist[pixellocation]
                sfcurvelist3.append([59 - x, y])
        BasicInputStaticProps[:, 2] = remapfaults(RawFaultData, Numberxpixels, Numberypixels, sfcurvelist3).astype(
            np.float32)
        sfcurvelist4 = []
        for yloc in range(0, Numberypixels):
            for xloc in range(0, Numberxpixels):
                pixellocation = yloc * Numberxpixels + xloc
                [x, y] = sfcurvelist[pixellocation]
                sfcurvelist4.append([59 - x, 39 - y])
        BasicInputStaticProps[:, 3] = remapfaults(RawFaultData, Numberxpixels, Numberypixels, sfcurvelist4).astype(
            np.float32)

    addRundleEMA = 1
    if Dailyunit != 14:
        addRundleEMA = 0
    RundleLambda = [2.5]
    RundleSteps = [144]
    NpropperTimeDynamicCalculated = 11 + addRundleEMA
    InputIndextogenerateEMA = 3
    FirstEMAIndex = 15
    NpropperTimeDynamic = NpropperTimeDynamicInput + NpropperTimeDynamicCalculated

    NpropperTimeStatic = NumFaultLabels
    #  NumpredbasicperTime = NpropperTimeDynamic
    NumpredbasicperTime = 1  # Can be 1 upto NpropperTimeDynamic
    NumpredFuturedperTime = NumpredbasicperTime

    # Setup Transformed Data
    # MagnitudeMethodTransform = 0 No Transform
    # MagnitudeMethodTransform = 1 E^0.25
    # MagnitudeMethodTransform = 2 E^0.5
    MagnitudeMethodTransform = 1
    TransformName = 'E^0.25'

    NpropperTime = NpropperTimeStatic + NpropperTimeDynamic
    InputPropertyNames = [' '] * NpropperTime

    DynamicNames = ['Magnitude Now', 'Depth Now', 'Multiplicity Now', 'Mult >3.29 Now', 'Mag 2/3 Month Back',
                    'Mag 1.5 Month Back', 'Mag 3 Months Back', 'Mag 6 Months Back',
                    'Mag Year Back', TransformName + ' Now', TransformName + ' 2/3 Month Back',
                    TransformName + ' 1.5 Month Back', TransformName + ' 3 Months Back',
                    TransformName + ' 6 Months Back', TransformName + ' Year Back']
    if Dailyunit == 14:
        DynamicNames = ['Magnitude 2 weeks Now', 'Depth 2 weeks Now', 'Multiplicity 2 weeks Now',
                        'Mult >3.29 2 weeks Now',
                        'Mag 4 Weeks Back', 'Mag 2 Months Back', 'Mag 3 Months Back', 'Mag 6 Months Back',
                        'Mag Year Back',
                        TransformName + ' 2 weeks Back', TransformName + ' 4 weeks Back',
                        TransformName + ' 2 Months Back', TransformName + ' 3 Months Back',
                        TransformName + ' 6 Months Back', TransformName + ' Year Back']
        if addRundleEMA != 0:
            for i in range(0, addRundleEMA):
                DynamicNames.append('MultEMA' + str(RundleSteps[i]) + ' L' + str(round(RundleLambda[i], 2)))
    Property_is_Intensive = np.full(NpropperTime, True, dtype=np.bool)
    for iprop in range(0, NpropperTimeStatic):
        InputPropertyNames[iprop] = 'Fault ' + str(iprop)
    for iprop in range(0, NpropperTimeDynamic):
        InputPropertyNames[iprop + NpropperTimeStatic] = DynamicNames[iprop]
    Num_Extensive = 0

    CDSpecial = False
    ScaleProperties = True
    GenerateFutures = False
    GenerateSequences = True
    PredictionsfromInputs = True
    ConvertDynamicPredictedQuantity = False
    AddSpecialstoSummedplots = True
    UseRealDatesonplots = True
    EarthquakeImagePlots = False
    UseFutures = False
    PopulationNorm = False
    OriginalNloc = Nloc
    MapLocation = False

    # Add summed magnitudes as properties to use in prediction and Calculated Properties for some
    # Calculated Properties are sums starting at given time and are set to NaN if necessary
    NumTimeSeriesCalculatedBasic = 9
    NumTimeSeriesCalculated = 2 * NumTimeSeriesCalculatedBasic + 1
    NamespredCalculated = ['Mag 2/3 Month Ahead', 'Mag 1.5 Month Ahead', 'Mag 3 Months Ahead', 'Mag 6 Months Ahead',
                           'Mag Year Ahead Ahead', 'Mag 2 Years Ahead', 'Mag 4 years Ahead', 'Mag Skip 1, Year ahead',
                           'Mag 2 years 2 ahead',
                           TransformName + ' Daily Now', TransformName + ' 2/3 Month Ahead',
                           TransformName + ' 1.5 Month Ahead', TransformName + ' 3 Months Ahead',
                           TransformName + ' 6 Months Ahead', TransformName + ' Year Ahead',
                           TransformName + ' 2 Years Ahead', TransformName + ' 4 years Ahead',
                           TransformName + ' Skip 1, Year ahead', TransformName + ' 2 years 2 ahead']
    Unitjumps = [23, 46, 92, 183, 365, 730, 1460, 365, 730]
    Unitdelays = [0, 0, 0, 0, 0, 0, 0, 365, 730]
    Plottingdelay = 1460
    if Dailyunit == 14:
        NumTimeSeriesCalculatedBasic = 9
        NumTimeSeriesCalculated = 2 * NumTimeSeriesCalculatedBasic + 1
        NamespredCalculated = ['Mag 4 Weeks Ahead', 'Mag 2 Month Ahead', 'Mag 3 Months Ahead', 'Mag 6 Months Ahead',
                               'Mag Year Ahead', 'Mag 2 Years Ahead', 'Mag 4 years Ahead', 'Mag Skip 1, Year ahead',
                               'Mag 2 years 2 ahead',
                               TransformName + ' 2 Weeks Now', TransformName + ' 4 Weeks Ahead',
                               TransformName + ' 2 Months Ahead', TransformName + ' 3 Months Ahead',
                               TransformName + ' 6 Months Ahead',
                               TransformName + ' Year Ahead', TransformName + ' 2 Years Ahead',
                               TransformName + ' 4 years Ahead', TransformName + ' Skip 1, Year ahead',
                               TransformName + ' 2 years 2 ahead']
        Unitjumps = [2, 4, 7, 13, 26, 52, 104, 26, 52]
        Unitdelays = [0, 0, 0, 0, 0, 0, 0, 26, 52]
        Plottingdelay = 104

    NumpredbasicperTime += NumTimeSeriesCalculated
    CalculatedTimeSeries = np.empty([Num_Time, Nloc, NumTimeSeriesCalculated], dtype=np.float32)
    for icalc in range(0, NumTimeSeriesCalculatedBasic):
        newicalc = icalc + 1 + NumTimeSeriesCalculatedBasic
        for itime in range(0, Num_Time):
            MagnitudeMethod = 0
            CalculatedTimeSeries[itime, :, icalc] = AggregateEarthquakes(itime, Unitdelays[icalc], Unitjumps[icalc],
                                                                         Nloc,
                                                                         BasicInputTimeSeries[:, :, 0], 0)
            MagnitudeMethod = MagnitudeMethodTransform
            CalculatedTimeSeries[itime, :, newicalc] = TransformMagnitude(CalculatedTimeSeries[itime, :, icalc])
            MagnitudeMethod = 0
        current_time = timenow()
        print(startbold + startred + 'Earthquake ' + str(icalc) + ' ' + NamespredCalculated[
            icalc] + ' ' + current_time + ' ' + RunName + resetfonts)
        print(startbold + startred + 'Earthquake ' + str(newicalc) + ' ' + NamespredCalculated[
            newicalc] + ' ' + current_time + ' ' + RunName + resetfonts)
    MagnitudeMethod = MagnitudeMethodTransform
    CalculatedTimeSeries[:, :, NumTimeSeriesCalculatedBasic] = TransformMagnitude(BasicInputTimeSeries[:, :, 0])
    MagnitudeMethod = 0
    print(startbold + startred + 'Earthquake ' + str(NumTimeSeriesCalculatedBasic) + ' ' + NamespredCalculated[
        NumTimeSeriesCalculatedBasic] + ' ' + current_time + ' ' + RunName + resetfonts)

    for iprop in range(0, NumTimeSeriesCalculated):
        InputPropertyNames.append(NamespredCalculated[iprop])

