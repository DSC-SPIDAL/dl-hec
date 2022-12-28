
# Read Hydrology
if Hydrology:
    PreparedDataFile = APPLDIR + '/data.tar.bz2'
    content = Shell.run("ls '/content/gdrive/My Drive/Colab Datasets/Hydrology'")
    print(content)
    content = Shell.run ("tar xjf '/content/gdrive/My Drive/Colab Datasets/Hydrology/data.tar.bz2' -C "
               "'/content/gdrive/My Drive/Colab Datasets/Hydrology'")
    print(content)
    import json

    RawInputStaticProps = np.load(APPLDIR + '/BasicInputStaticProps.npy', allow_pickle=True)
    RawInputTimeSeries = np.load(APPLDIR + '/BasicInputTimeSeries.npy', allow_pickle=True)
    NuminputSeries = RawInputTimeSeries.shape[1]
    NuminputProps = RawInputStaticProps.shape[1]
    print('Input Hydrology Shapes ' + str(RawInputTimeSeries.shape) + ' ' + str(RawInputStaticProps.shape))

    with open(APPLDIR + '/metadata.json', 'r') as f:
        metadata = json.load(f)
    Nloc = metadata['Nloc']
    TimeSeriesmetadata = metadata['BasicInputTimeSeries']
    InitialDate = datetime.strptime(TimeSeriesmetadata['initial_date'], '%Y-%m-%dT%H:%M:%S.%f000')
    FinalDate = datetime.strptime(TimeSeriesmetadata['end_date'], '%Y-%m-%dT%H:%M:%S.%f000')
    NumberofTimeunits = (FinalDate - InitialDate).days + 1
    print(InitialDate.strftime("%d/%m/%Y") + ' To ' + FinalDate.strftime("%d/%m/%Y") + ' days ' + str(
        NumberofTimeunits) + ' Locations ' + str(Nloc))
    TimeSeriesLabels = TimeSeriesmetadata['fields']
    print(TimeSeriesLabels)
    StaticPropsmetadata = metadata['BasicInputStaticProps']
    RawLabels = StaticPropsmetadata['fields']
    print(RawLabels)
    BasicInputTimeSeries = np.delete(RawInputTimeSeries, [0, 1], 1)
    BasicInputTimeSeries = np.reshape(BasicInputTimeSeries, [NumberofTimeunits, Nloc, NuminputSeries - 2])
    BasicInputStaticProps = np.delete(RawInputStaticProps, [0, 12, 21, 22], 1)
    StaticLabels = np.delete(RawLabels, [0, 12, 21, 22], 0)

    Num_Time = NumberofTimeunits
    NFIPS = Nloc
    Locationfips = np.empty(NFIPS, dtype=int)  # integer version of FIPs/gauge_id
    Locationcolumns = []  # String version of FIPS/gauge_id
    FIPSintegerlookup = {}
    FIPSstringlookup = {}
    Locationname = ['Empty'] * NFIPS
    Locationstate = [' '] * NFIPS
    Locationpopulation = np.ones(NFIPS, dtype=int)
    gauge_idvalues = metadata['locs']
    placenames = metadata['loc_names']
    for iloc in range(0, Nloc):
        fips = str(gauge_idvalues[iloc])
        Locationfips[iloc] = int(fips)
        Locationcolumns.append(fips)
        FIPSintegerlookup[int(fips)] = iloc
        FIPSstringlookup[fips] = iloc
        Locationname[iloc] = placenames[iloc]

    CDSpecial = False
    NpropperTimeDynamic = 6
    NpropperTimeStatic = 27
    NumpredbasicperTime = NpropperTimeDynamic
    NumpredFuturedperTime = NumpredbasicperTime
    NpropperTime = NpropperTimeStatic + NpropperTimeDynamic
    InputPropertyNames = [' '] * NpropperTime
    Property_is_Intensive = np.full(NpropperTime, True, dtype=np.bool)
    for iprop in range(0, NpropperTimeStatic):
        InputPropertyNames[iprop] = StaticLabels[iprop]
    for iprop in range(0, NpropperTimeDynamic):
        InputPropertyNames[iprop + NpropperTimeStatic] = TimeSeriesLabels[iprop + 2]
    Num_Extensive = 0

    ScaleProperties = True
    GenerateFutures = False
    GenerateSequences = True
    PredictionsfromInputs = True
    ConvertDynamicPredictedQuantity = False

    UseFutures = False
    PopulationNorm = False
    DynamicPropertyTimeSeries = np.empty_like(BasicInputTimeSeries, dtype=np.float32)
    CountNaN = np.zeros(NpropperTimeDynamic, dtype=np.int)
    for itime in range(0, NumberofTimeunits):
        for iloc in range(0, Nloc):
            for iprop in range(0, NpropperTimeDynamic):
                localval = BasicInputTimeSeries[itime, iloc, iprop]
                if np.math.isnan(localval):
                    localval = NaN
                    CountNaN[iprop] += 1
                else:
                    if (localval < 0.0) and (iprop == 5):
                        localval = NaN
                        CountNaN[iprop] += 1
                DynamicPropertyTimeSeries[itime, iloc, iprop] = localval
    print(startbold + startred + 'Input NaN values ' + resetfonts)
    for iprop in range(0, NpropperTimeDynamic):
        print(InputPropertyNames[iprop + NpropperTimeStatic] + ' ' + str(CountNaN[iprop]))

    BasicInputTimeSeries = None
    if GarbageCollect:
        gc.collect()

    # Overall Parameters set for Hydrology
    SymbolicWindows = True
    Tseq = 21
    Plotsplitsize = 6