"""##Earthquake Routines"""


def printeq():
    qsort = np.argsort(Specialdate)
    for jquake in range(0, numberspecialeqs):
        iquake = qsort[jquake]
        print(str(iquake) + ' ' + str(Specialdate[iquake]) + ' ' + str(round(Specialmags[iquake], 1)) + ' ' + Specialeqname[
            iquake])


def Addfixedearthquakes(plotpointer, graphmin, graphmax, ylogscale=False, quakecolor=None, Dateplot=True, vetoquake=None):
    if vetoquake is None:  # Vetoquake = True means do not plot this quake
        vetoquake = np.full(numberspecialeqs, False, dtype=np.bool)
    if quakecolor is None:  # Color of plot
        quakecolor = 'black'
    Place = np.arange(numberspecialeqs, dtype=np.int)
    Place[8] = 11
    Place[10] = 3
    Place[12] = 16
    Place[7] = 4
    Place[2] = 5
    Place[4] = 14
    Place[11] = 18

    ymin, ymax = plotpointer.get_ylim()  # Or work with transform=ax.transAxes
    qindex = 1
    fudge = 0.03
    qsort = np.argsort(Specialdate)
    for jquake in range(0, numberspecialeqs):
        iquake = qsort[jquake]

        if qindex == 16:
            continue
        # This is the x position for the vertical line
        if Dateplot:
            x_line_annotation = Specialdate[iquake]  # numpy date format
        else:
            x_line_annotation = Numericaldate[iquake]  # Float where each interval 1 and start is 0

        # This is the x position for the label
        if Dateplot:
            x_text_annotation = x_line_annotation - np.timedelta64(2 * Dailyunit, 'D')
        else:
            x_text_annotation = x_line_annotation - 2.0

        if Specialuse[iquake]:

            # Draw a text
            ascii = str(round(Specialmags[iquake], 1)) + '\n' + Specialeqname[iquake]
            ascii = 'EQ' + str(qindex)
            qindex += 1
            fudge = 0.09 - fudge

            if (x_line_annotation < graphmin) or (x_line_annotation > graphmax):
                continue

            if vetoquake[iquake]:
                continue

            # Draw a line at the position
            ydelta = 0.15
            plotpointer.axvline(x=x_line_annotation, linestyle='solid', alpha=1.0, linewidth=2.0, color='black',
                                ymin=1.0 - ydelta)

            acfudge = fudge
            if qindex == 13:
                acfudge = 0.09
            if ylogscale:
                yminl = max(0.01 * ymax, ymin)
                yminl = math.log(yminl, 10)
                ymaxl = math.log(ymax, 10)
                logyplot = yminl + (0.1 + 0.8 * (float(Place[iquake]) / float(numberspecialeqs - 1))) * (ymaxl - yminl)
                yplot = pow(10, logyplot)
                yplot = ymax - (ydelta + acfudge) * (ymax - ymin)
            else:
                yplot = ymax - (0.1 + 0.8 * (float(Place[iquake]) / float(numberspecialeqs - 1))) * (ymax - ymin)
                yplot = ymax - (ydelta + acfudge) * (ymax - ymin)
            if Dateplot:
                if x_text_annotation > graphmax - np.timedelta64(1200, 'D'):
                    x_text_annotation = graphmax - np.timedelta64(1200, 'D')
                x_text_annotation = max(x_text_annotation, graphmin + np.timedelta64(400, 'D'))
            else:
                if x_text_annotation > graphmax - 60:
                    x_text_annotation = graphmax - 60
                x_text_annotation = max(x_text_annotation, graphmin + 20)
            #      print(str(yplot) + " " + str(ymin) + " " + str(ymax) + " " + str(x_text_annotation) + " " + str(x_line_annotation)) + " " + ascii
            #      print(str(qindex-1) + ' ' + str(iquake) + ' ' +str(Specialdate[iquake]) + ' ' + str(x_line_annotation) + ' ' + str(x_text_annotation) + ' ' + ascii + ' ' +  str(round(Specialmags[iquake],1)) + ' ' + Specialeqname[iquake])
            plotpointer.text(x=x_text_annotation, y=yplot, s=wraptotext(ascii, size=10), alpha=1.0, color='black',
                             fontsize=10)


def quakesearch(iquake, iloc):
    # see if top earthquake iquake llies near location iloc
    # result = 0 NO; =1 YES Primary: locations match exactly; = -1 Secondary: locations near
    # iloc is location before mapping
    xloc = iloc % 60
    yloc = (iloc - xloc) / 60
    if (xloc == Specialxpos[iquake]) and (yloc == Specialypos[iquake]):
        return 1
    if (abs(xloc - Specialxpos[iquake]) <= 1) and (abs(yloc - Specialypos[iquake]) <= 1):
        return -1
    return 0


# Read Earthquake Data
def log_sum_exp10(ns, sumaxis=0):
    max_v = np.max(ns, axis=None)
    ds = ns - max_v
    sum_of_exp = np.power(10, ds).sum(axis=sumaxis)
    return max_v + np.log10(sum_of_exp)


def log_energyweightedsum(nvalue, ns, sumaxis=0):
    max_v = np.max(ns, axis=None)
    ds = ns - max_v
    ds = np.power(10, 1.5 * ds)
    dvalue = (np.multiply(nvalue, ds)).sum(axis=sumaxis)
    ds = ds.sum(axis=0)
    return np.divide(dvalue, ds)


# Set summed magnitude as log summed energy = 10^(1.5 magnitude)
def log_energy(mag, sumaxis=0):
    return log_sum_exp10(1.5 * mag, sumaxis=sumaxis) / 1.5


def AggregateEarthquakes(itime, DaysDelay, DaysinInterval, Nloc, Eqdata, Approach, weighting=None):
    if (itime + DaysinInterval + DaysDelay) > NumberofTimeunits:
        return np.full([Nloc], NaN, dtype=np.float32)
    if Approach == 0:  # Magnitudes
        if MagnitudeMethod == 0:
            TotalMagnitude = log_energy(Eqdata[itime + DaysDelay:itime + DaysinInterval + DaysDelay])
        else:
            TotalMagnitude = Eqdata[itime + DaysDelay:itime + DaysinInterval + DaysDelay, :].sum(axis=0)
        return TotalMagnitude
    if Approach == 1:  # Depth -- energy weighted
        WeightedResult = log_energyweightedsum(Eqdata[itime + DaysDelay:itime + DaysinInterval + DaysDelay],
                                               weighting[itime + DaysDelay:itime + DaysinInterval + DaysDelay])
        return WeightedResult
    if Approach == 2:  # Multiplicity -- summed
        SimpleSum = Eqdata[itime + DaysDelay:itime + DaysinInterval + DaysDelay, :].sum(axis=0)
        return SimpleSum


# MagnitudeMethodTransform = 0 No Transform
# MagnitudeMethodTransform = 1 E^0.25
# MagnitudeMethodTransform = 2 E^0.5
def TransformMagnitude(mag):
    if MagnitudeMethod == 0:
        return mag
    if MagnitudeMethod == 1:
        return np.power(10, 0.375 * (mag - 3.29))
    return np.power(10, 0.75 * (mag - 3.29))


# Change Daily Unit
# Accumulate data in Dailyunit chunks.
# This changes data so it looks like daily data bu really collections of chunked data.
# For earthquakes, the aggregations uses energy averaging for depth and magnitude. It just adds for multiplicity
def GatherUpData(OldInputTimeSeries):
    Skipped = NumberofTimeunits % Dailyunit
    NewInitialDate = InitialDate + timedelta(days=Skipped)
    NewNum_Time = int(Num_Time / Dailyunit)
    NewFinalDate = NewInitialDate + Dailyunit * timedelta(days=NewNum_Time - 1)
    print(' Daily Unit ' + str(Dailyunit) + ' number of ' + TimeIntervalUnitName + ' Units ' + str(NewNum_Time) + ' ' +
          NewInitialDate.strftime("%d/%m/%Y") + ' To ' + NewFinalDate.strftime("%d/%m/%Y"))
    NewInputTimeSeries = np.empty([NewNum_Time, Nloc, NpropperTimeDynamicInput], dtype=np.float32)
    for itime in range(0, NewNum_Time):
        NewInputTimeSeries[itime, :, 0] = AggregateEarthquakes(Skipped + itime * Dailyunit, 0, Dailyunit, Nloc,
                                                               BasicInputTimeSeries[:, :, 0], 0)
        NewInputTimeSeries[itime, :, 1] = AggregateEarthquakes(Skipped + itime * Dailyunit, 0, Dailyunit, Nloc,
                                                               BasicInputTimeSeries[:, :, 1], 1,
                                                               weighting=BasicInputTimeSeries[:, :, 0])
        NewInputTimeSeries[itime, :, 2] = AggregateEarthquakes(Skipped + itime * Dailyunit, 0, Dailyunit, Nloc,
                                                               BasicInputTimeSeries[:, :, 2], 2)
        NewInputTimeSeries[itime, :, 3] = AggregateEarthquakes(Skipped + itime * Dailyunit, 0, Dailyunit, Nloc,
                                                               BasicInputTimeSeries[:, :, 3], 2)
    return NewInputTimeSeries, NewNum_Time, NewNum_Time, NewInitialDate, NewFinalDate


# make numpy array not a list
# Return Exponential Moving Average
def MakeEMA(TimeSeries, Nsteps):
    datasize = len(TimeSeries)
    EMASeries = np.empty(datasize, dtype=np.float32)
    WeightedSum = 0.0
    WeightedCount = 0.0
    beta = (Nsteps - 1.0) / (Nsteps + 1.0)
    for i in range(0, datasize):
        WeightedSum = TimeSeries[i] + beta * WeightedSum
        WeightedCount = 1.0 + beta * WeightedCount
        EMASeries[i] = WeightedSum / WeightedCount
    return EMASeries


# Insist on a minimum count at each time value and then find Exponential Moving Average
def MakeEMAMinCT(TimeSeries, Nsteps, Lambda):
    Mincount = np.mean(TimeSeries) * Lambda
    Temporaryarray = np.maximum(TimeSeries, Mincount)
    EMASeries = MakeEMA(Temporaryarray, Nsteps)
    return EMASeries


"""###Space Filling Curves"""

from gilbert import cal_gilbert2d

"""###Read Earthquake Data"""


def ReadEarthquakeData():
    read1950 = True
    Eigenvectors = 2
    UseEarthquakeEigenSystems = False
    Dailyunit = 14
    addwobblingposition = False

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
    print(' Pixels ' + str(Nloc) + ' x dimension ' + str(Nlocaxislengths[0]) + ' y dimension ' + str(Nlocaxislengths[1]))

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
        line += str(round(Speciallong[iquake], 2)) + ' ' + str(round(Speciallong[iquake], 2)) + ' ' + np.datetime_as_string(
            Specialdate[iquake])
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
                    printexit(
                        'EXIT: Incorrect row length Fault Label Data ' + str(iloc) + ' ' + str(len(nextrow)) + ' ' + str(
                            Numberxpixels))
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
                    TransformName + ' 1.5 Month Back', TransformName + ' 3 Months Back', TransformName + ' 6 Months Back',
                    TransformName + ' Year Back']
    if Dailyunit == 14:
        DynamicNames = ['Magnitude 2 weeks Now', 'Depth 2 weeks Now', 'Multiplicity 2 weeks Now', 'Mult >3.29 2 weeks Now',
                        'Mag 4 Weeks Back', 'Mag 2 Months Back', 'Mag 3 Months Back', 'Mag 6 Months Back', 'Mag Year Back',
                        TransformName + ' 2 weeks Back', TransformName + ' 4 weeks Back', TransformName + ' 2 Months Back',
                        TransformName + ' 3 Months Back', TransformName + ' 6 Months Back', TransformName + ' Year Back']
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
            CalculatedTimeSeries[itime, :, icalc] = AggregateEarthquakes(itime, Unitdelays[icalc], Unitjumps[icalc], Nloc,
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


"""###Earthquake Eigensystems"""

if Earthquake:
    if UseEarthquakeEigenSystems:
        content = Shell.run("pip install scipy -U")
        print(content)
        import scipy as sc
        import scipy.linalg as solver

        version = sc.version.version
        print('SciPy version ' + str(version))
        # x = np.array([[1,2.0],[2.0,0]])
        # w, v = solver.eigh(x,  driver='evx')
        # print(w)
        # print(v)

"""###Multiplicity Data"""


def histogrammultiplicity(Type, numbins, Data):
    hitcounts = np.zeros(Nloc, dtype=np.int)
    rawcounts = np.zeros(Nloc, dtype=np.int)
    for iloc in range(0, Nloc):
        rawcounts[iloc] = np.int(0.1 + Data[:, iloc].sum(0))
        hitcounts[iloc] = np.int(min(numbins, rawcounts[iloc]))
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)
    plt.rcParams.update({'font.size': 9})
    plt.rcParams["figure.figsize"] = [8, 6]
    plt.hist(hitcounts, numbins, facecolor='b', alpha=0.75, log=True)
    plt.title('\n'.join(wrap(RunComment + ' ' + RunName + ' ' + Type + ' Earthquake Count per location ', 70)))
    plt.xlabel('Hit Counts')
    plt.ylabel('Occurrences')
    plt.grid(True)
    plt.show()
    return rawcounts


def threebythree(pixellocation, numxlocations, numylocations):
    indices = np.empty([3, 3], dtype=np.int)

    y = int(0.1 + pixellocation / numxlocations)
    x = pixellocation - y * numxlocations
    bottomx = max(0, x - 1)
    bottomx = min(bottomx, numxlocations - 3)
    bottomy = max(0, y - 1)
    bottomy = min(bottomy, numylocations - 3)
    for ix in range(0, 3):
        for iy in range(0, 3):
            x = bottomx + ix
            y = bottomy + iy
            pixellocation = y * numxlocations + x
            indices[ix, iy] = pixellocation
    return indices


if Earthquake:
    MappedLocations = np.arange(0, Nloc, dtype=np.int)
    LookupLocations = np.arange(0, Nloc, dtype=np.int)
    MappedNloc = Nloc
    histogrammultiplicity('Basic', 100, BasicInputTimeSeries[:, :, 2])
    nbins = 10
    if read1950:
        nbins = 50
    rawcounts1 = histogrammultiplicity('Rundle > 3.29', nbins, BasicInputTimeSeries[:, :, 3])
    TempTimeSeries = np.zeros([Num_Time, Nloc], dtype=np.float32)
    for iloc in range(0, Nloc):
        indices = threebythree(iloc, 60, 40)
        for itime in range(0, Num_Time):
            sum3by3 = 0.0
            for ix in range(0, 3):
                for iy in range(0, 3):
                    pixellocation = indices[ix, iy]
                    sum3by3 += BasicInputTimeSeries[itime, pixellocation, 3]
            TempTimeSeries[itime, iloc] = sum3by3
    nbins = 40
    if read1950:
        nbins = 150
    rawcounts2 = histogrammultiplicity('3x3 Rundle > 3.29', nbins, TempTimeSeries)
    #
    # Define "Interesting Locations"
    if read1950:
        singleloccut = 25
        groupedloccut = 110
        singleloccut = 7.1
        groupedloccut = 34.1

    #    groupedloccut = 1000000000
    else:
        singleloccut = 5.1
        groupedloccut = 24.9
    MappedLocations.fill(-1)
    MappedNloc = 0
    ct1 = 0
    ct2 = 0
    for iloc in range(0, Nloc):
        if rawcounts1[iloc] >= singleloccut:
            ct1 += 1
        if rawcounts2[iloc] >= groupedloccut:
            ct2 += 1
        if rawcounts1[iloc] < singleloccut and rawcounts2[iloc] < groupedloccut:
            continue
        MappedLocations[iloc] = MappedNloc
        MappedNloc += 1

    LookupLocations = None
    LookupLocations = np.empty(MappedNloc, dtype=np.int)
    for iloc in range(0, Nloc):
        jloc = MappedLocations[iloc]
        if jloc >= 0:
            LookupLocations[jloc] = iloc

    TempTimeSeries = None
    print('Total ' + str(MappedNloc) + ' Single location multiplicity cut ' + str(singleloccut) +
          ' ' + str(ct1) + ' 3x3 ' + str(groupedloccut) + ' ' + str(ct2))

    if UseEarthquakeEigenSystems:
        if Eigenvectors > 0:
            UseTopEigenTotal = 16
            UseTopEigenLocal = 0
            if Eigenvectors > 1:
                UseTopEigenLocal = 4
            Num_EigenProperties = UseTopEigenTotal + UseTopEigenLocal
            EigenTimeSeries = np.empty([Num_Time, MappedNloc], dtype=np.float32)
            PsiTimeSeries = np.empty([Num_Time, MappedNloc], dtype=np.float32)
            FiTimeSeries = np.empty([Num_Time, MappedNloc], dtype=np.float32)
            EigenTimeSeries[:, :] = BasicInputTimeSeries[:, LookupLocations, 3]
            StoreEigenvectors = np.zeros([Num_Time, MappedNloc, MappedNloc], dtype=np.float32)
            StoreEigencorrels = np.zeros([Num_Time, MappedNloc, MappedNloc], dtype=np.float32)
            StoreNormingfactor = np.zeros([Num_Time], dtype=np.float32)
            StoreNormingfactor1 = np.zeros([Num_Time], dtype=np.float32)
            StoreNormingfactor2 = np.zeros([Num_Time], dtype=np.float32)
            current_time = timenow()
            print(startbold + startred + 'Start Eigen Earthquake '
                  + current_time + ' ' + RunName + resetfonts)

            for itime in range(0, Num_Time):
                imax = itime
                imin = max(0, imax - 25)
                Result = np.zeros(MappedNloc, dtype=np.float64)
                Result = AggregateEarthquakes(imin, 0, imax - imin + 1, MappedNloc, EigenTimeSeries[:, :], 2)
                PsiTimeSeries[itime, :] = Result
                FiTimeSeries[itime, :] = EigenTimeSeries[itime, :]

            current_time = timenow()
            print(startbold + startred + 'End Eigen Earthquake 1 '
                  + current_time + ' ' + RunName + resetfonts)
            Eigenvals = np.zeros([Num_Time, MappedNloc], dtype=np.float32)
            Chi1 = np.zeros(Num_Time, dtype=np.float32)
            Chi2 = np.zeros(Num_Time, dtype=np.float32)
            Sumai = np.zeros(Num_Time, dtype=np.float32)
            Bestindex = np.zeros(Num_Time, dtype=np.int)
            Numbereigs = np.zeros(Num_Time, dtype=np.int)
            Besttrailingindex = np.zeros(Num_Time, dtype=np.int)
            Eig0coeff = np.zeros(Num_Time, dtype=np.float32)
            meanmethod = 0
            if meanmethod == 1:
                Meanovertime = np.empty(MappedNloc, dtype=np.float32)
                sigmaovertime = np.empty(MappedNloc, dtype=np.float32)
                Meanovertime = FiTimeSeries.mean(axis=0)
                Meanovertime = Meanovertime.reshape(1, MappedNloc)
                sigmaovertime = FiTimeSeries.std(axis=0)
                sigmaovertime = sigmaovertime.reshape(1, MappedNloc)
            countbad = 0
            OldActualNumberofLocationsUsed = -1
            for itime in range(25, Num_Time):
                LocationCounts = FiTimeSeries[0:itime, :].sum(axis=0)
                NumLocsToday = np.count_nonzero(LocationCounts)
                Nonzeromapping = np.zeros(NumLocsToday, dtype=np.int)
                ActualNumberofLocationsUsed = 0
                for ipos in range(0, MappedNloc):
                    if LocationCounts[ipos] == 0:
                        continue
                    Nonzeromapping[ActualNumberofLocationsUsed] = ipos
                    ActualNumberofLocationsUsed += 1
                if ActualNumberofLocationsUsed <= 1:
                    print(str(itime) + ' Abandoned ' + str(ActualNumberofLocationsUsed))
                    continue
                FiHatTimeSeries = np.empty([itime + 1, ActualNumberofLocationsUsed], dtype=np.float32)
                if meanmethod == 1:
                    FiHatTimeSeries[:, :] = np.divide(
                        np.subtract(FiTimeSeries[0:(itime + 1), Nonzeromapping], Meanovertime[0, Nonzeromapping]),
                        sigmaovertime[0, Nonzeromapping])
                else:
                    FiHatTimeSeries[:, :] = FiTimeSeries[0:(itime + 1), Nonzeromapping]
                #          FiHatTimeSeries[:,:] = PsiTimeSeries[0:(itime+1),Nonzeromapping]
                CorrelationMatrix = np.corrcoef(FiHatTimeSeries, rowvar=False)
                bad = np.count_nonzero(np.isnan(CorrelationMatrix))
                if bad > 0:
                    countbad += 1
                    continue
                evalues, evectors = solver.eigh(CorrelationMatrix)
                Newevector = evectors[:, ActualNumberofLocationsUsed - 1]
                Newevalue = evalues[ActualNumberofLocationsUsed - 1]
                debug = False
                if debug:
                    if OldActualNumberofLocationsUsed == ActualNumberofLocationsUsed:
                        Mapdiff = np.where(np.not_equal(OldNonzeromapping, Nonzeromapping), 1, 0.).sum()
                        if Mapdiff > 0:
                            print(str(itime) + ' Change in mapping ' + str(ActualNumberofLocationsUsed) + ' Change ' + str(
                                Mapdiff))
                        else:
                            Corrdiff = np.absolute(np.subtract(OldCorrelationMatrix, CorrelationMatrix)).sum()
                            Corrorg = np.absolute(CorrelationMatrix).sum()
                            yummy = CorrelationMatrix.dot(Oldevector)
                            vTMv = yummy.dot(Oldevector)
                            Doubleyummy = CorrelationMatrix.dot(Newevector)
                            newvTMv = Doubleyummy.dot(Newevector)
                            print(str(itime) + ' Change in correlation ' + str(ActualNumberofLocationsUsed) + ' Change '
                                  + str(Corrdiff) + ' original ' + str(Corrorg) + ' eval ' + str(Oldevalue) + ' new '
                                  + str(Newevalue) + ' vTMv ' + str(vTMv) + ' New ' + str(newvTMv))

                    else:
                        print(str(itime) + ' Change in size ' + str(OldActualNumberofLocationsUsed) + ' ' +
                              str(ActualNumberofLocationsUsed))

                OldActualNumberofLocationsUsed = ActualNumberofLocationsUsed
                OldNonzeromapping = Nonzeromapping
                OldCorrelationMatrix = CorrelationMatrix
                Oldevector = Newevector
                Oldevalue = Newevalue

                normcoeff = 100.0 / evalues.sum()
                evalues = np.multiply(evalues, normcoeff)
                Numbereigs[itime] = ActualNumberofLocationsUsed

                for ieig in range(0, ActualNumberofLocationsUsed):
                    Eigenvals[itime, ieig] = evalues[ActualNumberofLocationsUsed - ieig - 1]
                chival = 0.0
                sumaieig = 0.0
                Checkvector = np.zeros(ActualNumberofLocationsUsed, dtype=np.float32)
                largesteigcoeff = -1.0
                largestindex = -1

                Keepaisquared = np.zeros(ActualNumberofLocationsUsed, dtype=np.float32)
                for ieig in range(0, ActualNumberofLocationsUsed):
                    aieig = 0.0
                    backwards = ActualNumberofLocationsUsed - ieig - 1
                    for vectorindex in range(0, ActualNumberofLocationsUsed):
                        StoreEigenvectors[itime, backwards, Nonzeromapping[vectorindex]] = evectors[vectorindex, ieig]
                        aieig += evectors[vectorindex, ieig] * PsiTimeSeries[itime, Nonzeromapping[vectorindex]]
                    for vectorindex in range(0, ActualNumberofLocationsUsed):
                        Checkvector[vectorindex] += aieig * evectors[vectorindex, ieig]
                    aieig *= aieig
                    chival += aieig * evalues[ieig]
                    sumaieig += aieig
                    Keepaisquared[backwards] = aieig

                for ieig in range(0, ActualNumberofLocationsUsed):
                    backwards = ActualNumberofLocationsUsed - ieig - 1
                    aieig = Keepaisquared[backwards]
                    aieig = aieig / sumaieig
                    if backwards == 0:
                        Eig0coeff[itime] = aieig
                    test = evalues[ieig] * aieig
                    if test > largesteigcoeff:
                        largesteigcoeff = test
                        largestindex = backwards
                Bestindex[itime] = largestindex

                discrep = 0.0
                for vectorindex in range(0, ActualNumberofLocationsUsed):
                    discrep += pow(Checkvector[vectorindex] - PsiTimeSeries[itime, Nonzeromapping[vectorindex]], 2)
                if discrep > 0.01:
                    print('Eigendecomposition Failure ' + str(itime) + ' ' + str(discrep))
                Chi1[itime] = chival
                Chi2[itime] = chival / sumaieig
                Sumai[itime] = sumaieig

                largesteigcoeff = -1.0
                largestindex = -1
                sumaieig = 0.0
                Trailingtimeindex = itime - 3
                if itime > 40:
                    Trailinglimit = Numbereigs[Trailingtimeindex]
                    KeepTrailingaisquared = np.zeros(Trailinglimit, dtype=np.float32)
                    for ieig in range(0, Trailinglimit):
                        aieig = 0.0
                        for vectorindex in range(0, MappedNloc):
                            #              aieig += StoreEigenvectors[Trailingtimeindex,ieig,vectorindex]*PsiTimeSeries[itime,vectorindex]
                            aieig += StoreEigenvectors[Trailingtimeindex, ieig, vectorindex] * StoreEigenvectors[itime,
                                                                                                                 Bestindex[
                                                                                                                     itime], vectorindex]
                        aieig *= aieig
                        sumaieig += aieig
                        KeepTrailingaisquared[ieig] = aieig

                    for ieig in range(0, Trailinglimit):
                        aieig = KeepTrailingaisquared[ieig]
                        aieig = aieig / sumaieig
                        test = Eigenvals[Trailingtimeindex, ieig] * aieig
                        if test > largesteigcoeff:
                            largesteigcoeff = test
                            largestindex = ieig
                    Besttrailingindex[itime] = largestindex

                if itime > 40:  # Calculate eigenvector tracking
                    Leader = StoreEigenvectors[itime, :, :]
                    Trailer = StoreEigenvectors[itime - 3, :, :]
                    StoreEigencorrels[itime, :, :] = np.tensordot(Leader, Trailer, ((1), (1)))
                    StrippedDown = StoreEigencorrels[itime, Bestindex[itime], :]
                    Normingfactor = np.multiply(StrippedDown, StrippedDown).sum()
                    Normingfactor1 = np.multiply(StrippedDown[0:8], StrippedDown[0:8]).sum()
                    Normingfactor2 = np.multiply(StrippedDown[0:30], StrippedDown[0:30]).sum()
                    StoreNormingfactor[itime] = Normingfactor
                    StoreNormingfactor1[itime] = Normingfactor1
                    StoreNormingfactor2[itime] = Normingfactor2

            averagesumai = Sumai.mean()
            Chi1 = np.divide(Chi1, averagesumai)
            print('Bad Correlation Matrices ' + str(countbad))
            print(startbold + startred + 'End Eigen Earthquake 2 '
                  + current_time + ' ' + RunName + resetfonts)


def makeasmalldateplot(figure, ax, Dateaxis):
    plt.rcParams.update({'font.size': 9})
    months = mdates.MonthLocator(interval=2)  # every month
    datemin = np.datetime64(Dateaxis[0], 'M')
    datemax = np.datetime64(Dateaxis[-1], 'M') + np.timedelta64(1, 'M')
    ax.set_xlim(datemin, datemax)

    months_fmt = mdates.DateFormatter('%y-%b')
    locator = mdates.AutoDateLocator()
    locator.intervald['MONTHLY'] = [2]
    formatter = mdates.ConciseDateFormatter(locator)
    #  ax.xaxis.set_major_locator(locator)
    #  ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(months_fmt)

    figure.autofmt_xdate()
    return datemin, datemax


def plotquakeregions(HalfSize, xaxisdates, SetofPlots, Commontitle, ylabel, SetofColors, Startx, ncols):
    numplotted = SetofPlots.shape[1]
    totusedquakes = 0
    for iquake in range(0, numberspecialeqs):
        x_line_index = Numericaldate[iquake]
        if (x_line_index <= Startx) or (x_line_index >= Num_Time - 1):
            continue
        if Specialuse[iquake]:
            totusedquakes += 1
    nrows = math.ceil(totusedquakes / ncols)
    sortedquakes = np.argsort(Numericaldate)

    jplot = 0
    kplot = -1
    for jquake in range(0, numberspecialeqs):
        iquake = sortedquakes[jquake]
        if not Specialuse[iquake]:
            continue
        x_line_annotation = Specialdate[iquake]
        x_line_index = Numericaldate[iquake]
        if (x_line_index <= Startx) or (x_line_index >= Num_Time - 1):
            continue

        kplot += 1
        if kplot == ncols:
            plt.savefig(APPLDIR + '/Outputs/QRegions' + str(jplot) + RunName + '.png ', format='png')
            plt.show()
            kplot = 0
            jplot += 1
        if kplot == 0:
            plt.rcParams["figure.figsize"] = [16, 6]
            figure, axs = plt.subplots(nrows=1, ncols=ncols, squeeze=False)

        beginplotindex = x_line_index - HalfSize
        beginplotindex = max(beginplotindex, Startx)
        endplotindex = x_line_index + HalfSize
        endplotindex = min(endplotindex, Num_Time - 1)

        eachplt = axs[0, kplot]
        ascii = ''
        if Specialuse[iquake]:
            ascii = np.datetime_as_string(Specialdate[iquake]) + ' ' + str(round(Specialmags[iquake], 1)) + ' ' + \
                    Specialeqname[iquake]
        eachplt.set_title(str(iquake) + ' ' + RunName + ' Best Eigenvalue (Black) Trailing (Red) \n' + ascii)
        datemin, datemax = makeasmalldateplot(figure, eachplt, xaxisdates[beginplotindex:endplotindex + 1])
        for curves in range(0, numplotted):
            eachplt.plot(xaxisdates[beginplotindex:endplotindex + 1], SetofPlots[beginplotindex:endplotindex + 1, curves],
                         'o', color=SetofColors[curves], markersize=1)

        ymin, ymax = eachplt.get_ylim()
        if ymax >= 79.9:
            ymax = 82
        eachplt.set_ylim(bottom=-1.0, top=max(ymax, 20))
        eachplt.set_ylabel(ylabel)
        eachplt.set_xlabel('Time')
        eachplt.grid(True)
        eachplt.set_yscale("linear")
        eachplt.axvline(x=x_line_annotation, linestyle='dashed', alpha=1.0, linewidth=2.0, color='red')
        for kquake in range(0, numberspecialeqs):
            if not Specialuse[kquake]:
                continue
            if kquake == iquake:
                continue
            anotherx_line_index = Numericaldate[kquake]
            if (anotherx_line_index < beginplotindex) or (anotherx_line_index >= endplotindex):
                continue
            eachplt.axvline(x=Specialdate[kquake], linestyle='dashed', alpha=1.0, linewidth=1.0, color='purple')
        eachplt.tick_params('x', direction='in', length=15, width=2, which='major')

    plt.savefig(APPLDIR + '/Outputs/QRegions' + str(jplot) + RunName + '.png ', format='png')
    plt.show()


EigenAnalysis = False
if Earthquake and EigenAnalysis:

    UseTopEigenTotal = 40
    FirstTopEigenTotal = 10
    PLTlabels = []
    for ieig in range(0, UseTopEigenTotal):
        PLTlabels.append('Eig-' + str(ieig))

    plt.rcParams["figure.figsize"] = [12, 10]
    figure, ax = plt.subplots()
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    plt.rcParams["figure.figsize"] = [12, 10]
    for ieig in range(0, FirstTopEigenTotal):
        ax.plot(Dateaxis[26:], np.maximum(Eigenvals[26:, ieig], 0.1))

    ax.set_title(RunName + ' Multiplicity Eigenvalues')
    ax.set_ylabel('Eigenvalue')
    ax.set_xlabel('Time')
    ax.set_yscale("log")
    ax.grid(True)
    ax.legend(PLTlabels[0:FirstTopEigenTotal], loc='upper right')
    Addfixedearthquakes(ax, datemin, datemax, ylogscale=True)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    plt.show()

    plt.rcParams["figure.figsize"] = [12, 10]
    figure, ax = plt.subplots()
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    plt.rcParams["figure.figsize"] = [12, 10]
    for ieig in range(FirstTopEigenTotal, UseTopEigenTotal):
        ax.plot(Dateaxis[26:], np.maximum(Eigenvals[26:, ieig], 0.1))

    ax.set_title(RunName + ' Multiplicity Eigenvalues')
    ax.set_ylabel('Eigenvalue')
    ax.set_xlabel('Time')
    ax.set_yscale("linear")
    ax.grid(True)
    ax.legend(PLTlabels[FirstTopEigenTotal:], loc='upper right')
    Addfixedearthquakes(ax, datemin, datemax, ylogscale=False)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    plt.show()

    ShowEigencorrels = False
    if ShowEigencorrels:
        for mastereig in range(0, UseTopEigenTotal):
            figure, ax = plt.subplots()
            plt.rcParams["figure.figsize"] = [12, 8]
            datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
            for ieig in range(0, UseTopEigenTotal):
                alpha = 1.0
                width = 3
                if ieig == mastereig:
                    alpha = 0.5
                    width = 1
                ax.plot(Dateaxis[26:], np.power(StoreEigencorrels[26:, mastereig, ieig], 2), alpha=alpha, linewidth=width)
            ax.set_title(RunName + ' Eigenvalue ' + str(mastereig) + ' Current versus Past Total Correlation')
            ax.set_ylabel('Norm')
            ax.set_xlabel('Time')
            ax.grid(True)
            ax.legend(PLTlabels, loc='upper right')
            Addfixedearthquakes(ax, datemin, datemax, ylogscale=False)
            ax.tick_params('x', direction='in', length=15, width=2, which='major')
            ax.xaxis.set_minor_locator(mdates.YearLocator(1))
            ax.tick_params('x', direction='in', length=10, width=1, which='minor')
            plt.show()

    figure, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [12, 8]
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    alpha = 1.0
    width = 0.5
    ax.plot(Dateaxis[26:], StoreNormingfactor[26:], alpha=alpha, linewidth=width)
    ax.set_title(RunName + ' Eigenvalue Full Norming Factor with Past')
    ax.set_ylabel('Norming Factor')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax, ylogscale=False)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    plt.show()

    figure, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [12, 8]
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    alpha = 1.0
    width = 0.5
    ax.plot(Dateaxis[26:], StoreNormingfactor1[26:], alpha=alpha, linewidth=width)
    ax.set_title(RunName + ' Eigenvalue First 8 Norming Factor with Past')
    ax.set_ylabel('Norming Factor')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax, ylogscale=False)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    plt.show()

    figure, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [12, 8]
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    alpha = 1.0
    width = 0.5
    ax.plot(Dateaxis[26:], StoreNormingfactor2[26:], alpha=alpha, linewidth=width)
    ax.set_title(RunName + ' Eigenvalue First 30 Norming Factor with Past')
    ax.set_ylabel('Norming Factor')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax, ylogscale=False)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    plt.show()

    figure, ax = plt.subplots()
    plt.rcParams["figure.figsize"] = [12, 8]
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    plt.rcParams["figure.figsize"] = [12, 8]
    ax.plot(Dateaxis[26:], Chi1[26:])

    ax.set_title(RunName + ' Correlations Normalized on average over time')
    ax.set_ylabel('Chi1')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    ax.set_yscale("linear")
    plt.show()

    figure, ax = plt.subplots()
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    plt.rcParams["figure.figsize"] = [12, 8]
    ax.plot(Dateaxis[26:], Chi2[26:])

    ax.set_title(RunName + ' Correlations Normalized at each time')
    ax.set_ylabel('Chi2')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    ax.set_yscale("linear")
    plt.show()

    figure, ax = plt.subplots()
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    plt.rcParams["figure.figsize"] = [12, 8]
    norm = np.amax(Chi1[26:])
    Maxeig = 80
    # ax.plot(Dateaxis[26:],Chi1[26:]*Maxeig/norm)
    ax.plot(Dateaxis[26:], 0.5 + np.minimum(Maxeig, Bestindex[26:]), 'o', color='black', markersize=1)
    ax.plot(Dateaxis[26:], np.minimum(Maxeig, Besttrailingindex[26:]), 'o', color='red', markersize=1)

    ax.set_title(RunName + ' Best Eigenvalue (Black) Trailing (Red)')
    ax.set_ylabel('Eig#')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    ax.set_yscale("linear")
    plt.show()

    SetofPlots = np.empty([len(Bestindex), 2], dtype=np.float32)
    SetofPlots[:, 0] = 0.5 + np.minimum(Maxeig, Bestindex[:])
    SetofPlots[:, 1] = np.minimum(Maxeig, Besttrailingindex[:])
    SetofColors = ['black', 'red']
    plotquakeregions(25, Dateaxis, SetofPlots,
                     RunName + ' Best Eigenvalue (Black) Trailing (Red)', 'Eig#', SetofColors, 26, 2)

    plt.rcParams["figure.figsize"] = [12, 8]
    figure, ax = plt.subplots()
    datemin, datemax = makeadateplot(figure, ax, Dateaxis[26:])
    ax.plot(Dateaxis[26:], Eig0coeff[26:], 'o', color='black', markersize=2)
    ymin, ymax = ax.get_ylim()
    ax.plot(Dateaxis[26:], Chi1[26:] * ymax / norm)

    ax.set_title(RunName + ' Fraction Largest Eigenvalue')
    ax.set_ylabel('Eig 0')
    ax.set_xlabel('Time')
    ax.grid(True)
    Addfixedearthquakes(ax, datemin, datemax)
    ax.tick_params('x', direction='in', length=15, width=2, which='major')
    ax.xaxis.set_minor_locator(mdates.YearLocator(1))
    ax.tick_params('x', direction='in', length=10, width=1, which='minor')
    ax.set_yscale("linear")
    plt.show()

"""###End of Earthquake. Reset Timing"""

# Reset Start Date by a year so first entry has a 365 day sample ending at that day and so can be made an input as can all
# lower time intervals
# Do NOT include 2 year or 4 year in input stream
# So we reset start date by one year skipping first 364 daya except to calculate the first one year (and lower limit) observables
# Time indices go from 0 to NumberofTimeunits-1
# Sequence Indices go from Begin to Begin+Tseq-1 where Begin goes from 0 to NumberofTimeunits-1-Tseq
# So Num_Seq = Numberodays-Tseq and Begin has  Num_Seq values

if Earthquake:
    MagnitudeMethod = 0

    SkipTimeUnits = 364
    if Dailyunit == 14:
        SkipTimeUnits = 25
    Num_Time_old = NumberofTimeunits
    NumberofTimeunits = NumberofTimeunits - SkipTimeUnits
    Num_Time = NumberofTimeunits
    InitialDate = InitialDate + timedelta(days=SkipTimeUnits * Dailyunit)
    FinalDate = InitialDate + timedelta(days=(NumberofTimeunits - 1) * Dailyunit)
    print('Skip ' + str(SkipTimeUnits) + ' New dates: ' + InitialDate.strftime("%d/%m/%Y") + ' To '
          + FinalDate.strftime("%d/%m/%Y") + ' days ' + str(NumberofTimeunits * Dailyunit))

    DynamicPropertyTimeSeries = np.empty([Num_Time, Nloc, NpropperTimeDynamic], dtype=np.float32)
    CountNaN = np.zeros(NpropperTimeDynamic, dtype=np.int)
    # Skewtime makes certain propert ENDS at given cell and is the cell itself if size = DailyUnit
    # NpropperTimeDynamicInput is number of input quantities (4)
    # NpropperTimeDynamic adds calculated quantities: accumulated magnitudes (5 magnitude, 6 E^0.25)and Rundle EMA
    # Rundle EMA has same time cut as input time series i.e. no skew but all times are used to find EMA
    #   addRundleEMA = 4
    # RundleLambda = [0.75,1.5,0.75,1.5]
    # RundleSteps = [72,72,36,36]
    # NpropperTimeDynamicCalculated = 11 + addRundleEMA
    # InputIndextogenerateEMA = 3
    # FirstEMAIndex = 15
    SkewTime = [0] * NpropperTimeDynamicInput
    if Dailyunit == 1:
        SkewTime = SkewTime + [22, 45, 91, 182, 364, 0, 22, 45, 91, 182, 364]
    if Dailyunit == 14:
        SkewTime = SkewTime + [1, 3, 6, 12, 25, 0, 1, 3, 6, 12, 25]
    if addRundleEMA > 0:
        SkewTime = SkewTime + [0] * addRundleEMA

    for iprop in range(0, NpropperTimeDynamic):
        addtime = SkipTimeUnits - SkewTime[iprop]

        if iprop >= NpropperTimeDynamic - addRundleEMA:  # Rundle EMA 15-->18
            EMANumber = iprop - NpropperTimeDynamic + addRundleEMA
            for iloc in range(0, Nloc):
                localEMA = MakeEMAMinCT(BasicInputTimeSeries[:, iloc, InputIndextogenerateEMA], RundleSteps[EMANumber],
                                        RundleLambda[EMANumber])
                for itime in range(0, NumberofTimeunits):
                    localval = localEMA[itime + addtime]
                    if np.math.isnan(localval):
                        localval = NaN
                        CountNaN[iprop] += 1
                    DynamicPropertyTimeSeries[itime, iloc, iprop] = localval
            print('Dynamic ' + str(iprop) + ' Rundle EMA ' + str(EMANumber) + ' Prop ' + DynamicNames[
                iprop] + ' Time Shift ' + str(addtime) + ' NaN ' + str(CountNaN[iprop]) + ' Could be in Pred same name')

        else:
            if iprop < NpropperTimeDynamicInput:  # Input Data 0-->3
                for itime in range(0, NumberofTimeunits):
                    for iloc in range(0, Nloc):
                        localval = BasicInputTimeSeries[itime + addtime, iloc, iprop]
                        if np.math.isnan(localval):
                            localval = NaN
                            CountNaN[iprop] += 1
                        DynamicPropertyTimeSeries[itime, iloc, iprop] = localval
                print('Dynamic ' + str(iprop) + ' Input Time Series ' + str(iprop) + ' Prop ' + DynamicNames[
                    iprop] + ' Time Shift ' + str(addtime) + ' NaN ' + str(CountNaN[iprop]) + ' Could be in Pred same name')

            else:
                if iprop < (NpropperTimeDynamic - 6 - addRundleEMA):  # Transformed Magnitude Dynamic 4-8 from Calc 0 to 4
                    icalc = iprop - NpropperTimeDynamicInput
                else:  # Aggregated E^0.25 magnitude Dynamic 9-14 from Calc 9-14
                    icalc = iprop - NpropperTimeDynamicInput + 4
                for itime in range(0, NumberofTimeunits):
                    for iloc in range(0, Nloc):
                        localval = CalculatedTimeSeries[itime + addtime, iloc, icalc]
                        if np.math.isnan(localval):
                            localval = NaN
                            CountNaN[iprop] += 1
                        DynamicPropertyTimeSeries[itime, iloc, iprop] = localval
                print(
                    'Dynamic ' + str(iprop) + ' Calc ' + str(icalc) + ' Prop ' + DynamicNames[iprop] + ' Time Shift ' + str(
                        addtime) + ' NaN ' + str(CountNaN[iprop]) + ' Pred ' + NamespredCalculated[icalc])

    # Predictions
    NewNumTimeSeriesCalculated = NumTimeSeriesCalculated + addRundleEMA
    NewCalculatedTimeSeries = np.empty([Num_Time, Nloc, NewNumTimeSeriesCalculated], dtype=np.float32)
    for iprop in range(0, NumTimeSeriesCalculated):
        NewCalculatedTimeSeries[:, :, iprop] = CalculatedTimeSeries[SkipTimeUnits:Num_Time + SkipTimeUnits, :, iprop]
    if addRundleEMA > 0:
        for iEMA in range(0, addRundleEMA):
            NewCalculatedTimeSeries[:, :, iEMA + NumTimeSeriesCalculated] = DynamicPropertyTimeSeries[:, :,
                                                                            iEMA + FirstEMAIndex]
            NamespredCalculated.append(DynamicNames[iEMA + FirstEMAIndex])
            InputPropertyNames.append(DynamicNames[iEMA + FirstEMAIndex])

    NumTimeSeriesCalculated = NewNumTimeSeriesCalculated
    CalculatedTimeSeries = None
    CalculatedTimeSeries = NewCalculatedTimeSeries
    BasicInputTimeSeries = None
    if GarbageCollect:
        gc.collect()

    print(startbold + startred + 'Predicted NaN values ' + resetfonts)
    for icalc in range(0, NumTimeSeriesCalculated):
        CountofNaN = 0
        for itime in range(0, NumberofTimeunits):
            for iloc in range(0, Nloc):
                localval = CalculatedTimeSeries[itime, iloc, icalc]
                if np.math.isnan(localval):
                    localval = NaN
                    CountofNaN += 1
        print('Predictions(calc) ' + str(icalc) + ' NaN ' + str(CountofNaN) + ' Pred ' + NamespredCalculated[icalc])

    current_time = timenow()
    print(
        startbold + startred + 'Finish Basic Earthquake Setup ' + current_time + ' ' + RunName + ' ' + RunComment + resetfonts)

"""###Set Earthquake Execution Mode"""

if Earthquake:
    SymbolicWindows = True
    Tseq = 26
    if UseTFTModel == True:
        if Dailyunit == 14:
            GenerateFutures = True
            UseFutures = True
    else:
        GenerateFutures = False
        UseFutures = False

"""###Plot Earthquake Images"""

from matplotlib import colors


def plotimages(Array, Titles, nrows, ncols):
    usedcolormap = "YlGnBu"
    plt.rcParams["figure.figsize"] = [16, 6 * nrows]
    figure, axs = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False)
    iplot = 0
    images = []
    norm = colors.Normalize(vmin=fullmin, vmax=fullmax)
    for jplot in range(0, nrows):
        for kplot in range(0, ncols):
            eachplt = axs[jplot, kplot]
            if MapLocation:
                Plotit = np.zeros(OriginalNloc, dtype=np.float32)
                for jloc in range(0, Nloc):
                    Plotit[LookupLocations[jloc]] = Array[iplot][jloc]
                    TwoDArray = np.reshape(Plotit, (40, 60))
            else:
                TwoDArray = np.reshape(Array[iplot], (40, 60))
            extent = (-120, -114, 36, 32)
            images.append(eachplt.imshow(TwoDArray, cmap=usedcolormap, norm=norm, extent=extent))
            eachplt.label_outer()
            eachplt.set_title(Titles[iplot])
            iplot += 1
    figure.colorbar(images[0], ax=axs, orientation='vertical', fraction=.05)
    plt.show()


if Earthquake:
    # DynamicPropertyTimeSeries and CalculatedTimeSeries are dimensione by time 0 ...Num_Time-1
    # DynamicPropertyTimeSeries holds values upto and including that time
    # CalculatedTimeSeries holds values STARTING at that time
    # Plot magnitudes first and them chosen Energy power
    for transformedpointer in range(0, 2):
        localplot1 = 0
        localplot2 = NumTimeSeriesCalculatedBasic
        if transformedpointer == 0:
            fullmin = np.nanmin(CalculatedTimeSeries[:, :, 0:NumTimeSeriesCalculatedBasic])
            fullmax = np.nanmax(CalculatedTimeSeries[:, :, 0:NumTimeSeriesCalculatedBasic])
            fullmin = min(fullmin, np.nanmin(DynamicPropertyTimeSeries[:, :, 0]))
            fullmax = max(fullmax, np.nanmax(DynamicPropertyTimeSeries[:, :, 0]))
            print('Full Magnitude Ranges ' + str(fullmin) + ' ' + str(fullmax))
        else:
            localplot1 = NumTimeSeriesCalculatedBasic
            localplot2 = NumTimeSeriesCalculated
            fullmin = np.nanmin(CalculatedTimeSeries[:, :, localplot1:localplot2])
            fullmax = np.nanmax(CalculatedTimeSeries[:, :, localplot1:localplot2])
            print('Full Energy Transformed Ranges ' + str(fullmin) + ' ' + str(fullmax))
        Num_Seq = NumberofTimeunits - Tseq
        dayindexmax = Num_Seq - Plottingdelay
        Numdates = 4
        denom = 1.0 / np.float64(Numdates - 1)
        for plotdays in range(0, Numdates):
            dayindexvalue = math.floor(0.1 + (plotdays * dayindexmax) * denom)
            if dayindexvalue < 0:
                dayindexvalue = 0
            if dayindexvalue > dayindexmax:
                dayindexvalue = dayindexmax
            dayindexvalue += Tseq
            InputImages = []
            InputTitles = []
            InputImages.append(DynamicPropertyTimeSeries[dayindexvalue, :, 0])
            ActualDate = InitialDate + timedelta(days=dayindexvalue)
            if transformedpointer == 0:
                localmax1 = DynamicPropertyTimeSeries[dayindexvalue, :, 0].max()
                localmin1 = DynamicPropertyTimeSeries[dayindexvalue, :, 0].min()
                InputTitles.append(
                    'Day ' + str(dayindexvalue) + ' ' + ActualDate.strftime("%d/%m/%Y") + ' ' + InputPropertyNames[
                        NpropperTimeStatic] + ' max/min '
                    + str(round(localmax1, 3)) + ' ' + str(round(localmin1, 3)))

            for localplot in range(localplot1, localplot2):
                localmax1 = CalculatedTimeSeries[dayindexvalue, :, localplot].max()
                localmin1 = CalculatedTimeSeries[dayindexvalue, :, localplot].min()
                InputImages.append(CalculatedTimeSeries[dayindexvalue, :, localplot])
                InputTitles.append(
                    'Day ' + str(dayindexvalue) + ' ' + ActualDate.strftime("%d/%m/%Y") + ' ' + NamespredCalculated[
                        localplot] + ' max/min '
                    + str(round(localmax1, 3)) + ' ' + str(round(localmin1, 3)))
            plotimages(InputImages, InputTitles, 5, 2)  # Ten plots fixed
