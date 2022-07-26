import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as io
from scipy import interpolate
import glob
import datetime
from pyephem_sunpath.sunpath import sunpos

class uncertainty(object):
    def __init__(self,cruise):
        self.in_dir_main = "" # main directory of your input data 
        self.out_file_main = " " # main directory of your output data
        
        self.in_dir_CalFile = self.in_dir_main + cruise + '/HyperSAS_config/'+cruise
        self.in_dir_L0 = self.in_dir_main+ cruise + '/Processed/L0/'
        self.in_dir_L1 = self.in_dir_main+ cruise + '/Processed/L1/'
        self.in_dir_L2= self.in_dir_main+ cruise + '/Processed/L2/'
        self.out_file_H=self.out_file_main + cruise + '/'+cruise+'_HSAS_unc_v1.txt'
        self.rho_LUT_fn = 'rhoTable_AO2015.txt'
        self.rho_LUT = self.read_rho_LUT(self.rho_LUT_fn)
        self.Es550_Radtran = pd.read_excel('Es_Radtran.xlsx', sheet_name='Es550')

        self.Es_thermal = pd.read_excel('PML_8-010-20-thermal-0258.xlsx')
        self.Unc_Es_Temp0 = self.Es_thermal.iloc[:,[1,3]]
        self.Lt_thermal = pd.read_excel('PML_8-010-20-thermal-0222.xlsx')
        self.Unc_Lt_Temp0 = self.Lt_thermal.iloc[:, [1, 3]]
        self.Li_thermal = pd.read_excel('PML_8-010-20-thermal-0223.xlsx')
        self.Unc_Li_Temp0 = self.Li_thermal.iloc[:, [1, 3]]

        self.pdata = pd.read_excel('non-linearity coefficients.xlsx',sheet_name='258')
        self.Unc_Es_non_linearity_alfa = pdata['uncertainty'].as_matrix()
        self.non_linearity_alfa_wv_Es = pdata['Wavelength'].as_matrix()
        pdata = pd.read_excel('non-linearity coefficients.xlsx',sheetname='222')
        self.Unc_Li_non_linearity_alfa = pdata['uncertainty'].as_matrix()
        self.non_linearity_alfa_wv_Li = pdata['Wavelength'].as_matrix()
        pdata =pd.read_excel('non-linearity coefficients.xlsx',sheetname='223')
        self.Unc_Lt_non_linearity_alfa = pdata['uncertainty'].as_matrix()
        self.non_linearity_alfa_wv_Lt = pdata['Wavelength'].as_matrix()

        self.flag_sl = 1 # 1: include strayligth correctoin
        self.flag_temp = 1 # 1: include temperature correction



    def read_rho_LUT(self,rho_LUT_fn):
        f = open(rho_LUT_fn,'r')
        rho_LUT = {}
        wind_speed = []
        sun_angle = []
        for x in f:
            if 'WIND SPEED =' in x:
                wsp = x.split(';')[0].split('=')[1].strip()
                sa = x.split(';')[1].split('=')[1][0:-1].strip()
                wind_speed.append(float(wsp))
                sun_angle.append(float(sa))
                var_name = 'w'+wsp+'s'+sa
                rho = np.zeros([118,3])
                rho_LUT[var_name] = rho
                i = 0
                continue
            if rho_LUT == {}:
                continue
            rho_LUT[var_name][i] = np.array([float(j) for j in x.split()])
            i += 1
        rho_LUT['wsp'] = np.array(wind_speed)
        rho_LUT['sa'] = np.array(sun_angle)
        return rho_LUT

    def load_calfile(self,fn_CalFile):
        for fn in fn_CalFile:
            if 'LI' in fn:
                self.Cal_LI = io.loadmat(fn)
            if 'LT' in fn:
                self.Cal_LT = io.loadmat(fn)
            if 'ES' in fn:
                self.Cal_ES = io.loadmat(fn)


    def load(self,dir):
        fn_CalFile = glob.glob(self.in_dir_CalFile + '/' + "mean_*.mat")
        self.load_calfile(fn_CalFile)
        fn_L0 = glob.glob(self.in_dir_L0 + dir + '/' + "*Hsas_L0*.mat")
        fn_L1 = glob.glob(self.in_dir_L1 + dir + '/' + "*Hsas_L1*.mat")
        fn_L2 = glob.glob(self.in_dir_L2 + dir + '/' + "*Hsas_L2*.mat")
        self.L0 = io.loadmat(fn_L0[0], squeeze_me=True, struct_as_record=False)
        self.L1 = io.loadmat(fn_L1[0], squeeze_me=True, struct_as_record=False)
        self.L2 = io.loadmat(fn_L2[0], squeeze_me=True, struct_as_record=False)
        # L1 information processing:
        # get sun angle and convert time to datetime format
        Hsas_dt = []
        SA = []
        Hsas_time = self.L1['L1'].time
        lon = self.L1['L1'].gps.lon
        lat = self.L1['L1'].gps.lat
        for i in range(0, len(Hsas_time)):
            Hsas_datetime = (datetime.datetime.fromordinal(int(Hsas_time[i])) + datetime.timedelta(
                seconds=86400 * (Hsas_time[i] - float(int(Hsas_time[i])))) - datetime.timedelta(days=366))
            Hsas_dt.append(Hsas_datetime)
            alt, azm = sunpos(Hsas_datetime, lat[i], lon[i], 0, dst=False)
            SA.append(90 - alt)
        self.L1['L1'].SA = np.array(SA)
        self.L1['L1'].Hsas_dt = np.array(Hsas_dt)
        #L2 information processing:
        Hsas_dt = []
        Hsas_time = self.L2['L2'].time
        for i in range(0, len(Hsas_time)):
            Hsas_datetime = (datetime.datetime.fromordinal(int(Hsas_time[i])) + datetime.timedelta(
                seconds=86400 * (Hsas_time[i] - float(int(Hsas_time[i])))) - datetime.timedelta(days=366))
            Hsas_dt.append(Hsas_datetime)
        self.L2['L2'].Hsas_dt = np.array(Hsas_dt)


    def get_uncerOfrho_geometry(self,wsp,sza,vza,vaa,delta_vza,delta_vaa,Is_showFig=False):
        inx = (np.abs(self.rho_LUT['wsp'] - wsp) + np.abs(self.rho_LUT['sa'] - sza)).argmin()
        key_name = list(self.rho_LUT)[inx]
        print(key_name)
        LUT = self.rho_LUT[key_name]
        tp = np.zeros([13, 3])
        tp[:, 0] = np.zeros(1)
        tp[:, 1] = np.unique(LUT[:, 1])
        tp[:, 2] = LUT[0, 2] * np.ones(13)
        LUT = np.delete(LUT, 0, 0)
        LUT = np.concatenate([tp, LUT])

        # Get look-up table via nearest approach
        # inx = np.where(np.abs(LUT[:, 1] - vaa) == np.min(np.abs(LUT[:, 1] - vaa)))
        # LUT1 = LUT[inx[0], :]
        # get look-up table via interpret
        LUT1 = np.zeros([10,3])
        for i,j in zip(np.array(range(10)),np.unique(LUT[:,0])):
            id = LUT[:,0] == j
            LUT1[i,:]=np.array([j,vaa,np.interp(vaa,LUT[id,1],LUT[id,2])])
        if vza + delta_vza <= np.max(LUT1[:, 0]):
            U1 = np.interp([vza, vza + delta_vza], LUT1[:, 0], LUT1[:, 2])
        else:
            U1 = np.interp([vza - delta_vza, vza], LUT1[:, 0], LUT1[:, 2])

        # Get look-up table via nearest approach
        # inx = np.where(np.abs(LUT[:, 0] - vza) == np.min(np.abs(LUT[:, 0] - vza)))
        # LUT2 = LUT[inx[0], :]
        # get look-up table via interpret
        LUT2 = np.zeros([13,3])
        for i, j in zip(np.array(range(13)),np.unique(LUT[:,1])):
            id = LUT[:,1] == j
            LUT2[i, :] = np.array([vza, j, np.interp(vza, LUT[id, 0], LUT[id, 2])])
        if vaa + delta_vaa <= np.max(LUT2[:, 1]):
            U2 = np.interp([vaa, vaa + delta_vaa], LUT2[:, 1], LUT2[:, 2])
        else:
            U2 = np.interp([vaa - delta_vaa, vaa], LUT2[:, 1], LUT2[:, 2])

        U = np.sqrt((U1[1] - U1[0]) ** 2 + (U2[1] - U2[0]) ** 2)
        f = interpolate.interp2d(LUT[:, 0], LUT[:, 1], LUT[:, 2])
        rho = f(vza, vaa)
        U_per = U / rho

        if Is_showFig:
            fig = plt.figure()
            plt.subplot(2, 1, 1)
            plt.plot(LUT1[:, 0], LUT1[:, 2], '-or')
            plt.xlabel(r'$\theta_v (deg)$')
            plt.ylabel(r'$\rho$')
            plt.subplot(2, 1, 2)
            plt.plot(LUT2[:, 1], LUT2[:, 2], '-or')
            plt.xlabel(r'$\phi_v (deg)$')
            plt.ylabel(r'$\rho$')
            plt.show()
        return {'uncer':U, 'uncer_per':U_per, 'uncer_theta':U1, 'uncer_phi':U2, 'rho':rho}


    def detect_skycondition(self,measured, clearsky, overcast,SA,times, window_length,
                        mean_diff=75, max_diff=75, min_diff = 75,
                        lower_line_length=-5, upper_line_length=10,
                        var_diff=0.005, slope_dev=8,max_cv = 0.05,diff_ref = 10, max_iterations=20,
                        max_diff_detrend=5,low_limit=25,diff_overcast_reff=10,return_components=False):
        """
        Detects clear sky times according to the algorithm developed by Reno
        and Hansen for GHI measurements. The algorithm [1]_ was designed and
        validated for analyzing GHI time series only. Users may attempt to
        apply it to other types of time series data using different filter
        settings, but should be skeptical of the results.
        The algorithm detects clear sky times by comparing statistics for a
        measured time series and an expected clearsky time series.
        Statistics are calculated using a sliding time window (e.g., 10
        minutes). An iterative algorithm identifies clear periods, uses the
        identified periods to estimate bias in the clearsky data, scales the
        clearsky data and repeats.
        Clear times are identified by meeting 5 criteria. Default values for
        these thresholds are appropriate for 10 minute windows of 1 minute
        GHI data.
        Parameters
        ----------
        measured : array or Series
            Time series of measured values.
        clearsky : array or Series
            Time series of the expected clearsky values.
        times : DatetimeIndex
            Times of measured and clearsky values.
        window_length : int
            Length of sliding time window in minutes. Must be greater than 2
            periods.
        mean_diff : float, default 75
            Threshold value for agreement between mean values of measured
            and clearsky in each interval, see Eq. 6 in [1].
        max_diff : float, default 75
            Threshold value for agreement between maxima of measured and
            clearsky values in each interval, see Eq. 7 in [1].
        lower_line_length : float, default -5
            Lower limit of line length criterion from Eq. 8 in [1].
            Criterion satisfied when
            lower_line_length < line length difference < upper_line_length
        upper_line_length : float, default 10
            Upper limit of line length criterion from Eq. 8 in [1].
        var_diff : float, default 0.005
            Threshold value in Hz for the agreement between normalized
            standard deviations of rate of change in irradiance, see Eqs. 9
            through 11 in [1].
        slope_dev : float, default 8
            Threshold value for agreement between the largest magnitude of
            change in successive values, see Eqs. 12 through 14 in [1].
        max_iterations : int, default 20
            Maximum number of times to apply a different scaling factor to
            the clearsky and redetermine clear_samples. Must be 1 or larger.
        return_components : bool, default False
            Controls if additional output should be returned. See below.
        Returns
        -------
        clear_samples : array or Series
            Boolean array or Series of whether or not the given time is
            clear. Return type is the same as the input type.
        components : OrderedDict, optional
            Dict of arrays of whether or not the given time window is clear
            for each condition. Only provided if return_components is True.
        alpha : scalar, optional
            Scaling factor applied to the clearsky_ghi to obtain the
            detected clear_samples. Only provided if return_components is
            True.
        References
        ----------
        .. [1] Reno, M.J. and C.W. Hansen, "Identification of periods of clear
           sky irradiance in time series of GHI measurements" Renewable Energy,
           v90, p. 520-531, 2016.
        Notes
        -----
        Initial implementation in MATLAB by Matthew Reno. Modifications for
        computational efficiency by Joshua Patrick and Curtis Martin. Ported
        to Python by Will Holmgren, Tony Lorenzo, and Cliff Hansen.
        Differences from MATLAB version:
            * no support for unequal times
            * automatically determines sample_interval
            * requires a reference clear sky series instead calculating one
              from a user supplied location and UTCoffset
            * parameters are controllable via keyword arguments
            * option to return individual test components and clearsky scaling
              parameter
        """

        # calculate deltas in units of minutes (matches input window_length units)
        deltas = np.diff(times.values) / np.timedelta64(1, '60s')

        # determine the unique deltas and if we can proceed
        unique_deltas = np.unique(deltas)
        if len(unique_deltas) == 1:
            sample_interval = unique_deltas[0]
        else:
            raise NotImplementedError('algorithm does not yet support unequal '
                                      'times. consider resampling your data.')

        intervals_per_window = int(window_length / sample_interval)

        # generate matrix of integers for creating windows with indexing
        from scipy.linalg import hankel
        H = hankel(np.arange(intervals_per_window),  # noqa: N806
                   np.arange(intervals_per_window - 1, len(times)))

        # convert pandas input to numpy array, but save knowledge of input state
        # so we can return a series if that's what was originally provided
        ispandas = isinstance(measured, pd.Series)
        measured = np.asarray(measured)
        clearsky = np.asarray(clearsky)

        # calculate measurement statistics
        meas_mean = np.mean(measured[H], axis=0)
        meas_max = np.max(measured[H], axis=0)
        meas_min = np.min(measured[H], axis=0)
        meas_diff = np.diff(measured[H], n=1, axis=0)
        meas_slope = np.diff(measured[H], n=1, axis=0) / sample_interval
        # matlab std function normalizes by N-1, so set ddof=1 here
        meas_slope_nstd = np.std(meas_slope, axis=0, ddof=1) / meas_mean
        meas_line_length = np.sum(np.sqrt(
            meas_diff * meas_diff +
            sample_interval * sample_interval), axis=0)

        # calculate clear sky statistics
        clear_mean = np.mean(clearsky[H], axis=0)
        clear_max = np.max(clearsky[H], axis=0)
        clear_min = np.min(clearsky[H], axis=0)
        clear_diff = np.diff(clearsky[H], n=1, axis=0)
        clear_slope = np.diff(clearsky[H], n=1, axis=0) / sample_interval

        # calculate overcast statistics
        overcast_mean = np.mean(overcast[H],axis=0)


        from scipy.optimize import minimize_scalar

        alpha = 1
        for iteration in range(max_iterations):
            clear_line_length = np.sum(np.sqrt(
                alpha * alpha * clear_diff * clear_diff +
                sample_interval * sample_interval), axis=0)

            line_diff = meas_line_length - clear_line_length

            # evaluate comparison criteria
            c1 = np.abs(meas_mean - alpha * clear_mean) < mean_diff
            c2 = np.abs(meas_max - alpha * clear_max) < max_diff
            c22 = np.abs(meas_min - alpha * clear_min) < min_diff
            c3 = (line_diff > lower_line_length) & (line_diff < upper_line_length)
            c4 = meas_slope_nstd < var_diff
            c5 = np.max(np.abs(meas_slope -
                               alpha * clear_slope), axis=0) < slope_dev
            c6 = (clear_mean != 0) & ~np.isnan(clear_mean)
            c7 = np.abs((np.std(measured[H],axis=0, ddof=1)/np.mean(measured[H],axis=0))\
                 - np.std(clearsky[H]*alpha, axis=0, ddof=1) / np.mean(clearsky[H]*alpha, axis=0))< max_cv
            c8 = (np.abs(meas_mean - clear_mean)) < diff_ref
            c9 = meas_mean > low_limit

            diff = (measured[H] - clearsky[H] * alpha)
            diff_mean = np.mean(diff, axis=0)
            diff = diff - diff_mean[None, :]
            c10 = np.max(np.abs(diff), axis=0) < max_diff_detrend

            # corr = np.sum(measured[H]*clearsky[H],axis=0)/(np.sqrt(np.sum(measured[H]**2,axis=0))*np.sqrt(np.sum(clearsky[H]**2,axis=0)))
            # c11 = corr>0.99999999
            #clear_windows = c1 & c2 & c3 & c4 & c5 & c6
            clear_windows =  c1 & c2 & c22 & c7 & c8 & c9 & c6 & c10
            # create array to return
            clear_samples = np.full_like(measured, False, dtype='bool')
            # find the samples contained in any window classified as clear
            clear_samples[np.unique(H[:, clear_windows])] = True

            c70 = np.abs((np.std(measured[H], axis=0, ddof=1) / np.mean(measured[H], axis=0)) \
                   - np.std(overcast[H], axis=0, ddof=1) / np.mean(overcast[H], axis=0)) < max_cv
            c71 = np.abs(meas_mean - overcast_mean) < diff_overcast_reff
            #c72 = np.mean(SA[H],axis=0)<60
            c73 = meas_mean> low_limit * 0.2
            overcast_windows = c70 & c71 & c73 & c6 # & c72
            overcast_samples = np.full_like(measured, False, dtype='bool')
            # find the samples contained in any window classified as overcast
            overcast_samples[np.unique(H[:, overcast_windows])] = True

            # find a new alpha
            previous_alpha = alpha
            clear_meas = measured[clear_samples]
            clear_clear = clearsky[clear_samples]
            print(iteration)
            def rmse(alpha):
                return np.sqrt(np.mean((clear_meas - alpha * clear_clear) ** 2))

            alpha = minimize_scalar(rmse).x
            if round(alpha * 10000) == round(previous_alpha * 10000):
                print('alpha:',alpha,'previous_alpha:',previous_alpha)
                break
        else:
            import warnings
            warnings.warn('failed to converge after %s iterations'
                          % max_iterations, RuntimeWarning)

        # check the distance to reference of Radtran
        if previous_alpha < 0.8:
            clear_samples[:] = False

        std_slope_diff = np.std(meas_slope - clear_slope,axis=0)
        distance = abs(measured - clearsky * previous_alpha)/abs(clearsky * previous_alpha)
        dev_remove = np.mean(distance[clear_samples]) + 2 * np.std(distance[clear_samples])
        flag = np.where((distance > 0.05)&(clear_samples==True))
        for pos in flag[0]:
            try:
                pass #clear_samples[int(pos)-50:int(pos)+50] = False
            except Exception as ex:
                 print(ex)

        diff_std = np.std(measured[H]-clearsky[H]* previous_alpha,axis=0, ddof=1)
        diff = (measured[H]-clearsky[H]* previous_alpha)
        diff_mean = np.mean(diff,axis=0)
        diff = diff - diff_mean[None,:]
        diff = np.max(np.abs(diff),axis=0)
        # be polite about returning the same type as was input
        if ispandas:
            clear_samples = pd.Series(clear_samples, index=times)

        if return_components:
            alpha = previous_alpha
            components = {}
            components['mean_diff_flag'] = c1
            components['max_diff_flag'] = c2
            components['line_length_flag'] = c3
            components['slope_nstd_flag'] = c4
            components['slope_max_flag'] = c5
            components['mean_nan_flag'] = c6
            components['windows'] = clear_windows

            components['mean_diff'] = np.abs(meas_mean - alpha * clear_mean)
            components['max_diff'] = np.abs(meas_max - alpha * clear_max)
            components['min_diff'] = np.abs(meas_min - alpha * clear_min)
            components['line_length'] = meas_line_length - clear_line_length
            components['slope_nstd'] = meas_slope_nstd
            components['slope_max'] = (np.max(
                meas_slope - alpha * clear_slope, axis=0))
            components['relative_mean_diff'] = np.abs(meas_mean - alpha * clear_mean)/meas_mean
            components['relative_max_diff'] = np.abs(meas_max - alpha * clear_max)/meas_max
            components['cv'] = np.std(measured[H],axis=0, ddof=1)/np.mean(measured[H],axis=0)
            components['model_cv'] = np.std(clearsky[H]*alpha, axis=0, ddof=1) / np.mean(clearsky[H]*alpha, axis=0)
            components['diff_std'] = diff_std
            components['relative_mean_diff_ref'] = np.abs(meas_mean - clear_mean)
            components['std_slope_diff'] = std_slope_diff
            components['max_dis_ref'] = diff
            components['meas_mean'] = meas_mean
            tp = meas_mean * 1.
            tp[clear_windows == False] = np.nan
            components['clear_windows'] = tp

            return clear_samples,overcast_samples, pd.DataFrame.from_dict(components), alpha
        else:
            return clear_samples

    #-------Sky condition detection for L1 data
    #--check the sky is clear or no-clear
    def sky_detect(self,is_figure=False):
        wv = self.L1['L1'].wv
        tilt = self.L1['L1'].tilt
        Es = self.L1['L1'].instr.Es.data
        Lt = self.L1['L1'].instr.Lt.data
        Li = self.L1['L1'].instr.Li.data
        lon = self.L1['L1'].gps.lon
        lat = self.L1['L1'].gps.lat
        Hsas_time = self.L1['L1'].time
        Hsas_dt = self.L1['L1'].Hsas_dt
        SA = self.L1['L1'].SA

        tilt_correction = np.cos(np.deg2rad(tilt))
        Es_0 = Es / tilt_correction[:, None]

        Es_Radtran = self.Es550_Radtran.copy()
        Es_clear_model = np.interp(SA, Es_Radtran.iloc[:, 0].tolist(), Es_Radtran.iloc[:, 1].tolist())
        Es_overcast_model = np.interp(SA, Es_Radtran.iloc[:, 0].tolist(), Es_Radtran.iloc[:, 2].tolist())
        total_Es = Es_0[:, 100]  # 550
        total_Es800 = Es_0[:, 225]  # 800

        Es_480 = Es_0[:, 65]  # 480
        Es_600 = Es_0[:, 125]  # 600
        Es_412 = Es_0[:, 31]  # 412
        Es_820 = Es_0[:, 235]  # 820

        start = pd.Timestamp(Hsas_dt[0])
        end = pd.Timestamp(Hsas_dt[-1])
        delta = datetime.timedelta(seconds=5)
        steps = int((end - start) / delta)
        increments = range(0, steps) * np.array([delta] * steps)
        new_Hsas_dt = pd.to_datetime(start + increments)

        x0 = [a.total_seconds() for a in increments]
        x = [a.total_seconds() for a in Hsas_dt - Hsas_dt[0]]
        new_total_Es = np.interp(x0, x, total_Es)
        new_Es_clear_model = np.interp(x0, x, Es_clear_model)
        new_Es_overcast_model = np.interp(x0, x, Es_overcast_model)
        new_Es_480 = np.interp(x0, x, Es_480)
        new_Es_600 = np.interp(x0, x, Es_600)
        new_Es_412 = np.interp(x0, x, Es_412)
        new_Es_820 = np.interp(x0, x, Es_820)

        window_length = 15
        mean_diff = 8
        max_diff = 8
        min_diff = 8
        lower_line_length = -5
        upper_line_length = 200
        var_diff = 0.2
        slope_dev = 100
        max_cv = 0.02
        diff_ref = 12
        max_iterations = 20
        max_diff_detrend = 5
        low_limit = 25
        diff_overcast_reff = 3
        return_components = True
        clear_samples, overcast_samples, components, alph = self.detect_skycondition(new_total_Es, new_Es_clear_model,
                                                                                      new_Es_overcast_model, SA,
                                                                                      new_Hsas_dt, window_length,
                                                                                      mean_diff, max_diff, min_diff,
                                                                                      lower_line_length,
                                                                                      upper_line_length,
                                                                                      var_diff, slope_dev, max_cv,
                                                                                      diff_ref, max_iterations,
                                                                                      max_diff_detrend, low_limit,
                                                                                      diff_overcast_reff,
                                                                                      return_components)

        if is_figure:
            fig = plt.figure(figsize=(6, 8))
            plt.subplot(2, 1, 1)
            plt.plot(new_Hsas_dt, new_total_Es, '.r', label='Ed(550) measurement')
            plt.plot(new_Hsas_dt, new_Es_clear_model, '--b', label='Ed(550) clearsky ref')
            plt.plot(new_Hsas_dt, new_Es_clear_model * alph, '-b', label='clearsky fitting')
            plt.plot(new_Hsas_dt[clear_samples], new_total_Es[clear_samples], '.g', label='detected clearsky')
            plt.plot(new_Hsas_dt, new_Es_overcast_model, '-k', label='Ed(550) overcast ref')
            plt.plot(new_Hsas_dt[overcast_samples], new_total_Es[overcast_samples], '.k', label='detected overcast')
            plt.legend()
            plt.grid()
            plt.xlabel('datetime')
            plt.ylabel(r'$Es(0+) (mw/cm^2/nm)$')
            plt.subplot(2, 1, 2)
            plt.plot(new_Hsas_dt, new_Es_600 / new_Es_480, '-r', label='Es(600)/Es(480)')
            plt.plot(new_Hsas_dt, new_Es_600 / new_Es_412 - 0.5, '-g', label='Es(600)/Es(412) - 0.5')
            plt.plot(new_Hsas_dt, new_Es_820 / new_Es_412 - 0.3, '-b', label='Es(820)/Es(412) - 0.3')
            plt.xlabel('datetime')
            plt.ylabel('band-ratio')
            plt.ylim([0, 1.2])
            plt.grid()
            plt.legend()
            plt.show()

            fig = plt.figure(figsize=(6, 4))
            plt.plot(new_Hsas_dt, new_total_Es, '.r', label=r'$E_s(550)$ measurements')
            plt.plot(new_Hsas_dt[clear_samples], new_total_Es[clear_samples], '.b', label=r'detected clearsky')
            plt.plot(new_Hsas_dt, new_Es_clear_model, '--k', label=r'$E_s(550)$ clearsky reference')
            #plt.plot(new_Hsas_dt, new_Es_clear_model * alph, '-k', label=r'clearsky fitting')
            plt.legend()
            plt.grid()
            plt.xlabel('datetime (MM-DD HH)')
            plt.ylabel(r'$E_s(0^+) (mw/cm^2/nm)$')
            plt.savefig('figs/SkyDetection_288.png', dpi=450, bbox_inches='tight', pad_inches=0.1)

        self.clearsky = {'new_Hsas_dt':new_Hsas_dt,'clear_samples':clear_samples}


##########################----Uncertianty Calculation Equation------######
    def env_unc(self,es_meas_counts,li_meas_counts,lt_meas_counts):
        L1_Hsdt = self.L1['L1'].Hsas_dt
        L2_Hsdt = self.L2['L2'].Hsas_dt
        es_out_unc = np.zeros((len(L2_Hsdt),len(es_meas_counts[0,:])))
        li_out_unc = es_out_unc * 0
        lt_out_unc = es_out_unc * 0
        Rrs_CV = es_out_unc * 0
        rho_out_unc = np.zeros((len(L2_Hsdt)))
        delta_L_out_unc = rho_out_unc * 1.0
        mean_rho = rho_out_unc * np.nan
        mean_delta_L = rho_out_unc * np.nan
        for i in range(len(L2_Hsdt)):
            pos1 = np.where(np.abs(L1_Hsdt - L2_Hsdt[i] + datetime.timedelta(seconds=120)) == np.min(
                np.abs(L1_Hsdt - L2_Hsdt[i] + datetime.timedelta(seconds=120))))
            pos2 = np.where(np.abs(L1_Hsdt - L2_Hsdt[i] - datetime.timedelta(seconds=120)) == np.min(
                np.abs(L1_Hsdt - L2_Hsdt[i] - datetime.timedelta(seconds=120))))
            pos1, pos2 = pos1[0][0], pos2[0][0]
            es_meas_counts_chuck = es_meas_counts[pos1:pos2 + 1, :]
            li_meas_counts_chuck = li_meas_counts[pos1:pos2 + 1, :]
            lt_meas_counts_chuck = lt_meas_counts[pos1:pos2 + 1, :]
            #rho =lt_meas_counts_chuck[:,200:226] / li_meas_counts_chuck[:,200:226] # 750-800nm
            #rho = np.mean(rho,axis=1)
            rho=self.L1['L1'].rho_fitted[pos1:pos2 + 1, 0]
            delta_L = self.L1['L1'].rho_fitted[pos1:pos2 + 1, 1]
            Trs = lt_meas_counts_chuck / es_meas_counts_chuck
            Trs_mean=np.mean(Trs[:,200:226],axis=1) # 750 nm - 800 nm
            sort_id=np.argsort(Trs_mean)
            id=sort_id[0: int(np.floor(len(sort_id) * 0.2))] # 0.2
            es_out_unc[i,:] = np.std(es_meas_counts_chuck[id], axis=0)
            li_out_unc[i, :] = np.std(li_meas_counts_chuck[id], axis=0)
            lt_out_unc[i, :] = np.std(lt_meas_counts_chuck[id], axis=0)
            rho_out_unc[i] = np.std(rho[id],axis=0)
            delta_L_out_unc[i] = np.std(delta_L[id],axis=0)
            mean_rho[i] = np.nanmean(rho)
            mean_delta_L[i] = np.nanmean(delta_L)
            rho0= rho[id]
            delta_L0 = delta_L[id]
            Rrs = ((lt_meas_counts_chuck[id]-rho0[:,None]*li_meas_counts_chuck[id])-delta_L0[:,None]) / es_meas_counts_chuck[id]
            Rrs_CV[i,:] = np.std(Rrs,axis=0)/np.mean(Rrs,axis=0)
        return es_out_unc[:, 25:176]**2, li_out_unc[:, 25:176]**2, lt_out_unc[:, 25:176]**2,rho_out_unc**2,delta_L_out_unc**2,Rrs_CV,mean_rho,mean_delta_L

    def ins_meas_unc_es(self,sig_darkcount_sq,sig_calcoeff_sq, meas_counts, sig_env,unc_per_temp):
        count_term = (sig_darkcount_sq)
        coeff_term = (sig_calcoeff_sq)
        sig_ins_sq = count_term + coeff_term + sig_env

        if self.flag_sl == 0:
            sl_term = (meas_counts * 0) ** 2
        else:
            sl_term = (meas_counts * 0.0025) ** 2  # 0.01
        cos_term = (meas_counts * 0.02) ** 2
        pol_term = (meas_counts * 0.006) ** 2
        if self.flag_temp == 0:
            tempe_corr_term = (meas_counts * 0) ** 2
        else:
            tempe_corr_term = (meas_counts * unc_per_temp) ** 2
            # tempe_corr_term = (meas_counts * 0.002) ** 2
        # -----non-linearity correction uncertainty
        alfa = np.interp(self.L2['L2'].wv[25:176], self.non_linearity_alfa_wv_Es, self.Unc_Es_non_linearity_alfa )
        gain_ES = np.interp(self.L1['L1'].wv[25:176], self.Cal_ES['wv'][0], self.Cal_ES['gain']['mean'][0][0][0])
        non_linearity = (-(meas_counts/gain_ES )** 2 * alfa * gain_ES) ** 2
        non_linearity_pc = non_linearity ** 0.5 / meas_counts
        sig_ins_sq = sig_ins_sq + sl_term + cos_term + pol_term + tempe_corr_term + non_linearity
        pc_meas_var = count_term / sig_ins_sq
        pc_cal_coef = coeff_term / sig_ins_sq
        pc_ct = cos_term / sig_ins_sq
        pc_sl = sl_term / sig_ins_sq
        pc_pt = pol_term / sig_ins_sq
        pc_env = sig_env / sig_ins_sq

        return sig_ins_sq, pc_meas_var, pc_cal_coef, pc_ct, pc_sl, pc_pt, pc_env, non_linearity_pc

    def ins_meas_unc_lt(self,sig_darkcount_sq,sig_calcoeff_sq, meas_counts, sig_env,unc_per_temp):
        count_term = (sig_darkcount_sq)
        coeff_term = (sig_calcoeff_sq)
        sig_ins_sq = count_term + coeff_term + sig_env
        if self.flag_sl == 0:
            sl_term = (meas_counts * 0) ** 2
        else:
            sl_term = (meas_counts * 0.005) ** 2  # 0.04
        cos_term = 0
        pol_term = (meas_counts * 0.013) ** 2
        if self.flag_temp == 0:
            tempe_corr_term = (meas_counts * 0) ** 2
        else:
            tempe_corr_term = (meas_counts * unc_per_temp) ** 2
            # tempe_corr_term = (meas_counts * 0.002) ** 2
        # -----non-linearity correction uncertainty
        alfa = np.interp(self.L2['L2'].wv[25:176], self.non_linearity_alfa_wv_Lt, self.Unc_Lt_non_linearity_alfa)
        gain_LT = np.interp(self.L1['L1'].wv[25:176], self.Cal_LT['wv'][0], self.Cal_LT['gain']['mean'][0][0][0])
        non_linearity = (- (meas_counts/gain_LT) ** 2 * alfa * gain_LT) ** 2
        non_linearity_pc = non_linearity ** 0.5 /meas_counts
        #print(non_linearity_pc)
        sig_ins_sq = sig_ins_sq + sl_term + cos_term + pol_term + tempe_corr_term + non_linearity
        pc_meas_var = count_term / sig_ins_sq
        pc_cal_coef = coeff_term / sig_ins_sq
        pc_ct = cos_term / sig_ins_sq
        pc_sl = sl_term / sig_ins_sq
        pc_pt = pol_term / sig_ins_sq
        pc_env = sig_env / sig_ins_sq

        return sig_ins_sq, pc_meas_var, pc_cal_coef, pc_ct, pc_sl, pc_pt, pc_env, non_linearity_pc


    def ins_meas_unc_li(self,sig_darkcount_sq,sig_calcoeff_sq, meas_counts, sig_env,unc_per_temp):
        count_term = (sig_darkcount_sq)
        coeff_term = (sig_calcoeff_sq)
        sig_ins_sq = count_term + coeff_term + sig_env

        if self.flag_sl == 0:
            sl_term = (meas_counts * 0) ** 2
        else:
            sl_term = (meas_counts * 0.0025) ** 2  # 0.01
        cos_term = 0
        pol_term = (meas_counts * 0.013) ** 2
        if self.flag_temp == 0:
            tempe_corr_term = (meas_counts * 0) ** 2
        else:
            tempe_corr_term = (meas_counts * unc_per_temp) ** 2
            # tempe_corr_term = (meas_counts * 0.002) ** 2
        # -----non-linearity correction uncertainty
        alfa = np.interp(self.L2['L2'].wv[25:176], self.non_linearity_alfa_wv_Li, self.Unc_Li_non_linearity_alfa)
        gain_LI = np.interp(self.L1['L1'].wv[25:176], self.Cal_LI['wv'][0], self.Cal_LI['gain']['mean'][0][0][0])
        non_linearity = (-(meas_counts/gain_LI) ** 2 * alfa * gain_LI) ** 2
        non_linearity_pc = non_linearity ** 0.5 / meas_counts
        #print(non_linearity_pc)
        sig_ins_sq = sig_ins_sq + sl_term + cos_term + pol_term + tempe_corr_term + non_linearity
        pc_meas_var = count_term / sig_ins_sq
        pc_cal_coef = coeff_term / sig_ins_sq
        pc_ct = cos_term / sig_ins_sq
        pc_sl = sl_term / sig_ins_sq
        pc_pt = pol_term / sig_ins_sq
        pc_env = sig_env / sig_ins_sq

        return sig_ins_sq, pc_meas_var, pc_cal_coef, pc_ct, pc_sl, pc_pt, pc_env, non_linearity_pc

    def darkcount_unc(self):
        wv0 = self.L1['L1'].wv[25:176]
        wv = self.Cal_ES['wv'][0]
        offset_ES_unc = self.Cal_ES['offset']['unc'][0][0][0]
        offset_LT_unc = self.Cal_LT['offset']['unc'][0][0][0]
        offset_LI_unc = self.Cal_LI['offset']['unc'][0][0][0]
        gain_ES = self.Cal_ES['gain']['mean'][0][0][0]
        gain_LT = self.Cal_LT['gain']['mean'][0][0][0]
        gain_LI = self.Cal_LI['gain']['mean'][0][0][0]

        offset_ES_unc = np.interp(wv0,wv,offset_ES_unc)
        offset_LT_unc = np.interp(wv0,wv,offset_LT_unc)
        offset_LI_unc = np.interp(wv0,wv,offset_LI_unc)
        gain_ES = np.interp(wv0, wv, gain_ES)
        gain_LT = np.interp(wv0, wv, gain_LT)
        gain_LI = np.interp(wv0, wv, gain_LI)

        sig_darkcount_ES = offset_ES_unc * gain_ES
        sig_darkcount_LT = offset_LT_unc * gain_LT
        sig_darkcount_LI = offset_LI_unc * gain_LI

        return sig_darkcount_ES**2, sig_darkcount_LT**2, sig_darkcount_LI**2

    def coeff_unc(self,Es_meas_count,Lt_meas_count,Li_meas_count):
        wv0 = self.L1['L1'].wv[25:176]
        wv = self.Cal_ES['wv'][0]
        gain_ES_unc = self.Cal_ES['gain']['unc'][0][0][0]
        gain_LT_unc = self.Cal_LT['gain']['unc'][0][0][0]
        gain_LI_unc = self.Cal_LI['gain']['unc'][0][0][0]
        gain_ES = self.Cal_ES['gain']['mean'][0][0][0]
        gain_LT = self.Cal_LT['gain']['mean'][0][0][0]
        gain_LI = self.Cal_LI['gain']['mean'][0][0][0]

        gain_ES_unc = np.interp(wv0,wv,gain_ES_unc)
        gain_LT_unc = np.interp(wv0,wv,gain_LT_unc)
        gain_LI_unc = np.interp(wv0,wv,gain_LI_unc)
        gain_ES = np.interp(wv0, wv, gain_ES)
        gain_LT = np.interp(wv0, wv, gain_LT)
        gain_LI = np.interp(wv0, wv, gain_LI)

        sig_constant = 0.0 #--2% other calibration uncertainty
        sig_coeff_ES = Es_meas_count * ((gain_ES_unc/gain_ES)+sig_constant)
        sig_coeff_LT = Lt_meas_count * ((gain_LT_unc/gain_LT)+sig_constant)
        sig_coeff_LI = Li_meas_count * ((gain_LI_unc/gain_LI)+sig_constant)

        return sig_coeff_ES**2, sig_coeff_LT**2, sig_coeff_LI**2

    def calc_cc_terms(self,lt_in, li_in, es_in,rho_para):
        L1_Hsdt = self.L1['L1'].Hsas_dt
        L2_Hsdt = self.L2['L2'].Hsas_dt
        cc_lt_li = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_lt_es = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_li_es = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_li_rho = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_lt_rho = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_li_dlt = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_lt_dlt = np.zeros((len(L2_Hsdt),len(lt_in[0,:])))
        cc_rho_dlt = np.zeros(len(L2_Hsdt))
        for i in range(len(L2_Hsdt)):
            print(i)
            pos1 = np.where(np.abs(L1_Hsdt - L2_Hsdt[i] + datetime.timedelta(seconds=120)) == np.min(
                np.abs(L1_Hsdt - L2_Hsdt[i] + datetime.timedelta(seconds=120))))
            pos2 = np.where(np.abs(L1_Hsdt - L2_Hsdt[i] - datetime.timedelta(seconds=120)) == np.min(
                np.abs(L1_Hsdt - L2_Hsdt[i] - datetime.timedelta(seconds=120))))
            pos1, pos2 = pos1[0][0], pos2[0][0]
            for j in np.arange(0, np.shape(lt_in)[1]):
                cc_lt_li[i,j] = np.corrcoef(lt_in[pos1:pos2 + 1, j], li_in[pos1:pos2 + 1, j])[0, 1]
                cc_lt_es[i,j] = np.corrcoef(lt_in[pos1:pos2 + 1, j], es_in[pos1:pos2 + 1, j])[0, 1]
                cc_li_es[i,j] = np.corrcoef(li_in[pos1:pos2 + 1, j], es_in[pos1:pos2 + 1, j])[0, 1]
                cc_li_rho[i,j] = np.corrcoef(li_in[pos1:pos2 + 1, j], rho_para[pos1:pos2 + 1, 0])[0, 1]
                cc_lt_rho[i,j] = np.corrcoef(lt_in[pos1:pos2 + 1, j], rho_para[pos1:pos2 + 1, 0])[0, 1]
                cc_li_dlt[i,j] = np.corrcoef(li_in[pos1:pos2 + 1, j], rho_para[pos1:pos2 + 1, 1])[0, 1]
                cc_lt_dlt[i, j] = np.corrcoef(lt_in[pos1:pos2 + 1, j], rho_para[pos1:pos2 + 1, 1])[0, 1]
            cc_rho_dlt[i] = np.corrcoef(rho_para[pos1:pos2 + 1, 0], rho_para[pos1:pos2 + 1, 1])[0, 1]
        return cc_lt_li, cc_lt_es, cc_li_es, cc_li_rho, cc_lt_rho, cc_li_dlt, cc_lt_dlt, cc_rho_dlt



    def calc_combined(self,es_samp, sig_es, lt_samp, sig_lt, li_samp, sig_li, rho_samp,pc_sig_rho, sig_rho_sq,sig_delta_L_sq,cc_lt_li, cc_lt_es, cc_li_es,cc_li_rho, cc_lt_rho, cc_li_dlt, cc_lt_dlt, cc_rho_dlt, rrs_in):
        rrs_lt_term = ((1 / es_samp) * sig_lt) ** 2
        rrs_li_term = ((-rho_samp / es_samp) * sig_li) ** 2
        rrs_es_term = (((-lt_samp + rho_samp * li_samp) / es_samp ** 2) * sig_es) ** 2
        # rrs_rho_term = ((-li_samp / es_samp) * (rho_samp * pc_sig_rho)) ** 2
        rrs_rho_term = (-li_samp / es_samp)**2 * sig_rho_sq
        rrs_delta_L_term = (-1 / es_samp)**2 * sig_delta_L_sq
        c_lt_li_term = 2 * (1 / es_samp) * (-rho_samp / es_samp) * sig_lt * sig_li * cc_lt_li
        c_lt_es_term = 2 * (1 / es_samp) * ((-lt_samp + rho_samp * li_samp) / es_samp ** 2) * sig_lt * sig_es * cc_lt_es
        c_li_es_term = 2 * ((-rho_samp / es_samp) * (
                    (-lt_samp + rho_samp * li_samp) / es_samp ** 2)) * sig_li * sig_es * cc_li_es
        c_li_rho_term = 2 * (-rho_samp / es_samp) * (- li_samp / es_samp) * sig_li * (sig_rho_sq ** 0.5) * cc_li_rho
        c_lt_rho_term = 2 * (1 / es_samp) * (- li_samp / es_samp) * sig_lt * (sig_rho_sq ** 0.5) * cc_lt_rho
        c_li_dlt_term = 2 * (-rho_samp / es_samp) * (-1 / es_samp) * sig_li * (sig_delta_L_sq ** 0.5) * cc_li_dlt
        c_lt_dlt_term = 2 * (1 / es_samp) * (-1 / es_samp) * sig_lt * (sig_delta_L_sq ** 0.5) * cc_lt_dlt
        c_rho_dlt_term = 2 * (-li_samp / es_samp) * (-1 / es_samp) * (sig_rho_sq ** 0.5) * (sig_delta_L_sq ** 0.5) * cc_rho_dlt
        sig_rrs = rrs_lt_term + rrs_li_term + rrs_es_term + rrs_rho_term + rrs_delta_L_term + c_lt_li_term + c_lt_es_term + c_li_es_term \
        + c_li_rho_term + c_lt_rho_term + c_li_dlt_term + c_lt_dlt_term + c_rho_dlt_term



        return sig_rrs, rrs_lt_term, rrs_li_term, rrs_es_term, rrs_rho_term, rrs_delta_L_term,c_lt_li_term, c_lt_es_term , c_li_es_term, \
               c_li_rho_term, c_lt_rho_term, c_li_dlt_term, c_lt_dlt_term, c_rho_dlt_term

    def rho_unc(self,clearsky, clearsky_dt, L2_Hsas_dt):
        x0 = [a.total_seconds() for a in L2_Hsas_dt - clearsky_dt[0].to_pydatetime()]
        x = [a.total_seconds() for a in clearsky_dt - clearsky_dt[0].to_pydatetime()]
        L2_clearsky = np.interp(x0, x, clearsky)
        self.L2_clearsky = L2_clearsky
        rho = L2_clearsky * np.nan
        rho[L2_clearsky==1] = 0.3
        rho[L2_clearsky != 1] = 0.5
        return rho

    def cal_unc_nlw(self,t_sig_rrs):
        exLwn = self.L2['L2'].exLwn.data[:,25:176]
        Rrs = self.L2['L2'].Rrs.data[:,25:176]
        Lw = self.L2['L2'].Lw.data[:,25:176]
        Lw0 = self.L2['L2'].Lw0.data[:,25:176]

        sf= exLwn/Rrs
        #---ur is set as 0.5, note (Lw0-Lw)/Lw0 with a convinent way
        sig_nLw = (sf**2) * t_sig_rrs + (exLwn**2) * (0.5*(Lw0-Lw)/Lw0)**2
        return sig_nLw





    def perform_unc_summation(self):
        t_in = self.L2
        wv = self.L2['L2'].wv[25:176]

        #----Chla estimation with OC4
        Rrs = self.L2['L2'].Rrs.data
        nLw = self.L2['L2'].exLwn.data
        Rrs1=np.nanmax(Rrs[:,[47,70,80]],axis=1) # 443>490>510
        Rrs2= Rrs[:,100] #555
        r = np.log10(Rrs1/Rrs2)
        p = 0.3272 -2.994*r +2.7218 * r**2 -1.2259 * r**3 -0.5683 * r**4
        chla = 10**p
        ##### Run through each step of uncertainty derivation #####

        # 1. Calculate the uncertainty in the calibration coefficient

        t_lt_unc_sq = (t_in['L2'].instr.Lt.data[:, 25:176] * 0.01) ** 2
        t_li_unc_sq = (t_in['L2'].instr.Li.data[:, 25:176] * 0.01) ** 2
        t_es_unc_sq = (t_in['L2'].instr.Es.data[:, 25:176] * 0.01) ** 2

        # print np.shape

        # 2. Calculate the measurement uncertainty for each measurement (lt, li, es)

        # Load dummy variable to get size correct for initialising arrays

        shape_init = t_in['L2'].instr.Es.data[:, 25:176]

        t_sig_lt = np.zeros((np.shape(shape_init)))
        t_sig_li = np.zeros((np.shape(shape_init)))
        t_sig_es = np.zeros((np.shape(shape_init)))

        pc_meas_lt = np.zeros((np.shape(shape_init)))
        pc_dark_lt = np.zeros((np.shape(shape_init)))
        pc_cal_lt = np.zeros((np.shape(shape_init)))
        pc_ct_lt = np.zeros((np.shape(shape_init)))
        pc_sl_lt = np.zeros((np.shape(shape_init)))
        pc_pt_lt = np.zeros((np.shape(shape_init)))
        pc_env_lt = np.zeros((np.shape(shape_init)))

        pc_meas_li = np.zeros((np.shape(shape_init)))
        pc_dark_li = np.zeros((np.shape(shape_init)))
        pc_cal_li = np.zeros((np.shape(shape_init)))
        pc_ct_li = np.zeros((np.shape(shape_init)))
        pc_sl_li = np.zeros((np.shape(shape_init)))
        pc_pt_li = np.zeros((np.shape(shape_init)))
        pc_env_li = np.zeros((np.shape(shape_init)))

        pc_meas_es = np.zeros((np.shape(shape_init)))
        pc_dark_es = np.zeros((np.shape(shape_init)))
        pc_cal_es = np.zeros((np.shape(shape_init)))
        pc_ct_es = np.zeros((np.shape(shape_init)))
        pc_sl_es = np.zeros((np.shape(shape_init)))
        pc_pt_es = np.zeros((np.shape(shape_init)))
        pc_env_es = np.zeros((np.shape(shape_init)))
        pc_non_linearity_lt = np.zeros((np.shape(shape_init)))
        pc_non_linearity_li = np.zeros((np.shape(shape_init)))
        pc_non_linearity_es = np.zeros((np.shape(shape_init)))

        t_es = t_in['L2'].instr.Es.data[:, 25:176]
        t_lt = t_in['L2'].instr.Lt.data[:, 25:176]
        t_li = t_in['L2'].instr.Li.data[:, 25:176]
        t_time = t_in['L2'].time

        # Calculate uncertainty associated with environmental factors (based on moving std with time)
        env_term_es,env_term_li,env_term_lt,unc_rho, unc_delta_L,Rrs_CV,mean_rho, mean_delta_L = self.env_unc(
                            self.L1['L1'].instr.Es.data,
                            self.L1['L1'].instr.Li.data,
                            self.L1['L1'].instr.Lt.data)

        sig_darkcount_Es, sig_darkcount_Lt, sig_darkcount_Li=self.darkcount_unc()
        sig_coeff_Es, sig_coeff_Lt, sig_coeff_Li = self.coeff_unc(t_es, t_lt, t_li)

        #---Temperature Uncertainty:
        #(1)---Uncertainty factors
        self.Unc_Es_Temp = np.interp(wv,self.Unc_Es_Temp0.iloc[:,0].tolist(),self.Unc_Es_Temp0.iloc[:,1].tolist())
        self.Unc_Lt_Temp = np.interp(wv, self.Unc_Lt_Temp0.iloc[:, 0].tolist(), self.Unc_Lt_Temp0.iloc[:, 1].tolist())
        self.Unc_Li_Temp = np.interp(wv, self.Unc_Li_Temp0.iloc[:, 0].tolist(), self.Unc_Li_Temp0.iloc[:, 1].tolist())
        #-Temperature correction factor
        self.c_Es_Temp = np.interp(self.L2['L2'].wv, self.Es_thermal.wl.tolist(), self.Es_thermal.c.tolist())
        self.c_Lt_Temp = np.interp(self.L2['L2'].wv, self.Lt_thermal.wl.tolist(), self.Lt_thermal.c.tolist())
        self.c_Li_Temp = np.interp(self.L2['L2'].wv, self.Li_thermal.wl.tolist(), self.Li_thermal.c.tolist())
        #self.c_temp = (self.c_Es_Temp + self.c_Lt_Temp + self.c_Li_Temp) / 3
        #(2)---Air temperature interplation to L2 time
        try:
            temp_pcb_Es = np.interp(self.L2['L2'].time,self.L0['L0'].time,self.L0['L0'].instr.Es.temp_pcb)
            temp_pcb_Li = np.interp(self.L2['L2'].time,self.L0['L0'].time,self.L0['L0'].instr.Li.temp_pcb)
            temp_pcb_Lt = np.interp(self.L2['L2'].time,self.L0['L0'].time,self.L0['L0'].instr.Lt.temp_pcb)
            #(3)---Unc_Temp In Percentage
            self.unc_Es_Temp_inPer = np.abs(-(temp_pcb_Es[:,None]-21)*self.Unc_Es_Temp[None,:])*0.001
            self.unc_Lt_Temp_inPer = np.abs(-(temp_pcb_Lt[:, None] - 21) * self.Unc_Lt_Temp[None, :])*0.001
            self.unc_Li_Temp_inPer = np.abs(-(temp_pcb_Li[:, None] - 21) * self.Unc_Li_Temp[None, :])*0.001
        except Exception as e:
            print(e,'Probaly air temperture unavailable!')
            self.temp_pcb = np.ones(np.shape(self.L2['L2'].time)) * 25
            self.unc_Es_Temp_inPer = np.abs((self.temp_pcb[:,None] - 21) * self.Unc_Es_Temp[None, :]) * 0.001
            self.unc_Lt_Temp_inPer = np.abs((self.temp_pcb[:,None] - 21) * self.Unc_Lt_Temp[None, :]) * 0.001
            self.unc_Li_Temp_inPer = np.abs((self.temp_pcb[:,None] - 21) * self.Unc_Li_Temp[None, :]) * 0.001
       

        for i in range(0, len(t_time)):
            t_sig_lt[i, :], pc_meas_lt[i, :], pc_cal_lt[i, :], pc_ct_lt[i, :], pc_sl_lt[i, :], pc_pt_lt[i,:], pc_env_lt[i, :], pc_non_linearity_lt[i, :] = self.ins_meas_unc_lt(
                sig_darkcount_Lt,sig_coeff_Lt[i,:], t_in['L2'].instr.Lt.data[i, 25:176], env_term_lt[i, :],self.unc_Lt_Temp_inPer[i,:])
            t_sig_li[i, :], pc_meas_li[i, :], pc_cal_li[i, :], pc_ct_li[i, :], pc_sl_li[i, :], pc_pt_li[i,:], pc_env_li[i,:], pc_non_linearity_li[i, :] = self.ins_meas_unc_li(
                sig_darkcount_Li,sig_coeff_Li[i,:], t_in['L2'].instr.Li.data[i, 25:176], env_term_li[i, :],self.unc_Li_Temp_inPer[i,:])
            t_sig_es[i, :], pc_meas_es[i, :], pc_cal_es[i, :], pc_ct_es[i, :], pc_sl_es[i, :], pc_pt_es[i,:], pc_env_es[i,:], pc_non_linearity_es[i, :] = self.ins_meas_unc_es(
                sig_darkcount_Es,sig_coeff_Es[i,:], t_in['L2'].instr.Es.data[i, 25:176], env_term_es[i, :],self.unc_Es_Temp_inPer[i,:])
        try:
            pc_sig_rho = self.rho_unc(self.clearsky['clear_samples'], self.clearsky['new_Hsas_dt'], self.L2['L2'].Hsas_dt)
        except Exception as error:
            print('Cannot save clearsky info')
            print(error)

        t_sig_lt = np.sqrt(t_sig_lt)
        t_sig_li = np.sqrt(t_sig_li)
        t_sig_es = np.sqrt(t_sig_es)

        # 3. Calculate the uncertainty in the resulting Rrs

        # Determine correlation coefficients for covarying terms (i.e. all light based measurements!!)

        cc_lt_li, cc_lt_es, cc_li_es, cc_li_rho, cc_lt_rho, cc_li_dlt, cc_lt_dlt, cc_rho_dlt =\
            self.calc_cc_terms(self.L1['L1'].instr.Lt.data[:, 25:176],
            self.L1['L1'].instr.Li.data[:, 25:176],
            self.L1['L1'].instr.Es.data[:, 25:176],
            self.L1['L1'].rho_fitted)

        # Perform summation of uncertainties

        t_sig_rrs = np.zeros((np.shape(shape_init)))
        t_sig_rrs_lt = np.zeros((np.shape(shape_init)))
        t_sig_rrs_li = np.zeros((np.shape(shape_init)))
        t_sig_rrs_es = np.zeros((np.shape(shape_init)))
        t_sig_rrs_rho = np.zeros((np.shape(shape_init)))
        t_sig_rrs_delta_L = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_lt_li = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_lt_es = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_li_es = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_li_rho = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_lt_rho = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_li_dlt = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_lt_dlt = np.zeros((np.shape(shape_init)))
        t_sig_rrs_c_rho_dlt = np.zeros((np.shape(shape_init)))



        for i in range(0, len(t_time)):
            t_sig_rrs[i, :], t_sig_rrs_lt[i,:], t_sig_rrs_li[i,:], t_sig_rrs_es[i,:], t_sig_rrs_rho[i,:],\
            t_sig_rrs_delta_L[i,:],t_sig_rrs_c_lt_li[i,:], t_sig_rrs_c_lt_es[i,:] , t_sig_rrs_c_li_es[i,:],\
            t_sig_rrs_c_li_rho[i,:], t_sig_rrs_c_lt_rho[i,:],t_sig_rrs_c_li_dlt[i,:],t_sig_rrs_c_lt_dlt[i,:],t_sig_rrs_c_rho_dlt[i,:]= self.calc_combined(
                t_in['L2'].instr.Es.data[i, 25:176], t_sig_es[i, :], t_in['L2'].instr.Lt.data[i, 25:176],
                t_sig_lt[i, :], t_in['L2'].instr.Li.data[i, 25:176], t_sig_li[i, :], t_in['L2'].rho_fitted[i, 0],0,unc_rho[i],unc_delta_L[i],
                cc_lt_li[i,:], cc_lt_es[i,:], cc_li_es[i,:],
                cc_li_rho[i,:], cc_lt_rho[i,:], cc_li_dlt[i,:], cc_lt_dlt[i,:], cc_rho_dlt[i],
                t_in['L2'].Rrs.data[i, 25:176])

        # Calculate uncertainties as percent - OR IS THIS CORRECT? I THINK THIS IS (error comes out in units, so need to calculate as a percentage of each measurement?

        unc = (np.sqrt(t_sig_rrs) / t_in['L2'].Rrs.data[:, 25:176]) * 100

        # mean_perwv=np.mean(unc[:,0:101],0)
        mean_persamp = np.mean(unc, 1)

        #--Calculated uncertainty for nlw:
        t_sig_nLw=self.cal_unc_nlw(t_sig_rrs)

        self.unc_non_linearity_Lt_inPer = pc_non_linearity_lt
        self.unc_non_linearity_Li_inPer = pc_non_linearity_li
        self.unc_non_linearity_Es_inPer = pc_non_linearity_es
        #-----yield result
        self.uncer_result = {'wv':wv,'mean_persamp':mean_persamp,'unc':unc,'Hsas_dt':self.L2['L2'].Hsas_dt,
                             'Lw':self.L2['L2'].Lw.data, 'Lw0':self.L2['L2'].Lw0.data,
                             'Rrs':t_in['L2'].Rrs.data[:, 25:176],'sig_rrs':t_sig_rrs,
                             'nLw':self.L2['L2'].exLwn.data[:,25:176],'sig_nLw':t_sig_nLw,
                             't_sig_es':t_sig_es,'t_sig_lt':t_sig_lt,'t_sig_li':t_sig_li,
                             'sig_rrs_lt':t_sig_rrs_lt,'sig_rrs_li':t_sig_rrs_li,
                             'sig_rrs_es':t_sig_rrs_es,'sig_rrs_rho':t_sig_rrs_rho,'sig_rrs_delta_L':t_sig_rrs_delta_L,
                             'sig_rrs_c_lt_li':t_sig_rrs_c_lt_li,'sig_rrs_c_lt_es':t_sig_rrs_c_lt_es,
                             'sig_rrs_c_li_es':t_sig_rrs_c_li_es,'sig_rrs_c_li_rho':t_sig_rrs_c_li_rho,
                             'sig_rrs_c_lt_rho':t_sig_rrs_c_lt_rho,'sig_rrs_c_li_dlt': t_sig_rrs_c_li_dlt,
                             'sig_rrs_c_lt_dlt': t_sig_rrs_c_lt_dlt,'sig_rrs_c_rho_dlt':t_sig_rrs_c_rho_dlt,
                             'cc_lt_li':cc_lt_li,'cc_lt_es':cc_lt_es,'cc_li_es':cc_li_es,
                             'env_term_lt':env_term_lt,'env_term_li':env_term_li,'env_term_es':env_term_es,
                             't_es':t_es,'t_li':t_li,'t_lt':t_lt,'Rrs_CV':Rrs_CV,
                             'sig_darkcount_Es':sig_darkcount_Es,'sig_darkcount_Lt':sig_darkcount_Lt,'sig_darkcount_Li':sig_darkcount_Li,
                             'sig_coeff_Es':sig_coeff_Es, 'sig_coeff_Lt':sig_coeff_Lt,'sig_coeff_Li':sig_coeff_Li,'sig_rho':unc_rho ** 0.5,
                             'air_temp':self.air_temp,'mean_rho': mean_rho,'mean_delta_L': mean_delta_L,'time':t_time}


    def run(self,dir):
        self.load(dir)
        self.sky_detect(is_figure=True)