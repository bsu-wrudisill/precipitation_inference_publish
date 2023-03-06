import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import sys
import xarray as xr
from scipy import optimize

# ---THIS IS THE FILE OF STATIC FILE PATHS ; NEEDS TO BE EDITED ----
from model_paths import  base_path, data_path

# parameters
def DayLength_Calc(dayOfYear, lat=40.0):
    """Computes the length of the day (the time between sunrise and
    sunset) given the day of the year and latitude of the location.
    Function uses the Brock model for the computations.
    For more information see, for example,
    Forsythe et al., "A model comparison for daylength as a
    function of latitude and day of year", Ecological Modelling,
    1995.
    Parameters
    ----------
    dayOfYear : int
        The day of the year. 1 corresponds to 1st of January
        and 365 to 31st December (on a non-leap year).
    lat : float
        Latitude of the location in degrees. Positive values
        for north and negative for south.
    Returns
    -------
    d : float
        Daylength in hours.
    """
    latInRad = np.deg2rad(lat)
    declinationOfEarth = 23.45*np.sin(np.deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) <= -1.0:
        return 24.0
    elif -np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth)) >= 1.0:
        return 0.0
    else:
        hourAngle = np.rad2deg(np.arccos(-np.tan(latInRad) * np.tan(np.deg2rad(declinationOfEarth))))
        return 2.0*hourAngle/15.0

    # https://gist.github.com/anttilipp/ed3ab35258c7636d87de6499475301ce


class ForcingModel:
    def __init__(self):
        self.discharge_unit_mm = None
        self.watershed_area = None
        self.hypsometry = None
        self.snotel_elev = None
        self.DailyTempForcing = None
        self.DailyPrecipForcing = None
        self.DayLength = None
        self.days = None

    def read(self):
        self.PrepareDEM()
        self.PrepareDischarge()
        self.PrepareDEM()
        self.PrepareSnowObs()
        self.ApplyDayLength(self.days)

    def PrepareDEM(self):
        ## ---  Get the basin as among other things ---##

        # ds = xr.open_rasterio("/Volumes/Transcend/ASOdata/DEMs/3mdem_upsample_50m_clipped_to_east.tif").sel(band=1)
        dem = data_path + "/3mdem_upsample_50m_clipped_to_east.tif"
        ds = xr.open_rasterio(dem).sel(band=1)
        flat_ds = ds.values.flatten()
        flat_ds_east = flat_ds[np.where(flat_ds > 0)]
        flat_ds_east_sort = np.sort(flat_ds_east)
        east_river_area_m2 = 748983000.0 # m2
        self.watershed_area = east_river_area_m2
        self.hypsometry = flat_ds_east_sort
        ds.close()

    def PrepareDischarge(self):
        ## --- Gather the USGS stream observations ----#
#         columns = ['usgs', 'id', 'date', 'timezone', 'discharge', 'flag']
#         usgs_df = pd.read_csv('./data/USGS_almont.csv', skiprows=41, sep='\t', names = columns)
#         usgs_df['date'] = pd.to_datetime(usgs_df.date)
#         usgs_df.set_index('date', inplace=True)

#         # there might be missing data -- reinterpolate to correct time (adds timesteps)
#         # go from hourly --> daily
#         usgs_df = usgs_df.resample("D").mean()
#         usgs_df = usgs_df.interpolate()

#         del usgs_df['id']
#         usgs_df['discharge_m3s'] = usgs_df['discharge']* 0.0283  # convert to m3/s
#         usgs_df['discharge_unit_mm'] = usgs_df['discharge_m3s'] * (24*60*60)/self.watershed_area*1000
# #       usgs_df = usgs_df.loc[start:end]
#         self.discharge_unit_mm = usgs_df['discharge_unit_mm']

          # this one ahas already been processed...
          #usgs_df = pd.read_csv('/Users/williamrudisill/Documents/Conceptual_Runoff_Model/data/USGS_almont.csv')
          usgs_df = pd.read_csv(data_path + "/USGS_almont.csv")
          usgs_df = usgs_df.set_index('date')
          self.discharge_unit_mm = usgs_df['discharge_unit_mm']


    def PrepareSnowObs(self):
        ## ---  Gather the Snotel Forcings  ----##
        #database = "/Volumes/Transcend/EastRiverClimatePaper/Snotel/CO_snotel.db"
        database = data_path + "/CO_snotel.db"

        dbengine = "sqlite:///{}".format(database)

        # select the snotel station X
        df = pd.read_sql_table('380', dbengine)
        df_schofield = pd.read_sql_table('737', dbengine)

        # butte snotel
        df['Date'] = pd.to_datetime(df['Date'])
        df.set_index('Date', inplace=True)

        # scho snotel
        df_schofield['Date'] = pd.to_datetime(df_schofield['Date'])
        df_schofield.set_index('Date', inplace=True)


        # comptue/fix the forcing data
        DailyPrecipForcing = df.IncPrecip.fillna(0) * 25.4
        DailyTempForcing = ((df.TavgF-32.0)*5/9)
        DailyTempForcing = DailyTempForcing.fillna(DailyTempForcing.mean())

        # butte, schofield
        snotel_elev = 3096.768, 3261.36


        ####################################################################################
        # 0) Remove Bad dates
        ####################################################################################
        df = df.loc[pd.to_datetime("1986-01-01"):]
        df_schofield = df_schofield.loc[pd.to_datetime("1986-01-01"):]

        # Filter out the big jumps in temperature
        df.TavgF = df.TavgF.where(abs(df.TavgF.diff()) < 60, np.nan)
        df_schofield.TavgF = df_schofield.TavgF.where(abs(df_schofield.TavgF.diff()) < 50, np.nan)

        # SUSPECT DATES --- BUTTE
        bd1 = (pd.to_datetime("2002-07-01"), pd.to_datetime("2003-10-01"))
        bd2 = (pd.to_datetime("1989-08-20"), pd.to_datetime("1989-08-21"))


        # Remove the data from the bad date periods
        df.TavgF.loc[bd1[0]:bd1[1]] = np.nan
        df.TavgF.loc[bd2[0]:bd2[1]] = np.nan

        ####################################################################################
        # 1) Adjust the schof. snotel by the adiabatic lapse rate
        # The other snotel is at a different elevation. adjust temperature accordingly
        ####################################################################################

        # adjust the schofield temp by the adiabatic lapse rate
        lapse_rate_f = -.00356 # defgrees F per 1000m
        df_schofield['TavgF_adjusted_schof'] = df_schofield.TavgF + lapse_rate_f * (3261.36 - 3096.768)


        ####################################################################################
        # 2)  Replace bad Butte values with values from the adjusted Schofield !!!
        ####################################################################################

        replace_nan = lambda a, b: b if np.isnan(a) else a
        fixed = df.TavgF.combine(df_schofield['TavgF_adjusted_schof'], replace_nan)
        df['fixed_v1'] = fixed

        # Cut out the bad/missing data from the combined dataset (could do this eralier )

        ####################################################################################
        # 3) Fit a sinusoidal function to fill in data for regions where niether station is good. should be small ish
        ####################################################################################

        # interpolate bad values (must be done for curve fitting)
        fixed_int = fixed.interpolate()

        # take the weekly average
        fixed_roll = fixed_int.rolling(7).mean() # 7 day rolling mean average. smooths out the sinish wave


        # Remove the nans
        valid = ~(np.isnan(fixed_roll))
        ydata = fixed_roll.values[valid]

        # create x data
        xdata = np.linspace(0, len(ydata), len(ydata))

        # fucntion to be fitted.
        def sin_func(x, a, w, c, d):
            return a*np.sin(w*x + c) + d

        # initial guess of the parameters
        p0 = [30., 2*np.pi/365., 0., 30.]

        # fit the curve
        p, pcov = optimize.curve_fit(sin_func, xdata, ydata, p0=p0)

        ## make a plot to compare
        ## plt.plot(xdata, ydata)
        ## plt.plot(xdata, sin_func(xdata, *p))
        df['fitted_data'] = np.nan # ugh more nans.
        df.fitted_data.iloc[7:] = sin_func(xdata, *p) # we lose the first week when we do the rolling mean (and mask nans)
        ######################################
        # 4) Sub in the fitted values
        ######################################
        replace_nan = lambda a, b: b if np.isnan(a) else a
        fixed_v2 = df.fixed_v1.combine(df.fitted_data, replace_nan)

        #######################################
        # 5) Create output timeseries
        ######################################
        df = df.iloc[7:]
        df['final_T2'] = fixed_v2

        # get the number of times...
        ntimes = len(df.index)
        days = df.index.dayofyear.values

        # Now assign the final things to
        DailyPrecipForcing = df.IncPrecip.fillna(0) * 25.4
        DailyTempForcing = ((df.final_T2-32.0)*5/9)  # Convert to degC
        self.DailyTempForcing = DailyTempForcing
        self.DailyPrecipForcing = DailyPrecipForcing
        self.days = days
        self.SWE = df.SWE*25.4 # convert to mm

    def ApplyDayLength(self, days):
        LenOfDayHr = [DayLength_Calc(d, lat=38.0) for d in days]
        self.DayLength = LenOfDayHr
        return LenOfDayHr


    def __call__(self, start, end):
        sub_dtf = self.DailyTempForcing.loc[start:end]
        sub_pcp = self.DailyPrecipForcing.loc[start:end]
        sub_q = self.discharge_unit_mm.loc[start:end]
        sub_doy = sub_pcp.index.dayofyear.values
        sub_day_len = self.ApplyDayLength(sub_doy)
        sub_swe = self.SWE.loc[start:end]
        return sub_dtf, sub_pcp, sub_q, sub_day_len, sub_swe
        # sub_dlen = self.DailyTempForcing.loc[start:end]

if __name__ == "__main__":
    start = "2017-10-01"
    end = "2018-09-30"
    f = ForcingModel()
    f.read()
    daily_temp, daily_precip, q, day_len_hrs, sub_swe = f(start, end)

    # # print(len(daily_temp))
    # # print(len(daily_precip))
    # # print(len(q))
    # # print(len(day_len_hrs))

    # # np.save("./data/daily_precip.npy", daily_precip)
    # # np.save("./data/daily_q_observed.npy", q)
    # # np.save("./data/daily_temp.npy", daily_temp)
    # # np.save("./data/day_len_hrs.npy", day_len_hrs)
    # np.save("./data/daily_swe_mm.npy", sub_swe)





