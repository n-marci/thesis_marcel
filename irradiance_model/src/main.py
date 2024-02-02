import toml
# import requests
import json
import pathlib
import matplotlib.pyplot as plt
import pvlib as pv
import pandas as pd

def instantiate():
    # load the configuration file
    with open('eei-tower.toml', 'r') as file:
        inputs = toml.load(file)

    # setup location
    cauerstreet = pv.location.Location(
        inputs['lighting']['natural']['latitude'],
        inputs['lighting']['natural']['longitude'],
        inputs['lighting']['natural']['timezone'],
        inputs['lighting']['natural']['altitude'],
        'Cauerstraße 7'
    )

    tmy, tmy_meta = pv.iotools.read_epw(inputs['weather']['datafile'], coerce_year=2024)
    fig, ax = plt.subplots(nrows=2, ncols=1)

    # calculate the solar position for the DateTimeIndex of the TMY file
    timepoints = tmy.index - pd.Timedelta('30min')          # data has to be shifted because we need the middle of the hour for solar position
    # timepoints2 = timepoints.resample['D']
    solarpos = cauerstreet.get_solarposition(timepoints)
    sun_rise_set_transit = cauerstreet.get_sun_rise_set_transit(timepoints).resample('D').min()
    solarpos.index += pd.Timedelta('30min')                 # shift back

    # sun_rise_set_transit.map(lambda x: x.hour + x.minute / 60).plot()
    # print(sun_rise_set_transit.loc['2024-01-01':'2024-01-01'])

    # calculate the irradiance on the specific orientation of the farm provided by eei-tower.toml file
    poa = pv.irradiance.get_total_irradiance(               # plane of array (poa)
        surface_tilt=90,
        surface_azimuth=inputs['lighting']['natural']['orientation'],
        dni=tmy['dni'],
        ghi=tmy['ghi'],
        dhi=tmy['dhi'],
        solar_zenith=solarpos['apparent_zenith'],
        solar_azimuth=solarpos['azimuth'],
        model='isotropic'                                   # Perez or Hay Davies are other models but require additional weather inputs
    )
    # poa.loc['2024-01-01':'2024-01-01'].plot()
    # print(poa.loc['2024-01-01':'2024-01-01'])

    # calculate natural irradiance for plants
    par = poa['poa_global'] * inputs['lighting']['natural']['solar-to-par']       # photosynthetic active radiation (par)
    hppfd = par.map(lambda x: min(x, inputs['plant']['max-ppfd'])) * 60 * 60 / 1000000    # hourly ppfd - photosynthetic photon flux density (ppfd)
    dli = hppfd.resample('D').sum()                         # daily light integral (dli)
    dli.resample('M').mean().plot(ax=ax[0])

    # draw figure
    ax[0].set_title("Monthly mean of daily light integral by natural irradiance")
    ax[0].set_xlabel("TMY")
    ax[0].set_ylabel("DLI in mol/m2 per day")

    # heat = poa['poa_global'].map(lambda x: x * 60 * 60).resample('D').sum()
    # tempchange = heat.map(lambda x: x / 1210).map(lambda x: x / (24 * 60)).resample('M').mean().plot()

    # hourly['poa_global'].resample('D').mean().plot()
    # dli.resample('D').sum().resample('M').mean().plot.bar()
    # horizon.resample('M').mean()
    # print(horizon)
    # print(tmy['temp_air'])

    # calculate the needed supplemental lighting based on the natural irradiance
    # dli.resample('D').sum().resample('M').mean().plot.bar()
    missing_dli = dli.map(lambda x: inputs['plant']['optimal-dli'] - x)
    hli_opt = inputs['plant']['optimal-ppfd'] * 60 * 60 / 1000000            # hourly light integral (hli) at optimal ppfd
    hli = inputs['lighting']['artificial']['ppf'] * 60 * 60 / 1000000            # hourly light integral (hli) at actual ppfd
    hours_of_suppl_light = missing_dli / hli
    # hours_of_suppl_light.map(lambda x: x * inputs['lighting']['artificial']['ppf']).resample('M').mean().plot()
    led_energy_con = hours_of_suppl_light.map(lambda x: x * inputs['lighting']['artificial']['led-watts']).resample('M').mean().plot(ax=ax[1])

    # draw figure
    ax[1].set_title("Monthly mean of daily energy consumption of LEDs")
    ax[1].set_xlabel("TMY")
    ax[1].set_ylabel("Energy consumption in Wh/m2 per day")
    # dli.loc['2024-06-04':'2024-06-04'].plot()

    # total yearly energy consumption
    yearly_energy_con = hours_of_suppl_light.map(lambda x: x * inputs['lighting']['artificial']['led-watts']).agg("sum") 
    print(f"Yearly energy consumption: {yearly_energy_con} Wh")
    
    # timeframe = pd.date_range(start=inputs['timeframe']['start'], end=inputs['timeframe']['end'], freq=inputs['timeframe']['step'], tz=cauerstreet.tz)
    # get solar position
    # solarpos = cauerstreet.get_solarposition(timeframe)
    # solarpos[['elevation']].plot()
    # solarpos.plot()
    # airmass = cauerstreet.get_airmass(timeframe, solarpos)
    # airmass.plot()
    # cauerclearsky = cauerstreet.get_clearsky(timeframe)
    # cauerclearsky.plot();
    plt.show()

    # TODO AT SOME POINT calculate λ, U-value - see https://de.wikipedia.org/wiki/W%C3%A4rmed%C3%A4mmung#Geb%C3%A4ude

    # Calculate heat generation from natural lighting
    # 1. ???

    # Calculate heat generation from leds
    # 1. use efficiency to calculate heat generation

    # Calculate energy consumption for heat pump
    # 1. get heat pump efficiency
    # 2. get TMY data to get temps and stuff for the time of year

    # Calculate total energy consumption
    # 1. get led consumption
    # 2. get heat pump consumption
    # 3. get consumption of other consumers - how???
    
instantiate()

def archive():
    
    # # the hourly data from pvgis is strange...
    # # not using it for now
    # hourly, hourly_in, hourly_meta = pv.iotools.get_pvgis_hourly(
    #     cauerstreet.latitude,
    #     cauerstreet.longitude,
    #     components=False,
    #     surface_tilt=90,
    #     surface_azimuth=inputs['lighting']['natural']['orientation'],
    #     start=2016,
    #     end=2016
    # )
    # print(hourly.loc['2016-01-01':'2016-01-01'])            # sun only was there for 3 hours on this day?

    # get weather data
    tmy, months, input, metadata = pv.iotools.get_pvgis_tmy(
        cauerstreet.latitude,
        cauerstreet.longitude,
        startyear=2005,
        endyear=2016
    )
    print(tmy)
    tmy['ghi'].plot()

    pass

