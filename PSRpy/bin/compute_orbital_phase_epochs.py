#! /usr/bin/python

from PSRpy.parfile import Parfile
from astropy.time import Time
import PSRpy.orbit.elements as orbelem
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(description="A script that computes the epoch(s) for which a certain orbital phase is achieved, given a TEMPO/TEMPO2 solution of a binary pulsar.")
parser.add_argument("epoch", nargs="+", action="store", type=np.float, help="Epoch to refine for one or more orbital phases.")
parser.add_argument("-f", nargs=1, action="store", dest="parfile", required=True, type=str, help="Input TEMPO/TEMPO2 parameter file with binary parameters set.")
parser.add_argument("-o", nargs=1, action="store", dest="observatory", default=["GBT"], type=str, help="Name of observatory (e.g. GBT, Arecibo, CHIME, etc). Default: GBT.")
parser.add_argument("-p", nargs="+", action="store", dest="orbital_phases", default=[0.25], type=np.float, help="One or more orbital phases (0 <= phase < 1) for which to determine epoch(s) of occurence. Default is 0.25, commonly known as superior conjunction.")

args = parser.parse_args()
epoch_to_refine, = args.epoch
parfile, = args.parfile
observatory, = args.observatory
orbital_phases = args.orbital_phases
params = Parfile()
params.read(parfile)

# define parameters describing desired observatory.
observatories = {
    "CHIME"  : [49.32, -119.62],
    "GBT"    : [38.433121, -79.839835]
}

obs_coord_lat, obs_coord_lon = observatories[observatory]
obs_coord = coord.EarthLocation(lat=obs_coord_lat, lon=obs_coord_lon)

# if input parameter file has low-ecentricity model, compute eccentric terms.
if params.BINARY["value"] == "ELL1":
    eps1 = params.EPS1["value"]
    eps2 = params.EPS2["value"]
    print(eps1, eps2)
    print(np.sqrt(eps1**2 + eps2**2))
    params.E["value"] =  np.sqrt(eps1**2 + eps2**2)
    print(params.E)
    params.OM["value"] = np.arctan2(eps1, eps2) * 180 / np.pi % 360
    params.T0["value"] = params.TASC["value"] + params.OM["value"] / 360. * params.PB["value"]

#if hasattr(params, "ECC"):
#    params.E["value"] = params.ECC["value"]

# define desired orbital phases.
tolerance = 1e-10

# loop over phases and compute refined epochs.
for current_phase in orbital_phases:
    current_epoch = epoch_to_refine
    coord_psr = coord.SkyCoord(params.RAJ["value"], params.DECJ["value"], unit=(u.hourangle, u.deg), frame='icrs') 
    #coord_psr = coord.SkyCoord(params.LAMBDA, params.BETA, unit=(u.hourangle, u.deg), 
    #    frame=coord.BarycentricTrueEcliptic())

    for ii in range(100):
        ma = orbelem.anomaly_mean(params.PB["value"], current_epoch, params.T0["value"])
        ea = orbelem.anomaly_eccentric(np.array(ma), params.E["value"])
        ta = orbelem.anomaly_true(ea, params.E["value"])

        # compute function to minimize.
        total = (params.OM["value"] + ta) % 360
        f = total / 360. - current_phase
    
        # now compute the time derivative for the above function, which is just 
        # the time derivative in the true anomaly. 
        # TODO: add component in derivative equation for OMDOT term.
        dea_dt = 1. / params.PB["value"] / (1 - params.E["value"] * np.cos(ea * np.pi / 180))
        ecc_ratio = (1 + params.E["value"]) / (1 - params.E["value"])
        df_dt = np.sqrt(ecc_ratio) / (1 + ecc_ratio * np.tan(ea / 2 * np.pi / 180)**2) * \
                (1 / np.cos(ea / 2 * np.pi / 180)**2) * dea_dt

        dt = f / df_dt
        current_epoch -= dt

        if (np.fabs(dt) < tolerance):
            break

    # now compute toppcentric time.
    t = Time(current_epoch, format='mjd', scale='tdb', location=obs_coord)
    altaz_psr = coord_psr.transform_to(coord.AltAz(obstime=t, location=obs_coord))
    print("For orbital phase = {0:.3f}:".format(current_phase)) 
    print("... MJD (bary): {0:.10f}".format(t.tdb.value)) 
    print("... MJD (UTC):  {0:.10f}".format(t.utc.value))
    print("... Alt. (deg): {0:.10f}".format(altaz_psr.alt.deg))
    print("... Az.  (deg): {0:.10f}".format(altaz_psr.az.deg))
    t.format = 'fits'
    print("... TDB date:   {0}".format(t.tdb))
    print("... UTC date:   {0}".format(t.utc))
