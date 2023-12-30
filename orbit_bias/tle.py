import numpy as np
from orbit_predictor.sources import EtcTLESource # needed for tle -> orbit info

RE = 6371.2
d2r = np.pi/180

def get_satellite_state(dates, tlefile, id):
    """ get satellite position and velocity for a at given time from tle file

    Parameters
    ----------
    dates: datetime or list/array of datetimes
        The time(s) of interest. Can be iterable, in which case
        a corresponding number of return parameters are produced
    tlefile: string
        name and path of tle file
    id: string
        ID of spacecraft, must match information in tle file

    Returns
    -------
    r: array
        radius at satellite position
    theta: array
        colatitudes/polar angle at satellite position (degrees)
        geocentric coordinates
    phi: array
        longitude at satellite position (degrees) geocentric 
        coordinates
    v: array
        (N, 3) array of satellite velocity in geocentric ENU coordinates.
        N is the number of datetimes passed to the function
    """

    if not hasattr(dates, '__iter__'):
        assert type(dates) == type(dt.datetime(2000, 1, 1))
        dates = [dates]

    source = EtcTLESource(filename = tlefile)
    predictor = source.get_predictor(id)

    # get position in ECEF coords:
    r_ecef = np.array([predictor.get_position(date).position_ecef for date in dates]).T
    v_ecef = np.array([predictor.get_position(date).velocity_ecef for date in dates]).T

    # calculate spherical/geocentric coords:
    r     = np.linalg.norm(r_ecef, axis = 0)
    phi   = np.arctan2(r_ecef[1], r_ecef[0]) / d2r
    theta = np.arccos(r_ecef[2] / r) / d2r

    # calculate ENU components of velocity:
    R = np.stack((np.vstack((-                      np.sin(phi * d2r),                        np.cos(phi * d2r), np.zeros_like(phi) )),
                  np.vstack((-np.cos(theta * d2r) * np.cos(phi * d2r), -np.cos(theta * d2r) * np.sin(phi * d2r), np.sin(theta * d2r))),
                  np.vstack(( np.sin(theta * d2r) * np.cos(phi * d2r),  np.sin(theta * d2r) * np.sin(phi * d2r), np.cos(theta * d2r))) ), axis = 0)
    v_enu = np.einsum('jik, ik->kj', R, v_ecef)

    return tuple(map(np.squeeze, [r, theta, phi, v_enu]))
