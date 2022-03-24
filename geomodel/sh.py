""" tools that are useful for spherical harmonic analysis

    SHkeys       -- class to contain n and m - the indices of the spherical harmonic terms
    nterms       -- function which calculates the number of terms in a 
                    real expansion of a poloidal (internal + external) and toroidal expansion 
    get_legendre -- calculate associated legendre functions - with option for Schmidt semi-normalization
"""
import numpy as np
#from pyamps.mlt_utils import mlon_to_mlt
d2r = np.pi/180

class SHkeys(object):

    def __init__(self, Nmax, Mmax):
        """ container for n and m in spherical harmonics

            keys = SHkeys(Nmax, Mmax)

            keys will behave as a tuple of tuples, more or less
            keys['n'] will return a list of the n's
            keys['m'] will return a list of the m's
            keys[3] will return the fourth n,m tuple

            keys is also iterable

        """

        self.Nmax = Nmax
        self.Mmax = Mmax
        keys = []
        for n in range(self.Nmax + 1):
            for m in range(self.Mmax + 1):
                keys.append((n, m))

        self.keys = tuple(keys)
        self.make_arrays()

    def __getitem__(self, index):
        if index == 'n':
            return [key[0] for key in self.keys]
        if index == 'm':
            return [key[1] for key in self.keys]

        return self.keys[index]

    def __iter__(self):
        for key in self.keys:
            yield key

    def __len__(self):
        return len(self.keys)

    def __repr__(self):
        return ''.join(['n, m\n'] + [str(key)[1:-1] + '\n' for key in self.keys])[:-1]

    def __str__(self):
        return ''.join(['n, m\n'] + [str(key)[1:-1] + '\n' for key in self.keys])[:-1]

    def setNmin(self, nmin):
        """ set minimum n """
        self.keys = tuple([key for key in self.keys if key[0] >= nmin])
        self.make_arrays()
        return self

    def MleN(self):
        """ set m <= n """
        self.keys = tuple([key for key in self.keys if abs(key[1]) <= key[0]])
        self.make_arrays()
        return self

    def Mge(self, limit):
        """ set m <= n """
        self.keys = tuple([key for key in self.keys if abs(key[1]) >= limit])
        self.make_arrays()
        return self

    def NminusModd(self):
        """ remove keys if n - m is even """
        self.keys = tuple([key for key in self.keys if (key[0] - abs(key[1])) % 2 == 1])
        self.make_arrays()
        return self

    def NminusMeven(self):
        """ remove keys if n - m is odd """
        self.keys = tuple([key for key in self.keys if (key[0] - abs(key[1])) % 2 == 0])
        self.make_arrays()
        return self

    def negative_m(self):
        """ add negative m to the keys """
        keys = []
        for key in self.keys:
            keys.append(key)
            if key[1] != 0:
                keys.append((key[0], -key[1]))
        
        self.keys = tuple(keys)
        self.make_arrays()
        
        return self


    def make_arrays(self):
        """ prepare arrays with shape ( 1, len(keys) )
            these are used when making G matrices
        """

        if len(self) > 0:
            self.m = np.array(self)[:, 1][np.newaxis, :]
            self.n = np.array(self)[:, 0][np.newaxis, :]
        else:
            self.m = np.array([])[np.newaxis, :]
            self.n = np.array([])[np.newaxis, :]



def nterms(NT = 0, MT = 0, NVi = 0, MVi = 0, NVe = 0, MVe = 0):
    """ return number of coefficients in an expansion in real spherical harmonics of
        toroidal magnetic potential truncated at NT, MT
        poloidal magnetic potential truncated at NVi, MVi for internal sources
        poloidal magnetic potential truncated at NVe, MVe for external sources
    """

    return len(SHkeys(NT , MT ).setNmin(1).MleN().Mge(0)) + \
           len(SHkeys(NT , MT ).setNmin(1).MleN().Mge(1)) + \
           len(SHkeys(NVe, MVe).setNmin(1).MleN().Mge(0)) + \
           len(SHkeys(NVe, MVe).setNmin(1).MleN().Mge(1)) + \
           len(SHkeys(NVi, MVi).setNmin(1).MleN().Mge(0)) + \
           len(SHkeys(NVi, MVi).setNmin(1).MleN().Mge(1))



def get_legendre(nmax, mmax, theta, schmidtnormalize = True, negative_m = False, minlat = 0, keys = None):
    """ calculate associated Legendre functions 

        nmax             -- maximum total wavenumber
        mmax             -- maximum zonal wavenumber
        theta            -- colatitude in degrees (not latitude!), with N terms
        schmidtnormalize -- True if Schmidth seminormalization is wanted, False otherwise
        negative_m       -- True if you want the functions for negative m (complex expansion)
        keys             -- pass SHkeys object to return an array instead of a dict
                            default None


        returns:
          P, dP -- dicts of legendre functions, and derivatives, with wavenumber tuple as keys
        
        or, if keys != None:
          PdP   -- array of size N, 2*M, where M is the number of terms. The first half of the 
                   columns are P, and the second half are dP


        algorithm from "Spacecraft Attitude Determination and Control" by James Richard Wertz
        
        could be unstable for large nmax...

        KML 2016-04-22

    """

    theta = theta.flatten()[:, np.newaxis]

    P = {}
    dP = {}
    sinth = np.sin(d2r*theta)
    costh = np.cos(d2r*theta)

    if schmidtnormalize:
        S = {}
        S[0, 0] = 1.

    # initialize the functions:
    for n in range(nmax +1):
        for m in range(nmax + 1):
            P[n, m] = np.zeros_like(theta, dtype = np.float64)
            dP[n, m] = np.zeros_like(theta, dtype = np.float64)

    P[0, 0] = np.ones_like(theta, dtype = np.float64)
    P[0, 0][np.abs(90 - theta) < minlat] = 0
    for n in range(1, nmax +1):
        for m in range(0, min([n + 1, mmax + 1])):
            # do the legendre polynomials and derivatives
            if n == m:
                P[n, n]  = sinth * P[n - 1, m - 1]
                dP[n, n] = sinth * dP[n - 1, m - 1] + costh * P[n - 1, n - 1]
            else:

                if n == 1:
                    Knm = 0.
                    P[n, m]  = costh * P[n -1, m]
                    dP[n, m] = costh * dP[n - 1, m] - sinth * P[n - 1, m]

                elif n > 1:
                    Knm = ((n - 1)**2 - m**2) / ((2*n - 1)*(2*n - 3))
                    P[n, m]  = costh * P[n -1, m] - Knm*P[n - 2, m]
                    dP[n, m] = costh * dP[n - 1, m] - sinth * P[n - 1, m] - Knm * dP[n - 2, m]

            if schmidtnormalize:
                # compute Schmidt normalization
                if m == 0:
                    S[n, 0] = S[n - 1, 0] * (2.*n - 1)/n
                else:
                    S[n, m] = S[n, m - 1] * np.sqrt((n - m + 1)*(int(m == 1) + 1.)/(n + m))


    if schmidtnormalize:
        # now apply Schmidt normalization
        for n in range(1, nmax + 1):
            for m in range(0, min([n + 1, mmax + 1])):
                P[n, m]  *= S[n, m]
                dP[n, m] *= S[n, m]

    if negative_m:
        for n  in range(1, nmax + 1):
            for m in range(0, min([n + 1, mmax + 1])):
                P[n, -m]  = -1.**(-m) * factorial(n-m)/factorial(n+m) *  P[n, m]
                dP[n, -m] = -1.**(-m) * factorial(n-m)/factorial(n+m) * dP[n, m]


    if keys is None:
        return P, dP
    else:
        Pmat  = np.hstack(tuple(P[key] for key in keys))
        dPmat = np.hstack(tuple(dP[key] for key in keys)) 
    
        return np.hstack((Pmat, dPmat))


def get_G(r, latitude, longitude, N, M, R = 6371.2e3):
    """ Make matrix that relates magnetic field components to Gauss coefficients.
        The matrix will have shape (3J x K) where J is the number coordinates given, 
        and K is the number of Gauss coefficients implied by N and M (spherical harmonic
        degree and order). The 3J rows correspond to the North, East, and Center components,
        stacked vertically.

        The Gauss coefficients (matrix columns) are organized with cos terms first, and sin 
        terms second. A set of spherical harmonic keys are also returned, which indicates
        which term the coefficient corresponds to
    """


    th  = (90 - latitude[:, np.newaxis]) * np.pi / 180
    phi =       longitude[:, np.newaxis]  * np.pi / 180
    r   =       r[:, np.newaxis]
    sin_keys = SHkeys(N, M).setNmin(1).MleN().Mge(1)
    cos_keys = SHkeys(N, M).setNmin(1).MleN().Mge(0)

    P, dP = get_legendre(N, M, 90 - latitude)
    P_cos, dP_cos = np.hstack([P[key] for key in cos_keys]), np.hstack([dP[key] for key in cos_keys])
    P_sin, dP_sin = np.hstack([P[key] for key in sin_keys]), np.hstack([dP[key] for key in sin_keys])

    rr = R / r

    BN_cos_G =   rr ** (cos_keys.n + 2)                    * dP_cos              * np.cos(cos_keys.m * phi) 
    BE_cos_G =   rr ** (cos_keys.n + 2)                    *  P_cos * cos_keys.m * np.sin(cos_keys.m * phi) / np.sin(th)
    BC_cos_G = - rr ** (cos_keys.n + 2) * (cos_keys.n + 1) *  P_cos              * np.cos(cos_keys.m * phi)

    BN_sin_G =   rr ** (sin_keys.n + 2)                    * dP_sin              * np.sin(sin_keys.m * phi) 
    BE_sin_G = - rr ** (sin_keys.n + 2)                    *  P_sin * sin_keys.m * np.cos(sin_keys.m * phi) / np.sin(th)
    BC_sin_G = - rr ** (sin_keys.n + 2) * (sin_keys.n + 1) *  P_sin              * np.sin(sin_keys.m * phi)

    GN = np.hstack((BN_cos_G, BN_sin_G))
    GE = np.hstack((BE_cos_G, BE_sin_G))
    GC = np.hstack((BC_cos_G, BC_sin_G))
    G = np.vstack((GN, GE, GC))
    
    return G, cos_keys, sin_keys
