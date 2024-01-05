from astropy.table import Table
import numpy as np
from astropy import constants as const
from astropy import units as u

wave_eff = {'W1': 33526.0,
            'W2': 46028.0,
            'G': 6230.06}
w1_nu_effective = (const.c / (33526.0 * u.AA)).to(u.Hz).value
w2_nu_effective = (const.c / (46028.0 * u.AA)).to(u.Hz).value
G_nu_effective = (const.c / (6230.06 * u.AA)).to(u.Hz).value


def convert_ab_mag_to_fnu(_mag_ab):
    """
    abmag=-2.5log10(fnu)-48.6  where fnu in units of erg/s/cm2/Hzs
    :return:
    """
    return np.power(10., (48.6 + _mag_ab) / -2.5)


def compute_gaia_features(_t: Table) -> Table:
    """
    Compute features based on the X-ray and Gaia properties of the source.
    :param _t:
    :type _t:
    :return:
    :rtype:
    """
    # Parallax and motion features
    suffix_gaia = '_GAIA_DR3'
    _t['PLX_SIG'] = _t['parallax%s' % suffix_gaia] / _t['parallax_error%s' % suffix_gaia]
    _t['PM_SIG'] = np.sqrt(np.power(_t['pmra%s' % suffix_gaia] / _t['pmra_error%s' % suffix_gaia], 2.) +
                                 np.power(_t['pmdec%s' % suffix_gaia] / _t['pmdec_error%s' % suffix_gaia], 2.))

    # Distance and absolute mag
    _t['INV_PLX'] = 1000. / _t['parallax%s' % suffix_gaia]  # in parsecs
    _t['GAIA_ABS_MAG'] = _t['phot_g_mean_mag%s' % suffix_gaia] - 5 * np.log10(0.1 * _t['INV_PLX'])

    # Flux ratios
    _t['Bp_Rp'] = _t['phot_bp_mean_mag_GAIA_DR3'] - _t['phot_rp_mean_mag_GAIA_DR3']
    _t['G_AB'] = _t['phot_g_mean_mag%s' % suffix_gaia] - (25.6884 - 25.7934)
    _t['G_FNU'] = convert_ab_mag_to_fnu(_t['G_AB'])
    _t['G_NUFNU'] = G_nu_effective * _t['G_FNU']
    _t['Fx_over_FG'] = np.log10(_t['Fx'] / _t['G_NUFNU'])
    return _t


def compute_wise_features(_t: Table) -> Table:
    """
    Compute features based on the X-ray and WISE (mid-infrared) properties of the source.
    :param _t:
    :type _t:
    :return:
    :rtype:
    """
    suffix_allwise = '_ALLWISE'
    suffix_gaia = '_GAIA_DR3'
    _t['W1_W2'] = _t['W1mag%s' % suffix_allwise] - _t['W2mag%s' % suffix_allwise]

    # Compute flux ratios
    _t['W1_AB'] = _t['W1mag%s' % suffix_allwise] + 2.699
    _t['W1_FNU'] = convert_ab_mag_to_fnu(_t['W1_AB'])
    _t['W1_NUFNU'] = w1_nu_effective * _t['W1_FNU']
    _t['Fx_over_FW1'] = np.log10(_t['Fx'] / _t['W1_NUFNU'])

    # and colours
    _t['G_W1'] = _t['phot_g_mean_mag%s' % suffix_gaia] - _t['W1mag%s' % suffix_allwise]
    _t['G_W2'] = _t['phot_g_mean_mag%s' % suffix_gaia] - _t['W2mag%s' % suffix_allwise]
    return _t
