"""

"""
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import yaml
#from xclass.utils import test_fn
#from . import utils
from astropy import constants as const
from astropy import units as u

stilts = '../data/tools/stilts.jar'
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


def read_cat(_f: str, _cfg: dict, _label: str) -> Table:

    cols = ['RA_BEST', 'DEC_BEST', 'LABEL']
    t = Table.read(_f)
    t['RA_BEST'] = t[_cfg['RA_BEST']]
    t['DEC_BEST'] = t[_cfg['DEC_BEST']]
    t['LABEL'] = _label
    return t[cols]


def stilts_match(_cat1, _cat2, _values1, _values2, _suffix2, _output_cat):
    """
    Run stilts cross-match between primary and secondary catalogue
    :return:
    """
    try:
        cmds = ["java -jar %s tmatch2" % stilts,
                "in1=%s" % _cat1,
                "in2=%s" % _cat2,
                "ifmt1=fits",
                "ifmt2=fits",
                "matcher=skyerr",
                "find=best1",
                "params=10",
                "values1=%s" % str(_values1),
                "values2=%s" % str(_values2),
                "suffix1=''",
                "suffix2=%s" % str(_suffix2),
                "join=1and2",
                "out=%s" % _output_cat,
                "ofmt=fits"]
        cmd = ' '.join(cmds)
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        print(e)


def add_multiwavelength_counterparts(_cat: str, _cat_output: str):
    cat_tmp = '%s_allwise.fits' % _cat[:-5]
    cmds = ["java -jar %s cdsskymatch" % stilts,
            "cdstable=ALLWISE",
            "in=%s" % _cat,
            "radius=3",
            "ra=RA_BEST",
            "dec=DEC_BEST",
            "find=each",
            "omode=out",
            "fixcols=all",
            "suffixin=''",
            "suffixremote=_ALLWISE",
            "out=%s" % cat_tmp]
    cmd = ' '.join(cmds)
    subprocess.check_call(cmd, shell=True)

    cmds = ["java -jar %s cdsskymatch" % stilts,
            "in=%s" % cat_tmp,
            "radius=3",
            "cdstable='I/350/gaiaedr3'",
            "ra=RA_BEST",
            "dec=DEC_BEST",
            "find=each",
            "omode=out",
            "fixcols=all",
            "suffixin=''",
            "suffixremote=_GAIA_DR3",
            "out=%s" % _cat_output]
    cmd = ' '.join(cmds)
    subprocess.check_call(cmd, shell=True)


def compute_gaia_features(_t: Table) -> Table:
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


def compute_features(_t: Table) -> Table:
    """

    :return:
    :rtype:
    """
    _t['Fx'] = _t['SC_EP_8_FLUX']  # 0.2--12 keV. 0.2-2 keV is EP_6.

    _t = compute_gaia_features(_t)
    _t = compute_wise_features(_t)

    # Clean up features file
    features = ['RA_BEST', 'DEC_BEST', 'LABEL',
                'Fx', 'PLX_SIG', 'PM_SIG', 'Bp_Rp', 'W1_W2', 'G_W1', 'G_W2',
                'Fx_over_FG', 'Fx_over_FW1']
    return _t[features]


if __name__ == '__main__':
    cfg = yaml.safe_load(Path('../cfg/train.yml').read_text())
    ver = 'v0001'
    f_training = '../data/train_%s.fits' % ver
    f_xmm = '../data/cats/4XMM_slim_DR13cat_v1.0.fits'
    f_training_xray = '../data/train_%s_xray.fits' % ver
    f_training_multiwavelength = '../data/train_%s_xray_w_multiwavelength.fits' % ver
    f_training_features = '../data/train_%s_features.fits' % ver

    # Read in all catalogues and stack into single file.
    tables = []
    for label, cfg_data in cfg[ver].items():
        print(label, cfg_data)
        print(read_cat(cfg_data['PATH'], cfg_data, label))
        tables.append(read_cat(cfg_data['PATH'], cfg_data, label))
    tables = vstack(tables)
    tables.write(f_training, format='fits', overwrite=True)
    print(tables)

    match_params = "'RA_BEST DEC_BEST 1'"
    match_params_xmm = "'SC_RA SC_DEC 5'"
    suffix_xmm = 'XMM_'


    stilts_match(f_training, f_xmm, match_params, match_params_xmm, suffix_xmm, f_training_xray)
    add_multiwavelength_counterparts(f_training_xray, f_training_multiwavelength)
    print(cfg)

    t = Table.read(f_training_multiwavelength)
    t = compute_features(t)
    t.write(f_training_features, format='fits', overwrite=True)