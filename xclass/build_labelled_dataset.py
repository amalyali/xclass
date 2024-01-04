"""

"""
from astropy.io import fits
from astropy.table import Table, vstack
import pandas as pd
from pathlib import Path
import subprocess
import yaml

stilts = '../data/tools/stilts.jar'


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
                "out=%s" % _output_cat,  # todo: change this to a new directory
                "ofmt=fits"]
        cmd = ' '.join(cmds)
        subprocess.check_call(cmd, shell=True)
    except Exception as e:
        print(e)


def add_multiwavelength_counterparts(_cat, _cat_output):
    cat_tmp = '%s_allwise.fits' % _cat[:-5]
    cmds = ["java -jar %s cdsskymatch" % stilts,
            "cdstable=ALLWISE",
            "in=%s" % _cat,
            "radius=3",
            "ra=RA_BEST",
            "dec=DEC_BEST",
            "find=best",
            "omode=out",
            "out=%s" % cat_tmp]
    cmd = ' '.join(cmds)
    print(cmd)
    subprocess.check_call(cmd, shell=True)

    cmds = ["java -jar %s cdsskymatch" % stilts,
            "in=%s" % cat_tmp,
            "radius=3",
            "cdstable='I/350/gaiaedr3'",
            "ra=RA_BEST",
            "dec=DEC_BEST",
            "find=best",
            "omode=out",
            "out=%s" % _cat_output]
    cmd = ' '.join(cmds)
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def enrich_cat():
    """
    Cross-match the set of known objects with their X-ray and further multi-wavelength counterparts.
    :return:
    :rtype:
    """
    #vals_xmm =
    stilts_match()
    pass


def compute_features():
    pass


if __name__ == '__main__':

    cfg = yaml.safe_load(Path('../cfg/train.yml').read_text())
    ver = 'v0001'
    f_training = '../data/train_%s.fits' % ver
    f_xmm = '../data/cats/4XMM_slim_DR13cat_v1.0.fits'
    f_training_xray = '../data/train_%s_xray.fits' % ver
    f_training_multiwavelength = '../data/train_%s_xray_w_multiwavelength.fits' % ver

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


    #stilts_match(f_training, f_xmm, match_params, match_params_xmm, suffix_xmm, f_training_xray)
    add_multiwavelength_counterparts(f_training_xray, f_training_multiwavelength)
    print(cfg)
