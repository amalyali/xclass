from astropy.table import Table, vstack
from pathlib import Path
import subprocess
import yaml
from xclass.utils import *

stilts = '../data/tools/stilts.jar'


def read_cat(_f: str, _cfg: dict, _label: str) -> Table:
    """
    Read in the basic labelled source catalogue.
    """
    cols = ['RA_BEST', 'DEC_BEST', 'LABEL']
    t = Table.read(_f)
    t['RA_BEST'] = t[_cfg['RA_BEST']]
    t['DEC_BEST'] = t[_cfg['DEC_BEST']]
    t['LABEL'] = _label
    return t[cols]


def stilts_match(_cat1: str, _cat2: str, _values1: str, _values2: str, _suffix2: str, _output_cat: str) -> None:
    """
    Run stilts cross-match between primary and secondary catalogue
    :param _cat1:
    :type _cat1:
    :param _cat2:
    :type _cat2:
    :param _values1:
    :type _values1:
    :param _values2:
    :type _values2:
    :param _suffix2:
    :type _suffix2:
    :param _output_cat:
    :type _output_cat:
    :return:
    :rtype:
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


def add_multiwavelength_counterparts(_cat: str, _cat_output: str) -> None:
    """
    Cross-match the input catalogue with the AllWISE and Gaia catalogues.
    """
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


def compute_features(_t: Table) -> Table:
    """
    Add in additional feature columns to the table.
    """
    _t['Fx'] = _t['SC_EP_8_FLUX']  # 0.2--12 keV. 0.2-2 keV is EP_6.
    _t = compute_gaia_features(_t)
    _t = compute_wise_features(_t)

    # Clean up features file
    features = ['RA_BEST', 'DEC_BEST', 'LABEL',
                'Fx', 'PLX_SIG', 'PM_SIG', 'Bp_Rp', 'W1_W2', 'G_W1', 'G_W2',
                'Fx_over_FG', 'Fx_over_FW1']
    return _t[features]


def main(_cfg: dict, _ver: str):
    f_training = '../data/train_%s.fits' % ver
    f_xmm = '../data/cats/4XMM_slim_DR13cat_v1.0.fits'
    f_training_xray = '../data/train_%s_xray.fits' % ver
    f_training_multiwavelength = '../data/train_%s_xray_w_multiwavelength.fits' % ver
    f_training_features = '../data/train_%s_features.fits' % ver

    # Read in all catalogues and stack into single file.
    tables = []
    for label, cfg_data in cfg[ver].items():
        tables.append(read_cat(cfg_data['PATH'], cfg_data, label))
    tables = vstack(tables)
    tables.write(f_training, format='fits', overwrite=True)

    # Match to multi-wavelength counterparts
    match_params = "'RA_BEST DEC_BEST 1'"
    match_params_xmm = "'SC_RA SC_DEC 5'"
    suffix_xmm = 'XMM_'
    stilts_match(f_training, f_xmm, match_params, match_params_xmm, suffix_xmm, f_training_xray)
    add_multiwavelength_counterparts(f_training_xray, f_training_multiwavelength)

    # Compute features and store in file
    t = Table.read(f_training_multiwavelength)
    t = compute_features(t)
    t.write(f_training_features, format='fits', overwrite=True)


if __name__ == '__main__':
    cfg = yaml.safe_load(Path('../cfg/train.yml').read_text())
    ver = 'v0001'
    main(cfg, ver)
