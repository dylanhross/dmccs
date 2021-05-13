"""
    dmim_analysis/analysis.py
    Dylan H. Ross
    2019/10/30

    description:
        tools for analyzing data
"""

from numpy import array

from dmim_analysis.util import remove_counter_ions, get_and_fit_chrom
from dmim_analysis.scraping.pubchem import get_pchem_monoiso, get_pchem_smi
from dmim_analysis.scraping.cactus import get_cactus_img
from dmim_analysis.data import Compound, Adduct, Metabolite
from dhrmasslynxapi.reader import MassLynxReader
from biotrans.pybiotransformer.wrapper import PBTWrapper


def setup_compound(cmpd_id, name):
    """
setup_compound
    description:
        Gathers some metadata on a compound including compound name with counterions removed, monoisotopic mass and
        SMILES structure (from PubChem), and png image of the compound structure (from cactus), storing all of it in
        a Compound object and returning it. Returns None if monoisotopic mass or SMILES structure fail to be retrieved
        since both are needed for later steps in the data processing.
    parameters:
        cmpd_id (int) -- integer identifier for the compound
        name (str) -- compound name
    returns:
        (dmim_analysis.data.Compound or None) -- initial Compound instance (or None if anything goes wrong)
"""
    # remove unwanted counter ions from name
    name = remove_counter_ions(name)

    # get a monoisotopic mass from PubChem
    monoiso = get_pchem_monoiso(name)
    if not monoiso:
        # need monoisotopic mass to proceed ...
        return None

    # get a SMILES structure from PubChem
    smi = get_pchem_smi(name)
    if not smi:
        # need to get a SMILES structure to proceed ...
        return None

    # get an image of the compound's chemical structure
    png = get_cactus_img(smi)

    # actually set up the Compound object
    return Compound(cmpd_id, name, monoiso, smi, png=png)


def find_adducts(raw_file, neutral_monoiso, tolerance, esi_mode, ccs_calibration, lc_func, im_func,
                 is_metab=False, ctrl_raw_file=None, filter_by_dt=True):
    """
find_adducts
    description:
        Try extracting data for masses corresponding to multiple possible MS adducts ([M]+, [M+Na]+, ...)
    parameters:
        raw_file (str) -- file name for raw data file
        neutral_monoiso (float) -- monoisotopic mass of neutral species
        tolerance (float) -- tolerance to use when filtering on mass
        esi_mode (str) -- ionization mode, may be 'pos' or 'neg'
        ccs_calibration (dhrmasslynxapi.ccs_calibration.CCSCalibrationRaw) -- CCS calibration instance
        lc_func (int) -- function number for LC data
        im_func (int) -- function number for IM data
        [is_metab (bool)] -- flag for specifying looking for metabolite adducts instead of parent compound adducts,
                                in which case the control data file is also checked for cofactor dependence
                                [optional, default=False]
        [ctrl_raw_file (str)] -- if is_metab is set to True, file name of the control raw file [optional, default=None]
        [filter_by_dt (bool)] -- filter the LC chromatogram using the fitted drift time [optional, default=True]
    yields:
        (DMIM_data.Adduct) -- data for the MS adduct (if fit is acceptable)
"""
    # define mass shifts for various MS adducts
    adducts = {
        'pos': {
            '[M]+': 0.,
            '[M+H]+': 1.0078,
            '[M+Na]+': 22.9898,
            '[M+K]+': 38.9637,
            '[M+H-H2O]+': -17.0028
        },
        'neg': {
            '[M-H]-': -1.0078
        }
    }

    # initialize a MassLynxReader using the raw file
    rdr = MassLynxReader(raw_file)
    if is_metab:
        if not ctrl_raw_file:
            m = 'find_adducts: is_metab was set to True but no ctrl_raw_file was provided'
            raise ValueError(m)

        ctrl_rdr = MassLynxReader(ctrl_raw_file)

    for adduct in adducts[esi_mode]:
        # MS adduct m/z
        adduct_mz = neutral_monoiso + adducts[esi_mode][adduct]

        # attempt to find an observed m/z
        #mz_obs = get_mz_obs(rdr, adduct_mz, tolerance)
        mz_obs = adduct_mz

        # attempt a fit on ATD
        atd, atd_fit_params = get_and_fit_chrom(rdr, im_func, mz_obs, tolerance)

        # apply rough acceptance criteria
        if atd_fit_params is not None:
            low_intensity = atd_fit_params[0] < 1000 
            peak_too_broad = atd_fit_params[2] > 0.75
            peak_too_narrow = atd_fit_params[2] < 0.025
            if low_intensity or peak_too_broad or peak_too_narrow:
                atd_fit_params = None  # do not accept the fit

        # only proceed if ATD fit was successful
        if atd_fit_params is not None:
            # get the drift time and calibrated CCS
            dt = atd_fit_params[1]
            ccs = ccs_calibration.calibrated_ccs(mz_obs, dt)

            # attempt to fit LC chromatogram
            if filter_by_dt:
                dt_tol = atd_fit_params[2] 
                dt_min, dt_max = dt - dt_tol, dt + dt_tol
                lc, lc_fit_params = get_and_fit_chrom(rdr, lc_func, mz_obs, tolerance, dt_bounds=(dt_min, dt_max))
            else:
                lc, lc_fit_params = get_and_fit_chrom(rdr, lc_func, mz_obs, tolerance)
            rt = lc_fit_params[1] if lc_fit_params is not None else -1

            # finally, if it is a metabolite, get the no-cofactor ATD
            ctrl_atd = None
            if is_metab:
                # noinspection PyUnboundLocalVariable
                _, ctrl_dti = array(ctrl_rdr.get_chrom(im_func, mz_obs, tolerance))
                yield Adduct(adduct, adduct_mz, mz_obs, tolerance, atd, atd_fit_params, lc, lc_fit_params, dt, rt, ccs,
                             raw_file, ctrl_d_file=ctrl_raw_file, no_cofac_dti=ctrl_dti)
            else:
                yield Adduct(adduct, adduct_mz, mz_obs, tolerance, atd, atd_fit_params, lc, lc_fit_params, dt, rt, ccs,
                             raw_file)


def predict_metabolites(cmpd):
    """
predict_metabolites
    description:
        predicts metabolites from the parent compound and yields them (as long as the mass difference from the parent
        compound is > 1 Da)
    parameters:
        cmpd (dmim_analysis.data.Compound) -- parent Compound instance
    yields:
        (dmim_analysis.data.Metabolite) -- data structure for a predicted metabolite
"""
    # initialize the Biotransformer wrapper
    pbtw = PBTWrapper()

    i = 0
    for metab in pbtw.predict_metabolites(cmpd.SMILES, n_steps=2, unique=True):
        if abs(metab.mass - cmpd.monoiso) > 1.:
            i += 1
            # grab an image for the metabolite structure
            png = get_cactus_img(metab.InChI)

            yield Metabolite(i, metab.reaction, metab.mass, metab.InChI, png=png)


def analyze_compound(cmpd_id, name, raw_file, lc_func, im_func, tolerance, esi_mode, ccs_calibration, check_ctrl=False,
                     ctrl_raw_file=None):
    """
analyze_compound
    description:
        Performs a complete analysis on a single compound:
            - set up Compound instance with metadata
            - find MS adducts of the parent compound
            - predict metabolites of the parent compound
            - look for MS adducts of predicted metabolites
            - return the Compound instance
    parameters:
        cmpd_id (str) -- parent compound identifier
        name (str) -- parent compound name
        raw_file (str) -- file name for raw data file
        lc_func (int) -- function number for LC data
        im_func (int) -- function number for IM data
        tolerance (float) -- tolerance to use when filtering on mass
        esi_mode (str) -- ionization mode, may be 'pos' or 'neg'
        ccs_calibration (dhrmasslynxapi.ccs_calibration.CCSCalibrationRaw) -- CCS calibration instance
        [check_ctrl (bool)] -- check the control data file for metabolite signals [optional, default=False]
        [ctrl_raw_file (str)] -- file name of the control raw file [optional, default=None]
    returns:
        (dmim_analysis.data.Compound) -- fully analyzed Compound instance
"""
    # setup Compound instance
    cmpd = setup_compound(cmpd_id, name)
    if not cmpd:
        # compound setup did not work
        return None

    # find MS adducts
    cmpd.adducts = [_ for _ in find_adducts(raw_file, cmpd.monoiso, tolerance, esi_mode, ccs_calibration, lc_func,
                                            im_func)]
    if not cmpd.adducts:
        # no MS adducts were observed for the parent compound
        return None

    # predict metabolites
    cmpd.metabolites = [_ for _ in predict_metabolites(cmpd)]

    # look for MS adducts of the predicted metabolites
    for metab in cmpd.metabolites:
        if check_ctrl:
            metab.adducts = [_ for _ in find_adducts(raw_file, metab.monoiso, tolerance, esi_mode, ccs_calibration,
                                                     lc_func, im_func, is_metab=True, ctrl_raw_file=ctrl_raw_file)]
        else:
            metab.adducts = [_ for _ in find_adducts(raw_file, metab.monoiso, tolerance, esi_mode, ccs_calibration,
                                                     lc_func, im_func)]

    # filter the list of metabolites down to those that were actually observed in the MS data
    cmpd.metabolites = [metab for metab in cmpd.metabolites if metab.adducts]

    # return the analyzed compound data
    return cmpd


def analyze_compound_targeted(name_adduct, mz, raw_file, lc_func, im_func, tolerance, ccs_calibration, 
                              filter_by_dt=True):
    """
analyze_compound_targeted
    description:
        Attempts to obtain a CCS value for an alr
    parameters:
        name_adduct (str) -- parent compound name and MS adduct
        mz (float) -- m/z
        raw_file (str) -- file name for raw data file
        lc_func (int) -- function number for LC data
        im_func (int) -- function number for IM data
        tolerance (float) -- tolerance to use when filtering on mass
        ccs_calibration (dhrmasslynxapi.ccs_calibration.CCSCalibrationRaw) -- CCS calibration instance
        [filter_by_dt (bool)] -- filter the LC chromatogram using the fitted drift time [optional, default=True]
    returns:
        (dict(...) or None) -- compound data and associated metadata or None if unsuccessful
"""
    # initialize a MassLynxReader using the raw file
    rdr = MassLynxReader(raw_file)

    # attempt to find an observed m/z
    #mz_obs = get_mz_obs(rdr, mz, tolerance)
    mz_obs = mz

    # attempt a fit on ATD
    atd, atd_fit_params = get_and_fit_chrom(rdr, im_func, mz_obs, tolerance)

    # apply rough acceptance criteria
    if atd_fit_params is not None:
        low_intensity = atd_fit_params[0] < 1000 
        peak_too_broad = atd_fit_params[2] > 0.75
        peak_too_narrow = atd_fit_params[2] < 0.025
        if low_intensity or peak_too_broad or peak_too_narrow:
            atd_fit_params = None  # do not accept the fit

    # only proceed if ATD fit was successful
    if atd_fit_params is not None:
        # get the drift time and calibrated CCS
        dt = atd_fit_params[1]
        ccs = ccs_calibration.calibrated_ccs(mz_obs, dt)

        # attempt to fit LC chromatogram
        if filter_by_dt:
            dt_tol = atd_fit_params[2] 
            dt_min, dt_max = dt - dt_tol, dt + dt_tol
            lc, lc_fit_params = get_and_fit_chrom(rdr, lc_func, mz_obs, tolerance, dt_bounds=(dt_min, dt_max))
        else:
            lc, lc_fit_params = get_and_fit_chrom(rdr, lc_func, mz_obs, tolerance)
        rt = lc_fit_params[1] if lc_fit_params is not None else -1

        # assemble the data and metadata to return
        data = {
            'dt': dt, 
            'ccs': ccs, 
            'rt': rt, 
            'meta': {
                'mz_obs': mz_obs, 
                'atd': atd, 
                'atd_fit_params': atd_fit_params,
                'lc': lc,
                'lc_fit_params': lc_fit_params
            }
        }
        return data

    else:
        return None
