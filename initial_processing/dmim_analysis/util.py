"""
    dmim_analysis/util.py
    Dylan H. Ross
    2019/10/30

    description:

"""


from numpy import array

from dhrmasslynxapi.ccs_calibration import CCSCalibrationRaw


def remove_counter_ions(name, separator=' '):
    """
remove_counter_ions
    description:
        removes unwanted counter ions (and other extras) from a compound name
    parameters:
        name (str) -- compound name
    returns:
        (str) -- compound name with counter ions removed
"""
    # counter ions to remove from compound names
    counter_ions = [
        'MALEATE', 'SULFATE', 'HEMISULFATE', 'FUMARATE', 'POTASSIUM', 'SODIUM', 'CITRATE', 'MESYLATE', 'SALICYLATE',
        'PHOSPHATE', 'CALCIUM', 'BROMIDE', 'HYDROBROMIDE', 'HYDROCHLORIDE', 'TARTRATE', 'SUCCINATE', 'HYDRATE', 
        'NITRATE', 'IODIDE', 'CHLORIDE', 'SALT', '(+/-)', '(D[+])', '[t(-)]', '(dl)'
    ]
    # check for removable counter ions
    for counter_ion in counter_ions:
        if counter_ion in name.split(separator):
            trimmed = name.split(separator)
            trimmed.remove(counter_ion)
            name = separator.join(trimmed)
    return name


def get_and_fit_chrom(reader, f_n, mz, tolerance, dt_bounds=None, rt_bounds=None):
    """
get_and_fit_chrom
    description:
        Attempts to fit an LC chromatogram or arrival time distribution using a mass and tolerance, returns the
        chromatogram or arrival time distribution as a 2D numpy array as well as the fit parameters on success,
        or (atd/chromatogram, None) on failure
    parameters:
        reader (dhrmasslynxapi.reader.MassLynxReader) -- interface for reading raw data
        f_n (int) -- MS function number
        mz (float) -- MS adduct m/z
        tolerance (float) -- mass tolerance
        rt_bounds, dt_bounds (None or tuple(float)) -- bounds for filtering by retention/drift time [optional, 
                                                       default=None]
    returns:
        (numpy.ndarray, tuple(float)) -- arrival time distribution and fit parameters on success or None
"""
    # use for peak fitting
    ccr = CCSCalibrationRaw(None, None, None, no_init=True)

    if dt_bounds is None and rt_bounds is None:
        t, i = reader.get_chrom(f_n, mz, tolerance)
    elif dt_bounds is not None:
        dt_min, dt_max = dt_bounds
        t, i = reader.get_filtered_chrom(f_n, mz, tolerance, dt_min=dt_min, dt_max=dt_max)
    elif rt_bounds is not None:
        rt_min, rt_max = rt_bounds
        t, i = reader.get_filtered_chrom(f_n, mz, tolerance, rt_min=rt_min, rt_max=rt_max)

    try:
        fit_params = ccr.peak_fit(t, i)
    except:
        fit_params = None

    return array([t, i]), fit_params
