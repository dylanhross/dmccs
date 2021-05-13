"""
    dmim_analysis/plots.py
    Dylan H. Ross
    2019/10/30

    description:
        Functions for generating plots of various extracted data.
"""


import os
from matplotlib import pyplot as plt, rcParams
from PIL import Image
from io import BytesIO
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRaw


# set fonts and stuff
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['font.size'] = 8


def atd_with_fit_and_structure(png, name, adduct, monoiso, ccs, dt_fit_params, atd, rt_fit_params, chromatogram, 
                               mztol, prefix, no_cofac_dti=None, met_rxn=None, met_mass_diff=None):
    """
atd_with_fit_and_structure
    description:
        generates a plot of arrival time distribution with fitted gaussian peak and residuals, alongside the compound
        structure and the LC chromatogram 
    parameters:
        png (bytes) -- PNG image of compound/metabolite structure
        name (str) -- compound name
        adduct (str) -- type of MS adduct
        monoiso (float) -- monoisotopic mass of MS adduct
        ccs (float) -- CCS of the adduct
        dt_fit_params (tuple(float)) -- parameters from gaussian fit on ATD
        atd (numpy.ndarray(float)) -- ATD (2D array: [[drift time], [intensity]])
        rt_fit_params (tuple(float)) -- parameters from gaussian fit on chromatogram
        chromatogram (numpy.ndarray(float)) -- LC chromatogram (2D array: [[retention time], [intensity]])
        mztol (float) -- m/z tolerance
        prefix (str) -- path to save images under
        [no_cofac_dti (numpy.ndarray(float))] -- intensities from extracted ATD from non-cofactor reaction, for 
                                                    metabolites only [optional, default=None]
        [met_rxn (str)] -- metabolic reaction, for metabolites only [optional, default=None]
        [met_mass_diff (float)] -- mass difference of metabolite relative to parent compound [optional, default=None]
"""
    # going to use the gauss method from here
    ccr = CCSCalibrationRaw(None, None, None, no_init=True)

    gs_kw = {'height_ratios': [3, 1, 1]}
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, gridspec_kw=gs_kw, figsize=[3, 6])

    if png:
        ax1.imshow(Image.open(BytesIO(png)))
    ax1.set_axis_off()

    ttl = '{} {} \n(m/z: {:.4f} ± {:.4f}, CCS: {:.1f})'.format(name, adduct, monoiso, mztol, ccs)
    if met_rxn:
        ttl += '\n{:40.40s}'.format(met_rxn)
        ttl += '\n({:+.4f})'.format(met_mass_diff)

    ax1.set_title(ttl)

    A, dt, C = dt_fit_params
    t, i = atd
    gauss = ccr.gauss(t, A, dt, C)
    resids = i - gauss 
    ax2.plot(t, i, 'b-', lw=1.5, label='raw')
    ax2.plot(t, gauss, 'r--', lw=1, label='fit')
    ax2.plot(t, resids, 'k-', lw=0.5, label='residuals')
    if no_cofac_dti is not None:
        ax2.plot(t, no_cofac_dti, '-', c='lime', lw=0.75, label='no cofactors')
    ax2.set_title('ATD')
    ax2.set_ylabel('intensity')
    ax2.set_xlabel('drift time (ms) → {:.2f} ms (FWHM ~ {:.2f} ms)'.format(dt, C * 2.355))
    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax2.legend(fontsize=5)
    for d in ['top', 'right']:
        ax2.spines[d].set_visible(False)

    xlab = 'retention time (min)'
    t, i = chromatogram
    ax3.plot(t, i, 'b-', lw=1.5, label='raw')
    if rt_fit_params is not None:
        A, dt, C = rt_fit_params
        gauss = ccr.gauss(t, A, dt, C)
        resids = i - gauss
        ax3.plot(t, gauss, 'r--', lw=1, label='fit')
        ax3.plot(t, resids, 'k-', lw=0.5, label='residuals')
        xlab += '→ {:.2f} min (FWHM ~ {:.2f} min)'.format(dt, C * 2.355)
    ax3.set_xlabel(xlab)
    ax3.set_title('LC chromatogram')
    ax3.set_ylabel('intensity')

    ax3.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
    ax3.legend(fontsize=5)
    for d in ['top', 'right']:
        ax3.spines[d].set_visible(False)

    plt.tight_layout()

    fname = os.path.join(prefix, '{}_{}.png'.format(name.replace(' ', '_'), adduct))
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()


def parent_atd_with_fit_and_structure(cmpd, prefix):
    """
parent_atd_with_fit_and_structure
    description:
        generates a plot of arrival time distribution with fitted gaussian peak and residuals, alongside the compound
        structure and the LC chromatogram -- for all adducts of the parent compound
    parameters:
        cmpd (DMIM_data.Compound) -- data structure with all of the analyzed data on a particular compound
        prefix (str) -- path to save images under
"""
    for adduct in cmpd.adducts:

        atd_with_fit_and_structure(
            cmpd.png, 
            cmpd.name, 
            adduct.adduct, 
            adduct.monoiso,
            adduct.ccs,
            adduct.atd_fit_params,
            adduct.atd,
            adduct.rt_fit_params,
            adduct.chromatogram,
            adduct.tolerance,
            prefix
        )


def metabolite_atd_with_fit_and_structure(cmpd, prefix):
    """
metabolite_atd_with_fit_and_structure
    description:
        generates a plot of arrival time distribution with fitted gaussian peak and residuals, alongside the compound
        structure and the LC chromatogram -- for all adducts of all metabolites
    parameters:
        cmpd (DMIM_data.Compound) -- data structure with all of the analyzed data on a particular compound
        prefix (str) -- path to save images under
"""
    for metab in cmpd.metabolites:
        for adduct in metab.adducts:

            atd_with_fit_and_structure(
                metab.png, 
                '{} met{:03d}'.format(cmpd.name, metab.id), 
                adduct.adduct, 
                adduct.monoiso,
                adduct.ccs,
                adduct.atd_fit_params,
                adduct.atd,
                adduct.rt_fit_params,
                adduct.chromatogram,
                adduct.tolerance,
                prefix,
                no_cofac_dti=adduct.no_cofac_dti,
                met_rxn=metab.reaction,
                met_mass_diff=(metab.monoiso - cmpd.monoiso)
            )



def atd_with_fit(name_adduct, monoiso, ccs, dt_fit_params, atd, rt_fit_params, chromatogram, mztol, prefix):
    """
atd_with_fit_and_structure
    description:
        generates a plot of arrival time distribution with fitted gaussian peak and residuals, alongside the LC 
        chromatogram 
    parameters:
        name_adduct (str) -- compound name and adduct
        monoiso (float) -- monoisotopic mass of MS adduct
        ccs (float) -- CCS of the adduct
        dt_fit_params (tuple(float)) -- parameters from gaussian fit on ATD
        atd (numpy.ndarray(float)) -- ATD (2D array: [[drift time], [intensity]])
        rt_fit_params (tuple(float)) -- parameters from gaussian fit on chromatogram
        chromatogram (numpy.ndarray(float)) -- LC chromatogram (2D array: [[retention time], [intensity]])
        mztol (float) -- m/z tolerance
        prefix (str) -- path to save images under
"""
    # going to use the gauss method from here
    ccr = CCSCalibrationRaw(None, None, None, no_init=True)

    fig = plt.figure(figsize=(3, 4))
    ax2 = fig.add_subplot(211)
    ax3 = fig.add_subplot(212)

    ttl = '{} \n(m/z: {:.4f} ± {:.4f}, CCS: {:.1f})'.format(name_adduct, monoiso, mztol, ccs)

    ax2.set_title(ttl)

    A, dt, C = dt_fit_params
    t, i = atd
    gauss = ccr.gauss(t, A, dt, C)
    resids = i - gauss 
    ax2.plot(t, i, 'b-', lw=1.5, label='raw')
    ax2.plot(t, gauss, 'r--', lw=1, label='fit')
    ax2.plot(t, resids, 'k-', lw=0.5, label='residuals')
    ax2.set_ylabel('intensity')
    ax2.set_xlabel('drift time (ms) → {:.2f} ms (FWHM ~ {:.2f} ms)'.format(dt, C * 2.355))
    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
    ax2.legend(fontsize=5)
    for d in ['top', 'right']:
        ax2.spines[d].set_visible(False)

    xlab = 'retention time (min)'
    t, i = chromatogram
    ax3.plot(t, i, 'b-', lw=1.5, label='raw')
    if rt_fit_params is not None:
        A, dt, C = rt_fit_params
        gauss = ccr.gauss(t, A, dt, C)
        resids = i - gauss
        ax3.plot(t, gauss, 'r--', lw=1, label='fit')
        ax3.plot(t, resids, 'k-', lw=0.5, label='residuals')
        xlab += '→ {:.2f} min (FWHM ~ {:.2f} min)'.format(dt, C * 2.355)
    ax3.set_xlabel(xlab)
    ax3.set_title('LC chromatogram')
    ax3.set_ylabel('intensity')

    ax3.ticklabel_format(style='sci', scilimits=(0, 0), axis='y')
    ax3.legend(fontsize=5)
    for d in ['top', 'right']:
        ax3.spines[d].set_visible(False)

    plt.tight_layout()

    fname = os.path.join(prefix, '{}.png'.format(name_adduct))
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()


