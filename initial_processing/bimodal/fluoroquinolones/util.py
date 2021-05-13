"""  """


from numpy import array, exp, nan_to_num, savetxt
from scipy.optimize import curve_fit
from pickle import load as pload
from matplotlib import pyplot as plt, ticker as mtick
from matplotlib import rcParams

rcParams['font.size'] = 8


from dhrmasslynxapi.reader import MassLynxReader


def get_ccs_dist(rawf, func, mz, tol, calf):
    """ obtains  the ATD for a mass then converts drift time axis to CCS """
    # initialize reader
    rdr = MassLynxReader(rawf)

    # get the ATD and convert to intensity to an array
    t, i = rdr.get_chrom(func, mz, tol)
    i = array(i)

    # load and apply CCS calibration
    with open(calf, 'rb') as pf:
        cal = pload(pf)
    ccs = array([cal.calibrated_ccs(mz, dt) for dt in t])

    return t, ccs, i


def plot_replicate_ccs_dists(name, cd1, cd2, cd3, norm=False, x_range=None):
    """ generates a plot of the CCS distributions from all replicates """

    # split up the data, normalize if needed
    _, x1, y1 = cd1
    _, x2, y2 = cd2
    _, x3, y3 = cd3
    if norm:
        y1 = y1 / max(y1)
        y2 = y2 / max(y2)
        y3 = y3 / max(y3)

    # make the plot
    fig = plt.figure(figsize=(4, 2.5))
    ax = fig.add_subplot(111)
    ax.plot(x1, y1, 'b-', label='r1', lw=0.75)
    ax.plot(x2, y2, 'r-', label='r2', lw=0.75)
    ax.plot(x3, y3, '-', c='darkkhaki', label='r3', lw=0.75)
    ax.set_ylabel('intensity (normalized)' if norm else 'intensity')
    ax.set_xlabel(r'CCS ($\AA^2$)')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend()

    if x_range is not None:
        ax.set_xlim(x_range)

    # save and close
    plt.savefig('{}_ccsdists{}.png'.format(name, '_norm' if norm else ''), dpi=400, bbox_inches='tight')
    plt.close()


def plot_replicate_atds(name, cd1, cd2, cd3, norm=False):
    """ generates a plot of the CCS distributions from all replicates """

    # split up the data, normalize if needed
    x1, _, y1 = cd1
    x2, _, y2 = cd2
    x3, _, y3 = cd3
    if norm:
        y1 = y1 / max(y1)
        y2 = y2 / max(y2)
        y3 = y3 / max(y3)

    # make the plot
    fig = plt.figure(figsize=(4, 2.5))
    ax = fig.add_subplot(111)
    ax.plot(x1, y1, 'b-', label='r1', lw=0.75)
    ax.plot(x2, y2, 'r-', label='r2', lw=0.75)
    ax.plot(x3, y3, '-', c='darkkhaki', label='r3', lw=0.75)
    ax.set_ylabel('intensity (normalized)' if norm else 'intensity')
    ax.set_xlabel('drift time (ms)')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend()

    # save and close
    plt.savefig('{}_atds{}.png'.format(name, '_norm' if norm else ''), dpi=400, bbox_inches='tight')
    plt.close()


def two_peak_fit(t, i, t1, t2, sd):
    """ attempts to fit two peaks from ATD or CCS dist. with two guesses for peak centers """

    # if t has nans, convert those to 0
    t = nan_to_num(t)

    # sum of two gaussians for fitting
    def double_gauss(x, A1, B1, C1, A2, B2, C2):
        return A1 * exp(-(x - B1)**2 / (2. * C1**2)) + A2 * exp(-(x - B2)**2 / (2. * C2**2))

    # fit with double gaussian
    p0 = (max(i), t1, sd, max(i), t2, sd)
    opt, cov = curve_fit(double_gauss, t, i, maxfev=5000, p0=p0)

    return opt[1], opt[4]


def plot_fragment_spectra(label, rawf, peak1_dt, peak2_dt, dt_tol, mz_min, mz_max, dump_txt=False):
    """ plots mass spectra extracted for two dt peaks """

    # set up reader
    rdr = MassLynxReader(rawf)

    # get the spectra
    p1_m, p1_i = rdr.get_spectrum(3, peak1_dt - dt_tol, peak1_dt + dt_tol)
    p2_m, p2_i = rdr.get_spectrum(3, peak2_dt - dt_tol, peak2_dt + dt_tol)

    # save spectra to a .txt file
    if dump_txt:
        savetxt('{}_MS2_peak1.txt'.format(label), array([p1_m, p1_i]).T, '%.4f %d')
        savetxt('{}_MS2_peak2.txt'.format(label), array([p2_m, p2_i]).T, '%.4f %d')

    # find tha maximum within the viewing range
    ymax = 0
    for spect in [(p1_m, p1_i), (p2_m, p2_i)]:
        for m, i in zip(*spect):
            if m >= mz_min and m <= mz_max and i > ymax:
                ymax = i

    # generate the plot
    fig = plt.figure(figsize=(3.5, 2.5))
    ax = fig.add_subplot(111)
    ax.plot(p1_m, p1_i, 'b-', lw=0.75, alpha=0.5, label='peak 1')
    ax.plot(p2_m, p2_i, 'r-', lw=0.75, alpha=0.5, label='peak 2')
    ax.set_xlim([mz_min, mz_max])
    ax.set_ylabel('intensity')
    ax.set_xlabel('m/z')
    ax.yaxis.set_ticks([])
    ax.set_ylim([0, ymax])
    ax.xaxis.set_tick_params(rotation=45, size=6)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend()

    # save and close
    plt.savefig('{}_MS2_spectra.png'.format(label), dpi=400, bbox_inches='tight')
    plt.close()


def plot_fragment_atd(name, rawf, func, f_mz, tol, cd1, dt_scale=1.0):
    """ generates a plot of the ATD for a specified fragment """

    # split up the data, normalize if needed
    x1, _, y1 = cd1
    y1 = y1 / max(y1)


    # get the ATD and convert to intensity to an array
    rdr = MassLynxReader(rawf)
    x2, y2 = rdr.get_chrom(func, f_mz, tol)
    y2 = array(y2) / max(y2)
    x2 = array(x2) * dt_scale

    # make the plot
    fig = plt.figure(figsize=(4, 2.5))
    ax = fig.add_subplot(111)
    ax.plot(x1, y1, 'k-', label='precursor', lw=0.75)
    ax.plot(x2, y2, 'r-', label='fragment', lw=0.75)
    ax.set_ylabel('intensity')
    ax.yaxis.set_ticks([])
    ax.set_xlabel('drift time (ms)')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend()

    # save and close
    plt.savefig('{}_atd.png'.format(name), dpi=400, bbox_inches='tight')
    plt.close()

