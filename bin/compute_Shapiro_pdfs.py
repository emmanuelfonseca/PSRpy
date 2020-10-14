#! /usr/bin/python

from scipy.interpolate import interp1d, interp2d, griddata
from scipy.stats import chi2
from matplotlib.font_manager import FontProperties
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pickle
import sys

# set font parameters.
font = FontProperties()
font.set_name('serif')
plt.style.use('classic')

# parse command-line arguments.
parser = argparse.ArgumentParser(description="Uses output from the 'SDgrid.py' program and creates interpolated probability maps and marginalized PDFs of the Shapiro-delay parameters (M2, COSI) and the pulsar mass (M1).")
parser.add_argument('infile', nargs=1, action='store', help='input file that contain grid info.')
parser.add_argument('--COSI', nargs=2, type=float, metavar=('min','max'), action='store', dest='COSIlim', default=[0., 1.], help='limits for M2 to use when computing both maps, in units of solar masses. (Default = 0, 1)')
parser.add_argument('--COSIci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='COSIconfints', help='median and confidence interval for M2; used for plots only. ')
parser.add_argument('--M2', nargs=2, type=float, default=[0.01, 3.], metavar=('min','max'), action='store', dest='M2lim', help='limits for M2 to use when computing both maps, in units of solar masses. (Default = [min, max] from input M2 grid.)')
parser.add_argument('--M2ci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='M2confints', help='median and confidence interval for M2; used for plots only. ')
parser.add_argument('--M1', nargs=2, type=float, default=[0., 3.], metavar=('min','max'), action='store', dest='M1lim', help='limits for M1 to use when computing the M2-M1 map, in units of solar masses. (Default = [0.0, 3.0])')
parser.add_argument('--M1ci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='M1confints', help='median and confidence interval for M2; used for plots only. ')
parser.add_argument('--KOM', nargs=2, type=float, default=[0., 360.], metavar=('min','max'), action='store', dest='KOMlim', help='limits for KOM to use for 3D grids, in units of degrees. (Default = [0, 360])')
parser.add_argument('--KOMci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='KOMconfints', help='median and confidence interval for KOM; used for plots only. ')
parser.add_argument('--H3ci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='H3confints', help='median and confidence interval for H3; used for plots only. ')
parser.add_argument('--H4ci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='H4confints', help='median and confidence interval for H4; used for plots only. ')
parser.add_argument('--STIGci', nargs=3, type=float, default=[0., 0., 0.], metavar=('lower', 'median', 'upper'), action='store', dest='STIGconfints', help='median and confidence interval for STIG; used for plots only. ')
parser.add_argument('--M2prior', nargs=2, type=float, default=[0., 0.], metavar=('mean', 'sigma'), action='store', dest='M2priorint', help='Mean and standard deviation for prior on M2; assumed to be Gaussian.')
parser.add_argument('--overplotFile', nargs=1, type=str, default=[None], action='store', dest='overplotFile', help="a Python pickle file containing posterior PDFs to be overplotted with results derived in current execution.")
parser.add_argument('-i', action='store_true', dest='info_only', help='only print info from input grid files.')
parser.add_argument('-n', nargs=1, metavar='num', type=int, default=[100], action='store', dest='N_iter', help='number of M2-COSI/M1 grid points to generate from interpolated grid. (Default = 100)')

args = parser.parse_args()

infile = args.infile
m1_lo, m1_hi = (args.M1lim)[0], (args.M1lim)[1]
m2_lo, m2_hi = (args.M2lim)[0], (args.M2lim)[1]
cosi_lo, cosi_hi = (args.COSIlim)[0], (args.COSIlim)[1]
kom_lo, kom_hi = (args.KOMlim)[0], (args.KOMlim)[1]
info_only = args.info_only
n_iter = (args.N_iter)[0]
m1_confints = args.M1confints
m2_confints = args.M2confints
cosi_confints = args.COSIconfints
kom_confints = args.KOMconfints
h3_confints = args.H3confints
h4_confints = args.H4confints
stig_confints = args.STIGconfints
m2_priorint = args.M2priorint
overplotFile, = args.overplotFile

# read in grid data written using the output format from 
# the 'SDgrid.py' and 'SDgrid.ortho.py' programs.
m2_1D, cosi_1D = 0, 0
h3_1D, h4_1D, stig_1D = 0, 0, 0
kom_1D = 0
min_M2, max_M2 = 0, 0
min_COSI, max_COSI = 0, 0
px_1D, px_lo, px_hi = 0, 0, 0
dist_1D, dist_lo, dist_hi = 0, 0, 0
Tsun = 0

# define function to compute m2 as a function of 
def m2_massfunc(m2_in, m1_in, sini_in, mass_func):
    """
    Use Newton-Raphson's method to compute the companion mass, given a value of 
    the pulsar mass, the sine of the inclination angle, and the mass function.
    """
    m2b = m2_in
    for jj in range(100):
        g = (m2_in * sini_in)**3 / (m1_in + m2_in)**2 - mf
        dgdm2 = m2_in**2 * sini_in**3 * (m2_in + 3*m1_in) / (m1_in + m2_in)**3
        m2_in = m2_in - (g / dgdm2)
        if (np.fabs(m2_in - m2b) < 1e-7):
            break
        m2b = m2_in
    return m2_in

def compute_pdf2D_m1cosi(cosi_vals, m2_vals, m1_bins, mass_function, pdf2D_m2cosi):
    """
    Computes and interpolates a grid of M1-COSI from an input M2-COSI probability 
    density function.
    """

    # define some key variables and date structures.
    n_iter = len(cosi_vals)
    m1_bin_halfstep = (m1_bins[1] - m1_bins[0]) / 2
    pdf2D_m1cosi_mid = np.zeros((n_iter, n_iter))
    pdf2D_m1cosi_exp = np.zeros((n_iter, n_iter))

    # first, compute the "raw" M1-COSI grid.
    for x_cosi in range(n_iter):
        cosi = cosi_vals[x_cosi]
        sini = np.sqrt(1 - cosi**2)
        m1_row = np.sqrt((m2_vals * sini)**3 / mass_function) - m2_vals
        inds = np.digitize(m1_row, m1_bins) - 1
        inds_uniq = np.unique(inds)

        for idx in inds_uniq:
            if (idx >= 0 and idx < n_iter):
                m1_idx_good = np.where(idx == inds)
                m1_bin = m1_row[m1_idx_good]
                m2_bin = m2_vals[m1_idx_good]
                mtot_bin = m1_bin + m2_bin
                jacobian = np.fabs(1 / (3 * m2_bin**2 * sini**3 / 2 / mf / mtot_bin - 1))
                pdf2D_m1cosi_mid[idx, x_cosi] = np.mean(pdf2D_m2cosi[m1_idx_good, x_cosi] * jacobian)

    # next, interpolate the raw M1-COSI grid in order to smooth the PDF in case 
    # certain bins artificially have 0 probability density
    for l_x in range(n_iter):
        m1_slice = pdf2D_m1cosi_mid[:, l_x]
        idx_good = np.where(m1_slice != 0)
        m1_slice_good = m1_slice[idx_good]
        n_good = len(m1_slice_good)

        if (n_good > 1 and n_good < n_iter):
            m1_1D_good = m1_exp[idx_good]
            m1_slice_interpol = interp1d(m1_1D_good, m1_slice_good, bounds_error=False, fill_value=0)
            min_m1_slice, max_m1_slice = min(m1_exp[idx_good]), max(m1_exp[idx_good])
            m1_1D_new = np.arange(min_m1_slice, max_m1_slice, m1_bin_halfstep*2)
            pdf2D_m1cosi_exp[np.digitize(m1_1D_new, m1_bins)-1, l_x] = m1_slice_interpol(m1_1D_new)

    # delete intermediate, pre-interpolated grid from memory.
    del pdf2D_m1cosi_mid

    return pdf2D_m1cosi_exp

pin = open(infile[0])
GridDict = pickle.load(pin)

sini_bestfit = GridDict['SINI_bestfit']
m2_bestfit = GridDict['M2_bestfit']
chi2 = GridDict['chisq']
dof = GridDict['DOF']
mf = GridDict['massfunc']

# if the input grid is in H3 + H4/STIG, get those arrays;
# else, get the M2/COSI arrays.
if ('H3' in GridDict):
    h3_1D = GridDict['H3'] # FYI: the units are in micro-seconds.
    if ('STIG' in GridDict):
        stig_1D = GridDict['STIG']
    elif ('H4' in GridDict):
        h4_1D = GridDict['H4']
    Tsun = GridDict['Tsun']
else:
    cosi_1D = GridDict['COSI']
    m2_1D  = GridDict['M2']

# check if input grid is 3D, i.e. one axis is KOM.

if ('THETA' in GridDict):
    kom_1D = GridDict['THETA']

if ('DIST' in GridDict):
    dist_1D = GridDict['DIST'] 
    dist_lo = np.min(dist_1D)
    dist_hi = np.max(dist_1D)

if ('PX' in GridDict):
    px_1D = GridDict['PX'] 
    px_lo = np.min(px_1D)
    px_hi = np.max(px_1D)

# print out some info about the input grid.
cosi_bestfit = np.sqrt(1 - sini_bestfit**2)
m1_bestfit = np.sqrt((m2_bestfit * sini_bestfit)**3 / mf) - m2_bestfit
m_tot_bestfit = m1_bestfit + m2_bestfit

print "Input grid info for {}:".format(GridDict['PSR'])
print "    * best-fit chisq: {0:.3f}".format(np.min(chi2))
print "    * best-fit COSI:  {0:.3f}".format(cosi_bestfit)
print "    * best-fit M2:    {0:.3f}".format(m2_bestfit)
print "    * best-fit M1:    {0:.3f}".format(m1_bestfit)
print "    * mass function:  {0:.10f}".format(mf)
if ('THETA' in GridDict):
    if ('dist' in GridDict):
        print "    * grid is 3D (KOM as one axis)."
    else: 
        print "    * grid is 3D (KOM as one axis)."
gridtype = ''

# print out H3/STIG grid info, or M2/COSI info.

if ('H3' in GridDict):

    print "    * H3 min, max:   ({0:.3f}, {1:.3f})".format(min(h3_1D), max(h3_1D))

    if ('H4' in GridDict):
        print "    * H4 min, max:   ({0:.3f}, {1:.3f})".format(min(h4_1D), max(h4_1D))
        gridtype = 'h3h4'

    else:
        print "    * STIG min, max: ({0:.3f}, {1:.3f})".format(min(stig_1D), max(stig_1D))
        gridtype = 'h3stig'

else:

    print "    * M2 min, max:   ({0:.3f}, {1:.3f})".format(min(m2_1D), max(m2_1D))
    print "    * COSI min, max: ({0:.3f}, {1:.3f})".format(min(cosi_1D), max(cosi_1D))
    gridtype = 'rs'

if ('THETA' in GridDict):

    print "    * KOM min, max:  ({0:.3f}, {1:.3f})".format(min(kom_1D), max(kom_1D))
    gridtype = gridtype + '.fixXDOT'

# exit if only info within input-grid files are desired.

if (info_only):
    sys.exit()

print "Parameters for interpolated grids:"
print "    * COSI min, max = {0:.3f}, {1:.3f}".format(cosi_lo, cosi_hi)
print "    * M1 min, max   = {0:.3f}, {1:.3f}".format(m1_lo, m1_hi)
print "    * M2 min, max   = {0:.3f}, {1:.3f}".format(m2_lo, m2_hi)
if ('THETA' in GridDict):
    print "    * KOM min, max  = {0:.3f}, {1:.3f}".format(kom_lo, kom_hi)

# define common arrays and matrices that will be used for the contour 
# plots in the end, whether or not you start with the H3/STIG or M2/COSI maps.
m1_1D   = np.linspace(m1_lo, m1_hi, num=len(m2_1D))
m1_exp   = np.linspace(m1_lo, m1_hi, num=n_iter)
m2_exp   = np.linspace(m2_lo, m2_hi, num=n_iter)
px_exp   = np.linspace(px_lo, px_hi, num=n_iter)
kom_exp  = np.linspace(kom_lo, kom_hi, num=n_iter)
cosi_exp = np.linspace(cosi_lo, cosi_hi, num=n_iter)
px_exp = np.linspace(px_lo, px_hi, num=n_iter)
dist_exp = np.linspace(dist_lo, dist_hi, num=n_iter)
pdf2D_m2cosi_exp   = np.zeros((n_iter, n_iter))
pdf2D_m2cosi_orig  = np.zeros((n_iter, n_iter))
pdf2D_m2m1_mid     = np.zeros((n_iter, n_iter))
pdf2D_m2m1_exp     = np.zeros((n_iter, n_iter))
pdf2D_m1cosi_mid   = np.zeros((n_iter, n_iter))
pdf2D_m1cosi_exp   = np.zeros((n_iter, n_iter))
pdf2D_m1cosi_orig  = np.zeros((n_iter, n_iter))
pdf2D_komcosi_exp  = np.zeros((n_iter, n_iter))
pdf2D_komstig_exp  = np.zeros((n_iter, n_iter))
pdf2D_pxm1_exp     = np.zeros((n_iter, n_iter))
pdf2D_pxm2_exp     = np.zeros((n_iter, n_iter))
pdf2D_pxcosi_exp   = np.zeros((n_iter, n_iter))
pdf_px = 0

if ('THETA' in GridDict):
    pdf2D_komcosi_exp  = np.zeros((n_iter, n_iter))
    pdf2D_komstig_exp  = np.zeros((n_iter, n_iter))

pdf2D_m2cosi_interpol = 0

# set contour levels in chisq for two degrees of freedom.
# levels = 1, 2, and 3 sigma.
chisq = np.array([11.82915813, 6.1800743, 2.29574893])
pdf2D_chisq = np.exp(-0.5 * chisq)

# if grid is made using H3/STIG, treat this separately
# and interpolate this grid first; then convert to M2/COSI.
if ('H3' in GridDict and 'H4' in GridDict):

    deltachi2_h3h4 = chi2 - np.min(chi2)
    pdf2D_h3h4 = 0.5 * np.exp(-0.5 * deltachi2_h3h4)
    pdf2D_h3h4_area = np.sum(pdf2D_h3h4)
    pdf2D_h3h4 /= pdf2D_h3h4_area
    pdf2D_h3h4_interpol = interp2d(h4_1D, h3_1D, pdf2D_h3h4)

    h3_exp = np.linspace(min(h3_1D), max(h3_1D), num=n_iter)
    h4_exp = np.linspace(min(h4_1D), max(h4_1D), num=n_iter)
    pdf2D_h3h4_exp = np.zeros((n_iter, n_iter))

    # computed expanded grid from interpolated map.
    for l_h3 in range(n_iter):
        h3_new = h3_exp[l_h3]
        for l_h4 in range(n_iter):
            h4_new = h4_exp[l_h4]
            if (h4_new/h3_new > 1):
                break
            else: 
                pdf2D_h3h4_exp[l_h3, l_h4] = pdf2D_h3h4_interpol(h4_new, h3_new)

    pdf2D_h3h4_exp = pdf2D_h3h4_exp / np.max(pdf2D_h3h4_exp) * 0.5
    pdf_H3, pdf_H4 = np.zeros(n_iter), np.zeros(n_iter)

    # get PDFs for H3/H4 map.
    for m_pdf in range(n_iter):
        pdf_H3[m_pdf] = sum(pdf2D_h3h4_exp[m_pdf, :])
        pdf_H4[m_pdf] = sum(pdf2D_h3h4_exp[:, m_pdf])

    min_H3, max_H3 = np.min(h3_1D), np.max(h3_1D)
    min_STIG, max_STIG = np.min(h4_1D), np.max(h4_1D)

    # plot the H3/H4 map with the PDFs in slimmed panels.
    fig1 = plt.figure(1)
    fig1.set_figheight(8)
    fig1.set_figwidth(9)
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3]) 

    # configure H4 panel.
    ax0 = plt.subplot(gs[0,0])
    ax0.plot(h4_exp, pdf_H4 / np.max(pdf_H4), 'b-', lw=3)
    ax0.get_yaxis().set_visible(False)

    if (all(h4_confints)):
        h4_lowerlim, h4_median, h4_upperlim = h4_confints
        d1 = interp1d(h4_exp, pdf_H4 / np.max(pdf_H4))
        ax0.plot([h4_lowerlim, h4_lowerlim], [0., d1(h4_lowerlim)], 'r--', lw=2)
        ax0.plot([h4_median, h4_median], [0., d1(h4_median)], 'r-', lw=2)
        ax0.plot([h4_upperlim, h4_upperlim], [0., d1(h4_upperlim)], 'r--', lw=2)

    # coonfigure 2D-map panel.
    ax1 = plt.subplot(gs[1,0], sharex=ax0)
    ax1.pcolormesh(h4_exp, h3_exp, pdf2D_h3h4_exp, cmap='Blues', vmax=np.max(pdf2D_h3h4_exp), vmin=0.)
    ax1.contour(h4_exp, h3_exp, pdf2D_h3h4_exp, levels=pdf2D_chisq, colors='r', linewidths=1.5)
    ax1.fill_between(h4_exp, h4_exp, 0, color='green', alpha=0.5)
    ax1.set_xlabel(r'$h_4$ ($\mu$s)', fontproperties=font, fontsize=15)
    ax1.set_ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
    ax1.set_xlim([min(h4_exp), max(h4_exp)])
    ax1.set_ylim([min(h3_exp), max(h3_exp)])
    ax1.grid()
    ax2 = plt.subplot(gs[1,1], sharey=ax1)
    ax2.plot(pdf_H3 / np.max(pdf_H3), h3_exp, 'b-', lw=3)
    ax2.get_xaxis().set_visible(False)

    if (all(h3_confints)):
        h3_lowerlim, h3_median, h3_upperlim = h3_confints
        d1 = interp1d(h3_exp, pdf_H3 / np.max(pdf_H3))
        ax2.plot([0., d1(h3_lowerlim)], [h3_lowerlim, h3_lowerlim], 'r--', lw=2)
        ax2.plot([0., d1(h3_median)], [h3_median, h3_median], 'r-', lw=2)
        ax2.plot([0., d1(h3_upperlim)], [h3_upperlim, h3_upperlim], 'r--', lw=2)

    plt.figtext(0.74, 0.8, GridDict['PSR'], fontproperties=font, fontsize=15)
    plt.savefig('h3h4.'+GridDict['PSR']+'.png', format='png')

    # now convert the H3/H4 map to the M2/COSI map.
    # first, convert ortho map to H3-COSI space using the multi-bin method.
    pdf2D_h3cosi_mid = np.zeros((n_iter, n_iter))
    pdf2D_h3cosi_exp = np.zeros((n_iter, n_iter))
    cosi_bin_halfstep = (cosi_exp[1] - cosi_exp[0]) / 2 
    cosi_bin_loedge = cosi_exp[0] - cosi_bin_halfstep
    cosi_bin_hiedge = cosi_exp[len(cosi_exp)-1] + cosi_bin_halfstep

    cosi_exp_bins = np.linspace(cosi_bin_loedge, cosi_bin_hiedge, num=(n_iter+1))
    #cosi_row = (1 - h4_exp**2) / (1 + h4_exp**2)

    print "Computing the H3/COSI grid..."
   
    for n_h3 in range(n_iter):
        h3_row = h3_exp[n_h3]
        stig_row = h4_exp / h3_row
        cosi_row = (1. - stig_row**2) / (1. + stig_row**2)
        inds = np.digitize(cosi_row, cosi_exp_bins)-1
        inds_uniq = np.unique(inds)
        # loop over bin indices; this is needed if grid row covers more than one M1 bin.
        for idx in inds_uniq:
            if (idx >= 0 and idx < n_iter):
                cosi_idx_good = np.where(idx == inds)
                cosi_bin = cosi_row[cosi_idx_good]
                #if (len(cosi_bin) > 1):
                #    print cosi_bin
                stig_bin = np.sqrt((1. - cosi_bin) / (1. + cosi_bin))
                jacobian = 2 * h3_row * cosi_bin / stig_bin / (1. + cosi_bin**2)**2
                pdf2D_h3cosi_mid[n_h3, idx] = np.mean(pdf2D_h3h4_exp[n_h3, cosi_idx_good] * jacobian)

    pdf2D_h3cosi_mid = pdf2D_h3cosi_mid / np.max(pdf2D_h3cosi_mid) * 0.5

    #plt.figure(2)
    #plt.pcolormesh(cosi_exp, h3_exp, pdf2D_h3cosi_mid)
    #plt.colorbar()

    # no interpolate points where there should be probability.
    for o_h3 in range(n_iter):
        h3cosi_slice = pdf2D_h3cosi_mid[o_h3, :]
        idx_good = np.where(h3cosi_slice != 0)
        n_good = len(h3cosi_slice[idx_good])
        if (n_good != 0 and n_good > 1):
            min_cosi_slice, max_cosi_slice = min(cosi_exp[idx_good]), max(cosi_exp[idx_good])
            n_steps = int((max_cosi_slice - min_cosi_slice) / cosi_bin_halfstep / 2)
            cosi_exp_new = np.arange(min_cosi_slice, max_cosi_slice, cosi_bin_halfstep*2)
            h3cosi_slice_interpol = interp1d(cosi_exp[idx_good], h3cosi_slice[idx_good], bounds_error=False, fill_value=0.)
            pdf2D_h3cosi_exp[o_h3, np.digitize(cosi_exp_new, cosi_exp_bins)-1] = h3cosi_slice_interpol(cosi_exp_new)
        else:
            pdf2D_h3cosi_exp[o_h3, :] = h3cosi_slice
 
    pdf2D_h3cosi_exp = pdf2D_h3cosi_exp / np.max(pdf2D_h3cosi_exp) * 0.5

    #plt.figure(3)
    #plt.pcolormesh(cosi_exp, h3_exp, pdf2D_h3cosi_exp)
    #plt.colorbar()
    #plt.show()

    # now convert H3/COSI map to M2/COSI using the same method.
    pdf2D_m2cosi_mid = np.zeros((n_iter, n_iter))
    m2_bin_halfstep = (m2_exp[1] - m2_exp[0]) / 2
    m2_bin_loedge = m2_exp[0] - m2_bin_halfstep
    m2_bin_hiedge = m2_exp[len(m2_exp)-1] + m2_bin_halfstep
    m2_exp_bins = np.linspace(m2_bin_loedge, m2_bin_hiedge, num=(n_iter+1))

    print "Computing the M2/COSI grid..."

    for o_cosi in range(n_iter):
        cosi = cosi_exp[o_cosi]
        stig = np.sqrt((1 - cosi) / (1 + cosi))
        m2_row = h3_exp / (Tsun * 1e6) / stig**3
        inds = np.digitize(m2_row, m2_exp_bins)-1
        inds_uniq = np.unique(inds)
        # loop over bin indices; this is needed if grid row covers more than one M1 bin.
        for idx in inds_uniq:
            if (idx >= 0 and idx < n_iter):
                m2_idx_good = np.where(idx == inds)
                jacobian = (Tsun * 1e6) * stig**3
                pdf2D_m2cosi_mid[idx, o_cosi] = np.mean(pdf2D_h3cosi_exp[m2_idx_good, o_cosi] * jacobian)

    #plt.figure(3)
    #plt.pcolormesh(cosi_exp, m2_exp, pdf2D_m2cosi_mid)
    #plt.colorbar()
    #plt.show()

    #sys.exit()

    for o_cosi in range(n_iter):
        m2cosi_slice = pdf2D_m2cosi_mid[:, o_cosi]
        idx_good = np.where(m2cosi_slice != 0)
        if (len(m2cosi_slice[idx_good]) > 1):
            min_m2_slice, max_m2_slice = min(m2_exp[idx_good]), max(m2_exp[idx_good])
            n_steps = int((max_m2_slice - min_m2_slice) / m2_bin_halfstep / 2)
            m2_exp_new = np.arange(min_m2_slice, max_m2_slice, m2_bin_halfstep*2)
            m2cosi_slice_interpol = interp1d(m2_exp[idx_good], m2cosi_slice[idx_good], bounds_error=False, fill_value=0.)
            pdf2D_m2cosi_exp[np.digitize(m2_exp_new, m2_exp_bins)-1, o_cosi] = m2cosi_slice_interpol(m2_exp_new)
        else:
            pdf2D_m2cosi_exp[:, o_cosi] = pdf2D_m2cosi_mid[:, o_cosi]

    max_m2cosi = np.max(pdf2D_m2cosi_exp)
    pdf2D_m2cosi_exp = pdf2D_m2cosi_exp / max_m2cosi * 0.5


if ('H3' in GridDict and 'STIG' in GridDict):

    #deltachi2_h3stig = (chi2 - np.min(chi2)) * dof
    deltachi2_h3stig = chi2 - np.min(chi2)
    pdf2D_h3stig  = np.zeros((len(h3_1D), len(stig_1D)))
    pdf2D_komstig = 0

    if ('THETA' in GridDict):

        pdf2D_komstig = np.zeros((len(kom_1D), len(stig_1D)))
        pdf3D = np.exp(-0.5 * deltachi2_h3stig)
        
        # seperate 3D into two 2D grids: m2-cosi and kom-cosi.

        for ii in range(len(h3_1D)):
            for jj in range(len(stig_1D)):
                pdf2D_h3stig[ii, jj] = np.sum(pdf3D[ii, jj, :])

        for ii in range(len(kom_1D)):
            for jj in range(len(stig_1D)):
                pdf2D_komstig[ii, jj] = np.sum(pdf3D[:, jj, ii])

    else:

        pdf2D_h3stig = 0.5*np.exp(-0.5 * deltachi2_h3stig)

    pdf2D_h3stig_area = np.sum(pdf2D_h3stig)
    pdf2D_h3stig /= pdf2D_h3stig_area
    pdf2D_h3stig_interpol = interp2d(stig_1D, h3_1D, pdf2D_h3stig)

    h3_exp = np.linspace(min(h3_1D), max(h3_1D), num=n_iter)
    stig_exp = np.linspace(min(stig_1D), max(stig_1D), num=n_iter)
    pdf2D_h3stig_exp = np.zeros((n_iter, n_iter))

    # computed expanded grid from interpolated map.
    for l_h3 in range(n_iter):
        h3_new = h3_exp[l_h3]
        for l_stig in range(n_iter):
            stig_new = stig_exp[l_stig]
            pdf2D_h3stig_exp[l_h3, l_stig] = pdf2D_h3stig_interpol(stig_new, h3_new)

    pdf2D_h3stig_exp = pdf2D_h3stig_exp / np.max(pdf2D_h3stig_exp) * 0.5
    pdf_H3, pdf_STIG = np.zeros(n_iter), np.zeros(n_iter)

    # get PDFs for H3/STIG map.
    for m_pdf in range(n_iter):
        pdf_H3[m_pdf] = sum(pdf2D_h3stig_exp[m_pdf, :])
        pdf_STIG[m_pdf] = sum(pdf2D_h3stig_exp[:, m_pdf])

    min_H3, max_H3 = np.min(h3_1D), np.max(h3_1D)
    min_STIG, max_STIG = np.min(stig_1D), np.max(stig_1D)

    stig_1D_exp = np.linspace(min(stig_1D), max(stig_1D), num=n_iter)
    h3_1D_exp = np.linspace(min(h3_1D), max(h3_1D), num=n_iter)

    # plot the H3/H4 map with the PDFs as slimmed panels.
    fig1 = plt.figure(1)
    fig1.set_figheight(8)
    fig1.set_figwidth(9)
    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 3])
    ax0 = plt.subplot(gs[0,0])
    ax0.plot(stig_exp, pdf_STIG / np.max(pdf_STIG), 'b-', lw=3)
    ax0.get_yaxis().set_visible(False)

    if (all(stig_confints)):
        stig_lowerlim, stig_median, stig_upperlim = stig_confints
        d1 = interp1d(stig_exp, pdf_STIG / np.max(pdf_STIG))
        ax0.plot([stig_lowerlim, stig_lowerlim], [0., d1(stig_lowerlim)], 'r--', lw=2)
        ax0.plot([stig_median, stig_median], [0., d1(stig_median)], 'r-', lw=2)
        ax0.plot([stig_upperlim, stig_upperlim], [0., d1(stig_upperlim)], 'r--', lw=2)

    ax1 = plt.subplot(gs[1,0], sharex=ax0)
    ax1.pcolormesh(stig_1D_exp, h3_1D_exp, pdf2D_h3stig_exp, cmap='Blues')
    ax1.contour(stig_1D_exp, h3_1D_exp, pdf2D_h3stig_exp, levels=pdf2D_chisq, colors='r', lw=2)
    ax1.set_xlabel(r'$\zeta$', fontsize=15)
    ax1.set_ylabel(r'$h_3$ ($\mu$s)', fontproperties=font, fontsize=15)
    ax1.set_xlim([min(stig_exp), max(stig_exp)])
    ax1.set_ylim([min(h3_exp), max(h3_exp)])
    ax1.grid()
    ax2 = plt.subplot(gs[1,1], sharey=ax1)
    ax2.plot(pdf_H3 / np.max(pdf_H3), h3_exp, 'b-', lw=3)
    ax2.get_xaxis().set_visible(False)

    if (all(h3_confints)):
        h3_lowerlim, h3_median, h3_upperlim = h3_confints
        d1 = interp1d(h3_exp, pdf_H3 / np.max(pdf_H3))
        ax2.plot([0., d1(h3_lowerlim)], [h3_lowerlim, h3_lowerlim], 'r--', lw=2)
        ax2.plot([0., d1(h3_median)], [h3_median, h3_median], 'r-', lw=2)
        ax2.plot([0., d1(h3_upperlim)], [h3_upperlim, h3_upperlim], 'r--', lw=2)

    plt.savefig('h3stig.'+GridDict['PSR']+'.png', format='png')

    # now convert the H3/STIG map to the M2/COSI map.
    # first, convert ortho map to H3-COSI space using the multi-bin method.
    pdf2D_h3cosi_mid = np.zeros((n_iter, n_iter))
    pdf2D_h3cosi_exp = np.zeros((n_iter, n_iter))
    cosi_bin_halfstep = (cosi_exp[1] - cosi_exp[0]) / 2
    cosi_bin_loedge = cosi_exp[0] - cosi_bin_halfstep
    cosi_bin_hiedge = cosi_exp[len(cosi_exp)-1] + cosi_bin_halfstep

    cosi_exp_bins = np.linspace(cosi_bin_loedge, cosi_bin_hiedge, num=(n_iter+1))
    cosi_row = (1 - stig_exp**2) / (1 + stig_exp**2)

    print "Computing the H3/COSI grid..."
   
    inds = np.digitize(cosi_row, cosi_exp_bins)-1
    inds_uniq = np.unique(inds)

    for n_h3 in range(n_iter):
        # loop over bin indices; this is needed if grid row covers more than one M1 bin.
        for idx in inds_uniq:
            if (idx >= 0 and idx < n_iter):
                cosi_idx_good = np.where(idx == inds)
                cosi_bin = cosi_row[cosi_idx_good]
                stig_bin = np.sqrt((1 - cosi_bin) / (1 + cosi_bin))
                jacobian = cosi_bin / (stig_bin * (1 + cosi_bin**2)**2)
                pdf2D_h3cosi_mid[n_h3, idx] = np.mean(pdf2D_h3stig_exp[n_h3, cosi_idx_good] * jacobian)

    # interpolate empty points in the transformed/expanded grid.
    for o_h3 in range(n_iter):
        h3cosi_slice = pdf2D_h3cosi_mid[o_h3, :]
        idx_good = (np.where(h3cosi_slice != 0))[0]
        #if (len(h3cosi_slice[idx_good]) != 0):
        if (len(h3cosi_slice[idx_good]) > 1):
            min_cosi_slice, max_cosi_slice = min(cosi_exp[idx_good]), max(cosi_exp[idx_good])
            n_steps = int((max_cosi_slice - min_cosi_slice) / cosi_bin_halfstep / 2)
            cosi_exp_new = np.arange(min_cosi_slice, max_cosi_slice, cosi_bin_halfstep*2)
            h3cosi_slice_interpol = interp1d(cosi_exp[idx_good], h3cosi_slice[idx_good], bounds_error=False, fill_value=0.)
            pdf2D_h3cosi_exp[o_h3, np.digitize(cosi_exp_new, cosi_exp_bins)-1] = h3cosi_slice_interpol(cosi_exp_new)
        else:
            pdf2D_h3cosi_exp[o_h3, :] = h3cosi_slice

    #pdf2D_h3cosi_exp = pdf2D_h3cosi_mid
    max_h3cosi = np.max(pdf2D_h3cosi_exp)
    pdf2D_h3cosi_exp = pdf2D_h3cosi_exp / max_h3cosi * 0.5

    # if 3D, compute the transformed KOM/COSI grid.
    if ('THETA' in GridDict):
        
        min_KOM, max_KOM = np.min(kom_1D), np.max(kom_1D)
        min_STIG, max_STIG = np.min(stig_exp), np.max(stig_exp)

        pdf2D_komstig_int = interp2d(stig_1D, kom_1D, pdf2D_komstig)
        pdf2D_komstig_exp = np.zeros((n_iter, n_iter))
        pdf2D_komcosi_mid = np.zeros((n_iter, n_iter))

        for l_kom in range(n_iter):
            kom = kom_exp[l_kom]
            if (kom > 360. and kom < 720.): kom -= 360.
            for m_stig in range(n_iter):
                stig = stig_exp[m_stig]
                if (kom > min_KOM and kom < max_KOM and stig > min_STIG and stig < max_STIG):
                    pdf2D_komstig_exp[l_kom, m_stig] = pdf2D_komstig_int(stig, kom)
                else:
                    continue

        for n_kom in range(n_iter):
        # loop over bin indices; this is needed if grid row covers more than one M1 bin.
            for idx in inds_uniq:
                if (idx >= 0 and idx < n_iter):
                    cosi_idx_good = np.where(idx == inds)
                    cosi_bin = cosi_row[cosi_idx_good]
                    stig_bin = np.sqrt((1 - cosi_bin) / (1 + cosi_bin))
                    jacobian = cosi_bin / (stig_bin * (1 + cosi_bin**2)**2)
                    pdf2D_komcosi_mid[n_kom, idx] = np.mean(pdf2D_komstig_exp[n_kom, cosi_idx_good] * jacobian)

        # interpolate empty points in the transformed/expanded grid.
        for o_kom in range(n_iter):
            komcosi_slice = pdf2D_komcosi_mid[o_kom, :]
            idx_good = np.where(komcosi_slice != 0)
            if (len(komcosi_slice[idx_good]) != 0):
                min_cosi_slice, max_cosi_slice = min(cosi_exp[idx_good]), max(cosi_exp[idx_good])
                n_steps = int((max_cosi_slice - min_cosi_slice) / cosi_bin_halfstep / 2)
                cosi_exp_new = np.arange(min_cosi_slice, max_cosi_slice, cosi_bin_halfstep*2)
                komcosi_slice_interpol = interp1d(cosi_exp[idx_good], komcosi_slice[idx_good], bounds_error=False, fill_value=0.)
                pdf2D_komcosi_exp[o_kom, np.digitize(cosi_exp_new, cosi_exp_bins)-1] = komcosi_slice_interpol(cosi_exp_new)
            else:
                pdf2D_komcosi_exp[o_kom, :] = komcosi_slice

        pdf2D_komcosi_exp = pdf2D_komcosi_exp / np.max(pdf2D_komcosi_exp) * 0.5
        #pdf2D_komcosi_exp = pdf2D_komcosi_exp / np.sum(pdf2D_komcosi_exp)


    # now convert H3/COSI map to M2/COSI using the same method.
    pdf2D_m2cosi_mid = np.zeros((n_iter, n_iter))
    m2_bin_halfstep = (m2_exp[1] - m2_exp[0]) / 2
    m2_bin_loedge = m2_exp[0] - m2_bin_halfstep
    m2_bin_hiedge = m2_exp[len(m2_exp)-1] + m2_bin_halfstep
    m2_exp_bins = np.linspace(m2_bin_loedge, m2_bin_hiedge, num=(n_iter+1))

    print "Computing the M2/COSI grid..."

    for o_cosi in range(n_iter):
        cosi = cosi_exp[o_cosi]
        stig = np.sqrt((1 - cosi) / (1 + cosi))
        m2_row = h3_exp / (Tsun * 1e6) / stig**3
        inds = np.digitize(m2_row, m2_exp_bins)-1
        inds_uniq = np.unique(inds)
        # loop over bin indices; this is needed if grid row covers more than one M1 bin.
        for idx in inds_uniq:
            if (idx >= 0 and idx < n_iter):
                m2_idx_good = np.where(idx == inds)
                jacobian = (Tsun * 1e6) * stig**3
                pdf2D_m2cosi_mid[idx, o_cosi] = np.mean(pdf2D_h3cosi_exp[m2_idx_good, o_cosi] * jacobian)

    #plt.figure(3)
    #plt.pcolormesh(cosi_exp, m2_exp, pdf2D_m2cosi_mid)
    #plt.colorbar()
    #plt.show()
    #3sys.exit()    

    for o_cosi in range(n_iter):
        m2cosi_slice = pdf2D_m2cosi_mid[:, o_cosi]
        idx_good = np.where(m2cosi_slice != 0)
        if (len(m2cosi_slice[idx_good]) > 1):
            min_m2_slice, max_m2_slice = min(m2_exp[idx_good]), max(m2_exp[idx_good])
            n_steps = int((max_m2_slice - min_m2_slice) / m2_bin_halfstep / 2)
            m2_exp_new = np.arange(min_m2_slice, max_m2_slice, m2_bin_halfstep*2)
            m2cosi_slice_interpol = interp1d(m2_exp[idx_good], m2cosi_slice[idx_good], bounds_error=False, fill_value=0.)
            pdf2D_m2cosi_exp[np.digitize(m2_exp_new, m2_exp_bins)-1, o_cosi] = m2cosi_slice_interpol(m2_exp_new)
        else:
            pdf2D_m2cosi_exp[:, o_cosi] = pdf2D_m2cosi_mid[:, o_cosi]

    #pdf2D_m2cosi_exp = pdf2D_m2cosi_mid
    max_m2cosi = np.max(pdf2D_m2cosi_exp)
    pdf2D_m2cosi_exp = pdf2D_m2cosi_exp / max_m2cosi * 0.5
    
    #plt.figure(4)
    #plt.pcolormesh(cosi_exp, h3_exp, pdf2D_m2cosi_exp)
    #plt.colorbar()
    #plt.show()
    #sys.exit() 

    # for ease in computing M2/M1 grid, set the min/max values 
    # to those from the input/default values.
    
    min_M2, max_M2 = np.min(m2_exp), np.max(m2_exp)
    min_COSI, max_COSI = np.min(cosi_exp), np.max(cosi_exp)

elif ("M2" in GridDict and "COSI" in GridDict):

    min_M2, max_M2 = np.min(m2_1D), np.max(m2_1D)
    min_COSI, max_COSI = np.min(cosi_1D), np.max(cosi_1D)

    # get probability map in m2-cosi and marginalized PDFs.

    pdf2D_m2cosi = np.zeros((len(m2_1D), len(cosi_1D)))    
    #deltachi2_m2cosi_orig = (chi2 - np.min(chi2)) * dof
    deltachi2_m2cosi_orig = chi2 - np.min(chi2)

    # if grid is 3D, obtain the KOM-COSI map.
    if (len(deltachi2_m2cosi_orig.shape) == 3):

        # now compute 2D grid for M2 and COSI.
        pdf3D = np.exp(-0.5 * deltachi2_m2cosi_orig)
        pdf2D_m2cosi = np.sum(pdf3D, axis=2)
        pdf2D_komcosi = 0

        #if ('THETA' in GridDict):
        #    pdf3D = np.zeros((len(m2_1D), len(cosi_1D), len(kom_1D)))
        #    pdf2D_komcosi = np.zeros((len(kom_1D), len(cosi_1D)))
        #elif ('PX' in GridDict):
        #    pdf3D = np.zeros((len(m2_1D), len(cosi_1D), len(px_1D)))
        #elif ('DIST' in GridDict):
        #    pdf3D = np.zeros((len(m2_1D), len(cosi_1D), len(px_1D)))
           
        # if PX or DIST is the third dimension, compute the appropriate slices.
        # for the sake of reducing lines of code, i'm using the same variables 
        # to allocate either set of grids (i.e., PX or DICT) but they are 
        # differentiated by keyword name when saving to NumPy state file.
        if ('PX' in GridDict or 'DIST' in GridDict):

            if ('DIST' in GridDict):
                px_1D = dist_1D
                px_exp = dist_exp

            print "Computing 1D PDF for PX/DIST..."
            # compute marginalized PDF for PX.
            pdf_px_orig = np.sum(pdf3D, axis=(0, 1))
            pdf_px_orig /= np.sum(pdf_px_orig)
            pdf_px_interp = interp1d(px_1D, pdf_px_orig)
            pdf_px = pdf_px_interp(px_exp)
            pdf_px /= np.sum(pdf_px)

            # now compute interolated grids for M2-PX.
            pdf2D_current = np.sum(pdf3D, axis=1)
            pdf2D_current /= np.sum(pdf2D_current)
            pdf2D_interpol = interp2d(px_1D, m2_1D, pdf2D_current)

            print "Computing 2D PDF for PX/DIST-M2..."
            for idx_px in range(n_iter):
                for idx_m2 in range(n_iter):
                    pdf2D_pxm2_exp[idx_m2, idx_px] = pdf2D_interpol(px_exp[idx_px], m2_exp[idx_m2])

            pdf2D_pxm2_exp /= np.sum(pdf2D_pxm2_exp)

            # now compute interolated grids for COSI-PX.
            pdf2D_current = np.sum(pdf3D, axis=0)
            pdf2D_current /= np.sum(pdf2D_current)
            pdf2D_interpol = interp2d(px_1D, cosi_1D, pdf2D_current)

            print "Computing 2D PDF for PX/DIST-COSI..."
            for idx_px in range(n_iter):
                for idx_cosi in range(n_iter):
                    pdf2D_pxcosi_exp[idx_m2, idx_px] = pdf2D_interpol(px_exp[idx_px], cosi_exp[idx_cosi])

            pdf2D_pxcosi_exp /= np.sum(pdf2D_pxcosi_exp)

            # finally, compute interpolated grids for M1-PX.
            # currently, I first transform each M2-COSI grid at a given PX bin.
            pdf3D_m1cosipx = np.zeros(pdf3D.shape) 
            m1_bin_halfstep_orig = (m1_1D[1] - m1_1D[0]) / 2
            m1_bins_orig = np.linspace(m1_1D[0] - m1_bin_halfstep_orig, m1_1D[-1] + m1_bin_halfstep_orig, num=len(m1_1D)+1)

            print "Computing 3D PDF for M1-COSI-PX/DIST..."
            for idx_px in range(len(px_1D)):
                pdf2D_current = pdf3D[:, :, idx_px]
                pdf3D_m1cosipx[:, :, idx_px] = compute_pdf2D_m1cosi(cosi_1D, m2_1D, m1_bins_orig, mf, pdf2D_current)

            # then I extract the marginalized M1-PX grid and interpolate it.
            pdf2D_current = np.sum(pdf3D_m1cosipx, axis=1)
            pdf2D_current /= np.sum(pdf2D_current)
            pdf2D_interpol = interp2d(px_1D, cosi_1D, pdf2D_current)

            print "Computing 2D PDF for PX-M1..."
            for idx_px in range(n_iter):
                for idx_m1 in range(n_iter):
                    pdf2D_pxm1_exp[idx_m1, idx_px] = pdf2D_interpol(px_exp[idx_px], m1_exp[idx_m1])

            pdf2D_pxm1_exp /= np.sum(pdf2D_pxm2_exp)

        elif ('THETA' in GridDict):

            # seperate 3D into two 2D grids: m2-cosi and kom-cosi.
            pdf2D_m2kom = np.zeros((len(m2_1D), len(kom_1D)))

            for ii in range(len(kom_1D)):
                for jj in range(len(cosi_1D)):
                    pdf2D_komcosi[ii, jj] = np.sum(pdf3D[:, jj, ii])

            for ii in range(len(m2_1D)):
                for jj in range(len(kom_1D)):
                    pdf2D_m2kom[ii, jj] = np.sum(pdf3D[ii, :, jj])

            pdf2D_interpol = interp2d(cosi_1D, kom_1D, pdf2D_komcosi)
            min_KOM, max_KOM = np.min(kom_1D), np.max(kom_1D)

            for x_kom in range(n_iter):
                kom = kom_exp[x_kom]
                if (max_KOM <= 360. and kom > 360.):
                    kom -= 360.
                for x_cosi in range(n_iter):
                    cosi = cosi_exp[x_cosi]
                    if (kom > min_KOM and kom < max_KOM and cosi > min_COSI and cosi < max_COSI): 
                        pdf2D_komcosi_exp[x_kom, x_cosi] = pdf2D_interpol(cosi, kom)
                else:
                    continue

            pdf2D_komcosi_exp = pdf2D_komcosi_exp / np.max(pdf2D_komcosi_exp) * 0.5
    
    else:
        pdf2D_m2cosi = np.exp(-0.5 * deltachi2_m2cosi_orig)

    print "Interpolating probability map..."
    pdf2D_m2cosi_interpol = interp2d(cosi_1D, m2_1D, pdf2D_m2cosi)
    print "Computing probability map with higher resolution..."
    
    # now compute high-res m2-cosi map from interpolated grid.

    for k_m2 in range(len(m2_exp)):
        m2 = m2_exp[k_m2]
        for k_cosi in range(len(cosi_exp)):
            cosi = cosi_exp[k_cosi]
            if (m2 > min_M2 and m2 < max_M2 and cosi > min_COSI and cosi < max_COSI):
                pdf2D_m2cosi_exp[k_m2, k_cosi] = pdf2D_m2cosi_interpol(cosi, m2)
            else: 
                continue

pdf2D_m2cosi_orig = pdf2D_m2cosi_exp

# prior to transformation, apply bayes theorem using prior map if the latter isn't uniform.
if all(m2_priorint):
    prior_m2cosi_exp = np.zeros((n_iter, n_iter))
    m2_prior, cosi_prior = np.meshgrid(m2_exp, cosi_exp)
    prior2D_m2cosi_exp = np.exp(-0.5 * (m2_prior - m2_priorint[0])**2 / (m2_priorint[1])**2)
    pdf2D_m2cosi_exp *= np.transpose(prior2D_m2cosi_exp)

area_m2cosi = np.sum(pdf2D_m2cosi_exp)
max_m2cosi = np.max(pdf2D_m2cosi_exp)
pdf2D_m2cosi_exp = pdf2D_m2cosi_exp / max_m2cosi * 0.5 

print "Computing high-resolution grid in m2-m1..."
m1_bin_halfstep = (m1_exp[1] - m1_exp[0]) / 2
m1_bin_loedge = m1_exp[0] - m1_bin_halfstep
m1_bin_hiedge = m1_exp[len(m1_exp)-1] + m1_bin_halfstep

m1_exp_bins = np.linspace(m1_bin_loedge, m1_bin_hiedge, num=(n_iter+1))
sini_row = np.sqrt(1 - cosi_exp**2)
tani_row = sini_row / cosi_exp

# the following for-loops convert the M2-COSI map.

for j_m2 in range(len(m2_exp)):
    m2 = m2_exp[j_m2]
    m1_row = np.sqrt((m2 * sini_row)**3 / mf) - m2
    inds = np.digitize(m1_row, m1_exp_bins)-1
    inds_uniq = np.unique(inds)
    # loop over bin indices; this is needed if grid row covers more than one M1 bin.
    for idx in inds_uniq:
        if (idx >= 0 and idx < n_iter):
            m1_idx_good = np.where(idx == inds) 
            m1_bin = m1_row[m1_idx_good]
            tani_bin = tani_row[m1_idx_good]
            jacobian = np.fabs(2*mf/3 * tani_bin * (1 + m1_bin/m2) * (mf * (m1_bin + m2)**2)**(-2./3.))
            pdf2D_m2m1_mid[j_m2, idx] = np.mean(pdf2D_m2cosi_exp[j_m2, m1_idx_good] * jacobian)

# compute M1-COSI grid as well, since the M1 probability is less smeared out
# in this space.

print "Computing high-resolution grid in m1-cosi..."

for x_cosi in range(n_iter):
    cosi = cosi_exp[x_cosi]
    sini = np.sqrt(1 - cosi**2)
    m1_row = np.sqrt((m2_exp * sini)**3 / mf) - m2_exp
    inds = np.digitize(m1_row, m1_exp_bins) - 1
    inds_uniq = np.unique(inds)
    for idx in inds_uniq:
        if (idx >= 0 and idx < n_iter):
            m1_idx_good = np.where(idx == inds)
            m1_bin = m1_row[m1_idx_good]
            m2_bin = m2_exp[m1_idx_good]
            mtot_bin = m1_bin + m2_bin
            jacobian = np.fabs(1 / (3 * m2_bin**2 * sini**3 / 2 / mf / mtot_bin - 1))
            pdf2D_m1cosi_mid[idx, x_cosi] = np.mean(pdf2D_m2cosi_exp[m1_idx_good, x_cosi] * jacobian)

# interpolate both computed grids in case there are points that 
# should have probability but don't due to finite bin widths.

print "Interpolating both grids..."

for o_m2 in range(n_iter):
    m2m1_slice = pdf2D_m2m1_mid[o_m2, :]
    idx_good = np.where(m2m1_slice != 0)
    n_good = len(m2m1_slice[idx_good])
    if (n_good > 1 and n_good < n_iter):
        min_m1_slice, max_m1_slice = min(m1_exp[idx_good]), max(m1_exp[idx_good])
        n_steps = int((max_m1_slice - min_m1_slice) / m1_bin_halfstep / 2)
        m1_exp_new = np.arange(min_m1_slice, max_m1_slice, m1_bin_halfstep*2)
        m2m1_slice_interpol = interp1d(m1_exp[idx_good], m2m1_slice[idx_good], bounds_error=False, fill_value=0.)
        pdf2D_m2m1_exp[o_m2, np.digitize(m1_exp_new, m1_exp_bins)-1] = m2m1_slice_interpol(m1_exp_new)
    else:
        pdf2D_m2m1_exp[o_m2, :] = pdf2D_m2m1_mid[o_m2, :]

max_m2m1 = np.max(pdf2D_m2m1_exp)
pdf2D_m2m1_exp = pdf2D_m2m1_exp / max_m2m1 * 0.5

for l_x in range(n_iter):
    m1_slice = pdf2D_m1cosi_mid[:, l_x]
    idx_good = np.where(m1_slice != 0)
    m1_slice_good = m1_slice[idx_good]
    n_good = len(m1_slice_good)
    if (n_good > 1 and n_good < n_iter):
        m1_1D_good = m1_exp[idx_good]
        m1_slice_interpol = interp1d(m1_1D_good, m1_slice_good, bounds_error=False, fill_value=0)
        min_m1_slice, max_m1_slice = min(m1_exp[idx_good]), max(m1_exp[idx_good])
        m1_1D_new = np.arange(min_m1_slice, max_m1_slice, m1_bin_halfstep*2)
        pdf2D_m1cosi_exp[np.digitize(m1_1D_new, m1_exp_bins)-1, l_x] = m1_slice_interpol(m1_1D_new)

pdf2D_m1cosi_exp = compute_pdf2D_m1cosi(cosi_exp, m2_exp, m1_exp_bins, mf, pdf2D_m2cosi_exp)

max_m1cosi = np.max(pdf2D_m1cosi_exp)
pdf2D_m1cosi_exp = pdf2D_m1cosi_exp / max_m1cosi * 0.5

pdf_M1, pdf_M2, pdf_cosi = [], [], []
pdf_M1_m1cosi = []
pdf_KOM = []

# compute 1D PDFs and normalize to have unit area.
for elem in range(n_iter):
    pdf_M2.append(sum(pdf2D_m2cosi_exp[elem, :]))
    pdf_cosi.append(sum(pdf2D_m2cosi_exp[:, elem]))
    pdf_M1.append(sum(pdf2D_m2m1_exp[:, elem]))
    pdf_M1_m1cosi.append(sum(pdf2D_m1cosi_exp[elem, :]))
    if ('THETA' in GridDict):
        pdf_KOM.append(sum(pdf2D_komcosi_exp[elem, :]))

pdf_M1 = np.array(pdf_M1)
pdf_M1 /= np.sum(pdf_M1)
pdf_M1_m1cosi = np.array(pdf_M1_m1cosi)
pdf_M1_m1cosi /= np.sum(pdf_M1_m1cosi)
pdf_M2 = np.array(pdf_M2)
pdf_M2 /= np.sum(pdf_M2)
pdf_cosi = np.array(pdf_cosi)
pdf_cosi /= np.sum(pdf_cosi)

npanels = 2
width   = 12
ratios  = [3, 3, 1]

if ('THETA' in GridDict):
    npanels = 3
    width = 15
    ratios = [3, 3, 3, 1]
    pdf_KOM = np.array(pdf_KOM)
    pdf_KOM /= np.sum(pdf_KOM)

# the following commented bit was for creating tinted plots that 
# covered unphysical regions in the M1/M2 space.

#m2lim = np.zeros(len(m1_exp))
#for k_lim in range(len(m2lim)):
#    m2lim[k_lim] = m2_massfunc(0.3, m1_exp[k_lim], 1., mf)

#m1lim = np.sqrt(mf / np.sqrt(1 - cosi_exp**2)**3)

fig = plt.figure(2)
fig.set_figheight(6)
fig.set_figwidth(width)
gs = gridspec.GridSpec(2, npanels+1, width_ratios=ratios, height_ratios=[1, 3])

ax0 = plt.subplot(gs[0,0])
ax0.plot(m2_exp, pdf_M2 / np.max(pdf_M2), 'b-', lw=3)
if (all(m2_confints)):
    m2_lowerlim, m2_median, m2_upperlim = m2_confints
    d1 = interp1d(m2_exp, pdf_M2 / np.max(pdf_M2))
    ax0.plot([m2_lowerlim, m2_lowerlim], [0., d1(m2_lowerlim)], 'r--', lw=2)
    ax0.plot([m2_median, m2_median], [0., d1(m2_median)], 'r-', lw=2)
    ax0.plot([m2_upperlim, m2_upperlim], [0., d1(m2_upperlim)], 'r--', lw=2)
ax0.set_ylim([0., 1.5])
ax0.get_yaxis().set_visible(False)

ax1 = plt.subplot(gs[1,0], sharex=ax0)
ax1.contour(m2_exp, cosi_exp, np.transpose(pdf2D_m2cosi_exp), levels=pdf2D_chisq, colors='r', lw=3)
ax1.pcolormesh(m2_exp, cosi_exp, np.transpose(pdf2D_m2cosi_exp), cmap='Blues', vmax=np.max(pdf2D_m2cosi_exp), vmin=0.)
ax1.grid()
ax1.set_ylabel(r'$\cos i$', fontproperties=font, fontsize=15)
ax1.set_xlabel(r'Companion Mass (M$_{\odot}$)', fontproperties=font, fontsize=15)
ax1.set_ylim([min(cosi_exp), max(cosi_exp)])

ax2 = plt.subplot(gs[1,2], sharey=ax1)
if ('THETA' in GridDict):
    ax2 = plt.subplot(gs[1,3], sharey=ax1)
ax2.plot(pdf_cosi / np.max(pdf_cosi), cosi_exp, 'b-', lw=3)
if (all(cosi_confints)):
    cosi_lowerlim, cosi_median, cosi_upperlim = cosi_confints
    d1 = interp1d(cosi_exp, pdf_cosi / np.max(pdf_cosi))
    ax2.plot([0., d1(cosi_lowerlim)], [cosi_lowerlim, cosi_lowerlim], 'r--', lw=2)
    ax2.plot([0., d1(cosi_median)], [cosi_median, cosi_median], 'r-', lw=2)
    ax2.plot([0., d1(cosi_upperlim)], [cosi_upperlim, cosi_upperlim], 'r--', lw=2)
ax2.get_xaxis().set_visible(False)
ax2.set_ylim([min(cosi_exp), max(cosi_exp)])

ax3 = plt.subplot(gs[1,1])#, sharey=ax1)
ax3.contour(m1_exp, cosi_exp, np.transpose(pdf2D_m1cosi_exp), levels=pdf2D_chisq, colors='r', lw=3)
ax3.pcolormesh(m1_exp, cosi_exp, np.transpose(pdf2D_m1cosi_exp), cmap='Blues', vmax=np.max(pdf2D_m1cosi_exp), vmin=0.)
ax3.set_xlim([min(m1_exp), max(m1_exp)])
ax3.set_xlabel(r'Pulsar Mass (M$_{\odot}$)', fontproperties=font, fontsize=15)
ax3.grid()
ax3.set_ylim([min(cosi_exp), max(cosi_exp)])
ax3.set_yticklabels([])

ax4 = plt.subplot(gs[0,1], sharex=ax3, sharey=ax0)
ax4.plot(m1_exp, pdf_M1_m1cosi / np.max(pdf_M1_m1cosi), 'b-', lw=3)
if (all(m1_confints)):
    m1_lowerlim, m1_median, m1_upperlim = m1_confints
    d1 = interp1d(m1_exp, pdf_M1_m1cosi / np.max(pdf_M1_m1cosi))
    ax4.plot([m1_lowerlim, m1_lowerlim], [0., d1(m1_lowerlim)], 'r--', lw=2)
    ax4.plot([m1_median, m1_median], [0., d1(m1_median)], 'r-', lw=2)
    ax4.plot([m1_upperlim, m1_upperlim], [0., d1(m1_upperlim)], 'r--', lw=2)
#ax4.plot(m1_exp, pdf_M1, 'r-', lw=3)
ax4.get_yaxis().set_visible(False)

if ('THETA' in GridDict):   
 
    ax6 = plt.subplot(gs[1,2], sharey=ax1)
    ax6.set_xlim([min(kom_exp), max(kom_exp)])
    ax6.contour(kom_exp, cosi_exp, np.transpose(pdf2D_komcosi_exp), levels=pdf2D_chisq, colors='r', lw=3)
    ax6.pcolormesh(kom_exp, cosi_exp, np.transpose(pdf2D_komcosi_exp), cmap='Blues', vmax=np.max(pdf2D_komcosi_exp), vmin=0.)
    ax6.set_xlabel(r'$\Omega$ (deg)', fontproperties=font, fontsize=15)
    ax6.set_xlim([min(kom_exp), max(kom_exp)])
    #ax6.set_yticklabels([])
    ax6.grid()

    ax5 = plt.subplot(gs[0,2], sharex=ax6, sharey=ax0)
    ax5.plot(kom_exp, pdf_KOM / np.max(pdf_KOM), 'b-', lw=3)
    ax5.get_yaxis().set_visible(False)

    if (all(kom_confints)):
        kom_lowerlim, kom_median, kom_upperlim = kom_confints
        d1 = interp1d(kom_exp, pdf_KOM / np.max(pdf_KOM)) 
        ax5.plot([kom_lowerlim, kom_lowerlim], [0., d1(kom_lowerlim)], 'r--', lw=2)
        ax5.plot([kom_median, kom_median], [0., d1(kom_median)], 'r-', lw=2)
        ax5.plot([kom_upperlim, kom_upperlim], [0., d1(kom_upperlim)], 'r--', lw=2)

if (overplotFile is not None):
    overplot_GridDict = pickle.load(open(overplotFile, "rb"))
    op_M1 = overplot_GridDict['M1']
    op_pdf_M1 = overplot_GridDict['pdf_M1']
    op_M2 = overplot_GridDict['M2']
    op_pdf_M2 = overplot_GridDict['pdf_M2']
    op_COSI = overplot_GridDict['COSI']
    op_pdf_COSI = overplot_GridDict['pdf_COSI']
    ax0.plot(op_M2, op_pdf_M2 / np.max(op_pdf_M2), 'g-', alpha=0.7, lw=3)
    ax2.plot(op_pdf_COSI / np.max(op_pdf_COSI), op_COSI, 'g-', alpha=0.7, lw=3)
    ax4.plot(op_M1, op_pdf_M1 / np.max(op_pdf_M1), 'g-', alpha=0.7, lw=3)

plt.figtext(0.8, 0.8, GridDict['PSR'], fontproperties=font, fontsize=15)
plt.savefig(GridDict['PSR']+'.grids.'+gridtype+'.png',format='png')

# store the basic info for later processing.
SDpdfs = {}
SDpdfs['PSR'] = GridDict['PSR']
SDpdfs['pdf_M1'] = pdf_M1_m1cosi
SDpdfs['M1'] = m1_exp
SDpdfs['pdf_M2'] = pdf_M2
SDpdfs['M2'] = m2_exp
SDpdfs['pdf_COSI'] = pdf_cosi
SDpdfs['pdf2D_m2cosi'] = np.transpose(pdf2D_m2cosi_exp)
SDpdfs['pdf2D_m1cosi'] = np.transpose(pdf2D_m1cosi_exp)
SDpdfs['COSI'] = cosi_exp
SDpdfs['massfunc'] = mf
SDpdfs['bestfit_SINI'] = sini_bestfit
SDpdfs['bestfit_M1'] = m1_bestfit
SDpdfs['bestfit_M2'] = m2_bestfit

# if input grids are orthometric, store orthometric arrays/PDFs.
if ('H3' in GridDict):
    SDpdfs['H3'] = h3_exp
    SDpdfs['pdf_H3'] = pdf_H3
    if ('H4' in GridDict):
        SDpdfs['H4'] = h4_exp
        SDpdfs['pdf_H4'] = pdf_H4
    elif ('STIG' in GridDict):
        SDpdfs['STIG'] = stig_exp
        SDpdfs['pdf_STIG'] = pdf_STIG

# if input grid is 3D, store the KOM array/PDF.
if ('THETA' in GridDict):
    SDpdfs['pdf_KOM'] = pdf_KOM
    SDpdfs['KOM'] = kom_exp

elif ('DIST' in GridDict):
    SDpdfs['pdf_DIST'] = pdf_px
    SDpdfs['DIST'] = dist_exp
    SDpdfs["pdf2D_distm1"] = pdf2D_pxm1_exp
    SDpdfs["pdf2D_distm2"] = pdf2D_pxm2_exp
    SDpdfs["pdf2D_distcosi"] = pdf2D_pxcosi_exp

elif ('PX' in GridDict):
    SDpdfs['pdf_PX'] = pdf_px
    SDpdfs['PX'] = px_exp
    SDpdfs["pdf2D_pxm1"] = pdf2D_pxm1_exp
    SDpdfs["pdf2D_pxm2"] = pdf2D_pxm2_exp
    SDpdfs["pdf2D_pxcosi"] = pdf2D_pxcosi_exp

# write to a .npz file.
np.savez(
    "gridPDFs." + GridDict['PSR'] + "." + gridtype + ".npz", 
    grid_PDF_data = SDpdfs
)
