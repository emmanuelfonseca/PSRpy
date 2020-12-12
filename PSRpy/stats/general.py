import numpy as np
import sys

def compute_credible_interval(coordinate, pdf):
    """
    Computes the 68.3% and 95.4% equal-tailed credible interval, as well as the median, 
    of a supplied probability density function (PDF) of the given variable (or "coordinate").
    """

    cdf_at_interval = [0.0228, 0.1586, 0.5, 0.8414, 0.9772]
    intervals = []

    # first, ensure input PDF is normalized.
    pdf_normed = pdf / np.sum(pdf)

    # next, compute CDF.
    cdf = np.array([np.sum(pdf_normed[:idx]) for idx in range(1, len(pdf_normed))])

    for current_interval in cdf_at_interval:
        intervals.append(coordinate[(np.fabs(cdf - current_interval).argmin())])

    return intervals

def compute_pdf_from_chisq(chisq):
    """
    Computes a probability density function (PDF) from an input array or grid of 
    goodness-of-fit (chisq) values. This function uses the likelihood definition 
    given in Appendix A of Splaver et al. (2005, ApJ, 581, 509).
    """

    return 0.5 * np.exp(-0.5 * (chisq - np.min(chisq)))
