import scipy.special as ss
import sys

def compute_test_f(chisq_1, chisq_2, dof_1, dof_2):
    """
    Computes the p-value of an F-test between two models.
    """

    # derive some stats.
    delta_chisq = chisq_1 - chisq_2
    delta_dof = dof_1 - dof_2
    chisq_red_2 = chisq_2 / dof_2

    # now compute F statistic.
    stat_F = delta_chisq / delta_dof / chisq_red_2

    # finally, find p-value.
    prob = 1.

    if delta_chisq > 0.:
        prob -= ss.fdtr(delta_dof, dof_2, stat_F)

    return prob
