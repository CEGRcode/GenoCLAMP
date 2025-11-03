import numpy as np
from typing import Union
from scipy.stats import pearsonr

def highest_n_info_sum(pfm, n=4, w=3):
    '''
    Sum of the information content of the highest n windows of width w in a PFM
    '''
    pwm = pfm / np.sum(pfm, axis=1, keepdims=True)
    logpwm = np.array(pwm)
    logpwm[pwm > 0.] = np.log2(pwm[pwm > 0.])
    bits = np.sum(pwm * logpwm, axis=1) + 2.
    val = sum(bits[:w]) / w
    mean_bits = [val]
    for i in range(w, len(bits)):
        val += (bits[i] - bits[i - w]) / w
        mean_bits.append(val)
    return sum(sorted(mean_bits, reverse=True)[:n])

def check_periodicity(pfm, p=1):
    '''
    Periodicity score for a PFM based on the weighted Pearson correlation between
    positions separated by p bases, weighted by the information content of the bases
    '''
    pwm = pfm / np.sum(pfm, axis=1, keepdims=True)
    logpwm = np.array(pwm)
    logpwm[pwm > 0.] = np.log2(pwm[pwm > 0.])
    bits = np.sum(pwm * logpwm, axis=1) + 2.
    w = pfm.shape[0]
    corr_sum = 0
    total_bit_prod = 0
    for offset in range(p, w, p):
        for i in range(w - offset):
            corr = pearsonr(pfm[i, :], pfm[i + offset, :])[0]
            bit_prod = bits[i] * bits[i + offset]
            if not np.isnan(corr):
                corr_sum += corr * bit_prod
            total_bit_prod += bit_prod
    return corr_sum / total_bit_prod

def filter_motifs(items, nsites_thresh=10, evalue_thresh=.01, info_score_thresh=5., periodicity1_thresh=.6,
                  periodicity2_thresh=.75, periodicity3_thresh=.75):
    '''
    Filters a list of GreedyItems based on various QC thresholds
    '''
    new_items = []
    for item in items:
        nsites = item.source[1]
        evalue = item.source[2]
        pfm = item.pfm
        info_score = highest_n_info_sum(pfm)
        periodicity1 = check_periodicity(pfm, p=1)
        periodicity2 = check_periodicity(pfm, p=2)
        periodicity3 = check_periodicity(pfm, p=3)
        if nsites >= nsites_thresh and evalue <= evalue_thresh and info_score >= info_score_thresh \
           and periodicity1 <= periodicity1_thresh and periodicity2 <= periodicity2_thresh \
           and periodicity3 <= periodicity3_thresh:
            new_items.append(item)
    return new_items

def boltzmann(arr: np.ndarray, alpha: float, axis: Union[None, int, tuple] = None):
    '''
    Boltzmann operator as described by https://en.wikipedia.org/wiki/Smooth_maximum#Boltzmann_operator
    '''
    return np.sum(arr * np.exp(alpha * arr), axis=axis) / np.sum(np.exp(alpha * arr), axis=axis)
    
def trim_motif(aligned_pfms, info_thresh=.5, w=2):
    '''
    Trims bases from the start and end of a PFM stack with information content below info_thresh
    '''
    pfm = np.sum(aligned_pfms, axis=0)

    # If there is only one motif, return the original PFM
    if aligned_pfms.shape[0] == 1:
        return pfm, 0, 0, False

    # Calculate the posterior probability of each base (1 pseudocount)
    pwm = pfm / np.sum(pfm, axis=1, keepdims=True)
    pwm[pwm == 0.] = 1.
    # Calculate the information content of each base
    bits = np.sum(pwm * np.log2(pwm), axis=1) + np.log2(aligned_pfms.shape[2])
    val = sum(bits[:w]) / w
    mean_bits = np.zeros(aligned_pfms.shape[1] - w + 1, dtype=np.float64)
    mean_bits[0] = val
    for i in range(w, len(bits)):
        val += (bits[i] - bits[i - w]) / w
        mean_bits[i - w + 1] = val
    # Trim the PFM to only include bases with information content above the threshold
    informative_bits = np.flatnonzero(mean_bits > info_thresh)
    if len(informative_bits) <= 1:
        return pfm, 0, 0, False
    start = informative_bits[0]
    end = informative_bits[-1] + w
    return pfm[start:end, :], start, pfm.shape[0] - end, True
