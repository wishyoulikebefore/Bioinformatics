from scipy.stats import binom
import os
from decimal import getcontext
import numpy as np

def cal_diversity(input_array,type):
    count_sum = input_array.sum()
    p = input_array/count_sum
    if type == "shannon":
        p = input_array/count_sum
        shannon = -(p*np.log(p)).sum()
        value = np.e**shannon
    elif type == "simpson":
        value = 1/(p**2).sum()
    else:
        print("Type指定有误")
        return
    getcontext().prec = 2
    return value

def binom_test(m, n, fdr):
    cum_prob = binom.pmf(range(m), n, fdr).sum()
    if cum_prob >= 1:
        pValue = binom.pmf(range(m, n + 1), n, fdr).sum()
    else:
        pValue = 1 - cum_prob
    return pValue
