from scipy.stats import binom
import os
from decimal import getcontext

def cal_diversity(input_array,type):
    count_sum = input_array.sum()
    p = input_array/count_sum
    if type == "shannon":
        p = input_array/count_sum
        shannon = -(p*np.log(p)).sum()
        value = np.e**shannon
    elif type = "simpson":
        value = 1/(p**2).sum()
    getcontext().prec = 2
    return value    
