from time import perf_counter
from scipy.special import lambertw, wrightomega
import numpy as np
from random import random
import pyperf


def bench_lambertw():
    lambertw(random() * 100)


def bench_wrightomega():
    wrightomega(random() * 100)


runner = pyperf.Runner()
runner.bench_func("lambertw", bench_lambertw)
runner.bench_func("wrightomega", bench_wrightomega)
