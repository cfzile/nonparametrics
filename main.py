import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import *
from scipy.stats import norm
from tqdm import tqdm

from galton import GaltonTest


def galton(alternative="two-sided"):
    def test(sample1, sample2):
        return GaltonTest(sample1, sample2, alternative=alternative).p_value()
    return test

def wilcoxon(alternative="two-sided"):
    def test(sample1, sample2):
        return mannwhitneyu(sample1, sample2, alternative=alternative)[1]
    return test


def modeling(test, get_samples, sample_sizes, alpha_fun, num_iterations=100000):
    df = []

    for n, m in tqdm(sample_sizes):
        alpha = np.array(alpha_fun(n, m))
        rejected_H_0, rejected_H_1 = np.array([0] * len(alpha)), np.array([0] * len(alpha))
        for i in range(num_iterations):
            sample1, sample2 = get_samples(null_hypothesis=True, n=n, m=m)
            p_value_H_0 = test(sample1, sample2)
            sample1, sample2 = get_samples(null_hypothesis=False, n=n, m=m)
            p_value_H_1 = test(sample1, sample2)
            rejected_H_0 = rejected_H_0 + np.int_(p_value_H_0 <= alpha)
            rejected_H_1 = rejected_H_1 + np.int_(p_value_H_1 <= alpha)
        for i in range(len(alpha)):
            df.append([n + m, n, rejected_H_0[i] / num_iterations, alpha[i], rejected_H_1[i] / num_iterations])

    result = pd.DataFrame(df, columns=['N', 'n', 'type I error (generated)', 'type I error (real)', 'power'])

    result['error rate of type I error %'] = abs(result['type I error (generated)']
                                                 / result['type I error (real)'] - 1) * 100

    return result

def fun_samples(p):
    def get_samples_two_sided(null_hypothesis, n, m):
        smaple1, sample2 = [], []
        if null_hypothesis:
            sample1 = norm.rvs(loc=0, size=n)
            sample2 = norm.rvs(loc=0, size=m)
        else:
            sample1 = norm.rvs(loc=0, size=n)
            sample2 = norm.rvs(loc=p, size=m)
        return sample1, sample2
    return get_samples_two_sided

def f(n, m):
    return [0.025] # alpha

x = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8]
wy, gy = [], []
for n in [19, 39, 59, 79]:
    y = []
    ay = []
    for p in x:
        res = modeling(wilcoxon('less'), fun_samples(p), [(n, n)], f, num_iterations=10000)
        y.append(res['power'][0])
        ay.append(res['error rate of type I error %'][0])
    y2 = []
    ay2 = []
    for p in x:
        res = modeling(galton('less'), fun_samples(p), [(n, n)], f, num_iterations=10000)
        y2.append(res['power'][0])
        ay2.append(res['error rate of type I error %'][0])
    plt.plot(x, y2, label='galton')
    plt.plot(x, y, label='wilcoxon')
    plt.legend()
    plt.show()
    print(ay, ay2)
    wy = y
    gy = y2
    print(pd.DataFrame([wy, gy], columns=x, index=['wilcoxon', 'galton']).to_latex())