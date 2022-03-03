from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd

def generate_plottable_points(specificity_scores):
    """
    with a list of specificity scores in a list, generates a CDF
    :param specificity_scores: list of floating point values
    :return: pair of lists, x and y, plotting which will plot the CDF
    """
    vals = np.arange(0, 1.01, 0.025)
    #specificity_scores = [x for x in specificity_scores if x <= 1.0]
    s = np.array(specificity_scores)
    ys = []
    for v in vals:
        ys.append(np.sum(np.count_nonzero(s <= v)))
    yss = [y * 1.0 / len(specificity_scores) for y in ys]
    return vals, yss

human = pd.read_csv('krispmer-human.csv')
staph = pd.read_csv('krispmer-staph.csv')
yeast = pd.read_csv('krispmer-yeast.csv')
staph = staph[staph.krispmer_score < 2.5]

specificity_human_kr = [1.0/x for x in human['krispmer_score']]
specificity_human_ge = [1.0/x for x in human['genome_score']]
data1 = [generate_plottable_points(specificity_human_kr),generate_plottable_points(specificity_human_ge)]

specificity_staph_kr = [1.0/x for x in staph['krispmer_score']]
specificity_staph_ge = [1.0/x for x in staph['genome_score']]
data2 = [generate_plottable_points(specificity_staph_kr), generate_plottable_points(specificity_staph_ge)]

yeast_genome_scores = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.6933333336, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0767134141512873, 1.0, 1.0, 1.0, 1.0, 2.71860948564395, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.285185185137037, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.765625, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.4399999998, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.476190476, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2999999999820513, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
specificity_yeast_kr = [1.0/x for x in yeast['inverse_specificity']]
specificity_yeast_ge = [1.0/x for x in yeast_genome_scores]
data3 = [generate_plottable_points(specificity_yeast_kr), generate_plottable_points(specificity_yeast_ge)]

figure, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12,4))
labels = ['kRISP-meR', 'with reference']

ls = []
for i in range(2):
    l = ax1.plot(data2[i][0], data2[i][1])
    ls.append(l)
ax1.set_ylabel('CDF')

for i in range(2):
    ax2.plot(data3[i][0], data3[i][1])
ax2.set_xlabel('Specificity scores')

for i in range(2):
    ax3.plot(data1[i][0], data1[i][1])

ax1.set_title('S. aureus', fontsize=11)
ax2.set_title('S. cerevisiae', fontsize=11)
ax3.set_title('H. sapiens (Chr14)', fontsize=11)

figure.suptitle('Comparison of CDF of specificity scores', fontsize=13)
figure.legend(ls, labels=labels)

plt.savefig('cdf.pdf')
