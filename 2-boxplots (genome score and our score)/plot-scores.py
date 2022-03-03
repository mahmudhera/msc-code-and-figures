from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd

human = pd.read_csv('krispmer-human.csv')
staph = pd.read_csv('krispmer-staph.csv')
yeast = pd.read_csv('krispmer-yeast.csv')
staph = staph[staph.krispmer_score < 2.5]

data1 = [human['krispmer_score'], human['genome_score']]
data2 = [staph['krispmer_score'], staph['genome_score']]

yeast_genome_scores = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.6933333336, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0767134141512873, 1.0, 1.0, 1.0, 1.0, 2.71860948564395, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.285185185137037, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.765625, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.4399999998, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.476190476, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.2999999999820513, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
data3 = [yeast['inverse_specificity'], yeast_genome_scores]

figure, (ax2, ax3, ax1) = plt.subplots(1,3,constrained_layout=True)

for i in range(1, 3):
    y = data1[i-1]
    x = np.random.normal(i, 0.06, len(y))
    ax1.plot(x, y, 'x', alpha=0.2)
ax1.boxplot(data1, showfliers=False, widths=0.5) # Or you can use the boxplot from Pandas

for i in range(1, 3):
    y = data2[i-1]
    x = np.random.normal(i, 0.06, len(y))
    ax2.plot(x, y, '+', alpha=0.3)
ax2.boxplot(data2, showfliers=False, widths=0.5) # Or you can use the boxplot from Pandas

for i in range(1, 3):
    y = data3[i-1]
    x = np.random.normal(i, 0.06, len(y))
    ax3.plot(x, y, '*', alpha=0.3)
ax3.boxplot(data2, showfliers=False, widths=0.5) # Or you can use the boxplot from Pandas

ax1.set_xticklabels(['kRISP-meR', 'Genome'], rotation=40)
ax2.set_xticklabels(['kRISP-meR', 'Genome'], rotation=40)
ax3.set_xticklabels(['kRISP-meR', 'Genome'], rotation=40)

ax1.set_title('H. sapiens (Chr14)', fontsize=11)
ax2.set_title('S. aureus', fontsize=11)
ax3.set_title('S. cerevisiae', fontsize=11)

figure.suptitle('Comparison of kRISP-meR scores\nwith scores calculated using reference genome', fontsize=13)

plt.savefig('score-comparison-3.pdf')
