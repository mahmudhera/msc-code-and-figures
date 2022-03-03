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
    vals = np.arange(0, 200.0, 0.01)
    # specificity_scores = [x for x in specificity_scores if x <= 1.0]
    s = np.array(specificity_scores)
    ys = []
    for v in vals:
        ys.append(np.sum(np.count_nonzero(s <= v)))
    yss = [y * 1.0 / len(specificity_scores) for y in ys]
    return vals, yss


fig = plt.figure(figsize=(10,10))
fig.suptitle('upto 3 mismatches were used to determine the off-targets $in$ $silico$\n(Mutation matrix derived by $Doench$ $et.$ $al.$ was used)')

kr_scores = np.array(pd.read_csv('krispmer.csv').krispmer_score.tolist())
kr_predicted_ots = kr_scores - 1.0
actual_ot = np.array(pd.read_csv('krispmer.csv').genome_score.tolist()) - 1.0

all = zip(kr_scores.tolist(), kr_predicted_ots.tolist(), actual_ot.tolist())
labels = ['Predicted number of off-targets', 'Actual off-targets (in-silico)']

xs = []
ys1 = []
ys2 = []
ys3 = []
count = 0.0
for v in kr_scores:
    count += 1.0
    xs.append(1.0 / v)
    ys3.append(0)
    total_ot = 0.0
    total_predicted_ot = 0.0
    for a in all:
        if a[0] >= v:
            total_ot += a[2]
            total_predicted_ot += a[1]
    ys2.append(total_ot)
    ys1.append(total_predicted_ot)

ax = plt.subplot(221)
ax.set_xlabel('specificity')
ax.set_ylabel('number of off-targets')
ax.set_title('gRNAs with lower specificity')
ax.plot(xs, ys1)
ax.plot(xs, ys2)

xs = []
ys1 = []
ys2 = []
ys3 = []
count = 0.0
for v in kr_scores:
    count += 1.0
    xs.append(1.0 / v)
    ys3.append(0)
    total_ot = 0.0
    total_predicted_ot = 0.0
    for a in all:
        if a[0] <= v:
            total_ot += a[2]
            total_predicted_ot += a[1]
    ys2.append(total_ot)
    ys1.append(total_predicted_ot)

ax = plt.subplot(222)
ax.set_ylim(0,400)
ax.set_xlabel('specificity')
ax.set_ylabel('number of off-targets')
ax.set_title('gRNAs with higher specificity')
ax.plot(xs, ys1, label=labels[0])
ax.plot(xs, ys2, label=labels[1])

low, high = 0.3, 0.6

x1 = []
y1 = []
y2 = []
for (x, y) in zip(xs, ys1):
    if low <= x <= high:
        x1.append(x)
        y1.append(y)
for (x, y) in zip(xs, ys2):
    if low <= x <= high:
        y2.append(y)

ax = plt.subplot(212)
l1, = ax.plot(x1, y1)
l2, = ax.plot(x1, y2)

ax.set_xlabel('specificity')
ax.set_ylabel('number of off-targets')
ax.set_title('gRNAs with higher specificity (zoomed-in)')
fig.legend([l1, l2], labels, loc='lower center')
plt.savefig('cdf-with-gs.pdf')
