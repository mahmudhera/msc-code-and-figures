from matplotlib import pyplot as plt
import numpy as np
from collections import Counter
import pandas as pd

#figure, (ax2, ax1) = plt.subplots(2,1, constrained_layout=True)

actual_count_kr = 228
actual_count_gs = 241

ot_grna_count_kr = 3
ot_grna_count_gs = 21

genome=['Total detected gRNAs', 'Off-target producing gRNAs']
tools=['GuideScan','kRISP-meR']
pos = np.arange(2)
bar_width = 0.30

count_gs=[actual_count_gs, ot_grna_count_gs]
count_kr=[actual_count_kr, ot_grna_count_kr]

plt.ylim(0,270)
 
plt.bar(pos,count_gs,bar_width,color='red',edgecolor='black',alpha=0.4)
plt.bar(pos+bar_width,count_kr,bar_width,color='blue',edgecolor='black', alpha=0.4)

for i in range(2):
	plt.text(pos[i]-.03, count_gs[i]+2, str(count_gs[i]))
	plt.text(pos[i]+bar_width-.03, count_kr[i]+2, str(count_kr[i]))

plt.xticks(pos+bar_width/2, genome)
#plt.xlabel('Dataset', fontsize=13)
plt.ylabel('Number of gRNAs', fontsize=13)
plt.legend(tools,loc=1)
plt.savefig('variant-aware.pdf')
