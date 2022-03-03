from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from collections import Counter
import pandas as pd

#figure, (ax2, ax1) = plt.subplots(2,1, constrained_layout=True)

gs_hum = 109
kr_hum = 92
hum_com = 90

gs_staph = 236
kr_staph = 229
staph_com = 225

human_count_kr = 1
human_count_gs = 3

staph_count_kr = 1
staph_count_gs = 0

genome=['H. sapiens','S. aureus']
tools=['GuideScan','kRISP-meR']
pos = np.arange(2)
bar_width = 0.35
count_gs=[human_count_gs, staph_count_gs]
count_kr=[human_count_kr, staph_count_kr]

plt.ylim(0,4)
 
plt.bar(pos,count_gs,bar_width,color='red',edgecolor='black',alpha=0.4)
plt.bar(pos+bar_width,count_kr,bar_width,color='blue',edgecolor='black', alpha=0.4)

for i in range(2):
	plt.text(pos[i], count_gs[i]+0.1, str(count_gs[i]))
	plt.text(pos[i]+bar_width, count_kr[i]+0.1, str(count_kr[i]))

plt.xticks(pos+bar_width/2, genome)
plt.xlabel('Dataset', fontsize=13)
plt.ylabel('gRNAs yielding perfect off-targets', fontsize=13)
plt.legend(tools,loc=1)
plt.savefig('off-targets.pdf')
