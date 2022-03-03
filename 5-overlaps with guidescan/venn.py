from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from collections import Counter
import pandas as pd

figure, (ax2, ax1) = plt.subplots(2,1, constrained_layout=True)

gs_hum = 109
kr_hum = 92
hum_com = 90

gs_staph = 236
kr_staph = 229
staph_com = 225

gs_yeast = 187
kr_yeast = 185
yeast_com = 184

human_count_kr = 0
human_count_gs = 3

staph_count_kr = 1
staph_count_gs = 0

sets = Counter()
sets['01'] = 2
sets['11'] = 90
sets['10']  = 19

plt.subplot(1,3,3)
setLabels = ['', '']
v = venn2(subsets = sets, set_labels = setLabels)
h, l = [],[]
for i in sets:
    v.get_label_by_id(i).set_text('')
    h.append(v.get_patch_by_id(i))
    l.append(sets[i])

v.get_patch_by_id('10').set_alpha(0.35)
v.get_patch_by_id('10').set_color('red')

v.get_patch_by_id('01').set_alpha(0.35)
v.get_patch_by_id('01').set_color('blue')

v.get_patch_by_id('11').set_alpha(0.65)
v.get_patch_by_id('11').set_color('purple')

rot_angle = 90
plt.text(0.6,0,'kRISP-meR', rotation=rot_angle)
plt.text(-0.7,0,'GuideScan', rotation=rot_angle)

plt.legend(handles=h, labels=l, title="", loc=1)
plt.title('H. sapiens')

plt.subplot(1,3,1)
sets['01'] = kr_staph-staph_com
sets['11'] = staph_com
sets['10']  = gs_staph-staph_com
setLabels = ['', '']
v = venn2(subsets = sets, set_labels = setLabels)
h, l = [],[]
for i in sets:
    v.get_label_by_id(i).set_text('')
    h.append(v.get_patch_by_id(i))
    l.append(sets[i])

v.get_patch_by_id('10').set_alpha(0.35)
v.get_patch_by_id('10').set_color('red')

v.get_patch_by_id('01').set_alpha(0.35)
v.get_patch_by_id('01').set_color('blue')

v.get_patch_by_id('11').set_alpha(0.65)
v.get_patch_by_id('11').set_color('purple')

plt.text(0.6,0,'kRISP-meR', rotation=rot_angle)
plt.text(-0.7,0,'GuideScan', rotation=rot_angle)

plt.legend(handles=h, labels=l, title="", loc=1)
plt.title('S. aureus')

plt.subplot(1,3,2)
sets['01'] = kr_yeast-yeast_com
sets['11'] = yeast_com
sets['10']  = gs_yeast-yeast_com
setLabels = ['', '']
v = venn2(subsets = sets, set_labels = setLabels)
h, l = [],[]
for i in sets:
    v.get_label_by_id(i).set_text('')
    h.append(v.get_patch_by_id(i))
    l.append(sets[i])

v.get_patch_by_id('10').set_alpha(0.35)
v.get_patch_by_id('10').set_color('red')

v.get_patch_by_id('01').set_alpha(0.35)
v.get_patch_by_id('01').set_color('blue')

v.get_patch_by_id('11').set_alpha(0.65)
v.get_patch_by_id('11').set_color('purple')

plt.text(0.6,0,'kRISP-meR', rotation=rot_angle)
plt.text(-0.7,0,'GuideScan', rotation=rot_angle)

plt.legend(handles=h, labels=l, title="", loc=1)
plt.title('S. cerevisiae')

#figure.suptitle('Overlap of the guides identified by\nkRISPmeR and Guidescan')

plt.savefig('overlaps.pdf')
