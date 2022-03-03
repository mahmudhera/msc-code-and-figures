from matplotlib import pyplot as plt
import numpy as np

plt.figure(figsize=(10,6))

vals_rhodo = [(1, 100.0), (2, 100.0), (3, 99.96000000000001), (4, 99.94), (5, 99.88999999999999), (10, 99.10999999999999), (15, 95.62)]
vals_staph = [(1, 99.89999999999999), (2, 99.87999999999998), (3, 99.78999999999999), (4, 99.66000000000001), (5, 99.55), (10, 97.75999999999999), (15, 93.81)]

bar_width = 0.35
count_a = [x[1] for x in vals_rhodo]
pos = np.arange(7)
a = plt.bar(pos,count_a,bar_width,color='red',edgecolor='black',alpha=0.4)

count_b = [x[1] for x in vals_staph]
b = plt.bar(pos+bar_width,count_b,bar_width,color='blue',edgecolor='black', alpha=0.4)

legend_names = ['Rhodobacter sphaeroides', 'Staphylococcus aureus']

x_offset = 0.15
y_offset = 0.1
pos = pos - x_offset
#for i in range(7):
#    plt.text(pos[i], count_a[i]+y_offset, str(round(count_a[i],1)))
#    plt.text(pos[i]+bar_width, count_b[i]+y_offset, str(round(count_b[i],1)))

names = ['1% SV', '2% SV', '3% SV', '4% SV', '5% SV', '10% SV', '15% SV']
plt.xticks(pos+bar_width/2.0, names)
#plt.xlabel('Genome', fontsize=13)
plt.ylabel('Needleman-Wunch Similarity Score', fontsize=12)
plt.title('Similarity Score of reconstructed target sequence with\nthe original for various levels of introduced structural variations')
plt.legend([a, b], legend_names, loc='lower right')
plt.savefig('personalzed.pdf')
