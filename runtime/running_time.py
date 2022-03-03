import pandas as pd

genomes = ['worm', 'yeast', 'fruitfly']
times = ['jellyfish_time', 'alignment_time', 'preprocess_time']

df = pd.read_csv('rt.csv', delimiter=',')

for genome in genomes:
    print (genome)
    for time in times:
        print (time + ': ' + str(df[ df.genome.str.contains(genome) ].mean(axis=0)[time]))
