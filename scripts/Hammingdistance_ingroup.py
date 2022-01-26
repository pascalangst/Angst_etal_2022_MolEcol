# Compute Hamming distance between all sequences in finalconcatenation.fasta'

import numpy as np

with open ('finalconcatenation.fasta') as input_data:
        lines = input_data.readlines()
        data = [line.strip() for line in lines]

def HammingDistance(p, q):
	score = 0
	if len(p) == len(q):
		for i in range (len(p)):
			if (p[i]) != q [i]:
				score = score+1
	return (score)

for i in np.arange(1, 36, 2):
   for j in np.arange(1, 36, 2):
      print(1-(HammingDistance(data[i], data[j])/len(data[i])))


