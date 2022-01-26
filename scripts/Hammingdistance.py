# Compute Hamming distance of two sequences in alignement.fasta

with open ('alignment.fasta') as input_data:
	lines = input_data.readlines()
	data = [line.strip() for line in lines]
	BEOM = data[1]
	FIOER = data[3]

def HammingDistance(p, q):
	score = 0
	if len(p) == len(q):
		for i in range (len(p)):
			if (p[i]) != q [i]:
				score = score+1
	return (score)

print(HammingDistance(BEOM,FIOER)/len(FIOER))
