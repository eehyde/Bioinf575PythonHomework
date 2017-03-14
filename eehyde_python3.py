import sys
# The function should take the arguments "seq" and "kmer length" and return a list of kmers. 
# The k-mer window should be advanced by one position such that you obtain overlapping K-mers. 
# Please note to take the "seq" and "kmer length" as user input in the command line. 
def get_kmer(seq,kmer_length):
	kmers = [] #new empty list to add kmers to
	a = 0 #initialize an accumulator
	b = kmer_length
	while a < len(seq):
		a_kmer = seq[a:b]
		kmers.append(a_kmer)
		#advance the kmer position:
		a += 1
		b += 1
	return kmers

# The function should take the arguments "seq" and return a list of codons within the user input sequence. 
# Please note that codons are non-overlapping. 
def get_codon(seq):
	codons = [] #new empty list to add kmers to
	a = 0 # initilize an accumulator
	b = 3 # initilize the ending position of a codon
	while a < len(seq):
		a_codon = seq[a:b]
		codons.append(a_codon)
		# advance to the next three letters:
		a += 3
		b += 3
	return codons


seq = str(sys.argv[1]) #making sure the type is a string
kmer_length = int(sys.argv[2]) #making sure the type is an integer
print("The kmers for the given sequence are:",get_kmer(seq,kmer_length))
print("The codons for the given sequence are:",get_codon(seq))