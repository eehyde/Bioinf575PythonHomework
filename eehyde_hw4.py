#!/usr/bin/env python
'''Module containing functions and data for the winter 2017 BIOINF 575 class.
'''


################### import any need libraries here ##########################
from __future__ import print_function, division


######################## End of import section ##############################

#This is the one letter universal translation table
#It handles cases of DNA ambiguity where the encoded amino acid is unambiguous.
#You need to deal with the missing cases where ambiguity codes would result in
#an ambiguous amino acid assignment. It is suggested that you use 'X' in these
#cases as this is the standard character for an unknown amino acid.
#Only Y (pyrimidine), R (purine) and N (any) degeneracy symbols are handled at
#this time. (need to add M,K,W,S,B,D,H,V where appropirate)
#Stop codons are symbolized as X
#Reassign TAA, TAG, TAR and TGA to change the stop codon sybmol if desired.
transTab1L = {
'TTT': 'F', 'TTC': 'F', 'TTY': 'F', 'TTA': 'L', 'TTG': 'L', 'TTR': 'L', 
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TCN': 'S', 'TCY': 'S', 'TCR': 'S', 
'TAT': 'Y', 'TAC': 'Y', 'TAY': 'Y', 'TAA': 'X', 'TAG': 'X', 'TAR': 'X', 
'TGT': 'C', 'TGC': 'C', 'TGY': 'C', 'TGA': 'X', 'TGG': 'W', 
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CTY': 'L', 'CTR': 'L', 'CTN': 'L',
    'YTG': 'L', 'YTA': 'L', 
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CCY': 'P', 'CCR': 'P', 'CCN': 'P', 
'CAT': 'H', 'CAC': 'H', 'CAY': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CAR': 'Q', 
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'CGY': 'R', 'CGR': 'R', 'CGN': 'R', 
'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATY': 'I', 'ATG': 'M', 
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ACY': 'T', 'ACR': 'T', 'ACN': 'T', 
'AAT': 'N', 'AAC': 'N', 'AAY': 'N', 'AAA': 'K', 'AAG': 'K', 'AAR': 'K', 
'AGT': 'S', 'AGC': 'S', 'AGY': 'S', 'AGA': 'R', 'AGG': 'R', 'AGR': 'R', 
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GTY': 'V', 'GTR': 'V', 'GTN': 'V', 
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GCY': 'A', 'GCR': 'A', 'GCN': 'A', 
'GAT': 'D', 'GAC': 'D', 'GAY': 'D', 'GAA': 'E', 'GAG': 'E', 'GAR': 'E', 
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', 'GGY': 'G', 'GGR': 'G', 'GGN': 'G'
}
'''transTab1L is a dictionary that maps the codons to the corresponding
single-letter amino acid code. All codons are in upper case only. Stop
codons are translated to X.
'''

########################### Assigned Functions ###############################

def reverse_seq(seq):
    '''Accepts a string representing a nucleotide or protein sequence.
Returns the sequence in reverse order as a single string. Case of the input
string is preserved.

Arguments:
    seq: str
Returns:
    str: sequence in reverse order
'''
    return seq[::-1]


def complement_DNA(seq):
    '''Expects a DNA sequence as a single string. Returns a string
representing the complementary sequence, but does not reverse the order.
That is, if the input sequence is in the 5' -> 3' order, the complement
will be in the 3' -> 5' order. See the string module for a clue on how
to do this simply.

Arguments:
    seq (str)
Returns:
    str
    '''
    complement = ""
    for x in seq:
        if x == 'A' or 'a':
            complement.append('T')
        if x == "T" or 't':
            complement.append('A')
        if x == "G" or "g":
            complement.append("C")
        if x == "C" or "c":
            complement.append("G")
        else:
            complement.append("N")
    return complement



def reverse_complement_DNA(seq):
    '''Your function should not change the case of the characters, even if
there is a mix of upper and lower case characters in the sequence. See the
string module of Python for an efficient way to do this.

    Please write your own doc string here.
'''
    reverse_complement = complement_DNA(seq)
    return reverse_complement[::-1]


def percent_gc(seq):
    '''I'm leaving it to you to define exactly what this does - it should be
able to handle a sequence with a mix of upper and lower case characters,
and the DNA sequence may have N's in it.

Argurments:
    seq: str represting DNA sequence
Returns:
    foat
'''
    num_g = 0.0
    num_c = 0.0
    for x in seq:
        if x == "G" or "g":
            num_g += 1
        if x == "C" or "c":
            num_c += 1
    total_len = float(len(seq))
    g_and_c = num_c + num_g
    percent_gc = g_and_c / total_len
    return percent_gc * 100.0


def count_kmers(seq, kmerLength, DNA = True):
    '''Given a sequence, it will construct a dictionary that contains all
of the kmers of length k in the sequence and their counts - number of times
observed in the sequence. If the sequence is DNA, it will do this for the
input sequence and its reverse complement.

Arguments:
    seq (str), the sequence with all extraneous charcters removed
    kmerLength (int), the length of the kmers for the analysis
    DNA (bool), whether the sequence is DNA or not, default True
Returns:
    dict, dictionary with kmers as keys and the counts of the kmers as values
'''
    kmers = []
    a = 0
    b = kmerLength
    while a < len(seq):
        a_kmer = seq[a:b]
        kmers.append(a_kmer)
        a += 1
        b += 1
    kmer_dict = {}
    for x in kmers:
        if x in kmer_dict.keys():
            kmer_dict[x] += 1
        else:
            kmer_dict[x] = 1
    if DNA == True:
        reverse_seq = reverse_complement_DNA(seq)
        reverse_kmers = []
        a = 0
        b = kmerLength
        while a < len(reverse_seq):
            a_rev_kmer = seq[a:b]
            reverse_kmers.append(a_rev_kmer)
            a += 1
            b += 1
        reverse_kmer_dict = {}
        for x in reverse_kmers:
            if x in reverse_kmer_dict.keys():
                reverse_kmer_dict[x] += 1
            else:
                reverse_kmer_dict[x] = 1
    return kmer_dict,reverse_kmer_dict



def translate_DNA(seq):
    '''Translates a DNA sequence. Translation begins with at the first position
of the sequence and continues as long as there are complete codons (3-mers)
available.

Arguments:
    seq (str), the DNA sequence with all extraneous charcters removed
Returns:
    str, putative amino acid sequenc encoded by the DNA
'''
    codons = [] #new empty list to add kmers to
    a = 0 # initilize an accumulator
    b = 3 # initilize the ending position of a codon
    while a < len(seq):
        a_codon = seq[a:b]
        codons.append(a_codon)
        # advance to the next three letters:
        a += 3
        b += 3
    translation = ''
    for x in codons:
        if len(x) == 3:
            if x in transTab1L.keys():
                translation.append(transTab1L[x])
    return translation



######################## End of Assigned Functions ###########################

############################# FORMATTING #######################################

def fasta_format (header, seq, linelength = 60):
    """fasta_format will convert a sequence to a Fasta formated sequence.
The name (definition line) doesn't need to include the > symbol or
return characters. The function will add them. The header and seq (sequence)
should be strings.

The linelength controls the length of the seqence lines only.
The sequence should be a raw sequence, ie, contain no return characters.

Arguments,:
    header (str)
    seq (str)
    linelength (int), default 60
Returns:
    str, the formatted header and sequence as a single string
"""
    
    tempLines = []
    header = header.strip()
    if not header.startswith('>'):
        header = '>' + header
    tempLines.append(header)
    for x in xrange(0, len(seq), linelength):
        tempLines.append(seq[x:x + linelength])
    #add an empty string so that a \n will be joined at the end
    tempLines.append('')
    return '\n'.join(tempLines)


############################# END FORMATTING ###################################

######################## File Operations ###############################

def get_next_fasta (fileObject):
    '''usage: for header, seq in get_next_fasta(fileObject):
    
This is a generator that returns one fasta record's header and
sequence at a time from a multiple fasta file. Return character is removed
from the header. The sequence is returned as one continuous string
with no returns. The returned value is a tuple (header, sequence)
If their is no sequence associated with a header, seq will be an
empty string
Code simplification contributed by Dattatreya Mellacheruvu
01/16/2009, Jeffrey R. de Wet
08/02/2010 refactored to put lines into a temporary list

Arguments:
    fileObject: an open file object to a fasta file
Returns:
    str: the header
    str: the sequence
'''
    
    header = ''
    seq = ''

    lineList = []
    #The following for loop gets the header of the first fasta
    #record. Skips any leading junk in the file
    for line in fileObject:
        if line.startswith('>'):
            header = line.strip()
            break
    
    for line in fileObject:
        if line.startswith('>'):
            seq = ''.join(lineList)
            lineList = []
            yield header, seq
            header = line.strip()
            seq = ''
        else:
            #seq += line.strip()
            lineList.append(line.strip())
    #yield the last entry
    if header:
        seq = ''.join(lineList)
        yield header, seq


####################### End of File Functions #########################
