# -*- coding: utf-8 -*-
"""
Gene_Finder

@author: Minju Kang

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
# test = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if(nucleotide == 'A'):
        return 'T'
    elif(nucleotide == 'T'):
        return 'A'
    elif(nucleotide == 'C'):
        return 'G'
    elif(nucleotide == 'G'):
        return 'C'
    else:
        return False
    # or use dictionary?
    # how do i replace the list's element after setting up the dictionary?


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    complement = ''
    i = len(dna)-1
    while i >= 0:
        complement = complement + get_complement(dna[i])
        i = i-1
    return complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    i = 3
    while i < len(dna):
        if dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA':
            stop_codon = i
            # aliasing
            return dna[:stop_codon]
        i = i + 3
        # this function returns first nucleotide to "stop codon"
    return dna
    # returns dna sequence up to stop codon


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # ORFs should start on multiples of 3 from the start of the string
    # the function should not return ORFs that are nested within another ORFs
    # Use while loop
    i = 0
    ORFS = []
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            rest = rest_of_ORF(dna[i:])
            ORFS.append(rest)
            i = i + len(rest)
        else:
            i = i + 3
            # ORFs start on multiples of 3
    return ORFS


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    rest_orfs = []
    i = 0
    for i in range(0, 3):
        rest_orfs = rest_orfs + find_all_ORFs_oneframe(dna[i:])
        # update the value of rest_orfs
    return rest_orfs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence
        on both strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    result = find_all_ORFs(dna)
    result1 = find_all_ORFs(get_reverse_complement(dna))
    return result + result1


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # finding longest string of "find all ORFs both strands(dna)"
    return max(find_all_ORFs_both_strands(dna))


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    maxi = 0
    for i in range(0, num_trials):
        if len(longest_ORF(shuffle_string(dna))) > max:
            maxi = len(longest_ORF(shuffle_string(dna)))
        # update the max value, length of longest ORF off shuffled string
    return maxi


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # recall string of amino acids from aa_table
    # return a string of codon
    codon = ''
    for i in range(0, len(dna)-(len(dna) % 3), 3):
        # range(start,stop,step)
        codon = codon + aa_table[dna[i:i+3]]
        i += 3
    return codon


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    dna_list = find_all_ORFs_both_strands(dna)
    # threshold = longest_ORF_noncoding(dna, 1500)
    threshold = 200
    amino_list = []
    for i in dna_list:
        amino = coding_strand_to_AA(i)
        if len(amino) > threshold:
            amino_list.append(amino)
    return amino_list
# i have to put elements of the dna_list into coding_strand_to_AA

if __name__ == "__main__":
    # gene_finder(test)

    import doctest
    doctest.testmod(verbose=True)

    # doctest.run_docstring_examples(coding_strand_to_AA, globals())
