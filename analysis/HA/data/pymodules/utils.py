"""
Python utilities for divergence timing alignment creation
Written by Sarah Hilton.
"""
import re
from phydmslib import constants


def createDateDictionary(sequences):
    """
    This function organizes the sequences by year of collection.

    This function assumes that for a sequence collected in 2000 there is a string
            '2000_2000' somewhere in the header.

    A list of headers which cannot be parsed will be printed at the end.

    input: list of BioPython records
    output: dictionary keyed by date with a
            list of sequence records from that year.
    """
    datesDictionary = {}

    for seq in sequences: #go through each sequence
        regex = r"(/\d+_+\d+/)" #look for dates
        match = re.search(regex, seq)
        if match: # if a date is found, extract the exact date
            match = re.search(r"(\d\d\d\d)", (match.group(0)))
            if match:
                year =  int(match.group(0))
                if year not in range(1918,2018): #make sure the extracted date makes sense
                    pass
                elif year in datesDictionary.keys(): #add sequence to the dictionary
                    datesDictionary[year].extend([seq])
                else:
                    datesDictionary[year] = [seq]
    return datesDictionary

def translate_with_gaps(seq):
    """
    This function uses `phydmslib.constants` to translate a nucleotide sequence.
    input: `str`
    output: `str`
    """
    prot_seq = []
    for i in range(0, len(seq),3):
        codon = seq[i:i+3]
        codon = constants.CODON_TO_INDEX[codon]
        AA = constants.CODON_TO_AA[codon]
        AA = constants.INDEX_TO_AA[AA]
        prot_seq.append(AA)
    return "".join(prot_seq)
