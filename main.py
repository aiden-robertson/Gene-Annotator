# Stuff to import
import math
import os
import sqlite3
from asyncore import file_dispatcher

import matplotlib.pyplot as plt
import numpy as np

# Genome information & dataset location
genome = "saccharomyces_cerevisiae-genome"
bacterial = False

path = os.path.join("Genomes", f"{genome}.db")

codon_table = { # Copied from chatgpt
    "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
    "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
    "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
    "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
    "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
    "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "TAT":"Tyr", "TAC":"Tyr", "TAA":"*", "TAG":"*",
    "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "TGT":"Cys", "TGC":"Cys", "TGA":"*", "TGG":"Trp",
    "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"
}

connection = sqlite3.connect(path)
cursor = connection.cursor()

# Gets table in dataset
cursor.execute("SELECT chromosome_number FROM genome")
tables = cursor.fetchall()

def checks(bacterial, potential_genes, bases):
    # Reverses gene list
    potential_genes.sort(key=lambda x: x[1] - x[0], reverse=True)

    filtered = []

    # Looks through all the start and stops
    for start, stop in potential_genes:
        contained = False # Marks it as not contained by default
        for fstart, fstop in filtered:

            # Flags it if it is nested inside another gene
            if start >= fstart and stop <= fstop:
                contained = True
                break

        # Filters out any not nested genes
        if not contained:
            filtered.append((start, stop))

    potential_genes = filtered

    # Confidence ranking to give each gene a score
    confidence = [0 for i in range(len(potential_genes))]

    # Check for shine-dalgarno sequence

    # Partials
    partials = {"GAGG", "GGAGGT", "AGGA", "GGAG", "GAGGT", "TAAGG"}

    if bacterial:
        for i, (start, stop) in enumerate(potential_genes):
            upstream = ''.join(bases[max(0, start - 10):start])

            if "AGGAGG" in upstream:
                confidence[i] += .3
            elif True in [partial in upstream for partial in partials]:
                confidence[i] += .3

    # Checks for gc content
    for i, (start, stop) in enumerate(potential_genes):
        seq = bases[start:stop]
        gc = (seq.count("G") + seq.count("C")) / len(seq)

        # Genes generally have more GC contents
        if .35 < gc < .45:
            confidence[i] += .3

    # Codon usage bias checks

    # Gets total codon usage in organism
    codon_freq = {}

    for start, stop in potential_genes:
        if stop - start > 1000:  # Only check longer genes(just to be safer)
            for j in range(start, stop, 3):
                codon = bases[j:j + 3]
                codon_freq[codon] = codon_freq.get(codon, 0) + 1

    # Normalizes the values
    total = sum(codon_freq.values())
    codon_freq = {k: v / total for k, v in codon_freq.items()}

    # Increases confidence on genes with more preffered codons and vise versa
    for i, (start, stop) in enumerate(potential_genes):
        score = 0
        codon_count = 0

        for j in range(start, stop, 3):
            codon = bases[j:j + 3]

            if codon in codon_freq:
                score += codon_freq[codon]
                codon_count += 1

        if codon_count > 0:
            confidence[i] += score / codon_count

    filtered = []
    filtered_confidence = []

    # Filters out anything with too low of a score
    for i, val in enumerate(confidence):
        if val < .1:
            pass
        else:
            filtered.append(potential_genes[i])
            filtered_confidence.append(confidence[i])

    potential_genes = filtered
    confidence = filtered_confidence

    return potential_genes, confidence

# Primary loop for accessing data
while True:
    # Chromosome choice
    access = input(f"Enter chromosome index to access(1 - {len(tables)}): ") if len(tables) > 1 else 1

    # Checks for valid choice
    passed = True
    try:
        access = int(access)
    except:
        passed = False

    if (access > len(tables) or access < 1) and passed:
        print("Invalid index.")

    # Checks passed; gives info on chromosome
    else:
        # Gets chromosome base sequence
        chromosome = access
        bases = cursor.execute("SELECT base_sequence FROM genome WHERE chromosome_number = ?", (chromosome,)).fetchall()[0][0]
        print(f"chromosome: {access}")
        print(f"bases: {len(bases)}")

        # Creates number line
        length = len(bases)

        plt.figure()
        plt.hlines(y=0, xmin=0, xmax=length, color="black")
        plt.xlim(0, length)

        plt.yticks([])
        plt.xlabel("Base Position")
        plt.title(f"{genome} chromosome {access}")

        plt.scatter([0, 0, len(bases)], [-1, 1, 1], color="black", s=0)

        # Searches for potential genes
        stops = {"TAA", "TAG", "TGA"} # Which codons correspond to stop codons

        potential_genes = []

        # Seperates into three reading frames
        for frame in range(3):
            # Loops through each position in each reading frame
            for pos in range(frame, length - 2, 3):
                codon = ''.join(bases[pos:pos+3])

                # Checks for start codon
                if codon == "ATG":
                    # Looks for stop codon
                    for stop in range(pos+3, length-2, 3):
                        stop_codon = ''.join(bases[stop:stop+3])

                        if stop_codon in stops:
                            potential_genes.append((pos, stop+3))
                            break

        potential_genes, confidence = checks(bacterial, potential_genes, bases)

        # Plots gene map
        i = 1
        for i, (start, stop) in enumerate(potential_genes):
            plt.scatter(start, 0, c="red", s=0)
            plt.scatter(stop, 0, c="blue", s=0)

            # Labels genes
            if start % 3 == 0:
                plt.hlines(y=-.1, xmin=start, xmax=stop, color="red", alpha=min(confidence[i] *5, 1))
            elif start % 3 == 1:
                plt.hlines(y=-.2, xmin=start, xmax=stop, color="blue", alpha=min(confidence[i] *5, 1))
            elif start % 3 == 2:
                plt.hlines(y=-.3, xmin=start, xmax=stop, color="brown", alpha=min(confidence[i] *5, 1))
            #plt.annotate(text=str(i),
                         #xy=((start+stop)/2, 0),
                         #xytext=((start+stop)/2, (.1 if i % 2 == 1 else -.1)),
                         #ha="center")

            print(f"{i} of {len(potential_genes)} genes mapped." + ("  Confidence = " + str(confidence[i]) if confidence[i] > 0 else ""))

            # Loops through all codons to assign amino acids
            amino_acids = []

            # for j in range(int(start/3), int(stop/3)):
                # amino_acids.append(codon_table[codons[int(j)]])


            # print(f"total bps in gene {i}: {len(amino_acids*3)}")
            # print(f"amino acids in gene {i}: {amino_acids}\n")

            # Updates current letter
            i += 1

        # Plots TATA boxes
        #for i in range(len(bases)-3):
        #    if bases[i:i+6] == "TATAAA":
        #        plt.scatter(i, 0, color="green", s=10)

        plt.show()

        print(confidence)

        if input("Press ENTER to exit, or any key then enter to view other chromosome\n") == '':
            break