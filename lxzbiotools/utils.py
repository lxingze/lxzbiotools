import os
import re
import sys
import numpy as np
from io import StringIO
from pathlib import Path
from pandas.io.formats.format import GenericArrayFormatter
from Bio import SeqIO
import sys
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gzip
from mimetypes import guess_type
from functools import partial

codon_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'CGT':'R', 'CGC':'R',   
    'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', 'TCT':'S', 'TCC':'S',
    'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ATT':'I', 'ATC':'I',
    'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L',
    'CTG':'L', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'GTT':'V',
    'GTC':'V', 'GTA':'V', 'GTG':'V', 'ACT':'T', 'ACC':'T', 'ACA':'T',
    'ACG':'T', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'AAT':'N',
    'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'C', 'CAA':'Q',
    'CAG':'Q', 'GAA':'E', 'GAG':'E', 'CAT':'H', 'CAC':'H', 'AAA':'K',
    'AAG':'K', 'TTT':'F', 'TTC':'F', 'TAT':'Y', 'TAC':'Y', 'ATG':'M',
    'TGG':'W',
    'TAG':'STOP', 'TGA':'STOP', 'TAA':'STOP'
    }

def fasta_dict(file: Path):
    """
    Reads a file of FASTA entries and returns a dict where each key is a sequence ID.
    The value is another dict with two keys 'h' for header and 's' for sequence.  The
    header is all the other text after the id in the original FASTA header.  The
    sequence has all whitespace removed.  Obviously this should only be used on files
    where memory to load them isn't an issue.
    """
    seqs = dict()
    current_seq = ''
    current_id = None
    current_header = None
    
    for line in open(file):
        line = line.rstrip()
        m = re.search('>(\S+)\s*(.*)', line)
        if m:
            ## new residue line matched, purge the existing one, if not the first
            if current_id is not None:
                ## warn if it has already been found
                if current_id in seqs:
                    sys.stderr.write(f"WARN: Duplicate ID ({current_id}) found.  Only last one kept.\n")

                ## remove all whitespace and save
                current_seq = ''.join(current_seq.split())
                seqs[current_id] = {'h':current_header, 's':current_seq}
                    
            current_seq = ''
            current_id = m.group(1)
            current_header = m.group(2)
        else:
            current_seq += str(line)

    ## don't forget the last one
    current_seq = ''.join(current_seq.split())
    seqs[current_id] = {'h':current_header, 's':current_seq}

    return seqs

def fasta_sizes(file):
    """
    Reads a file of FASTA entries and returns a dict where each key is a sequence ID and
    the value is just the size of each sequence.  This is as almost as computationally
    intensive as fasta_dict_from_file() but takes less memory and is appropriate when you
    only care about the residue lengths.
    """
    seqs = dict()
    current_id = None
    
    for line in open(file):
        line = line.rstrip()
        m = re.search('>(\S+)\s*(.*)', line)
        if m:
            ## new residue line matched, set the new seq ID
            current_id = m.group(1)
            seqs[current_id] = 0
        else:
            seqs[current_id] += len(line)

    return seqs


def splitFileContents(f, delimiter, BLOCKSIZE=8192):
    """
    Same semantics as f.read().split(delimiter), but with memory usage
    determined by largest chunk rather than entire file size
    """
    remainder = StringIO()
    while True:
        block = f.read(BLOCKSIZE)
        if not isinstance(block, str):
            block = block.decode("utf-8")
        if not block:
            break
        parts = block.split(delimiter)
        remainder.write(parts[0])
        for part in parts[1:]:
            yield remainder.getvalue()
            remainder = StringIO()
            remainder.write(part)
    yield remainder.getvalue()


def get_genome_stats(infile, genome_size):
  """Calculate the statistics for an input fasta formatted genome file.
        infile should be a fasta formatted file
        genome_size should be in MBp
        Returns a dictionary of values for that genome.
  """
  
  # If infile ends in .gz us gzip.open, otherwise use standard open.
  encoding = guess_type(infile)[1]  # uses file extension
  _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
  
  with _open(infile) as f:
    records = list(SeqIO.parse(f, "fasta"))
  
  genome_dict={}
  
  genome_dict['number_of_scaffolds'] = len(records)
  genome_dict['len_seq'] = [len(rec) for rec in records]
  genome_dict['scaffold_lengths'] = pd.DataFrame(genome_dict['len_seq'])
  genome_dict['total_size_scaffolds'] = sum(genome_dict['len_seq'])

  genome_dict['total_scaffold_length_percentage_genome_size'] = ((genome_dict['total_size_scaffolds']/(genome_size*1000000))*100)
    
  genome_dict['seq_greater_1k'] = sorted(i for i in genome_dict['len_seq'] if i>1000)
  genome_dict['seq_greater_10k'] = sorted(i for i in genome_dict['seq_greater_1k'] if i>10000)
  genome_dict['seq_greater_25k'] = sorted(i for i in genome_dict['seq_greater_10k'] if i>25000)
  genome_dict['seq_greater_100k'] = sorted(i for i in genome_dict['seq_greater_25k'] if i>100000)
  genome_dict['seq_greater_1M'] = sorted(i for i in genome_dict['seq_greater_100k'] if i>1000000)
  genome_dict['seq_greater_10M'] = sorted(i for i in genome_dict['seq_greater_1M'] if i>10000000)

  genome_dict['sorted_len'] = sorted(genome_dict['len_seq'], reverse=True)
  genome_dict['sum_sorted_length'] = sum(genome_dict['sorted_len'])
  genome_dict['half_length'] = genome_dict['sum_sorted_length']/2.0

  #calculates N50 and L50 values
  testSum = 0
  genome_dict['N50'] = 0
  genome_dict['L50'] = 0
  i = 0
  half_length = genome_dict['sum_sorted_length']/2.0

  for con in genome_dict['sorted_len']:
      testSum += con
      i += 1
      if half_length < testSum:
        genome_dict['N50'] = con
        genome_dict['L50'] = i
        break

  #calculates NG50 and LG50 values
  half_genome = (genome_size*1000000)/2.0
  testSumNG50 = 0
  genome_dict['NG50'] = 0
  genome_dict['LG50'] = 0
  i = 0
  for conNG50 in genome_dict['sorted_len']:
      testSumNG50 += conNG50
      i += 1
      if  half_genome < testSumNG50:
        genome_dict['NG50'] = conNG50
        genome_dict['LG50'] = i
        break

  #calculates A,C,G,T,N percentages
  genome_dict['counterAT'] = 0
  genome_dict['counterGC'] = 0

  genome_dict['counterN'] = 0
  for record in records:
      genome_dict['counterAT'] += record.seq.count('A') 
      genome_dict['counterAT'] += record.seq.count('T')
      genome_dict['counterGC'] += record.seq.count('G') 
      genome_dict['counterGC'] += record.seq.count('C') 
      genome_dict['counterN'] += record.seq.count('N') 

  return genome_dict

def write_output_stats(genome_stats, genome_size, num_longest, outputfile):
  """ Writes the genome stats to output file
        genome_stats is a dictionary of dictionaries produced by get_genome_stats
        genome_size should be in MBp
        num_longest is int of number of longest scaffold lengths to print
        outputfile is the file to print to
  """
  OUT = open(outputfile, 'w')

  OUT.write('Genome:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(genome)
  OUT.write('\n')

  OUT.write('Number of scaffolds:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['number_of_scaffolds']))
  OUT.write('\n')

  OUT.write('Total sum of scaffold lengths:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['total_size_scaffolds']))
  OUT.write('\n')

  OUT.write('Percent of genome size:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['total_scaffold_length_percentage_genome_size']))
  OUT.write('\n')

  OUT.write('Longest scaffold:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(max(genome_stats[genome]['len_seq'])))
  OUT.write('\n')

  OUT.write('Shortest scaffold:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(min(genome_stats[genome]['len_seq'])))
  OUT.write('\n\n')

  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')

  OUT.write('Total no. scaffolds over 1KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_1k'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 1KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_1k'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 1KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_1k'])/(genome_size*1000000)*100))
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 10KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_10k'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 10KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_10k'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 10KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_10k'])/(genome_size*1000000)*100))
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 25KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_25k'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 25KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_25k'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 25KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_25k'])/(genome_size*1000000)*100))
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 100KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_100k'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 100KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_100k'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 100KBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_100k'])/(genome_size*1000000)*100))
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 1MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_1M'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 1MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_1M'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 1MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_1M'])/(genome_size*1000000)*100))
  OUT.write('\n')

  OUT.write('\n')
  OUT.write('Total no. scaffolds over 10MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(len(genome_stats[genome]['seq_greater_10M'])))
  OUT.write('\n')

  OUT.write(f'Sum of scaffold lengths over 10MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_10M'])))
  OUT.write('\n')

  OUT.write(f'Percent genome over 10MBp:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(sum(genome_stats[genome]['seq_greater_10M'])/(genome_size*1000000)*100))
  OUT.write('\n\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')
  
  OUT.write('N50:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['N50']))
  OUT.write('\n')

  OUT.write('L50:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['L50']))
  OUT.write('\n')

  OUT.write('NG50:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['NG50']))
  OUT.write('\n')

  OUT.write('LG50:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str(genome_stats[genome]['LG50']))
  OUT.write('\n\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')

  OUT.write('%AT:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str((genome_stats[genome]['counterAT']/genome_stats[genome]['total_size_scaffolds'])*100))
  OUT.write('\n')

  OUT.write('%GC:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str((genome_stats[genome]['counterGC']/genome_stats[genome]['total_size_scaffolds'])*100))
  OUT.write('\n')

  OUT.write('%N:')
  for genome in genome_stats:
    OUT.write(',')
    OUT.write(str((genome_stats[genome]['counterN']/genome_stats[genome]['total_size_scaffolds'])*100))
  OUT.write('\n')

  OUT.write(f'\n')
  OUT.write('--------------------------------------------------------------------\n')
  OUT.write('\n')
  # Print longest n scaffold lengths
  OUT.write(f'Longest {num_longest} scaffolds:')
  OUT.write(f'')
  for genome in genome_stats:
    OUT.write(',')
    nlongest=(genome_stats[genome]['scaffold_lengths'].sort_values(0,ascending=False).head(n=num_longest))
    for index in nlongest.index:
      OUT.write(f'{nlongest[0][index]} ')

  OUT.write('\n')

  OUT.close()


def plot_scaffold_dist(genome_stats, outputfile):
  """Print a histogram of the scaffold length distributions
      genome_stats is a dictionary of dictionaries produced by get_genome_stats
      outputfile is the file to print graph to
  """
  for genome in genome_stats:
    plt.hist(genome_stats[genome]['len_seq'], label=genome, bins=30, log=True, alpha=0.5)

  plt.legend()

  plt.xlabel("Scaffold Length", fontsize=16)  
  plt.ylabel("Log count", fontsize=16)

  plt.savefig(outputfile)



