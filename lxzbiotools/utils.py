import os
import re
import sys
import numpy as np
from io import StringIO
from pathlib import Path
from pandas.io.formats.format import GenericArrayFormatter
from Bio import SeqIO
import sys
import copy
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gzip
from mimetypes import guess_type
from functools import partial
from os.path import abspath, expanduser
from lxzbiotools.plot import plotDistribution

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


def stat2str(stats):
    if stats["gene_num"]:
        mean_gene_length = sum(stats["gene_length"])/1.0/len(stats["gene_length"])
    else:
        mean_gene_length = 0

    if stats["CDS_num"]:
        mean_cds_num = sum(stats["CDS_num"])/1.0/len(stats["CDS_num"])
        mean_cds_length = sum(stats["CDS_length"])/1.0/len(stats["CDS_length"])
    else:
        mean_cds_num = mean_cds_length = 0

    if stats["intron_num"]:
        mean_intron_num = sum(stats["intron_num"])/1.0/len(stats["intron_num"])
        mean_intron_length =  sum(stats["intron_length"])/1.0/len(stats["intron_length"])
    else:
        mean_intron_num = mean_intron_length = 0
    if stats["exon_num"]:
       mean_exon_num = sum(stats["exon_num"])/1.0/len(stats["exon_num"])
       mean_exon_length =  sum(stats["exon_length"])/1.0/len(stats["exon_length"])
    else:
       mean_exon_num = mean_exon_length = 0

    return "{0:<10,}\t{1:<10,}\t{2:<10,}\t{3:<10,}\t{4:<10,}\t{5:<10,}\t{6:<10,}".format(stats["gene_num"], float("%.2f" % mean_gene_length), float("%.2f" % mean_cds_length), float("%.2f" % mean_exon_num), float("%.2f" % mean_exon_length), float("%.2f" % mean_intron_num), float("%.2f" % mean_intron_length))

class GffRecord(object):
    def __init__(self, seqid, source, type, start, end, score, strand, phase, attrs):
        try:
            assert "\n" not in seqid
            assert "\n" not in source
            assert "\n" not in type
            assert "\n" not in start
            assert "\n" not in end
            assert "\n" not in score
            assert "\n" not in strand
            assert "\n" not in phase
            assert "\n" not in attrs
        except AssertionError:
            raise ValueError("Invalid GFF record data")

        self.seqid = seqid
        self.source = source
        self._type = type
        self.start = int(start)
        self.end = int(end)
        self.length = self.end-self.start+1
        self.score = score
        self.strand = strand
        self.phase = phase
        self._attrs = attrs
        self.attributes = self._split_attr(attrs)

    def _split_attr(self, attributes):
        r = {}
        contents = attributes.split(";")
        for content in contents:
            if not content:
                continue
            if "=" not in content:
                print("%r is not a good formated attribute: no tag!") 
                continue
            tag, value = content.split("=", 1)
            r[tag] = value

        return r
    
    def to_string(self):
        attr = []
        for key, value in self.attributes.items():
            if key in "ID":
                attr.insert(0, "%s=%s" % (key, value))
            else:
                attr.append("%s=%s" % (key, value))
        r = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.seqid, self.source, self._type, self.start, self.end, self.score, self.strand, self.phase, ";".join(attr))
        return r

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, s):
        self._type = s    

    @classmethod
    def from_string(cls,s):
        try:
            assert "\n" not in s
            parts = s.split("\t")
            assert len(parts) == 9
            seqid, source, type, start, end, score, strand, phase, attributes = parts
           # assert strand in "+-"
        except AssertionError:
            raise ValueError("%r not recognized as a valid GFF record" % s)

        return GffRecord(seqid, source, type, start, end, score, strand, phase, attributes)


def open_gff(fn):
    filename = abspath(expanduser(fn))
    if filename.endswith(".gz"):
        ofs = gzip.open(filename, mode)
    elif filename.endswith(".dexta"):
        ofs = stream_stdout("undexta -vkU -w60 -i", filename)
    else:
        ofs = open(filename)

    for line in ofs.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        if len(line) > 1:
            yield GffRecord.from_string(line)

    ofs.close()
    
    

def gff2dict(fn):
    """
    :param fn:
    :return:
    """

    # link gene to mRNA, mRNA to CDS, exon by parent
    r = {}
    gene_level = {}
    mRNA_level = {}

    for gff in open_gff(fn):
        _type = gff.type
        _attrs = gff.attributes
        _id = None
        _parents = None
        locus_tag = None
        

        if "ID" in _attrs:
            _id = _attrs["ID"]
        if "Parent" in _attrs:
            _parents = _attrs["Parent"]

        if _type == "gene":
            if not _id:
                print("gene %r has no ID, pass" % gff.to_string())
                continue
            if _id in r:  # gene named after mRNA....
                r[_id]["self"] = gff
                continue

            gene_level[_id] = {"self": gff, "mRNA": []}  # construct gene

        if _type in "mRNA":
            if not _id:
                print("mRNA %r has no ID, pass" % gff.to_string())
                continue

            if not _parents:  # mRNA has no parent, reconstruct gene
                # print("mRNA %r has no parent, reconstruct" % _id)
                gene = copy.deepcopy(gff)
                gid = _id+"-gene"
                gene.type = "gene"
                gene.attributes = {"ID": gid}
                gene_level[gid] = {"self": gene, "mRNA": [_id]}

                _attrs["Parent"] = gid
                mRNA_level[_id] = {"self": gff, "CDS": [], "exon": []}
                continue

            if "," in _parents:
                print("mRNA %r belongs to muti gene! pass it" % _id)
                continue

            if _parents not in gene_level:  # name gene by mRNA parent
                gene_level[_parents] = {"self": None, "mRNA": [_id]}
            else:
                gene_level[_parents]["mRNA"].append(_id)

            if _id in mRNA_level:
                mRNA_level[_id]["self"] = gff
            else:
                mRNA_level[_id] = {"self": gff, "CDS": [], "exon": []}

        if _type in "CDS|exon":
            if not _parents:
                print("%s '%s' has no parent, pass" % (_type, gff.to_string()))
                continue

            for parent in _parents.split(","):
                if parent not in mRNA_level:  # name mRNA by CDS parent
                    mRNA_level[parent] = {"self": None, "CDS": [], "exon": []}
                    mRNA_level[parent][_type].append(gff)
                    continue

                mRNA_level[parent][_type].append(gff)

    return gene_level, mRNA_level


def default_filter(gene_level, mRNA_level):

    filtered_gene_level = {}
    filtered_mRNA_level = mRNA_level

    for gid in gene_level:
        mids = gene_level[gid]["mRNA"]

        if not mids:
            print("%s has no mRNA, pass" % gid)
            continue

        mrna_status = 0
        for mid in mids:

            if mid not in mRNA_level:
                print(("%s has no mRNA, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue

            if not mRNA_level[mid]["CDS"]:
                print(("%s has no CDS, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue
            if not mRNA_level[mid]["CDS"]:
                print(("%s has no record, pass" % gid))
                gene_level[gid]["mRNA"].remove(mid)
                continue

            mrna_status = 1

        if mrna_status != 0:
            filtered_gene_level[gid] = gene_level[gid]

    return filtered_gene_level, filtered_mRNA_level


def longest_mRNA_filter(gene_level, mRNA_level):
    """
    :param gene_level:
    :param mRNA_level:
    :return:
    """
    print("Filter gff3 by longest mRNA")

    filtered_gene_level = {}
    filtered_mRNA_level = mRNA_level

    for gid in gene_level:
        mids = gene_level[gid]["mRNA"]

        filtered_gene_level[gid] = gene_level[gid]
        longest_mRNA = None

        for mid in mids:

            if not longest_mRNA:
                longest_mRNA = mid
            elif mRNA_level[mid]["self"].length > mRNA_level[longest_mRNA]["self"].length:
                longest_mRNA = mid
            else:
                pass

        filtered_gene_level[gid]["mRNA"] = [longest_mRNA]
        try:
            filtered_gene_level[gid]["self"].start = mRNA_level[longest_mRNA]["self"].start
            filtered_gene_level[gid]["self"].end = mRNA_level[longest_mRNA]["self"].end
        except:
            print("%s ?" % gid)


    return filtered_gene_level, filtered_mRNA_level


def no_UTR_filter(gene_level, mRNA_level):
    """
    :param gene_level:
    :param mRNA_level:
    :return:
    """
    print("Filter gff3 by no UTR")

    filtered_gene_level = {}
    filtered_mRNA_level = {}

    for mid, mv in mRNA_level.items():

        cds = sorted(mv["CDS"], key=lambda i: i.start)
        try:
            mv["self"].start = cds[0].start
            mv["self"].end = cds[-1].end
        except:
            print("%r is not a mRNA, maybe a ncRNA." % mid)
            continue
        mv["CDS"] = cds
        mv["exon"] = []
        filtered_mRNA_level[mid] = mv

    for gid, gv in gene_level.items():
        mids = gene_level[gid]["mRNA"]
        if not mids:
            print("%s has no mRNA, pass" % gid)
            continue

        start = end = 0
        for mid in mids:

            if start == 0 or end == 0:
                start = filtered_mRNA_level[mid]["self"].start
                end = filtered_mRNA_level[mid]["self"].end
                continue

            if filtered_mRNA_level[mid]["self"].start < start:
                start = filtered_mRNA_level[mid]["self"].start
            if filtered_mRNA_level[mid]["self"].end > end:
                end = filtered_mRNA_level[mid]["self"].end

        if hasattr(gv["self"], "start"):
            gv["self"].start = start
            gv["self"].end = end
            filtered_gene_level[gid] = gv
        else:
            logging.warning("%s format error" % gid)
            

    return filtered_gene_level, filtered_mRNA_level


def gff_format(gene_level, mRNA_level):
    r = {}

    # reconstruct exon and intron if necessary
    new_mRNA_level = {}
    for mid, mv in mRNA_level.items():
        if not mv["self"]: # mRNA has no named, pass
            continue

        if not mv["CDS"]:  # mRNA has no CDS, pass
            continue

        new_mRNA_level[mid] = mv
        cds = sorted(mv["CDS"], key=lambda i: i.start)
        new_mRNA_level[mid]["CDS"] = cds

        n = len(cds)
        if not mv["exon"]:  # construct exon
            # print("mRNA %s has no exon, reconstruct from CDS and mRNA" % mid)

            for i in range(n):
                exon = copy.deepcopy(cds[i])
                exon.type = "exon"
                exon.phase = "."
                if i == 0:
                    exon.start = mv["self"].start
                if i == n - 1:
                    exon.end = mv["self"].end
                new_mRNA_level[mid]["exon"].append(exon)
        else:
            new_mRNA_level[mid]["exon"] = sorted(mv["exon"], key=lambda i: i.start)

        # construct intron
        new_mRNA_level[mid]["intron"] = []
        for i in range(n):
            if i == 0:
                continue
            start = cds[i - 1].end + 1
            end = cds[i].start - 1
            if start >= end:
                continue

            intron = copy.deepcopy(cds[i])
            intron.type = "intron"
            intron.phase = "."
            intron.start = start
            intron.end = end
            if intron.length <= 0:
                print("intron %s <= 0" % intron.to_string())
                continue
            new_mRNA_level[mid]["intron"].append(intron)

    # add mRNA to gene
    for gid, gv in gene_level.items():
        if not gv["mRNA"]:
            print("%s has no mRNA, pass" % gid)
            continue

        if not gv["self"]:
            print("%s has no gff, pass" % gid)
            continue

        r[gid] = {"self": gv["self"], "mRNA": {}}
        for mid in gv["mRNA"]:
            r[gid]["mRNA"][mid] = new_mRNA_level[mid]
   # print(r)
    return r


def gff2str(gffdict, fn):

    out = open(fn, "w")
    for gid, gv in sorted(gffdict.items(), key=lambda d:(d[1]["self"].seqid, d[1]["self"].start)):
        out.write(gv["self"].to_string()+"\n")
        for mid, mv in gv["mRNA"].items():
            out.write(mv["self"].to_string()+"\n")
            for exon in mv["exon"]:
                out.write(exon.to_string()+"\n")
            for cds in mv["CDS"]:
                out.write(cds.to_string()+"\n")
    out.close()

    return fn


def stat(gffdict):
    r = {"gene_num": len(gffdict), "gene_length": [],
         "CDS_num": [], "CDS_length": [],
         "intron_num": [], "intron_length": [],
         "exon_num": [], "exon_length": [], }

    for gid, gv in sorted(gffdict.items(), key=lambda d:(d[1]["self"].seqid, d[1]["self"].start)):
        r["gene_length"].append(gv["self"].end-gv["self"].start+1)

        for mid, mv in gv["mRNA"].items():
            r["CDS_num"].append(len(mv["CDS"]))
            CDS_length = sum([i.end-i.start+1 for i in mv["CDS"]])
            r["CDS_length"].append(CDS_length)
            if mv["intron"]:
                r["intron_num"].append(len(mv["intron"]))
                intron_length = [i.end-i.start+1 for i in mv["intron"]]
                r["intron_length"] += intron_length
            else:
                r["intron_num"].append(0)

            r["exon_num"].append(len(mv["exon"]))
            exon_length = [i.end-i.start+1 for i in mv["exon"]]
            r["exon_length"] += exon_length

    return r


def gff_stat(fn, outdir):

    print("Loading gff3 information from %s" % fn)
    gene, mRNA = gff2dict(fn)

    filtered_gene_level, filtered_mRNA_level = default_filter(gene, mRNA)

    filtered_gene_level, filtered_mRNA_level = longest_mRNA_filter(filtered_gene_level, filtered_mRNA_level)
    filtered_gene_level, filtered_mRNA_level = no_UTR_filter(filtered_gene_level, filtered_mRNA_level)
    print("Format gff3 loaded")
    filtered_gff = gff_format(filtered_gene_level, filtered_mRNA_level)
    print("Statistics on gff3")
    stats = stat(filtered_gff)
    print("Output filted gff3")
    gff2str(filtered_gff, "%s/%s.longest.gff3" % (outdir, os.path.basename(fn)))
    return stats


def stat2str(stats):
    if stats["gene_num"]:
        mean_gene_length = sum(stats["gene_length"])/1.0/len(stats["gene_length"])
    else:
        mean_gene_length = 0

    if stats["CDS_num"]:
        mean_cds_num = sum(stats["CDS_num"])/1.0/len(stats["CDS_num"])
        mean_cds_length = sum(stats["CDS_length"])/1.0/len(stats["CDS_length"])
    else:
        mean_cds_num = mean_cds_length = 0

    if stats["intron_num"]:
        mean_intron_num = sum(stats["intron_num"])/1.0/len(stats["intron_num"])
        mean_intron_length =  sum(stats["intron_length"])/1.0/len(stats["intron_length"])
    else:
        mean_intron_num = mean_intron_length = 0
    if stats["exon_num"]:
       mean_exon_num = sum(stats["exon_num"])/1.0/len(stats["exon_num"])
       mean_exon_length =  sum(stats["exon_length"])/1.0/len(stats["exon_length"])
    else:
       mean_exon_num = mean_exon_length = 0

    return "{0:<10,}\t{1:<10,}\t{2:<10,}\t{3:<10,}\t{4:<10,}\t{5:<10,}\t{6:<10,}".format(stats["gene_num"], float("%.2f" % mean_gene_length), float("%.2f" % mean_cds_length), float("%.2f" % mean_exon_num), float("%.2f" % mean_exon_length), float("%.2f" % mean_intron_num), float("%.2f" % mean_intron_length))


def plot_image(labels, data_dict, outdir):
    plot_cfg = {
        "gene_length": {"window": 200, "x_max": 20000, "x_min": 0, "title": "Distribution of gene length", "x_label": "Gene length (bp)", "y_label": "Percentage (%)"},
        "CDS_length": {"window": 200, "x_max": 8000, "x_min": 0, "title": "Distribution of CDS length", "x_label": "CDS length per gene (bp)", "y_label": "Percentage (%)"},
        "exon_length": {"window": 50, "x_max": 2000, "x_min": 0, "title": "Distribution of exon length", "x_label": "Exon length (bp)", "y_label": "Number"},
        "intron_length": {"window": 50, "x_max": 2000, "x_min": 0, "title": "Distribution of intron length", "x_label": "Intron length (bp)", "y_label": "Percentage (%)"},
        "exon_num": {"window": 1, "x_max": 30, "x_min": 0, "title": "Distribution of exon number", "x_label": "Exon number per gene", "y_label": "Percentage (%)"},
        "intron_num": {"window": 1, "x_max": 30, "x_min": 0, "title": "Distribution of intron number", "x_label": "Intron number per gene", "y_label": "Percentage (%)"}
    }

    for p in plot_cfg:

        data = []
        for label in labels:
            d = data_dict[label]
            with open (os.path.join(outdir, "%s.%s.len" % (label, p)), "w") as out:
                 out.write("\n".join(map(str, sorted(d[p], reverse=True))))
            data.append(d[p])
        print("plot %s %s" % (p, labels))
        plotDistribution(data, label=labels, prefix=str(outdir)+"/"+str(p), **plot_cfg[p])


def read_cfg(fn):
    r = []
    with open(fn) as fh:
         for line in fh.readlines():
             line = line.strip()
             if len(line) <=1:
                 continue

             lines = line.split()
             r.append(lines)

    return r
