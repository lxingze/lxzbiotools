import os, re, sys
import logging
import time
import typer
import math
import hashlib
import shutil
import pandas as pd
import subprocess
from rich import print
from rich.progress import track
from typing import List, Union, Optional, NamedTuple
from openpyxl import load_workbook
from openpyxl.workbook.workbook import Workbook
from openpyxl.worksheet.worksheet import Worksheet
from collections import defaultdict
from collections import OrderedDict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from lxzbiotools.utils import *


# start=time.perf_counter()
app=typer.Typer(help="Xingze Li's bioinformatics analysis scripts. \n\n\n emali: lixingzee@gmail.com")

@app.command("m2ofa")
def to_one_line_fasta(inf: Path = typer.Argument(...,
                                                 help="input fasta file"),
                      outf: Path = typer.Argument(...,
                                                 help="output fasta file")):
    'Convert multi-line fasta to one-line fasta'
    inf = open(inf, 'r')
    outf = open(outf, 'w')
    db = {}
    for line in track(inf.readlines()):
        if line.startswith('>'):
            keys = line.strip()
            db[keys] = []
        else:
            db[keys].append(line.strip())
    for name, seq in db.items():
        outf.write(name + '\n')
        outf.write(''.join(seq) + '\n')
    inf.close()
    outf.close()

@app.command('extseq')
def extract_target_seq(inf: Path = typer.Argument(...,
                                                  help="input fasta file"), 
                       outf: Path = typer.Argument(...,
                                                   help="output result file"), 
                       ext: Path = typer.Argument(...,
                                                                    help="extract sequences name file")):
    'Extract sequences by sequence name or keyword'
    outf = open(outf,'w')
    dict = {}
    with open(inf, 'r') as fastaf:
        for line in fastaf:
            if line.startswith('>'):
                name = line
                dict[name] = ''
            else:
                dict[name] += line.replace('\n','') 
         
    with open(ext,'r') as listf:
        for row in listf:
            row = row.strip()
            for key in dict.keys(): 
                if row in key:
                    outf.write(key)
                    outf.write(dict[key] + '\n')
    outf.close()

@app.command()
def fq2fa(fq: Path = typer.Argument(...,
                                    help="input fastq file"),
          fa: Path = typer.Argument(...,
                                    help="output fasta file")):
    'Convert a fastq file to a fasta file'
    out=[]
    with open(fq,'r') as fin:
        c=0
        for line in fin:
            c+=1
            if (c%4) == 1:
                out.append('>' + line[1:])
            elif (c%4) == 2:
                out.append(line)
    with open(fa,'w') as fo:
        fo.write(''.join(out))
    return

@app.command()
def fa2fq(fa: Path = typer.Argument(...,
                                    help="input fasta file"),
          fq: Path = typer.Argument(...,
                                    help="output fastaq file")):
    'Convert a fasta file to a fastq file'
    seq=''
    out=[]
    n_reads=0
    with open(fa,'r') as fin:
        for line in track(fin):
            if line.startswith('>'):
                n_reads+=1
                if seq:
                    score='I'*len(seq)
                    out.append(seq+'\n+\n'+score+'\n')
                    seq=''
                out.append('@'+line[1:])
            else:
                seq+=line.strip()
    score='I'*len(seq)
    out.append(seq+'\n+\n'+score+'\n')
    with open(fq,'w') as fo:
        fo.write(''.join(out))
    return n_reads



@app.command()
def cds2pep(cds_input: Path = typer.Argument(...,
                                              help="input coding sequence(CDS) file"),
            cds_output: Path = typer.Argument(...),
            pep_output: Path = typer.Argument(...,
                                              help="output pep file")):
    'Convert cds file to pep file'
    inf = open(cds_input, 'r')
    cds_outf = open(cds_output, 'w')        
    pep_outf = open(pep_output, 'w')
    db = {}
    pep = {}
    for line in track(inf.readlines()):
        if line.startswith('>'):
            keys = line.strip()
            db[keys] = []
        else:
            db[keys].append(line.strip())

    for name, seq in track(db.items()):
        cds_outf.write(name + '\n')
        cds_outf.write(''.join(seq) + '\n')
    cds_outf.close()

    with open(cds_output, "r") as co:
        for line in co.readlines():
            if line.startswith('>'):
                key = line.strip()
                pep[key] = []
            else:
                seq = line.strip()
                # translate one frame at a time
                prot = '' 
                for i in range(0, len(seq), 3):
                    codon = seq[i:i + 3]
                    if codon in codon_table:
                        if codon_table[codon] == 'STOP':
                            prot = prot + '*'
                        else: 
                            prot = prot + codon_table[codon]
                pep[key].append(prot)
        for name, seq in track(pep.items()):
            pep_outf.write(name + '\n')
            pep_outf.write(''.join(seq) + '\n')
    inf.close()
    pep_outf.close()



@app.command()
def excel2txt(input: Path,
              output: Path):
    'Convert excel file to txt file'
    df = pd.read_excel(input, header=None)		
    df.to_csv(output, header=None, sep='\t', index=False)	

@app.command('length')
def get_fasta_len(fasta_file: Path = typer.Argument(...,
                                              help="input fasta file"),
                  len_file: Path = typer.Argument(...,
                                              help="ouput lens file")):
    'Get the length of each sequence'
    fasta_dict = OrderedDict()
    handle = open(fasta_file, "r")
    lenf = open(len_file, "w")
    active_sequence_name = ""
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"): 
            active_sequence_name = line[1:]
            active_sequence_name = active_sequence_name.split(" ")[0]
        if active_sequence_name not in fasta_dict:
            fasta_dict[active_sequence_name] = 0
            continue
        sequence = line
        fasta_dict[active_sequence_name] += len(sequence)
    print("seqs\tlens", file=lenf)
    for chrom,lens in track(fasta_dict.items()):
        print(f"{chrom}\t{lens}", file=lenf)
    lenf.close()    
    handle.close()
    return fasta_dict


@app.command('gff')
def parse_gff(gff: Path,
              sample_gff_file: Path):
    'Simplify gff3 file for WGD event analysis'
    gene_count = {}
    gene_dict = OrderedDict()
    tx_pos_dict = defaultdict(list)
    CDS_dict = defaultdict(list)
    handle = open(gff, "r")
    gff_file = open(sample_gff_file, "w")    
    for line in track(handle):
        if line.startswith("#"):
            continue
        content = line.split("\t")
        if len(content) <= 8:
            continue
        #print(content)
        if content[2] == "transcript" or content[2] == "mRNA":

            # use regual expression to extract the gene ID
            # match the pattern ID=xxxxx; or ID=xxxxx
            
            tx_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
            tx_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
            tx_parent = tx_parent.strip() # remove the 'r' and '\n'

            # if the parent of transcript is not in the gene_dict, create it rather than append
            if tx_parent in gene_dict:
                gene_dict[tx_parent].append(tx_id)
            else:
                gene_dict[tx_parent] = [tx_id]
            tx_pos_dict[tx_id] = [content[0],content[3], content[4], content[6] ]
        # GFF must have CDS feature
        if content[2] == 'CDS':
            width = int(content[4]) - int(content[3])
            CDS_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
            CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
            CDS_dict[CDS_parent].append(width)
    
    for gene, txs in gene_dict.items():
        tmp = 0
        for tx in txs:
            tx_len = sum(CDS_dict[tx])
            if tx_len > tmp:
                lst_tx = tx
                tmp = tx_len
        tx_chrom = tx_pos_dict[lst_tx][0]
        if tx_chrom not in gene_count:
            gene_count[tx_chrom] = 1
        else:
            gene_count[tx_chrom] += 1
        tx_start = tx_pos_dict[lst_tx][1]
        tx_end   = tx_pos_dict[lst_tx][2]
        tx_strand = tx_pos_dict[lst_tx][3]

        print(f"{tx_chrom}\t{gene}\t{tx_start}\t{tx_end}\t{tx_strand}\t{gene_count[tx_chrom]}\t{lst_tx}", file=gff_file )
    gff_file.close()
    handle.close()
    # return [gene_dict, tx_pos_dict, CDS_dict]




@app.command('parallel')
def parallel_run(cmd_path, 
                 thread, 
                 stdout=None, 
                 stderr=None, 
                 parallel="parallel", 
                 fg=False):
    'Parallelized running tasks'
    thread=int(thread)
    if stdout == None:
        stdout = "%s.log" % cmd_path
    if stderr == None:
        stderr = "%s.err" % cmd_path
    if thread >1:
        cmd = "nohup %s -j %s ::: < %s > %s 2> %s " % (parallel, thread, cmd_path, stdout, stderr)
    else:
        cmd = "nohup bash %s > %s 2> %s " %(cmd_path, stdout, stderr)
    if fg==False:
        cmd +=" &"
    sys.stdout.write("CMD: %s \n" % cmd)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')



@app.command()
def gfa2fa(gfa_file: Path,
           fa_file: Path):
    'Convert gfa file to fasta file'
    input = open(gfa_file, 'r').readlines()
    output = open(fa_file, 'w')
    print(f"Converting GFA {input} --> FASTA {output}...")
    num_seqs = 0
    for i in track(len(input)):
        line = input[i].strip('\n').split()
        if line[0] == 'S':
            num_seqs += 1
            simple_seq = Seq(line[2])
            simple_seq_r = SeqRecord(simple_seq, id=line[1])
            SeqIO.write(simple_seq_r, output, "fasta")
    output.close()
    print(f"FASTA of {num_seqs} sequences created at {output}.")

@app.command('rds')
def remove_duplicate_sequences(fasta_file:Path,
                               output_file:Path):
    'Read a multi-FASTA file sequence and remove duplicates (by MD5 hash)'
    ## output will either be a file or STDOUT
    '''
    Input: 
        - A (multi-)FASTA file.  
    Output:
        - Output will be multi-FASTA data printed to STDOUT or a file
    '''
    fout = sys.stdout
    if output_file is not None:
        fout = open(output_file, 'wt')
    
    seqs = fasta_dict(fasta_file)
    found = dict()

    m = hashlib.md5()

    for seq_id in track(seqs):
        seq = seqs[seq_id]['s']
        m.update(seq.encode())
        md5sum = m.hexdigest()
        
        ## write this sequence, 60bp per line
        if md5sum not in found:
            fout.write(f">{seq_id} {seqs[seq_id]['h']}\n")
            for i in range(0, len(seq), 60):
                fout.write(seq[i : i + 60] + "\n")
            found[md5sum] = 1


@app.command('genstats')
def genome_stats(infile: List[str] = typer.Option(...,"-i","--infile",
                                                help="input genome fasta"),
        genome_size: int = typer.Option(...,"-s","--size",
                                                help="input your genome size (MBp)"),
        outprefix: str = typer.Option("results","-o","--outprefix",
                                            help="specify output file prefix"),
        num_longest: int = typer.Option("10","-n","--nums",
                                          help="print the length of the longest scaffolds")):
    'single or multiple genome information statistics'
    
    print(f"Got {len(infile)} genome(s) to compare.")
        
    # Genome stats are stored in a dict of dicts, one for each genome.
    genomes_dict = {}
    append = 1

    # For each genome get the stats.
    for genome in infile:
        genome_name = os.path.split(genome)[1]

        if genome_name in genomes_dict:
            print(f"Multiple genome assemblies found with name {genome_name}, appending numbers to subsequent genomes")
            genome_name = genome_name + "_" + str(append)
            append += 1

        print(f"Getting statistics for {genome_name}.")
        genomes_dict[genome_name] = get_genome_stats(genome, genome_size)
        print(f"Found {(genomes_dict[genome_name]['number_of_scaffolds']):,} contigs for {genome_name}.")

    print(f"Done getting stats for {len(infile)} genomes. \nSummarizing data.")
    
    out_fasta = outprefix + '.csv'
    write_output_stats(genomes_dict, genome_size, num_longest, out_fasta)

    out_plots = outprefix + '.pdf'
    plot_scaffold_dist(genomes_dict, out_plots)

@app.command()
def movefile(file_dir: Path = typer.Argument(...,
                                    help="input your target files directory"),
             output_dir: Path = typer.Argument(...,
                                    help="output directory name"),
             dir_nums: int = typer.Argument(...,
                                    help="number of folders you want to create "),
             prefix: str = typer.Argument("",
                                    help="Folder prefix name, defaults to 1,2,3... 100")):
    'Randomly allocate files to a specified number of folders'
    file_dir = os.path.abspath(file_dir)
    folderPath = os.path.abspath(output_dir)
    sort_folder_number = [x for x in range(1,dir_nums+1)]
    for num in sort_folder_number:
        addfolder = prefix + str(num)
        new_folder_path = os.path.join(folderPath,addfolder) #new_folder_path is ‘folderPath\number'
        if not os.path.exists(new_folder_path):
            os.makedirs(new_folder_path)
            print("new a floder named "+str(num)+'at the path of '+ new_folder_path)

    #the taget file list
    file_list = os.listdir(file_dir)
    folderNumber = 0
    print('there are '+str(len(file_list))+' files at the path of '+file_dir)
    for i in range(0,len(file_list)):
        old_file_path = os.path.join(file_dir,file_list[i])
        if os.path.isdir(old_file_path):
            '''if the path is a folder,program will pass it'''
            print('file does not exist ,path=' + old_file_path+' it is a dir' )
            pass
        elif not os.path.exists(old_file_path):
            '''if the path does not exist,program will pass it'''
            print('file does not exist ,path='+old_file_path)
            pass
        else:
            '''define the number,it decides how many files each directory process'''
            number = math.ceil(len(file_list)/dir_nums)
            if(i%number == 0):
                folderNumber +=1
                addfile = prefix + str(folderNumber)
            new_file_path = os.path.join(folderPath,'%s'%(addfile))
            print(new_file_path)
            if not os.path.exists(new_file_path):
                print('not exist path:'+new_file_path)
                break
            shutil.move(old_file_path,new_file_path)
            print('success move file from '+ old_file_path +' to '+new_file_path)
    print('If you have any questions, please contact [bold magenta]lixingzee@gmail.com[/bold magenta]')
    print('[bold cyan]Thank you for using.[/bold cyan]')



@app.command()
def gffstat(gff3_file: Path = typer.Argument(...,
                                    help="input your target gff3 file"),
            output_dir: Path = typer.Argument(...,
                                    help="output directory")):
    'Format conversion and various information statistics of genome annotation file gff'
    fn = gff3_file
    outdir = output_dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out = open(f"{outdir}/gff.stat", "w")
    if gff3_file:
        stat_dict = gff_stat(fn, outdir)
        name = os.path.basename(fn)
        data = {name: stat_dict}
        names = [name]
        #plot_image(names, data, outdir)
        out.write("""
Number of gene	\tAverage gene length(bp)	\tAverage CDS length(bp)	\tAverage exons per gene	\tAverage exon length(bp)	\tAverage introns per gene	\tAverage intro length(bp)
""")
        out.write(stat2str(stat_dict))
        plot_image(names, data, outdir)
        out.close()


# end=time.perf_counter()
# print('Running time: %s Seconds'%(end-start))   



if __name__ == "__main__":
    app()


