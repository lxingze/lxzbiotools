## Introduction

**lxzbiotools** contains some bioinformation data processing scripts.

Obviously, it is not perfect, so it will be updated later until it is strong enough

## Install

```
pip3 install lxzbiotools
```

## Use

```bash
$ python3 ~/lxzbiotools/lxzbiotools/lxzbiotools.py  --help

 Usage: lxzbiotools.py [OPTIONS] COMMAND [ARGS]...

 Xingze Li's bioinformatics analysis scripts.
 emali: lixingzee@gmail.com

╭─ Options ───────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ --install-completion        [bash|zsh|fish|powershell|pwsh]  Install completion for the specified shell.        │
│                                                              [default: None]                                    │
│ --show-completion           [bash|zsh|fish|powershell|pwsh]  Show completion for the specified shell, to copy   │
│                                                              it or customize the installation.                  │
│                                                              [default: None]                                    │
│ --help                                                       Show this message and exit.                        │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ──────────────────────────────────────────────────────────────────────────────────────────────────────╮
│ cds2pep         Convert cds file to pep file                                                                    │
│ excel2txt       Convert excel file to txt file                                                                  │
│ fa2fq           Convert a fasta file to a fastq file                                                            │
│ fq2fa           Convert a fastq file to a fasta file                                                            │
│ genstats        single or multiple genome information statistics                                                │
│ gfa2fa          Convert gfa file to fasta file                                                                  │
│ gff             Simplify gff3 file for WGD event analysis                                                       │
│ len             Get the length of each sequence                                                                 │
│ rds             Read a multi-FASTA file sequence and remove duplicates (by MD5 hash)                            │
│ run             Parallelized running tasks                                                                      │
│ seq             Extract sequences by sequence name or keyword                                                   │
│ tolf            Convert multi-line fasta to one-line fasta                                                      │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

```
