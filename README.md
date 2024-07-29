[![Anaconda](https://anaconda.org/bioconda/carpedeam/badges/version.svg)](https://anaconda.org/bioconda/carpedeam) [![BioConda Install](https://anaconda.org/bioconda/carpedeam/badges/downloads.svg)](https://anaconda.org/bioconda/carpedeam)

# CarpeDeam - A metagenomic assembler for heavily damaged ancient datasets

CarpeDeam is a damage-aware metagenome assembler for ancient metagenomic DNA datasets. It takes (merged) reads and a damage matrix as input and prooved to work best for heavily damaged datasets.
CarpeDeam is built upon PenguiN (https://github.com/soedinglab/plass).

[Insert Citation Here](https://www.google.com).
 
# Install CarpeDeam

CarpeDeam can be install via [conda](https://github.com/conda/conda) or as statically compiled Linux version.

## Option 1 (easy): Static Binary

Step 1: Download the static binary. The list of releases is here: https://github.com/LouisPwr/CarpeDeam/tags Each release comes with a static binary. Find the URL of the static binary by right-clicking and selecting "copy link"

```
wget [URL TO RELEASE BINARY]
```

Where you paste the URL of the binary. If you have root access, simply install the executable by running

```
sudo cp CarpeDeam /usr/bin
```

Otherwise, just leave the executable where it is or copy it in a bin/ directory in your home directory:

```
mkdir -P $HOME/bin/
cp CarpeDeam $HOME/bin/
```

Step 2: Mark the binary executable:

```
chmod +x CarpeDeam
```

## Option 2 (easy): Conda

Install Conda (https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and Bioconda (https://bioconda.github.io/) and type the following:

```
 conda  install -c bioconda carpedeam
```
Thanks to martin-g for the bioconda support!

## Option 3 (cool): Compile from source

Compiling CarpeDeam from source optimizes it for your specific system, potentially boosting performance. 

Prerequisites:
- Git
- G++ (version 4.6 or higher)
- CMake (version 3.0 or higher)

After compilation, find the PLASS binary in the `build/bin` directory.


      git clone https://github.com/LouisPwr/CarpeDeam
      cd CarpeDeam
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4 && make install
      export PATH="$(pwd)/bin/:$PATH"


## How to assemble
We recommend merging reads before assembly as merging improves read quality. You need to supply a damage matrix (format see FAQ) to actually make it work.

      # assemble merged reads 
      carpedeam ancient_assemble examples/reads.fastq.gz assembly.fasta tmp --ancient-damage /path/to/damage/prefix

Important parameters: 

     --ancient-damage         Path to damage matrix. Must be tab separated file (format see FAQ)
     --num-iter-reads-only    Number of iterations which use raw-reads only (=sequences that have not been extended yet); damage aware; results in short (<500bp) but precise contigs; we suggest 3 to 5 iterations
     --num-iterations         Number of total iterations; somewhat damage aware; results in long contigs; we suggest 9 to 12 iterations
     --min-merge-seq-id       Minimum sequence identity threshold of overlapping sequences in contig-merging-step (=after correction). Lower than 99% will result in even longer contigs, but increases the risk of misassemblies
     --unsafe                 Turn On  to maximize contig lengths and sensitivity. Comes with a higher risk of assembling chimeric contigs. 
      





#### Dependencies

When compiling from source, PLASS requires `zlib` and `bzip`.

### Use the docker image

## Hardware requirements
CarpeDeam needs roughly 1 byte of memory per residue to work efficiently. CarpeDeam will scale its memory consumption based on the available main memory of the machine. CarpeDeam needs a CPU with at least the SSE4.1 instruction set to run. 

# FAQ 
* **Which format does the damage matrix need to have?** (Template see example dir in this repository)
  
  It must be a tab seperated file with a header line for the individual substitutions and one line per position. The C-to-T substitutions start at position 1 of the five-prime-end, the G-to-A substitutions start at position 1 of the three-prime-end.
  
  **Important: The damage matrices need to have the suffix 3p.prof and 5p.prof**
  
  myDamage_5p.prof could look like:
  
      A>C     A>G     A>T     C>A     C>G     C>T     G>A     G>C     G>T     T>A     T>C     T>G
      0       0       0       0       0       0.329405        0       0       0       0       0       0
      0       0       0       0       0       0.221745        0       0       0       0       0       0
      0       0       0       0       0       0.187678        0       0       0       0       0       0
      0       0       0       0       0       0.161196        0       0       0       0       0       0
      0       0       0       0       0       0.144011        0       0       0       0       0       0
  
   myDamage_3p.prof could look like:
  
       A>C     A>G     A>T     C>A     C>G     C>T     G>A     G>C     G>T     T>A     T>C     T>G
      0       0       0       0       0       0       0.32891 0       0       0       0       0
      0       0       0       0       0       0       0.223405        0       0       0       0       0
      0       0       0       0       0       0       0.188599        0       0       0       0       0
      0       0       0       0       0       0       0.164419        0       0       0       0       0
      0       0       0       0       0       0       0.146352        0       0       0       0       0

* **Do all other substitution rates have to be 0?**
  
  No, but focusing on C-to-T and G-to-A substitutions is the best approach (for now):
  
  Typical ancient DNA damage:
  C-to-T and G-to-A are the most common substitutions in ancient DNA. These are primarily caused by chemical processes over time, not evolution.
  
  Other substitution types:
  Non-zero rates for other substitutions usually indicate evolutionary changes or sequencing errors (among others). Including these could increase the risk of assembly errors.
  
  Focusing on C-to-T and G-to-A helps avoid incorrect "corrections" and ensures more accurate species-specific assembly.
  
* **How do I correctly porvide the damage matrix file path?**
  Assume our damage matrices are located at:

      /path/to/emotional/damage/sampleDamage_3p.prof
      /path/to/emotional/damage/sampleDamage_5p.prof

  The you can call CarpeDeam as follows:

      carpedeam ancient_assemble myInput.fasta myOutput.fasta tmpDir --ancient-damage /path/to/emotional/damage/sampleDamage_


