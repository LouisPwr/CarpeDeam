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
We recommend merging reads before assembly as merging improves read quality. You need to supply a damage matrix (format see below) to actually make it work.

      # assemble merged reads 
      carpedeam ancient_assemble examples/reads_1.fastq.gz assembly.fasta tmp --ancient-damage /path/to/damage/prefix

Important parameters: 

     --ancient-damage         Path to damage matrix. Must be tab separated file (format see below)
     --num-iter-reads-only    Number of iterations which use raw-reads only (=sequences that have not been extended yet); damage aware; results in short (<500bp) but precise contigs; we suggest 3 to 5 iterations
     --num-iterations         Number of total iterations; somewhat damage aware; results in long contigs; we suggest 9 to 12 iterations
     --min-merge-seq-id       Minimum sequence identity threshold of overlapping sequences in contig-merging-step (=after correction). Lower than 99% will result in even longer contigs, but increases the risk of misassemblies
      










#### Dependencies

When compiling from source, PLASS requires `zlib` and `bzip`.

### Use the docker image

## Hardware requirements
CarpeDeam needs roughly 1 byte of memory per residue to work efficiently. CarpeDeam will scale its memory consumption based on the available main memory of the machine. CarpeDeam needs a CPU with at least the SSE4.1 instruction set to run. 

## FAQ 
* One
* Two
