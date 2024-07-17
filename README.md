# CarpeDeam - A metagenomic assembler for heavily damaged ancient datasets


CarpeDeam - short description

[Insert Citation Here](https://www.google.com).
 
### Install CarpeDeam
CarpeDeam can be install via [conda](https://github.com/conda/conda) or as statically compiled Linux version. CarpeDeam requires a 64-bit Linux/MacOS system (check with `uname -a | grep x86_64`) with at least the SSE2 instruction set.

     # install from bioconda
     conda install -c conda-forge -c bioconda carpedeam 

     provide static builds!!!
     # static build with AVX2 (fastest)
     wget https://mmseqs.com/plass/plass-linux-avx2.tar.gz; tar xvfz plass-linux-avx2.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
     # static build with SSE4.1
     wget https://mmseqs.com/plass/plass-linux-sse41.tar.gz; tar xvfz plass-linux-sse41.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
     # static build with SSE2 (slowest, for very old systems)
     wget https://mmseqs.com/plass/plass-linux-sse2.tar.gz; tar xvfz plass-linux-sse2.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
 

## How to assemble
Plass can assemble both paired-end reads (FASTQ) and single reads (FASTA or FASTQ):

      # assemble paired-end reads 
      plass assemble examples/reads_1.fastq.gz examples/reads_2.fastq.gz assembly.fas tmp

      # assemble single-end reads 
      plass assemble examples/reads_1.fastq.gz assembly.fas tmp

      # assemble single-end reads using stdin
      cat examples/reads_1.fastq.gz | plass assemble stdin assembly.fas tmp


Important parameters: 

     --min-seq-id         Adjusts the overlap sequence identity threshold
     --min-length         minimum codon length for ORF prediction (default: 40)
     -e                   E-value threshold for overlaps 
     --num-iterations     Number of iterations of assembly
     --filter-proteins    Switches the neural network protein filter off/on

Modules: 

      plass assemble      Assembles proteins (i:Nucleotides -> o:Proteins)
      plass nuclassemble  Assembles nucleotides *experimental* (i:Nucleotides -> o:Nucleotides)
      


### Compile from source
Compiling CarpeDeam from source has the advantage that it will be optimized to the specific system, which should improve its performance. To compile PLASS `git`, `g++` (4.6 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the PLASS binary will be located in the `build/bin` directory.

      git clone https://github.com/LouisPwr/CarpeDeam
      cd plass
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4 && make install
      export PATH="$(pwd)/bin/:$PATH"
        

#### Dependencies

When compiling from source, PLASS requires `zlib` and `bzip`.

### Use the docker image

## Hardware requirements
CarpeDeam needs roughly 1 byte of memory per residue to work efficiently. CarpeDeam will scale its memory consumption based on the available main memory of the machine. CarpeDeam needs a CPU with at least the SSE4.1 instruction set to run. 

## FAQ 
* One
* Two
