## Building

### Build requirements

One needs `GNU make`, `Cmake 3+`, `Boost` libraries.
If you don't know if you have them or
not just type `make packages` and it will prompt you for install.

Also, need to have SSE4.1 or higher. If AVX2 is available
the second step will consume two times more sequences and will therefore work two times faster. In the second step
we use https://github.com/Martinsos/opal

### How to build?

Build by writing `make` in your command line. If everything went without errors
executables should be present in the directory called `build`. There is executable file
tachyon.

## Running
There are three modes in which you can run this application. Use `makedb` if you want to set up reduced database. If
you run

`$ tachyon makedb -d nr.fa -r reduced_nr`

it will create reduced database (required by algorithm) in the current directory.

You can run alignment task with the following commands:

`$ tachyon blastp -d nr.fa -i reduced_nr -q queries.fa -o results`
if you want to align protein sequences against a protein reference database, or

`$ tachyon blastx -d nr.fa -i reduced_nr -q queries.fa -o results`

if you want to align DNA sequences against a protein reference database.

##Commands

This commands determine the mode in which you want to run this app

| Command       | Description                                                             |
| ------------- |:----------------------------------------------------------------------- |
| makedb        | Create reduced database from FASTA reference file.                      |
| blastp        | Align protein query sequences against a protein reference database      |
| blastx        | Align DNA query against a protein reference database                    |

###General options
This options is necessary in all modes.

| Option        |  Short | Default | Description       |
| --------------|--------| --------| ------------------|
|--threads      | -p     |   8     | number of threads |
|--db           | -d     |         | path to reference database |
|--in           | -i     |         | path to reduced database |
|--kmer-len     | -l     |    5    | kmer length |

###Makedb options
Use this options when you want to create reduced database

| Option              |  Short | Default | Description       |
| --------------------|--------| --------| ------------------|
|--kmol-high          |        |  5,0    | most common k-kmer over which length, e.g. 5,0 means 5 k-mers over full length; 3,60 means 3 k-mers for every 60 AA |
|--kmol-low          |         |  7,0    | least common k-kmer over which length, e.g. 7,0 means 7 k-mers over full length; 3,60 means 3 k-mers for every 60 AA |
|--seg-window        |        |  12    | windows size of seg masker |
|--seg-low-cut          |        |  2.2    | low cut of seg masker |
|--seg-high-cut         |        |  2.5    | high cut of seg masker |

###Aligner options
Use this options when you want to align query against reduced database

| Option        |  Short | Default      | Description       |
| --------      |--------| --------     | ------------------|
|--query        | -q     |              | path to query file |
|--high-match           | 3     |              | minimum number of common kmer match |
|--low-match           | 2     |              | minimum number of least common kmer match |
|--gapopen      | -g     |     10       | gap open penalty |
|--gapext       | -e     |     1        | gap extend penalty |
|--matrix       | -m     | BLOSUM_62    | score matrix |
|--out          | -o     |              | path to output file |
|--out-format   |        |      bm9     | path to reduced database |
|--evalue       |   -v   |   0.001      | maximum e-value |
|--max-aligns   | -a     |   10         | maximum number of outputs per query |
|--algorithm    |   -A   |      SW        | algorithm used for alignment|




