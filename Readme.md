## Building

### Build requirements

One needs `GNU make`, `Cmake 3+`, `Boost` libraries. If you don't know if you have them or
not just type `make packages` and it will prompt you for install. If you do not have permissions
to install packages or want to put them in other folder (this requires installed `boost` in the default include directory)
you can install it and than uncomment first four lines in `deps/CMakeLists.txt`. In lines 3 and 4 you need
to set path to installed boost directories.


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

| Option  |  Short | Default | Description       |
| --------|--------| --------| ------------------|
|--threads| -p     |   8     | number of threads |
|--db     | -d     |         | path to reference database |

###Makedb options
Use this options when you want to create reduced database

| Option  |  Short | Default | Description       |
| --------|--------| --------| ------------------|
|--read   | -r     |         | path to reduced database |

###Aligner options
Use this options when you want to align query against reduced database

| Option        |  Short | Default      | Description       |
| --------      |--------| --------     | ------------------|
|--in           | -i     |              | path to reduced database |
|--query        | -q     |              | path to query file |
|--gapopen      | -g     |     10       | gap open penalty |
|--gapext       | -e     |     1        | gap extend penalty |
|--matrix       | -m     | BLOSUM_62    | score matrix |
|--out          | -o     |              | path to output file |
|--out-format   |        |      bm9     | path to reduced database |
|--evalue       |   -v   |   0.001      | maximum e-value |
|--max-aligns   | -a     |   10         | maximum number of outputs per query |
|--algorithm    |   -A   |      SW        | algorithm used for alignment|




