## Building

### Build requirements

One needs `GNU make`, `Cmake 3+`, `Boost` libraries. If you don't know if you have them or
not just type `make packages` and it will prompt you for install.

Built it with Boost 1.60 but it should work with older versions if interface
didn't change.

### How to build?

Build by writing `make` in your command line. If everything went without errors
executables should be present in the directory called `build`. There is executable file
tachyon.

## Running
There are three options in which you can run your application. Use `makedb` if you want to set up reduced database. If
you run

`$ tachyon makedb -d nr.fa -r reduced_nr`

it will create reduced database (required by algorithm) in the current directory.

You can run alignment task with the following commands:

`$ tachyon blastp -d nr.fa -i reduced_nr -q queries.fa -o results`
if you want to align protein sequences against nr database, or

`$ tachyon blastx -d nr.fa -i reduced_nr -q queries.fa -o results`

if you want to align DNA sequences against nr database.

##Commands

This commands determine the mode in which you want to run this app

| Command       | Description   |
| ------------- |:-------------:|
| makedb        | Create reduced database from FASTA reference file.                      |
| blastp        | Align protein query sequences against a protein reference database      |
| blastx        | Align DNA query against a protein reference database                    |
