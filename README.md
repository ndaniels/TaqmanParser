# TaqmanParser

This program uses a state machine approach to parsing Taqman files, with custom SNP-calling rules specific to Plasmodium falciparum 24-position barcodes.

## Requirements

Python 3.5+ must be installed; this will not run under the deprecated python 2.

## Usage

```bash
python3 parse_taqman.py infile [outfile]
```

or

```bash
./parse_taqman.py infile [outfile]
```

You must specify exactly one input file, and optionally one output file. If the output file is omitted, the program writes to STDOUT.

If you wish to parse an entire directory of input files, for now it is best to do so with a simple bash loop. At a bash prompt (e.g. MacOS terminal). Suppose you are sitting at the bash `$` prompt and have your input files in some_directory:

```bash
$ for file in some_directory/* ; do python3 parse_taqman.py $file $file.out ; done
```

Now some_directory will contain, for each parseable file it originally contained, a file with .out at the end of its name containing the SNP calls. For instance:

```bash
$ ls some_directory
```
```
sen1.txt sen2.txt
sen3.txt sen4.txt
```
```bash
$ for file in some_directory/* ; do python3 parse_taqman.py $file $file.out ; done
$ ls
```
```
sen1.txt sen1.txt.out
sen2.txt sen2.txt.out
sen3.txt sen3.txt.out
sen4.txt sen4.txt.out
```