# BIOSEQ - Bio-Sequence Decode Utilities in Perl

## Introduction

Perl scripts in this project, collectively called **bioseq utilities**, are designed to decode a bio-sequence such as nucleotides in RNA or DNA, or amino acids representing a protein.

In most cases, a bio-sequence is represented in text format for ease of sharing.
UNIX (or Linux) text tools are conveniently used to handle these sequence data.
The bioseq utilities are meant to augment Linux command-line tools to assist bio-sequence-specific decoding.

Currently the following three utilities exist (in written order):

| Utility | Name | Description |
|--|--|--|
| motifgrep | pattern matcher | searches a bio-sequence for occurences that satisfy a *motif* pattern (described later) |
| nucleoamino | nucleotide to amino acid converter | converts a nucleotide sequence to an amino acid sequence |
| modshunt | sequence diff | compares two bio-sequences and print differences |

You can run the scripts in any environment where you can use **Perl** version 5. 
I developed and tested the scripts with Perl v5.28.1 under Debian GNU/Linux 10 (buster) running under Windows Subsystem for Linux 2.

Before running the scripts, set the shell environment variable `BIOSEQ` to the directory which contains the files in this project:

```
$ sudo cp bioseq/* /usr/local/bin
$ export BIOSEQ=/usr/local/bin
```

## FASTA format

Bioseq utilities support a sequence data format called FASTA format.
A FASTA file begins with a header line prefixed with '>' character followed by one or more lines of sequence data.

Example:

```
>NC_045512.2 |Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
 ...
TTTAGTAGTGCTATCCCCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAA
```

Every bioseq utility strips the header line and newlines to form a single-line sequence data before processing.
A raw data without a header may be safely passed to it.

## motifgrep

MotifGrep searches an input sequence for a motif pattern.

A **motif** (or sequence motif) pattern is a biological way of representing a pattern in a bio-sequence.
It corresponds to a regular expresion in computing.

>https://en.wikipedia.org/wiki/Sequence_motif

### Usage

- **motifgrep** [-n|-nu] *motif* [*path*]

A *path* is a filename of a FASTA file to search.
If you omit *path*, MotifGrep reads from the standard input.

A *motif* is a pattern of a partial sequence of nucleotides or amino acids, such as 'AATTAT', 'AAUUAU' or 'NY'.

The `-n` or `-nu` option specifes conversion of an amino acid motif you specify on the command line to a corresponding nucleotide motif.
For example, 'NY' is converted to 'AA[TC]TA[TC]' and 'NxY' to 'AA[TC]xxxTA[TC]' where 'x' means any pattern.

If you specify `-nu` option instead of `-n` option, the converted pattern will use `U` (Uracil) instead of `T` (Thymine).
Thus 'NY' will be converted to 'AA[UC]UA[UC]'.

MotifGrep prints a match on each line with its position.
The position is one-based; the first character in a sequence is 1.

## Motif Pattern

A motif can contain:

- choice of characters within `[..]`, such as '[TG]' meaning T or G,
- exclusion with `{..}`, such as '{TG}' meaning any character except T or G.
- `x` for any character,

MotifGrep supports the PROSITE extensions:

- hyphen (-) delimiter for separating consecutive pattterns,
- range with *e*(*n*) and *e*(*n*,*m*), such as 'T(2)' for 'TT' and 'T(2,3)' for 'TT' or 'TTT'.

>Limitation: When you specify an amino acid pattern to be converted to nucleotide pattern (ie., `-n` or `-nu` option is present), the only available motif syntax is 'x'; no other ones can be used. For example, with `-n` or `-nu` option present, 'NY' or 'NxY' are OK but 'N[YA]' or 'N{YA}' are not.

### Example

```
$ cat test.fasta
>test data
ATATATCAGAG
AGCAGAGACC
$ motifgrep ATC test.fasta
 5:ATC
$ motifgrep A[TG]C test.fasta
 5:ATC
12:AGC
$ motifgrep AxC test.fasta
 5:ATC
12:AGC
19:ACC
$ motifgrep A{C}C test.fasta
 5:ATC
12:AGC
$ motifgrep 'AC(2)' test.fasta
19:ACC
```

Note the last example where quotes are added to the pattern. Otherwise a Linux shell tries to handle the parentheses in a special way.

## nucleoamino

NucleoAmino converts an input nucleotide sequence to an amino acid sequence.
A triplet of nucleotide codes such as "AAU" and "UAU" corresponds to 
a type of amino acid that comprises a protein, such as Asparagine and
Tyrosine, respectively, whose codes are N and Y.

NucleoAmino uses two mapping tables, `nuc2ami.txt` and `amino.txt`, to convert
a triplet to a short amino acid name and then to a single-letter code.
These are editable, tab-separated text files. 

| File | Table name | Description |
|--|--|--|
| nuc2ami.txt | nucleotide to amino acid | mapping table from nucleotide triplet (such as AAU) to amino acid short name (Asn)
| amino.txt | amino acid | lists amino acid names: code (such as N), short name (Asn) and long name (Asparagine) |

For example, a nucleotide triplet of `AAU` or `AAC` maps to amino acid short name of 'Asn' in nuc2ami.txt and 'Asn' leads to amino acid code of `N` by looking at amino.txt.

>https://www.bioinformatics.org/sms/iupac.html  
https://www.genome.gov/genetics-glossary/Genetic-Code 

You don't have to care about whether your source nucleotide sequence contains Uracil (U) or Thymine (T).
NecleoAmino automatically converts T's to U's before conversion.

### Usage

- **nucleoamino** [-**f**] [**-s***span*] [*fastafile*]

A *path* is a filename of a FASTA file to convert.
If you omit *path*, NucleoAmino reads from the standard input.

The `-f` option specifies the output should be in FASTA format.
Otherwise the sequence data is output in raw format (no header, no newlines).

Using the -**s***span* option, you can specify where to start and end conversion instead of converting an entire sequence.

- A *span* of *N-M* form takes characters from N'th to M'th positions, while
- *span* of *N+M* form takes M characters from N'th position.

The character position is one-based.
Take an example of string "hello".
Span '1-5' means the entire string (same as not specifying the `-s` option).
Both '2-4' and '2+3' extracts "ell".

If you omit N (that is '-M' or '+M'), the first position (1) is assumed.
If you omit M (just say 'N'), conversion procedes to the end of the sequence.

### Example

```
$ cat test.fasta
>test data
ATATATCAGAG
AGCAGAGACC
$ nucleoamino test.fasta
IYQRAES 
$ nucleoamino -s2 test.fasta
YIREQR
```

The first example starts to read the first character (ATA...) while the second example starts from the second character (TAT...).

## modshunt

ModsHunt searches two input sequences and shows differences in them.
ModsHunt produces comparison result with greatest identity and smallest difference.

### Usage

- **modshunt** [**-f***fmt*] *reference_fasta_file* *target_fasta_file*

You are required to pass two FASTA files.

- reference sequence (or source or original sequence)
- comparison target sequence

The **-f***fmt* option determines the output format:

| Option | Format | Description | Remark |
|--|--|--|--|
| -fs | Single view | Source sequence shown intermixed with differences | Default |
| -fd | Double view | Source and target sequences shown in parallel | |
| -ft | Tab-seperated ranges | list of paired regions (identical or different) in tab-separated format | |
| -fx | Tab-separated ranges | Same as `-ft` but easier to read for human | Debug | 
| -fj | Ranges in Json | list of paired regions in Json format | |

The single view format is used if the format options is not passed.

The next section describes ModsHunt output formats.

### Output Formats

Currently ModsHunt supports four types of output formats: two for human viewing and two for input to another program. 
You can switch the output format with the format option (-f*fmt*) as described above.

The `-fs` or `-fd` option specifies a human-readable output format: **single view** or **double view** respectively.
In the single view, the source sequence is shown intermixed with differences (insertion, deletion and replacement).
The double view shows the source and target sequences in parallel.

The `-ft` or `-fj` option specifies a **range list** formats: tab-separated text or Json respectively.
These formats may sometimes be useful for human but the main purpose is for use by another program used to process the comparison result, such as a graphical viewer. 

The raw output from the comparison function (split3) is an ordered list of *ranges*.
A range is a pair of corresponding portions of the source and target sequences called *regions*.
Each pair of source and target regions is classified as identical, inserted, deleted or replaced range.

| Range Type | Abbreviations | Description |
|--|--|--|
| Identity | sync or '=' | The corresponding regions in the source and target are exactly the same |
| Insertion | ins or '+' | A new region was inserted into the target |
| Deletion | del or '-' | An old region was deleted in the source |
| Replacement | diff or '^' or '~' | An old region in the source was replaced with a new region in the target |

The range list format allows you to get this raw data in text format (tab-separated text or Json).

Detailed description of the output formats will follow each example below.

### Example

The first example shows the contents of the source and target files and the single-view format output.

```
$ cat reference.fasta
>reference sequence
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCddddddAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT

$ cat variant.fasta
>variant sequence
ATTAAArrrrrrTACCTTCCCAGGTAACAAACCAACCAACTTTCGATCAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGiiiiiiTGCACT

$ modshunt reference.fasta variant.fasta
ATTAAA{7^GGTTTA/rrrrrr}TACCTTCCCAGGTAACAAACCAACCAACTTTCGATC{49-dddddd}AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAG{115+iiiiii}TGCACT
```

In the single-view format, the differences are enclosed in braces and shown intermixed with the source sequence.

| Difference | Syntax | Example |
|--|--|--|
| Replacement | { *source position* ^ *source region* / *target region* } | {7^GGTTTA/rrrrrr} |
| Deletion | { *source position* - *source region* } | {49-dddddd} |
| Insertion | { *source position* + *target region* } | {115+iiiiii} |

The next example shows the double-view format:

```
$ modshunt -fd reference.fasta variant.fasta
1     7     13                                  49    55
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCddddddAGATCTGTTCTCTAAA
......rrrrrr....................................      ................
1     7     13                                  49    49

                                            109   115
CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAG      TGCACT
............................................iiiiii......
                                            109   115
```

In the double-view format, the source (upper) and target (lower) sequences are shown in parallel.
The position of a source or target region is shown above or below the sequence data respectively. 
Identical target regions are shown in dots.

The next example shows the tab-separated range list format:

```
$ modshunt -ft reference.fasta variant.fasta
=       1+6     1+6     "ATTAAA"        "ATTAAA"
^       7+6     7+6     "GGTTTA"        "rrrrrr"
=       13+36   13+36   "TACCTTCCCAGGTAACAAACCAACCAACTTTCGATC"  "TACCTTCCCAGGTAACAAACCAACCAACTTTCGATC"
-       49+6    49+0    "dddddd"        ""
=       55+60   49+60   "AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAG"     "AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAG"
+       115+0   109+6   ""      "iiiiii"
=       115+6   115+6   "TGCACT"        "TGCACT"
```

In the tab-separated range list format, each line corresponds to a range that represents an identity or difference.

The elements in a range line are separated by a tab:

| Column | Description | Example |
|--|--|--|
| *symbol* | Range type: `=` for identical region, `+` for insertion, `-` for deletion or `^` for replacement | ^ |
| *pos*+*len* | position and length of the source region | 7+6 |
| *pos*+*len* | position and length of the target region | 7+6 |
| "*string*" | sequence of the source region | "GGTTTA" |
| "*string*" | sequence of the target regiion | "rrrrrr" |

Format option of `-fx` works the same way as the tab-separated option (-ft) excpet that a range type is shown in a simple name and sequences are displayed shorter.

```
$ modshunt -fx reference.fasta variant.fasta
sync    1+6     1+6     "ATTAAA"        "ATTAAA"
diff    7+6     7+6     "GGTTTA"        "rrrrrr"
sync    13+36   13+36   "TACCTTC..."    "TACCTTC..."
del     49+6    49+0    "dddddd"        ""
sync    55+60   49+60   "AGATCTG..."    "AGATCTG..."
ins     115+0   109+6   ""              "iiiiii"
sync    115+6   115+6   "TGCACT"        "TGCACT"
```

The next example shows the Json range list format:

```
$ modshunt -fj reference.fasta variant.fasta
[
 {
  type: "=",
  src: {
   str: "ATTAAA",
   pos: 1,
   len: 6
  }
  dst: {
   str: "ATTAAA",
   pos: 1,
   len: 6
  }
 }
  ...
 {
  type: "=",
  src: {
   str: "TGCACT",
   pos: 115,
   len: 6
  }
  dst: {
   str: "TGCACT",
   pos: 115,
   len: 6
  }
 }
]
```

In the Json range list format, a Json array contains one or more range objects.
A range object contains the following three attributes:

| Attribute | Description | Example |
|--|--|--|
| type | Region type: `=` for identical region, `+` for insertion, `-` for deletion or `^` for replacement | '^' |
| src | Source region object | See below |
| dst | Target region object | ditto |

The `src` and `dst` region objects contain the following three attributes:

| Attribute | Description | Example |
|--|--|--|
| pos | Position | 7 |
| len | Length  | 6 |
| str | Sequence data of the region | "GGTTTA" |

## Additional Information

A brief introduction page describes what and how you can do with these tools:

- [English](https://kobu.com/bioseq/index-en.html)
- [Japansese](https://kobu.com/bioseq/index.html)

Documentation of the tools:

- English - this page
- [Japansese](https://kobu.com/bioseq/github-README-ja.html) - now working

## License

Copyright (c) 2021 Kobu.Com. Some Rights Reserved.
Distributed under GNU General Public License 3.0.
Contact Kobu.Com for other types of licenses.

## History

2021-feb-09 first edition  
2021-feb-19 modshunt: exec time reduction of 40%  
2021-feb-19 modshunt: format option (-fx) for debug     
2021-feb-20 necuoamino: fasta option (-f) 
