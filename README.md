# Genomic-Quality-Control-Pipeline
A CLI tool for analyzing FASTA files (or variants like .fa, .fna, .fas, ect). The output consists of a table containing the lenth GC%, N-50 score, and L-50 score of each file. Results can be written to the console (default behavior) or written to a .txt or .md file.

## Usage

**Example**
```
python3 GQCP.py <flags> <path>
```

**Flags**
|Flag|Description|
|----|-----------|
|`-m`|Writes the output to a .md file|
|`-t`|Writes the output to a .txt file|

**Path**
The path does not need a flag (statement prefaced by a `-` char), you just input the path right after the flags if you use flags at all. The path can be either absolute or relativae. Path can also point to a single file or a directory, if it points to a directory then all FASTA files in the directory will be analyzed.
