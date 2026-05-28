# Genomic-Quality-Control-Pipeline
A CLI tool for analyzing FASTA files (or variants like .fa, .fna, .fas, ect). The output consists of a table containing the lenth GC%, N-50 score, and L-50 score of each file. Results can be written to the console (default behavior) or written to a .txt or .md file.

## Usage

```
python3 GQCP.py <flags> <path>
```

|Flag|Description|
|----|-----------|
|`-m`|Writes the output to a .md file|
|`-t`|Writes the output to a .txt file|
