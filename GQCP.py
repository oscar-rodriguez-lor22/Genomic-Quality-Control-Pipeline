import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import isfile
from rich.text import Text
from rich.table import Table
from rich.panel import Panel
from rich.console import Console
from Bio.SeqUtils import gc_fraction

console = Console()
entries = []
md = False
txt = False

## Files/Directories
def FileExistsAndIsFasta(filePath):
    valid_extensions = ['.fasta', '.fa', '.fna', '.fas']
    _, extension = os.path.splitext(filePath)

    if not os.path.isfile(filePath):
        return False
    if not extension.lower() in valid_extensions:
        return False
    else:
        return True

def ReturnFileStats(filePath):
    file_name = os.path.basename(filePath)
    n_occurences = 0
    gc_occurences = 0
    seq_lengths = []

    seq_length = 0
    seq_scaffolds = 0
    seq_n_percentage = 0
    seq_gc_percentage = 0
    seq_n50 = 0
    seq_l50 = 0

    for record in SeqIO.parse(filePath, "fasta"):
        seq_scaffolds += 1
        seq_length += len(record.seq)
        gc_occurences += record.count("G") + record.count("C")
        n_occurences += record.count("N")
        seq_lengths.append(len(record.seq))
        parts = [p for p in str(record.seq).split('N') if p]
    
    seq_lengths.sort(reverse=True)
    s = 0
    for length in seq_lengths:
        s += length
        seq_l50 += 1
        if s >= (seq_length / 2):
            seq_n50 = length
            break    

    seq_n_percentage = "{:.3f}".format((n_occurences/seq_length) * 100)
    seq_gc_percentage = "{:.3f}".format((gc_occurences/seq_length) * 100)
    return seq_length, seq_scaffolds, seq_gc_percentage, seq_n_percentage, seq_n50, seq_l50
def WriteToTxt():
    if (os.path.exists("results.txt")):
        error_message = Text("Error: ", style="Bold red")
        error_message.append("Cant write to results.txt, file already exists")
        console.print(error_message)
        return
    with open("results.txt", 'a') as f:
        f.write("Filename       Length      Scaffold's      %GC     %N      N-50 (Scaffolds)        L-50\n")
        for entry in entries:
            f.write(f"{entry[0]}       {entry[1]}     {entry[2]}       {entry[3]}      {entry[4]}      {entry[5]}      {entry[6]}\n")
    console.print(f"[bold green]Success:[/bold green] Results written to results.txt")
def WriteToMd():
    if os.path.exists("results.md"):
        error_message = Text("Error: ", style="Bold red")
        error_message.append(f"Cant write to results.md, file already exists")
        console.print(error_message)
        return

    with open("results.md", 'a') as f:
        f.write("| Filename | Length | Scaffolds | %GC | %N | N-50 (Scaffolds) | L-50 |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")
        
        for entry in entries:
            # Join entry elements with pipe characters
            row = " | ".join(map(str, entry))
            f.write(f"| {row} |\n")
            
    console.print(f"[bold green]Success:[/bold green] Results written to results.md")
# Output Formating and Displaying (Rich)
def DisplayFileStats(filePath):
    file_name = os.path.basename(filePath)
    # Table Formating
    table = Table(title="FASTA File Contents Summary", show_header=True, header_style="bold magenta")
    table.add_column("Length (bp/nt)", style="green")
    table.add_column("Scaffold Count", style="green")
    table.add_column("%GC", style="green")
    table.add_column("%N", style="green")
    table.add_column("N-50 (Scaffold)", style="green")
    table.add_column("L-50 (Scaffold)", style="green")
    
    ## Table Population and Formating
    try:
        seq_length, seq_scaffolds, seq_gc_percentage, seq_n_percentage, seq_n50, seq_l50 = ReturnFileStats(filePath)
        table.add_row(str(seq_length), str(seq_scaffolds), f"{seq_gc_percentage}%", f"{seq_n_percentage}%", str(seq_n50), str(seq_l50))
        entries.append([file_name, seq_length, seq_scaffolds, f"{seq_gc_percentage}%", f"{seq_n_percentage}%", seq_n50, seq_l50])
                
    except Exception as e:
        console.print(f"[bold red]FATAL PARSING ERROR:[/bold red] {e}")
        return

    console.print("\n")
    console.rule(f"[bold cyan]Analyzing ({file_name})[/bold cyan]")
    console.print(table)
def DisplayStats(path):
    if (FileExistsAndIsFasta(path)):
        DisplayFileStats(path)
    elif (os.path.isdir(path)):
        for filename in os.listdir(path):
            filepath = f"{path}{filename}"
            if (FileExistsAndIsFasta(filepath)):
                DisplayFileStats(filepath)
def DisplayHelpMenu():
    console.rule("[bold cyan]Help Menu[/bold cyan]")
    
    console.print("A tool for quick genomic data summaries.", justify="center", style="italic yellow")
    console.print("\n")

    usage_text = "python FASTA-Stats.py <flags> <path>"
    console.print(
        Panel(
            usage_text, 
            title="[bold green]Usage[/bold green]", 
            border_style="green", 
        ),
        justify="center"
    )
    console.print("\n")

    console.print("[bold magenta]Description: [/bold magenta]")
    console.print("A CLI tool for analyzing .fasta files (or variants like .fa, .fas, etc). The output consists of a table containing the length, GC%, N%, N-50 score, and L-50 score of each file. Results can be written to the console (default behavior) or written to a .txt or .md file using flags. The <path> argument can either be a path to a file or a path to a directory. If the <path> points to a file then the file will be analyzed, if it points to a directory then every fasta file in that directory will be analyzed.")
    console.print("\n")

    console.print("[bold magenta]Flags:[/bold magenta]")
    console.print("[cyan]-h / --help[/cyan]: Display this help menu") 
    console.print("[cyan]-m / --md[/cyan]: Writes the output to a .md file (results.md)")
    console.print("[cyan]-t / --txt[/cyan]: Writes the output to a .txt file (results.txt)")

    console.rule("")

## Argument Checking
if (all(s.startswith('-') for s in sys.argv[1:]) and ("-h" or "--help") not in sys.argv): # Check if all args start with a dash, invalid arguents list if true as a path is required and paths dont start with dashes
    error_message = Text("Error: ", style="bold red")
    error_message.append("No file path provided")
    console.print(error_message)

    console.print("[bold yellow]Usage:[/bold yellow] python FASTA-Stats.py <flags> <path>")
    console.print("[bold yellow]Help:[/bold yellow] for help use -h or --help")
    sys.exit(1)
if (len(sys.argv) < 2):
    error_message = Text("Error: ", style="bold red")
    error_message.append("At least 2 arguments required, fewer given")
    console.print(error_message)

    console.print("[bold yellow]Usage:[/bold yellow] python FASTA-Stats.py <flags> <path>")
    console.print("[bold yellow]Help:[/bold yellow] for help use -h or --help")
    sys.exit(1)
if (len(sys.argv) >= 2): 
    for path_idx in range(len(sys.argv)-1):
        arg = sys.argv[path_idx+1]
        if arg == "-h" or arg == "--help":
            DisplayHelpMenu()
            continue
        if arg == "-m" or arg == "--md":
            md = True
            continue
        if arg == "-t" or arg == "--txt":
            txt = True
        if (arg[0] != "-"): 
            console.print(f"[bold green]Analyzing file/dir:[/bold green] [italic blue] {arg} [/italic blue]")
            DisplayStats(arg)
    if (md):
        WriteToMd()
    if (txt):
        WriteToTxt()
