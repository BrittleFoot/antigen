# Oantigen


## Installation

```bash
python -m pip install -r requirements.txt
```

## Usage

```bash
# get gff with operons-specific features
# the result will be saved in b.gff, b.gff will contain the same features as a.gff,
# but with additional qualifiers for operons and antigen's genes
python merge-operons.py -t O a.gff list_of_operons orfs_coordinates -o b.gff


# You can also search for multiple operons at once
python merge-operons.py -t O -t H -t K a.gff list_of_operons orfs_coordinates -o b.gff

# draw operons
python draw-operons.py b.gff

```


# Deprecated

```bash
python oantigen.py -h
```