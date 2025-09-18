# HLA-PepClust  
![CI Status Ubuntu](https://github.com/Sanpme66/HLA-PepClust/actions/workflows/python-package.yml/badge.svg)  
![MultiOS](https://github.com/Sanpme66/HLA-PepClust/actions/workflows/matrix.yml/badge.svg)  
![test-hla-pepclust](https://github.com/Sanpme66/HLA-PepClust/actions/workflows/test-hla-pepclust.yml/badge.svg)  

# Demo Results
[Click here to view the results](https://sanpme66.github.io/HLA-PepClust/)

## Introduction  

HLA-PepClust is a `CLI` tool designed for clustering peptide sequences based on their HLA binding motifs.  

## Prerequisites  

Ensure your system meets the following requirements:  
- Python **3.9 or higher**  
- `pip` (Python package manager)  

## Download or Clone the Repository  

```bash
git clone https://github.com/Sanpme66/HLA-PepClust.git
cd HLA-PepClust/
```

## Setting Up the Python Environment  

1. **Create a virtual environment**  
    ```bash
    python3 -m venv hlapepclust-env
    ```

2. **Activate the virtual environment**  
    - **macOS / Linux**:  
      ```bash
      source hlapepclust-env/bin/activate
      ```
    - **Windows**:  
      ```bash
      .\hlapepclust-env\Scripts\activate
      ```

## Installing Dependencies  

1. **Navigate to the project directory** (if not already in it)  
    ```bash
    cd HLA-PepClust/
    ```

2. **Install the package and dependencies**  
    ```bash
    pip install -e .
    ```

## Running HLA-PepClust  

### Display Help Message  

To see available options and usage details:  
```bash
clust-search -h      
```

Example of running the help `clust-search -h` command:  

![Example Output](assets/img/clust-search-help.png)  

### Command Structure  

```bash
clust-search <input_data_path> <reference_data_path> \
  --hla_types <hla_types> \
  --n_clusters <number_of_clusters> \
  --output <output_path> \
  --species <Human or Mouse> \
  --threshold <similarity_threshold (default: 0.5)> \
  --log \
  --processes <number_of_threads>
```

### Example Usage  

```bash
clust-search data/D90_HLA_3844874 data/ref_data/Gibbs_motifs_human/output_matrices_human \
  --hla_types A0201,A0101,B1302,B3503,C0401 \   # Specify list of HLA alleles to search
  --n_clusters 6 \                               # Restrict analysis to 6 Gibbs clusters
  --species human \                              # Species to evaluate [human, mouse]
  --output My_results_directory \                # Directory where results will be saved
  --processes 4 \                                # Number of parallel processes to use
  --threshold 0.7 \                              # Correlation threshold for motif matching
  --topNHits 3                                   # Report top-N HLA matches for each Gibbs motif
```

## Command Line Arguments  
| Argument         | Type    | Description                                                                 | Default                  |
|------------------|---------|-----------------------------------------------------------------------------|--------------------------|
| `gibbs_folder`   | `str`   | Path to test folder containing matrices.                                    | *Required*               |
| `reference_folder` | `str` | Path to reference folder containing matrices.                               | *Required*               |
| `-o, --output`   | `str`   | Path to output folder.                                                      | `"output"`               |
| `-hla, --hla_types` | `list` | List of HLA types to search.                                               | *All*                    |
| `-p, --processes` | `int`  | Number of parallel processes to use.                                        | `4`                      |
| `-n, --n_clusters` | `str` | Number of clusters to search for.                                           | `"all"`                  |
| `-t, --threshold` | `float` | Motif similarity threshold.                                                 | `0.70`                   |
| `-s, --species`   | `str`   | Species to search [Human, Mouse].                                          | `"human"`                |
| `-db, --database` | `str`   | Generate a motif database from a configuration file.                        | `"data/config.json"`     |
| `-st, --Searchtype` | `str` | Type of search to perform [Numba, IO].                                      | `"Numba"`                |
| `-k, --best_KL`   | `bool`  | Find the best KL divergence only.                                          | `False`                  |
| `--topNHits`      | `int`   | Number of top hits to retain per Gibbs matrix.                             | `3`                      |
| `-l, --log`       | `bool`  | Enable logging.                                                            | `False`                  |
| `-im, --immunolyser` | `bool` | Enable immunolyser output.                                               | `False`                  |
| `-npDB, --NumbaDB` | `str`  | Path to the Array database folder.                                          | `"data/ref_data/human_db"` |
| `-c, --credits`   | `bool`  | Show credits for the motif database pipeline.                              | `False`                  |
| `-v, --version`   | `bool`  | Show the version of the pipeline.                                          | `False`                  |


## Example Output  

Example of running the `clust-search` command:  

![Example Output](assets/img/search-results.png)  

## Deactivating the Virtual Environment  

After finishing, deactivate the virtual environment with:  
```bash
deactivate
```


## Running HLA-PepClust in Google Colab  

You can try out **HLA-PepClust** directly in **Google Colab** without installing anything on your local system!  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Sanpme66/HLA-PepClust/blob/main/example/HLA_PepClust_testing.ipynb)  

### Steps to Run in Colab  

1. Click the **"Open In Colab"** button above.  
2. Once the notebook opens in Colab, go to **Runtime → Run All** to execute all cells.  
3. Modify input parameters if needed and run the pipeline in Colab.  

Example of `input` folder path:  

![Example Output](assets/img/google-colab.png)  

### Citation

If you use **MHC-TP** in your research, please cite:

**Immunolyser 2.0: an advanced computational pipeline for comprehensive analysis of immunopeptidomic data**
Prithvi Raj Munday¹,†, Sanjay S.G. Krishna¹,†, Joshua Fehring¹, Nathan P. Croft¹, Anthony W. Purcell¹, Chen Li¹,², and Asolina Braun¹
¹Department of Biochemistry and Molecular Biology and Biomedicine Discovery Institute, Monash University, Clayton, VIC, 3800, Australia
²Department of Medicine, School of Clinical Sciences at Monash Health, Monash University, Clayton, VIC 3168, Australia

*Computational and Structural Biotechnology Journal*

More detailed instructions coming soon... 
