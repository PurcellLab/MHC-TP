"""Verify that the database is generated correctly."""
import os
import sys
import json
import pandas as pd 
from cli.logger import CONSOLE
import time

def _prase_config_file(config_file):
    if not os.path.exists(config_file):
        # sys.exit(f"Config file {config_file} does not exist")
        CONSOLE.log(f"Config file {config_file} does not exist", style="red")
        CONSOLE.log("Using default configuration", style="blue")
        default_config ={
            "human":{
                "path": "Gibbs_motifs_human",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/hla_data.csv",
                "ref_data": "data/ref_data"
            },

            "mouse" :{
                "path": "Gibbs_motifs_mouse",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/mhc_data.csv",
                "ref_data": "data/ref_data"
            },

            "human_classii": {
                "path": "Gibbs_motifs_human_classII",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/hla_classII_data.csv",
                "ref_data": "data/ref_data"
            },

            "mouse_classii": {
                "path": "Gibbs_motifs_mouse_classII",
                "matrix": "matrices",
                "motif": "motif",
                "allotypes": "allotypes/mhc_classII_data.csv",
                "ref_data": "data/ref_data"
            }
        }
        print(default_config)
        return default_config 
           
    else:
        with open(config_file, "r") as f:
            config = json.load(f)
        return config

def _check_ref_files(config):
    for species in config:
        if not os.path.exists(os.path.join(config[species]["ref_data"], config[species]["path"])):
            # sys.exit(f"Path {os.path.join(config[species]['ref_data'], config[species]['path'])} does not exist")
            # print("Path", config[species]["path"])
            # print("Ref data", config[species]["ref_data"])
            # print("Join", os.path.join(config[species]["ref_data"], config[species]["path"].lstrip('/')))
            CONSOLE.log(f"Path {os.path.join(config[species]['ref_data'],config[species]['path'])} does not exist. Check the path.", style="red")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]["ref_data"],config[species]["path"],config[species]["matrix"])):
            # sys.exit(f"Matrix path {config[species]['path'] + config[species]['matrix']} does not exist")
            CONSOLE.log(f"Matrix path {os.path.join(config[species]['ref_data'],config[species]['path'],config[species]['matrix'])} does not exist")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]["motif"])):
            # sys.exit(f"Motif path {config[species]['path'] + config[species]['motif']} does not exist")
            CONSOLE.log(f"Motif path {os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]['motif'])} does not exist")
            sys.exit(1)
        if not os.path.exists(os.path.join(config[species]["ref_data"], config[species]["allotypes"])):
            # sys.exit(f"Allotype path {config[species]['path'] + config[species]['allotypes']} does not exist")
            CONSOLE.log(f"Allotype path {os.path.join(config[species]['ref_data'], config[species]['allotypes'])} does not exist")
            sys.exit(1)

def _HLA_liist(config):
    for species in config:
        allotype_file = os.path.join(config[species]["ref_data"], config[species]["allotypes"])
        allotypes = pd.read_csv(allotype_file)
        if species == "human":
            allotypes.rename(columns={ "formatted_HLA":"formatted_allotypes","HLA":"allotypes"}, inplace=True)
            allotypes['motif'] = "HLA_" + allotypes['formatted_allotypes'] + ".png"
            allotypes['species'] = species
        elif species == "mouse":
            allotypes.rename(columns={"formatted_MHC":"formatted_allotypes","MHC":"allotypes"}, inplace=True)
            allotypes['motif'] = allotypes['formatted_allotypes'] + ".png"
            allotypes['species'] = species
        elif species in ("human_classii", "mouse_classii"):
            # Class II allotype files: columns formatted_HLA/HLA (human) or formatted_MHC/MHC (mouse)
            # Motif filenames use the formatted allotype name directly (no HLA_ prefix)
            if "formatted_HLA" in allotypes.columns:
                allotypes.rename(columns={"formatted_HLA": "formatted_allotypes", "HLA": "allotypes"}, inplace=True)
            elif "formatted_MHC" in allotypes.columns:
                allotypes.rename(columns={"formatted_MHC": "formatted_allotypes", "MHC": "allotypes"}, inplace=True)
            allotypes['motif'] = allotypes['formatted_allotypes'] + ".png"
            allotypes['species'] = species
        else:
            # sys.exit(f"Species {species} not supported")
            CONSOLE.log(f"Species {species} not supported")
            sys.exit(1)
        allotypes['motif_path'] = ""
        allotypes['matrices_path'] = ""
        for motifs in allotypes['motif']:
            CONSOLE.log(f"[yellow]{species}[/yellow] {config[species]['path']}/motif/{motifs}", style="blue")
            if not os.path.exists(f"{config[species]['ref_data']}/{config[species]['path']}/motif/{motifs}"):
                # sys.exit(f"Motif {motifs} does not exist")
                CONSOLE.log(f"Motif {motifs} does not exist")
                sys.exit(1)
            elif not os.path.exists(f"{os.path.join(config[species]['ref_data'],config[species]['path'] ,config[species]['matrix'])}/{motifs.replace('.png', '.txt')}"):
                # sys.exit(f"Matrix {motifs.replace('.png', '.txt')} does not exist")
                CONSOLE.log(f"Matrix {motifs.replace('.png', '.txt')} does not exist")
                sys.exit(1)
            else:
                allotypes.loc[allotypes['motif'] == motifs, 'motif_path'] = f"{config[species]['path']}/motif/{motifs}"
                allotypes.loc[allotypes['motif'] == motifs, 'matrices_path'] = f"{config[species]['path']}/matrices/{motifs.replace('.png', '.txt')}"
                # CONSOLE.log(f"Motif {motifs} exists")
        # allotype_file_db = allotype_file.replace('.csv', '.db')        
        allotypes.to_csv(f"{config[species]['ref_data']}/{species}.db", index=False)
        CONSOLE.log(f"Database {species}.db created successfully ({config[species]['ref_data']}/{species}.db)", style="green")
        
def Database_gen(config_file):
    # CONSOLE.print("Database generation started", style="green")
    with CONSOLE.status("[bold green] Database generation started [bold green]") as status:
        config = _prase_config_file(config_file)
        # CONSOLE.log("Configuration file parsed and validated successfully.", style="green")
        status.update(
                        status=f"[bold green] Configuration file parsed and validated successfully. [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        
        _check_ref_files(config)
        time.sleep(3)
        
        # CONSOLE.log("Reference files checked and verified.", style="green")
        status.update(
                        status=f"[bold green] Reference files checked and verified [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        _HLA_liist(config)
        time.sleep(3)
        # CONSOLE.log("HLA list generated and database output created successfully.", style="green")
        status.update(
                        status=f"[bold green] HLA list generated and database output created successfully [/bold green]",
                        spinner="squish",
                        spinner_style="yellow",
                            )
        time.sleep(3)
        # CONSOLE.print("Database generation completed", style="green")
        status.update(
                        status=f"[bold green] Database generation completed [/bold green]",
                        spinner="squish",
                        spinner_style="green",
                            )
        sys.exit(0)
    
    
def generate_HLA_freq_database(out_put_ditectory=None):
    """Generate HLA frequency database."""
    if out_put_ditectory is None:
        CONSOLE.log("Output directory not specified. Please provide a valid output directory.", style="red")
        sys.exit(1)
    try:
        from .HLAfreq import get_list, makeURL, getAFdata, url_encode_name
    except ImportError:
        CONSOLE.log("HLAfreq function Not avilable please check....", style="red")
        sys.exit(1)

    # from .HLAfreq import HLAfreq_pymc as HLAhdi
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    #country_list = ["Australia", "Thailand", "United States of America", "United Kingdom", "Germany", "France", "Italy", "Spain", "Netherlands", "Sweden", "Norway", "Finland", "Denmark", "Belgium", "Switzerland", "Austria", "Poland", "Czech Republic", "Hungary", "Portugal"]
    # country_list = ["Australia"]
    CONSOLE.log("Processing HLA frequency data for countries...", style="blue")
    country_list = get_list("country")
    # print(country_list)
    # for country in country_list:
    #     print(f"Processing {url_encode_name(country)}...")
    CONSOLE.log(f"Total countries to process: {len(country_list)}", style="yellow")
    for country in country_list:
        CONSOLE.log(f"Processing {url_encode_name(country)}...", style="blue")
        base_url = makeURL(country)
        try:
            CONSOLE.log(f"Base URL for {country}: {base_url}", style="green")
            aftab = getAFdata(base_url)
            aftab.to_csv(f"{out_put_ditectory}/{country}_raw.csv", index=False)
            CONSOLE.log(f"Raw data for {country} saved to example/{country}_raw.csv", style="green")
        except Exception as e:
            CONSOLE.log(f"Error processing {country}: {e}", style="red")
            continue
    
    
def add_allefreq_to_db(db_pat, allefreq_path,allefreq_output,loci="A",freq=0.01):
    try:
        from .HLAfreq import only_complete, combineAF, decrease_resolution,check_resolution
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        CONSOLE.log("HLAfreq function Not available please check....", style="red")
        sys.exit(1)
    """Add allele frequency data to the database."""
    if not os.path.exists(db_pat):
        CONSOLE.log(f"Database file {db_pat} does not exist.", style="red")
        sys.exit(1)
    
    if not os.path.exists(allefreq_path):
        CONSOLE.log(f"Allele frequency file {allefreq_path} does not exist.", style="red")
        sys.exit(1)
    
    # Load the database and allele frequency data
    db = pd.read_csv(db_pat)
    # print(db)
    # allefreq = pd.read_csv(allefreq_path) 
    country_files = [f for f in os.listdir(allefreq_path) if f.endswith('_raw.csv')]
    if not country_files:
        CONSOLE.log(f"No allele frequency files found in {allefreq_path}.", style="red")
        sys.exit(1) 
    country_count = 0  
    cafs = []
    for country in country_files:
        country_count += 1
        country_path = os.path.join(allefreq_path, country)
        allefreq = pd.read_csv(country_path)
        CONSOLE.log(f"Processing allele frequency data for {country}... {country_count}/{len(country_files)}", style="blue")
        # print(allefreq.allele)
            # print(allefreq)
        # Drop any incomplete studies
        aftab = only_complete(allefreq)
        check = check_resolution(aftab)
        if check:
            aftab = decrease_resolution(aftab, 2)
            afloc = aftab[aftab.loci==loci]
            if afloc.empty:
                CONSOLE.log(f"No A locus data found for {country}. Skipping...", style="yellow")
                continue
            caf = combineAF(afloc)
            caf['country'] = country.split('_')[0]  # Extract country name from filename
            cafs.append(caf)
        else:
            continue
        # Ensure all alleles have the same resolution
        # aftab = decrease_resolution(aftab, 2)
    cafs = pd.concat(cafs, ignore_index=True)
    cafs.to_csv(f"{allefreq_output}/01HLA_freq_by_country_all_HLA-{loci}.csv", index=False)
    international = combineAF(cafs, datasetID='country')
    mask = international.allele_freq > freq
    international[mask].plot.barh('allele', 'allele_freq')
    plt.savefig(f"{allefreq_output}/01HLA_freq_by_country_Top20_HLA-{loci}_freq>{freq}.png")
    # plt.show()
# # # Database_gen("config.json")
import os

import os
import pandas as pd

def _formate_mouse_DB(ref_db_path, db_path, out_db_path=None):
    # Load DB
    db = pd.read_csv(db_path)

    # Walk through motif + matrices folders and rename files
    for folder in os.listdir(ref_db_path):
        if folder in ["motif", "matrices"]:
            folder_path = os.path.join(ref_db_path, folder)
            for fname in os.listdir(folder_path):
                if fname.startswith("MHC_"):
                    old_path = os.path.join(folder_path, fname)
                    new_name = fname.replace("MHC_", "", 1)  # replace prefix
                    new_path = os.path.join(folder_path, new_name)
                    os.rename(old_path, new_path)
                    print(f"Renamed: {old_path} → {new_path}")

                    # Update DB entries that reference this filename
                    db = db.replace(fname, new_name)

    # Save updated DB
    if out_db_path is None:
        out_db_path = db_path  # overwrite input DB
    db.to_csv(out_db_path, index=False)
    print(f"Database updated and saved to {out_db_path}")

# if __name__ == "__main__":
    # config_file = "config.json"
    # config = _prase_config_file(config_file)
    # _check_ref_files(config)
    # path="/home/sson0030/xy86_scratch2/SANJAY/MHC-TP/data/ref_data/Gibbs_motifs_mouse"
    # path="data/ref_data/Gibbs_motifs_mouse"
    # db_path="data/ref_data/mouse.db"
    # _formate_mouse_DB(path,db_path)
    # hla_list = _HLA_liist(config)
    # print(hla_list)
    # print(config)
    # CONSOLE.log("Config file parsed successfully")
#     sys.exit(0)
    
#     ## HLA frequency database generation
 

#         # Combine studies within country
#         caf = combineAF(aftab)
#         # Add country name to dataset, this is used as `datasetID` going forward
#         caf['country'] = country
#         cafs.append(caf)

#     cafs = pd.concat(cafs, ignore_index=True)
#     international = combineAF(cafs, datasetID='country')
#     print(international)
#     db_pat = "/home/sson0030/xy86_scratch2/SANJAY/MHC-TP/data/ref_data/human.db"
#     allefreq_path = "/home/sson0030/xy86_scratch2/SANJAY/MHC-TP/data/ref_data/HLAfreq/byCountry/"
#     allefreq_output = "/home/sson0030/xy86_scratch2/SANJAY/MHC-TP/data/ref_data/HLAfreq/output"
#     add_allefreq_to_db(db_pat, allefreq_path,allefreq_output,loci="C",freq=0.01)