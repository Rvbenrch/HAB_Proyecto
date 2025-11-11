#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================
WORKFLOW OF THE SCRIPT

1. STRING NETWORK FILTERING AND CONVERSION:
   - Reads the STRING network file (protein links) with combined scores.
   - Filters only interactions with combined_score >= 800.
   - Converts ENSP IDs to Entrez Gene IDs using MyGene.

   NOTE: The STRING network file (9606.protein.links.v12.0.txt) is **not included** due to large size.
   You need to download it manually from the STRING database:
   https://string-db.org/cgi/download?species_text=Homo+sapiens
   Once downloaded, provide its path via the --string_input argument.

2. SEED GENES READING:
   - Reads seed genes from the file.
   - Converts HUGO symbols to Entrez Gene IDs.

3. GRAPH CONSTRUCTION:
   - Builds a NetworkX graph with nodes representing genes and edges representing STRING interactions.
   - Edge weights correspond to the combined_score in the filtered network.

4. CHECK SEED GENE PRESENCE:
   - Checks which seed genes are present in the network.

5. ALGORITHMS EXECUTION:
   a) DIAMOnD:
      - Adds nodes to the network based on proximity to seed genes using hypergeometric test.
   b) GUILD NetScore:
      - Calculates a score for all nodes propagating information from seed genes.

6. SELECT TOP NODES:
   - Selects top nodes from DIAMOnD and GUILD for functional analysis.

7. FUNCTIONAL ENRICHMENT:
   - Converts Entrez IDs to HUGO symbols.
   - Uses Enrichr API to perform ORA on top nodes.
   - Generates bar plots for top terms.

8. OUTPUT:
   - All output files and plots are saved in the output folder.

Run the script from terminal using:

    python script.py --string_input <path_to_STRING_network_file> --seed_genes <path_to_seed_genes_file> --output_dir <output_folder> --top_n <number_of_top_genes>
============================================================
"""

import os
import argparse
import pandas as pd
import mygene
import networkx as nx
from scipy.stats import hypergeom
from tqdm import tqdm
import requests
import matplotlib.pyplot as plt
import seaborn as sns

# ==============================
# CLI ARGUMENTS (Command-line interface)
# ==============================
parser = argparse.ArgumentParser(description="STRING network + DIAMOnD + GUILD + Functional Analysis")
parser.add_argument("--string_input", type=str, required=True, help="Path to raw STRING network file (protein links)")
parser.add_argument("--seed_genes", type=str, required=True, help="Path to seed genes file (HUGO symbols)")
parser.add_argument("--output_dir", type=str, default="../results", help="Folder to save results")
parser.add_argument("--top_n", type=int, default=10, help="Number of top nodes to select for functional analysis")
args = parser.parse_args()

string_raw_file = args.string_input
genes_seed_file = args.seed_genes
results_dir = args.output_dir
TOP_N = args.top_n
os.makedirs(results_dir, exist_ok=True)

string_entrez_file = os.path.join(results_dir, 'string_filtered_entrez.txt')

# ==============================
# STEP 1: FILTER STRING NETWORK + CONVERT ENSP → ENTREZ
# ==============================
# Read the STRING network file (download manually as explained above)
mg = mygene.MyGeneInfo()
chunk_size = 1000000
unique_ids = set()
filtered_chunks = []

print("Filtering STRING network (score >= 800)...")
for chunk in tqdm(pd.read_csv(string_raw_file, sep=' ', header=0, chunksize=chunk_size)):
    filtered = chunk[chunk['combined_score'] >= 800]
    filtered_chunks.append(filtered)
    unique_ids.update(filtered['protein1'].unique())
    unique_ids.update(filtered['protein2'].unique())

filtered_data = pd.concat(filtered_chunks, ignore_index=True)
del filtered_chunks

# Convert STRING format '9606.ENSP...' → ENSP only
ensp_ids = [id_.split('.')[1] for id_ in unique_ids]

print("Mapping ENSP to ENTREZ...")
entrez_mapping = {}
batch_size = 1000
for i in tqdm(range(0, len(ensp_ids), batch_size)):
    batch = ensp_ids[i:i+batch_size]
    results = mg.querymany(batch, scopes='ensembl.protein', fields='entrezgene', species='human')
    for res in results:
        if 'entrezgene' in res:
            entrez_mapping[f'9606.{res["query"]}'] = res['entrezgene']

def map_to_entrez(id_):
    return entrez_mapping.get(id_, None)

filtered_data['protein1_entrez'] = filtered_data['protein1'].apply(map_to_entrez)
filtered_data['protein2_entrez'] = filtered_data['protein2'].apply(map_to_entrez)

# Remove interactions where mapping failed
filtered_data.dropna(subset=['protein1_entrez', 'protein2_entrez'], inplace=True)

final_data = filtered_data[['protein1_entrez', 'protein2_entrez', 'combined_score']]
final_data.to_csv(string_entrez_file, sep='\t', index=False)
print(f"Filtered STRING network saved to '{string_entrez_file}'")

# ==============================
# STEP 2: READ SEED GENES + CONVERT SYMBOL → ENTREZ
# ==============================
with open(genes_seed_file, 'r') as f:
    seed_genes = [g.strip() for g in f.read().splitlines()]
print("Seed genes read:", seed_genes)

res = mg.querymany(seed_genes, scopes='symbol', fields='entrezgene', species='human')
seed_entrez = [str(r['entrezgene']) for r in res if 'entrezgene' in r]
print("Seed genes converted to Entrez IDs:", seed_entrez)

# ==============================
# STEP 3: BUILD NETWORK (NetworkX graph)
# ==============================
def read_string_network(file):
    df = pd.read_csv(file, sep='\t')
    G = nx.Graph()
    for _, row in df.iterrows():
        # Weight = STRING combined score
        G.add_edge(str(row['protein1_entrez']), str(row['protein2_entrez']),
                   weight=row['combined_score'])
    return G

G_string = read_string_network(string_entrez_file)
print(f"STRING network loaded with {G_string.number_of_nodes()} nodes and {G_string.number_of_edges()} edges.")

# ==============================
# STEP 4: ALGORITHMS
# DIAMOnD: Hypergeometric interaction enrichment
# GUILD: Score propagation through network neighbors
# ==============================
def diamond_algorithm(G, seed_genes, num_nodes=10):
    """
    DIAMOnD:
    At each iteration, select the node most significantly enriched
    in interactions with the current seed set (hypergeometric test).
    """
    added_nodes = []
    seed_set = set(seed_genes)
    candidate_nodes = set(G.nodes()) - seed_set

    for _ in range(num_nodes):
        best_node, best_p = None, 1
        for node in candidate_nodes:
            neighbors = set(G.neighbors(node))
            k_s = len(neighbors & seed_set)  # links to seed genes
            k = len(neighbors)              # total degree
            K = len(seed_set)               # number of seeds
            N = len(G)                      # graph size

            # Hypergeometric p-value (lower = more enriched)
            p = hypergeom.sf(k_s - 1, N, K, k)

            if p < best_p:
                best_node, best_p = node, p
        if best_node is None:
            break
        added_nodes.append((best_node, best_p))
        seed_set.add(best_node)
        candidate_nodes.remove(best_node)
    return added_nodes

def guild_netscore(G, seed_genes, max_iter=5):
    """
    GUILD NetScore:
    Seed genes start with score=1 and propagate influence
    through neighboring nodes iteratively.
    """
    score = {n: 0.0 for n in G.nodes()}
    for g in seed_genes:
        if g in score:
            score[g] = 1.0
    for _ in range(max_iter):
        new_score = {}
        for node in G.nodes():
            neighbors = list(G.neighbors(node))
            new_score[node] = sum(score[n] for n in neighbors) / len(neighbors) if neighbors else score[node]
        score = new_score
    return score

# ==============================
# STEP 5: CHECK SEED PRESENCE
# ==============================
present = [g for g in seed_entrez if g in G_string.nodes()]
missing = [g for g in seed_entrez if g not in G_string.nodes()]
print(f"\nSeed genes present: {len(present)}, missing: {len(missing)}")
if missing:
    print("Missing:", missing)
if len(present) == 0:
    raise ValueError("None of the seed genes are present in the network!")

# ==============================
# STEP 6: RUN ALGORITHMS
# ==============================
diamond_results = diamond_algorithm(G_string, present, num_nodes=TOP_N)
pd.DataFrame(diamond_results, columns=['node', 'p_value']).to_csv(os.path.join(results_dir, 'diamond_results.csv'), index=False)

guild_scores = guild_netscore(G_string, present)
guild_df = pd.DataFrame(sorted(guild_scores.items(), key=lambda x: x[1], reverse=True),
                        columns=['node', 'score'])
guild_df.nlargest(TOP_N, 'score').to_csv(os.path.join(results_dir, 'guild_results.csv'), index=False)

# ==============================
# STEP 7: SELECT GENES FOR FUNCTIONAL ANALYSIS
# Combine top DIAMOnD + GUILD
# ==============================
diamond_df = pd.read_csv(os.path.join(results_dir, 'diamond_results.csv'))
top_diamond_nodes = diamond_df.nsmallest(TOP_N, 'p_value')['node'].astype(str).tolist()

guild_df = pd.read_csv(os.path.join(results_dir, 'guild_results.csv'))
top_guild_nodes = guild_df.nlargest(TOP_N, 'score')['node'].astype(str).tolist()

nodes_for_functional = list(set(top_diamond_nodes + top_guild_nodes))

# Convert Entrez → HUGO
def entrez_to_symbol(entrez_ids: list[str]) -> list[str]:
    symbols = []
    for chunk in range(0, len(entrez_ids), 1000):
        results = mg.getgenes(entrez_ids[chunk:chunk+1000], fields='symbol', species='human')
        for res in results:
            if 'symbol' in res:
                symbols.append(res['symbol'])
    return symbols

nodes_for_functional_symbols = entrez_to_symbol(nodes_for_functional)

with open(os.path.join(results_dir, "nodes_for_functional.txt"), "w") as f:
    for gene in nodes_for_functional_symbols:
        f.write(f"{gene}\n")

# ==============================
# STEP 8: FUNCTIONAL ENRICHMENT (Enrichr)
# ==============================
def enrichr_analysis(genes, libraries):
    add_list_url = "https://maayanlab.cloud/Enrichr/addList"
    enrich_url = "https://maayanlab.cloud/Enrichr/enrich"
    payload = {'list': (None, "\n".join(genes)), 'description': (None, 'Functional analysis')}
    response = requests.post(add_list_url, files=payload)
    if response.status_code != 200:
        raise RuntimeError("Error submitting gene list to Enrichr.")
    user_list_id = response.json()['userListId']
    results = {}
    for lib in libraries:
        enrich_response = requests.get(enrich_url, params={'userListId': user_list_id, 'backgroundType': lib})
        if enrich_response.status_code == 200:
            data = enrich_response.json()
            if lib in data:
                results[lib] = data[lib][:5]  # Top 5 terms
    return results

def plot_enrichment(enrichment_results, out_dir):
    sns.set(style="whitegrid")
    for lib, terms in enrichment_results.items():
        if not terms:
            continue
        term_names = [entry[1] for entry in terms]
        scores = [entry[4] for entry in terms]
        plt.figure(figsize=(8, 5))
        sns.barplot(x=scores, y=term_names)
        plt.xlabel("Combined Score")
        plt.ylabel("Term")
        plt.title(f"Top 5 Enriched Terms: {lib}")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"top5_{lib.replace(' ', '_')}.png"), dpi=150)
        plt.close()

print("\nPerforming functional enrichment...")
libraries = [
    'KEGG_2021_Human',
    'GO_Biological_Process_2023',
    'GO_Molecular_Function_2023',
    'GO_Cellular_Component_2023'
]
result = enrichr_analysis(nodes_for_functional_symbols, libraries)
plot_enrichment(result, results_dir)

print("Functional analysis completed successfully.")