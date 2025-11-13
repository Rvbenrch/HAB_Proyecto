#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HAB CLI - Command Line Interface for HAB Project (Network Analysis)
"""

import os
import sys
import argparse
import subprocess#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
funcnet_pipeline.py — Análisis funcional con propagación en redes (RWR) para listas de genes.

Flujo:
1) Carga de genes de entrada (símbolos).
2) Conversión a IDs de STRING.
3) Descarga/red de interacciones desde STRING (añadiendo vecinos).
4) Propagación (Random Walk with Restart) desde los genes semilla.
5) Ranking de genes por puntuación de propagación.
6) Enriquecimiento funcional (GO BP y KEGG) usando Enrichr.
7) Visualizaciones y exportación de resultados.
8) (Opcional) DIAMOnD y GUILD + enriquecimiento combinado para bonus.

Ejemplo:
    python scripts/funcnet_pipeline.py \
        --genes-file data/genes_input.txt \
        --outdir results \
        --species 9606 \
        --string-score 700 \
        --neighbors 50 \
        --alpha 0.5 \
        --top-k 200 \
        --plot-top 15 --plot-width 16 --plot-font 9 \
        --with-diamond --with-guild --combo-top 50

Requisitos (requirements.txt):
    requests
    pandas
    numpy
    networkx
    matplotlib
    scipy
"""

import argparse
import os
import sys
import json
import time
import logging
from typing import Dict, List, Set

import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import hypergeom


# ----------------------------- Utilidades generales -----------------------------

def setup_logging(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("Logger inicializado. Resultados en: %s", outdir)


def save_params(outdir: str, args: argparse.Namespace) -> None:
    d = vars(args).copy()
    d["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir, "run_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)


def read_gene_list(path: str) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            g = line.strip()
            if g:
                genes.append(g)
    genes = sorted(list(dict.fromkeys(genes)))  # únicos y ordenados
    if not genes:
        raise ValueError("La lista de genes está vacía.")
    return genes


# ----------------------------- STRING helpers -----------------------------------

STRING_API = "https://string-db.org/api"
STRING_FORMAT = "json"


def string_map_genes_to_ids(genes: List[str], species: int = 9606) -> pd.DataFrame:
    """
    Usa el endpoint 'get_string_ids' para mapear símbolos a IDs de STRING.
    Devuelve DataFrame con: queryItem, stringId, preferredName, taxonId.
    """
    url = f"{STRING_API}/{STRING_FORMAT}/get_string_ids"
    if not genes:
        raise ValueError("Lista de genes vacía para mapear.")

    rows: List[dict] = []
    batch_size = 2000
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]
        attempts = [("GET", "\n"), ("GET", "\r"), ("POST", "\n"), ("POST", "\r")]
        got = False
        last_text = None
        for method, sep in attempts:
            params = {
                "identifiers": sep.join(batch),
                "species": species,
                "limit": 1,
                "caller_identity": "HAB_Proyecto_FuncNet"
            }
            try:
                r = requests.get(url, params=params, timeout=60) if method == "GET" else \
                    requests.post(url, data=params, timeout=60)
                last_text = r.text
                r.raise_for_status()
                payload = r.json()
                if isinstance(payload, list) and payload:
                    rows.extend(payload)
                    got = True
                    break
            except Exception as e:
                logging.warning("Intento %s con sep %r falló: %s", method, repr(sep), e)

        if not got:
            logging.error("STRING no devolvió resultados en el mapeo para el lote %s. Última respuesta: %s",
                          i//batch_size + 1, (last_text[:300] + "...") if last_text else "<sin cuerpo>")

    if not rows:
        raise RuntimeError("STRING no devolvió resultados para el mapeo de IDs.")

    df = pd.DataFrame(rows)
    for col in ["queryItem", "stringId", "preferredName", "taxonId"]:
        if col not in df.columns:
            df[col] = np.nan

    if "taxonId" in df.columns and df["taxonId"].notna().any():
        df = df[df["taxonId"] == species]

    logging.info("Mapeados %d/%d genes a IDs de STRING.", df["stringId"].notna().sum(), len(genes))
    if df["stringId"].notna().sum() < len(genes):
        faltan = sorted(set(genes) - set(df["queryItem"].dropna()))
        if faltan:
            logging.warning("No se encontraron en STRING: %s", ", ".join(faltan))
    return df


def string_fetch_network(
    string_ids: List[str],
    species: int = 9606,
    required_score: int = 700,
    add_neighbors: int = 50
) -> pd.DataFrame:
    """
    Descarga red PPI desde STRING usando endpoint 'network'.
    Devuelve DataFrame con: preferredName_A, preferredName_B, score (normalizado a [0,1] si aplica)
    """
    if not string_ids:
        raise ValueError("No hay IDs de STRING para descargar red.")

    url = f"{STRING_API}/{STRING_FORMAT}/network"
    attempts = [("POST", "\n"), ("GET", "\n")]
    data = None
    last_text = None
    for method, sep in attempts:
        params = {
            "identifiers": sep.join(string_ids),
            "species": species,
            "required_score": required_score,
            "add_nodes": add_neighbors,
            "caller_identity": "funcnet_pipeline"
        }
        try:
            r = requests.get(url, params=params, timeout=120) if method == "GET" else \
                requests.post(url, data=params, timeout=120)
            last_text = r.text
            r.raise_for_status()
            data = r.json()
            if data:
                break
        except Exception as e:
            logging.warning("Descarga de red STRING intento %s falló: %s", method, e)

    if not data:
        raise RuntimeError(f"STRING devolvió red vacía o error. Última respuesta: {(last_text[:300] + '...') if last_text else '<sin cuerpo>'}")

    df = pd.DataFrame(data)
    needed = ["preferredName_A", "preferredName_B", "score"]
    for col in needed:
        if col not in df.columns:
            raise RuntimeError(f"Falta la columna '{col}' en la red descargada.")

    if df["score"].max() > 1.0:
        df["score"] = df["score"] / 1000.0

    df = df[df["preferredName_A"] != df["preferredName_B"]]
    df = df.drop_duplicates(subset=["preferredName_A", "preferredName_B"])
    return df


def build_graph_from_edges(df_edges: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for _, row in df_edges.iterrows():
        u = str(row["preferredName_A"])
        v = str(row["preferredName_B"])
        w = float(row["score"])
        if u and v:
            G.add_edge(u, v, weight=w)
    return G


# ----------------------------- Propagación (RWR) --------------------------------

def random_walk_with_restart(
    G: nx.Graph,
    seeds: Set[str],
    alpha: float = 0.5,
    tol: float = 1e-8,
    max_iter: int = 10_000
) -> pd.DataFrame:
    """
    p_{t+1} = (1 - alpha) * W * p_t + alpha * p0, con W columna-estocástica.
    Devuelve DataFrame con columnas: gene, score, is_seed.
    """
    if len(G) == 0:
        raise ValueError("El grafo está vacío.")

    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    deg = np.array([G.degree(n, weight="weight") for n in nodes], dtype=float)

    seeds = {s for s in seeds if s in idx}
    if not seeds:
        raise ValueError("Ninguna semilla está presente en la red.")
    p0 = np.zeros(n, dtype=float)
    for s in seeds:
        p0[idx[s]] = 1.0 / len(seeds)

    p = p0.copy()
    A = {u: {v: d.get("weight", 1.0) for v, d in G[u].items()} for u in nodes}

    for it in range(max_iter):
        p_next = np.zeros_like(p)
        for j, nj in enumerate(nodes):
            if deg[j] == 0 or p[j] == 0:
                continue
            for i, w in A[nj].items():
                ii = idx[i]
                p_next[ii] += (w / deg[j]) * p[j]
        p_next = (1 - alpha) * p_next + alpha * p0
        if np.linalg.norm(p_next - p, 1) < tol:
            p = p_next
            logging.info("RWR convergió en %d iteraciones.", it + 1)
            break
        p = p_next
    else:
        logging.warning("RWR alcanzó el máximo de iteraciones (%d) sin converger.", max_iter)

    out = pd.DataFrame({
        "gene": nodes,
        "score": p,
        "is_seed": [n in seeds for n in nodes]
    }).sort_values("score", ascending=False).reset_index(drop=True)
    return out


# ----------------------------- Enrichr helpers ----------------------------------

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"


def enrichr_submit_list(genes: List[str], description: str = "funcnet_pipeline") -> str:
    url = f"{ENRICHR_BASE}/addList"
    payload = {"list": "\n".join(genes), "description": description}
    r = requests.post(url, files=payload, timeout=60)
    r.raise_for_status()
    j = r.json()
    if "userListId" not in j:
        raise RuntimeError("No se obtuvo userListId de Enrichr.")
    return str(j["userListId"])


def enrichr_get_results(list_id: str, library: str) -> pd.DataFrame:
    url = f"{ENRICHR_BASE}/enrich"
    params = {"userListId": list_id, "backgroundType": library}
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    j = r.json()
    if library not in j or not j[library]:
        return pd.DataFrame()
    cols = [
        "rank", "term_name", "pvalue", "zscore", "combined_score",
        "overlapping_genes", "adjusted_pvalue", "old_pvalue", "old_adjusted_pvalue"
    ]
    rows = []
    for row in j[library]:
        if not isinstance(row, list) or len(row) < 9:
            continue
        rows.append({
            "rank": row[0],
            "term_name": row[1],
            "pvalue": row[2],
            "zscore": row[3],
            "combined_score": row[4],
            "overlapping_genes": row[5],
            "adjusted_pvalue": row[6],
            "old_pvalue": row[7],
            "old_adjusted_pvalue": row[8],
        })
    df = pd.DataFrame(rows, columns=cols)
    df = df.sort_values(["adjusted_pvalue", "pvalue", "combined_score"], ascending=[True, True, False])
    return df


def run_enrichment(genes: List[str], out_csv_dir: str, prefix: str) -> Dict[str, str]:
    """
    Ejecuta enriquecimiento en GO BP y KEGG.
    Devuelve dict {librería: ruta_csv}
    """
    os.makedirs(out_csv_dir, exist_ok=True)
    list_id = enrichr_submit_list(genes, description=f"{prefix}_enrichr")
    libraries = {
        "GO_Biological_Process_2021": f"{prefix}_enrichr_GO_BP.csv",
        "KEGG_2019_Human": f"{prefix}_enrichr_KEGG.csv",
    }
    saved: Dict[str, str] = {}
    for lib, fname in libraries.items():
        try:
            df = enrichr_get_results(list_id, lib)
            path = os.path.join(out_csv_dir, fname)
            df.to_csv(path, index=False)
            saved[lib] = path
            logging.info("Guardado enriquecimiento %s en %s", lib, path)
        except Exception as e:
            logging.exception("Fallo en enriquecimiento %s: %s", lib, e)
    return saved


# ----------------------------- Visualizaciones ----------------------------------

def plot_go_bar(
    enrich_csv: str,
    out_png: str,
    out_svg: str,
    top_n: int = 15,
    wrap_width: int = 35,
    fig_width: float = 16.0,
    font_size: int = 9,
    dpi: int = 300
) -> None:
    """
    Barplot legible para GO BP con etiquetas envueltas y altura dinámica. Guarda PNG+SVG.
    """
    import re, textwrap

    if not os.path.exists(enrich_csv):
        logging.warning("No existe %s para graficar GO.", enrich_csv)
        return
    df = pd.read_csv(enrich_csv)
    if df.empty:
        logging.warning("Enriquecimiento GO vacío, no se genera gráfico.")
        return

    df = df.head(top_n).copy()
    df["neglog10_adjP"] = -np.log10(np.clip(df["adjusted_pvalue"].replace(0, 1e-300), 1e-300, 1.0))
    df = df.iloc[::-1]

    def _wrap_label(s: str, width: int = 35) -> str:
        s = str(s)
        m = re.search(r"(.*)\s+\((GO:\d+)\)$", s)
        if m:
            name, go = m.groups()
            name_wrapped = textwrap.fill(name, width=width)
            return f"{name_wrapped}\n{go}"
        return textwrap.fill(s, width=width)

    labels = [_wrap_label(t, width=wrap_width) for t in df["term_name"]]
    n = len(labels)
    fig_height = max(3.5, min(0.8 * n + 1.5, 24.0))

    plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
    plt.barh(range(n), df["neglog10_adjP"])
    plt.yticks(range(n), labels, fontsize=font_size)
    plt.xlabel("-log10(FDR)", fontsize=font_size)
    plt.title("Top términos GO: Biological Process (Enrichr)", fontsize=font_size + 2, pad=12)
    plt.gcf().subplots_adjust(left=0.45, right=0.98, top=0.92, bottom=0.08)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()
    logging.info("Gráficos GO guardados en %s y %s", out_png, out_svg)


def plot_network_scores(G: nx.Graph, scores_df: pd.DataFrame, out_png: str) -> None:
    if len(G) == 0:
        logging.warning("Grafo vacío, no se grafica red.")
        return

    scores = dict(zip(scores_df["gene"], scores_df["score"]))
    nodes = list(G.nodes())
    vals = np.array([scores.get(n, 0.0) for n in nodes], dtype=float)
    sizes = 200 * (vals / vals.max()) + 50 if vals.max() > 0 else np.full_like(vals, 50) + 0.0

    pos = nx.spring_layout(G, seed=42, weight="weight", k=None)
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=vals, cmap=plt.cm.viridis)
    top_nodes = list(scores_df.head(30)["gene"])
    labels = {n: n for n in top_nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    plt.title("Red STRING con puntuación de propagación (RWR)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close()
    logging.info("Gráfico de red guardado en %s", out_png)

# ----------------------------- DIAMOnD / GUILD (opcionales) ---------------------

def diamond_algorithm(G: nx.Graph, seed_genes: Set[str], num_nodes: int = 10):
    added_nodes, seed_set = [], set(seed_genes)
    candidate_nodes = set(G.nodes()) - seed_set
    for _ in range(num_nodes):
        best_node, best_p = None, 1
        for node in candidate_nodes:
            neighbors = set(G.neighbors(node))
            k_s = len(neighbors & seed_set)
            k = len(neighbors)
            K = len(seed_set)
            N = len(G)
            p = hypergeom.sf(k_s - 1, N, K, k)
            if p < best_p:
                best_node, best_p = node, p
        if best_node is None:
            break
        added_nodes.append((best_node, best_p))
        seed_set.add(best_node)
        candidate_nodes.remove(best_node)
    return added_nodes


def guild_netscore(G: nx.Graph, seed_genes: Set[str], max_iter: int = 5):
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


# ----------------------------- CLI principal ------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Análisis funcional con propagación en redes (STRING + RWR + Enrichr)."
    )
    parser.add_argument("--genes-file", type=str, default="data/genes_input.txt",
                        help="Ruta al archivo de genes (una ID por línea, p. ej. símbolos HGNC).")
    parser.add_argument("--outdir", type=str, default="results",
                        help="Directorio de salida para resultados.")
    parser.add_argument("--species", type=int, default=9606,
                        help="TaxID NCBI para STRING (9606=Human).")
    parser.add_argument("--string-score", type=int, default=700,
                        help="Puntuación mínima de STRING (400/700/900 ~ med/alta/máxima).")
    parser.add_argument("--neighbors", type=int, default=50,
                        help="Nº de vecinos a añadir desde STRING (first shell).")
    parser.add_argument("--alpha", type=float, default=0.5,
                        help="Parámetro de reinicio (restart) para RWR (0<alpha<1).")
    parser.add_argument("--top-k", type=int, default=200,
                        help="Top genes por score para enriquecimiento/tabla.")
    parser.add_argument("--no-enrich", action="store_true",
                        help="Si se indica, omite enriquecimiento funcional.")
    parser.add_argument("--strict-symbols", action="store_true",
                        help="Mantiene solo matches cuyo preferredName coincide con el símbolo de entrada (case-insensitive).")

    # Opciones de visualización (mini-mejora)
    parser.add_argument("--plot-top", type=int, default=15,
                        help="Nº de términos a mostrar en el barplot GO.")
    parser.add_argument("--plot-width", type=float, default=16.0,
                        help="Ancho de figura para barplots.")
    parser.add_argument("--plot-font", type=int, default=9,
                        help="Tamaño de fuente en barplots.")

    # Opcionales extra (bonus)
    parser.add_argument("--with-diamond", action="store_true",
                        help="Ejecuta DIAMOnD y guarda diamond_results.csv")
    parser.add_argument("--with-guild", action="store_true",
                        help="Ejecuta NetScore (GUILD) y guarda guild_results.csv")
    parser.add_argument("--combo-top", type=int, default=50,
                        help="Top para combinar DIAMOnD+GUILD con RWR (para enriquecimiento combinado)")

    args = parser.parse_args()

    # --- preparar carpetas de salida organizadas ---
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    DIR_CSV     = os.path.join(outdir, "csv")
    DIR_FIG     = os.path.join(outdir, "figures")
    DIR_RWR     = os.path.join(outdir, "rwr")
    DIR_GRAPHS  = os.path.join(outdir, "graphs")
    for d in (DIR_CSV, DIR_FIG, DIR_RWR, DIR_GRAPHS):
        os.makedirs(d, exist_ok=True)

    setup_logging(outdir)
    save_params(outdir, args)

    # 1) Cargar genes
    genes_in = read_gene_list(args.genes_file)
    logging.info("Genes de entrada (%d): %s", len(genes_in), ", ".join(genes_in))

    # 2) Mapear a IDs de STRING
    logging.info("Mapeando a IDs de STRING...")
    map_df = string_map_genes_to_ids(genes_in, species=args.species)

    # Guardar mapeos
    map_raw_csv  = os.path.join(DIR_CSV, "mapping_string_ids_raw.csv")
    map_clean_csv = os.path.join(DIR_CSV, "mapping_string_ids.csv")
    map_df.to_csv(map_raw_csv, index=False)

    # Avisar de desajustes símbolo→preferredName
    try:
        cmp_df = map_df.dropna(subset=["stringId"]).copy()
        _q = cmp_df["queryItem"].astype(str).str.upper().str.strip()
        _p = cmp_df["preferredName"].astype(str).str.upper().str.strip()
        mismatch = cmp_df[_q != _p]
        if not mismatch.empty:
            logging.warning("Se detectaron mapeos no idénticos símbolo→preferredName. Revisa %s", map_clean_csv)
    except Exception:
        pass

    mapped = map_df.dropna(subset=["stringId"]).copy()
    if args.strict_symbols:
        mapped["__q"] = mapped["queryItem"].astype(str).str.upper().str.strip()
        mapped["__p"] = mapped["preferredName"].astype(str).str.upper().str.strip()
        before = len(mapped)
        mapped = mapped[mapped["__q"] == mapped["__p"]].drop(columns=["__q", "__p"])
        after = len(mapped)
        if after < before:
            logging.warning("STRICT: %d mapeos descartados por no coincidir símbolo-preferredName. Conservados: %d",
                            before - after, after)

    if mapped.empty:
        raise RuntimeError("No se obtuvieron IDs de STRING para las semillas tras el mapeo/filtrado.")

    mapped.to_csv(map_clean_csv, index=False)
    logging.info("Guardado mapeo filtrado en %s", map_clean_csv)

    seed_names = set(mapped["preferredName"].astype(str).tolist())
    string_ids = mapped["stringId"].astype(str).tolist()
    logging.info("Semillas mapeadas en STRING (%d): %s", len(seed_names), ", ".join(sorted(seed_names)))

    # 3) Descargar red de STRING
    logging.info("Descargando red de STRING con score>=%d y %d vecinos...", args.string_score, args.neighbors)
    edges_df = string_fetch_network(
        string_ids=string_ids,
        species=args.species,
        required_score=args.string_score,
        add_neighbors=args.neighbors
    )
    edges_csv = os.path.join(DIR_CSV, "string_network_edges.csv")
    edges_tsv = os.path.join(DIR_CSV, "string_network_edges.tsv")
    edges_df.to_csv(edges_csv, index=False)
    edges_df.to_csv(edges_tsv, sep="\t", index=False)  # copia TSV para traza
    logging.info("Red: %d aristas entre %d genes. Guardada en %s",
                 len(edges_df), len(pd.unique(edges_df[["preferredName_A", "preferredName_B"]].values.ravel())), edges_csv)

    # 3b) Guardar edgelist con pesos
    G = build_graph_from_edges(edges_df)
    edgelist_path = os.path.join(DIR_GRAPHS, "network_weighted.edgelist")
    nx.write_weighted_edgelist(G, edgelist_path)
    logging.info("Edgelist (ponderada) guardada en %s", edgelist_path)

    # 4) RWR
    logging.info("Grafo con %d nodos y %d aristas.", G.number_of_nodes(), G.number_of_edges())
    scores_df = random_walk_with_restart(G, seeds=seed_names, alpha=args.alpha)
    scores_csv = os.path.join(DIR_RWR, "propagation_scores.csv")
    scores_df.to_csv(scores_csv, index=False)
    logging.info("Puntuaciones de propagación guardadas en %s", scores_csv)

    # 5) Top genes para enriquecimiento
    top_k = min(args.top_k, len(scores_df))
    top_df = scores_df.head(top_k).copy()
    top_csv = os.path.join(DIR_RWR, f"top_genes_{top_k}.csv")
    top_df.to_csv(top_csv, index=False)
    logging.info("Top %d genes guardados en %s", top_k, top_csv)

    # 6) Enriquecimiento base (RWR)
    saved_enrich: Dict[str, str] = {}
    if not args.no_enrich:
        try:
            saved_enrich = run_enrichment(
                genes=top_df["gene"].tolist(),
                out_csv_dir=DIR_CSV,
                prefix=f"top{top_k}"
            )
        except Exception as e:
            logging.exception("Error en enriquecimiento funcional: %s", e)

    # 7) Visualizaciones
    # 7a) GO BP (RWR)
    if "GO_Biological_Process_2021" in saved_enrich:
        go_png = os.path.join(DIR_FIG, "go_bp_barplot.png")
        go_svg = os.path.join(DIR_FIG, "go_bp_barplot.svg")
        plot_go_bar(
            saved_enrich["GO_Biological_Process_2021"],
            out_png=go_png,
            out_svg=go_svg,
            top_n=args.plot_top,
            fig_width=args.plot_width,
            font_size=args.plot_font
        )

    # 7b) Grafo con scores
    net_png = os.path.join(DIR_FIG, "network_scores.png")
    plot_network_scores(G, scores_df, net_png)

    # 8) Opcional: DIAMOnD/GUILD + enriquecimiento combinado
    if args.with_diamond or args.with_guild:
        combo_genes = set(top_df["gene"].tolist())

        if args.with_diamond:
            dres = diamond_algorithm(G, seed_names, num_nodes=args.combo_top)
            d_csv = os.path.join(DIR_CSV, "diamond_results.csv")
            pd.DataFrame(dres, columns=["node", "p_value"]).to_csv(d_csv, index=False)
            logging.info("DIAMOnD guardado en %s", d_csv)
            try:
                d = pd.read_csv(d_csv)
                combo_genes |= set(d.nsmallest(args.combo_top, "p_value")["node"].astype(str))
            except Exception:
                pass

        if args.with_guild:
            gscores = guild_netscore(G, seed_names)
            gdf = pd.DataFrame(sorted(gscores.items(), key=lambda x: x[1], reverse=True),
                               columns=["node", "score"])
            g_csv = os.path.join(DIR_CSV, "guild_results.csv")
            gdf.head(args.combo_top).to_csv(g_csv, index=False)
            logging.info("GUILD guardado en %s", g_csv)
            try:
                g = pd.read_csv(g_csv)
                combo_genes |= set(g.nlargest(args.combo_top, "score")["node"].astype(str))
            except Exception:
                pass

        if not args.no_enrich and len(combo_genes) > 0:
            saved_combo = run_enrichment(sorted(combo_genes), DIR_CSV, prefix=f"combo_top{len(combo_genes)}")
            if "GO_Biological_Process_2021" in saved_combo:
                go_combo_png = os.path.join(DIR_FIG, "go_bp_barplot_combo.png")
                go_combo_svg = os.path.join(DIR_FIG, "go_bp_barplot_combo.svg")
                plot_go_bar(
                    saved_combo["GO_Biological_Process_2021"],
                    out_png=go_combo_png,
                    out_svg=go_combo_svg,
                    top_n=args.plot_top,
                    fig_width=args.plot_width,
                    font_size=args.plot_font
                )

    # 9) README de resultados (referenciando subcarpetas)
    results_md = os.path.join(outdir, "README_results.md")
    with open(results_md, "w", encoding="utf-8") as f:
        f.write(f"# Resultados – FuncNet Pipeline\n\n")
        f.write(f"**Fecha:** {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Estructura de salidas\n")
        f.write("- `csv/`: mapeos, red, enriquecimientos, etc.\n")
        f.write("- `figures/`: gráficos en PNG y SVG (GO, red).\n")
        f.write("- `rwr/`: puntuaciones de propagación y top genes.\n")
        f.write("- `graphs/`: representaciones de la red (edgelist).\n\n")
        f.write("## Archivos clave\n")
        f.write(f"- `csv/mapping_string_ids_raw.csv`, `csv/mapping_string_ids.csv`\n")
        f.write(f"- `csv/string_network_edges.csv` (+ `string_network_edges.tsv`)\n")
        f.write(f"- `graphs/network_weighted.edgelist`\n")
        f.write(f"- `rwr/propagation_scores.csv`, `rwr/top_genes_{top_k}.csv`\n")
        f.write(f"- `figures/go_bp_barplot.png|.svg`, `figures/network_scores.png`\n")
        f.write(f"- `pipeline.log`, `run_metadata.json`\n")
        f.write("(si activa bonus) `csv/diamond_results.csv`, `csv/guild_results.csv`, "
                "`figures/go_bp_barplot_combo.png|.svg`\n\n")
        f.write("## Parámetros\n```json\n")
        json.dump(vars(args), f, indent=2)
        f.write("\n```\n")

    logging.info("Pipeline completado. Revisa %s", outdir)


if __name__ == "__main__":
    main()


# Get the project root directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, 'scripts')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
NETWORK_SCRIPT = os.path.join(SCRIPTS_DIR, 'network_propagation_and_functional_analysis.py')

# Get virtual environment python
VENV_PYTHON = os.path.join(PROJECT_ROOT, 'venv', 'bin', 'python')
if not os.path.exists(VENV_PYTHON):
    VENV_PYTHON = sys.executable

def get_demo_args():
    """Returns demo arguments for testing"""
    # Look for STRING network file in data directory
    string_files = [f for f in os.listdir(DATA_DIR) if f.endswith('.tsv') and 'contrast' in f.lower()]
    
    if not string_files:
        print(f"Error: No network file found in {DATA_DIR}")
        print("Expected a file with pattern '*contrast*.tsv' or similar")
        return None
    
    string_input = os.path.join(DATA_DIR, string_files[0])
    seed_genes = os.path.join(DATA_DIR, 'genes_input.txt')
    
    if not os.path.exists(seed_genes):
        print(f"Error: Seed genes file not found at {seed_genes}")
        return None
    
    return {
        'string_input': string_input,
        'seed_genes': seed_genes,
        'output_dir': RESULTS_DIR,
        'top_n': 10
    }

def run_analysis(args):
    """Run the network analysis script with given arguments"""
    cmd = [
        VENV_PYTHON,
        NETWORK_SCRIPT,
        '--string_input', args['string_input'],
        '--seed_genes', args['seed_genes'],
        '--output_dir', args['output_dir'],
        '--top_n', str(args['top_n'])
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    print("-" * 80)
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error running analysis: {e}")
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HAB CLI - Command Line Interface for HAB Project (Network Analysis)
"""

import os
import sys
import argparse
import subprocess#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
funcnet_pipeline.py — Análisis funcional con propagación en redes (RWR) para listas de genes.

Flujo:
1) Carga de genes de entrada (símbolos).
2) Conversión a IDs de STRING.
3) Descarga/red de interacciones desde STRING (añadiendo vecinos).
4) Propagación (Random Walk with Restart) desde los genes semilla.
5) Ranking de genes por puntuación de propagación.
6) Enriquecimiento funcional (GO BP y KEGG) usando Enrichr.
7) Visualizaciones y exportación de resultados.
8) (Opcional) DIAMOnD y GUILD + enriquecimiento combinado para bonus.

Ejemplo:
    python scripts/funcnet_pipeline.py \
        --genes-file data/genes_input.txt \
        --outdir results \
        --species 9606 \
        --string-score 700 \
        --neighbors 50 \
        --alpha 0.5 \
        --top-k 200 \
        --plot-top 15 --plot-width 16 --plot-font 9 \
        --with-diamond --with-guild --combo-top 50

Requisitos (requirements.txt):
    requests
    pandas
    numpy
    networkx
    matplotlib
    scipy
"""

import argparse
import os
import sys
import json
import time
import logging
from typing import Dict, List, Set

import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import hypergeom


# ----------------------------- Utilidades generales -----------------------------

def setup_logging(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    log_path = os.path.join(outdir, "pipeline.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("Logger inicializado. Resultados en: %s", outdir)


def save_params(outdir: str, args: argparse.Namespace) -> None:
    d = vars(args).copy()
    d["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(os.path.join(outdir, "run_metadata.json"), "w", encoding="utf-8") as f:
        json.dump(d, f, indent=2, ensure_ascii=False)


def read_gene_list(path: str) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            g = line.strip()
            if g:
                genes.append(g)
    genes = sorted(list(dict.fromkeys(genes)))  # únicos y ordenados
    if not genes:
        raise ValueError("La lista de genes está vacía.")
    return genes


# ----------------------------- STRING helpers -----------------------------------

STRING_API = "https://string-db.org/api"
STRING_FORMAT = "json"


def string_map_genes_to_ids(genes: List[str], species: int = 9606) -> pd.DataFrame:
    """
    Usa el endpoint 'get_string_ids' para mapear símbolos a IDs de STRING.
    Devuelve DataFrame con: queryItem, stringId, preferredName, taxonId.
    """
    url = f"{STRING_API}/{STRING_FORMAT}/get_string_ids"
    if not genes:
        raise ValueError("Lista de genes vacía para mapear.")

    rows: List[dict] = []
    batch_size = 2000
    for i in range(0, len(genes), batch_size):
        batch = genes[i:i + batch_size]
        attempts = [("GET", "\n"), ("GET", "\r"), ("POST", "\n"), ("POST", "\r")]
        got = False
        last_text = None
        for method, sep in attempts:
            params = {
                "identifiers": sep.join(batch),
                "species": species,
                "limit": 1,
                "caller_identity": "HAB_Proyecto_FuncNet"
            }
            try:
                r = requests.get(url, params=params, timeout=60) if method == "GET" else \
                    requests.post(url, data=params, timeout=60)
                last_text = r.text
                r.raise_for_status()
                payload = r.json()
                if isinstance(payload, list) and payload:
                    rows.extend(payload)
                    got = True
                    break
            except Exception as e:
                logging.warning("Intento %s con sep %r falló: %s", method, repr(sep), e)

        if not got:
            logging.error("STRING no devolvió resultados en el mapeo para el lote %s. Última respuesta: %s",
                          i//batch_size + 1, (last_text[:300] + "...") if last_text else "<sin cuerpo>")

    if not rows:
        raise RuntimeError("STRING no devolvió resultados para el mapeo de IDs.")

    df = pd.DataFrame(rows)
    for col in ["queryItem", "stringId", "preferredName", "taxonId"]:
        if col not in df.columns:
            df[col] = np.nan

    if "taxonId" in df.columns and df["taxonId"].notna().any():
        df = df[df["taxonId"] == species]

    logging.info("Mapeados %d/%d genes a IDs de STRING.", df["stringId"].notna().sum(), len(genes))
    if df["stringId"].notna().sum() < len(genes):
        faltan = sorted(set(genes) - set(df["queryItem"].dropna()))
        if faltan:
            logging.warning("No se encontraron en STRING: %s", ", ".join(faltan))
    return df


def string_fetch_network(
    string_ids: List[str],
    species: int = 9606,
    required_score: int = 700,
    add_neighbors: int = 50
) -> pd.DataFrame:
    """
    Descarga red PPI desde STRING usando endpoint 'network'.
    Devuelve DataFrame con: preferredName_A, preferredName_B, score (normalizado a [0,1] si aplica)
    """
    if not string_ids:
        raise ValueError("No hay IDs de STRING para descargar red.")

    url = f"{STRING_API}/{STRING_FORMAT}/network"
    attempts = [("POST", "\n"), ("GET", "\n")]
    data = None
    last_text = None
    for method, sep in attempts:
        params = {
            "identifiers": sep.join(string_ids),
            "species": species,
            "required_score": required_score,
            "add_nodes": add_neighbors,
            "caller_identity": "funcnet_pipeline"
        }
        try:
            r = requests.get(url, params=params, timeout=120) if method == "GET" else \
                requests.post(url, data=params, timeout=120)
            last_text = r.text
            r.raise_for_status()
            data = r.json()
            if data:
                break
        except Exception as e:
            logging.warning("Descarga de red STRING intento %s falló: %s", method, e)

    if not data:
        raise RuntimeError(f"STRING devolvió red vacía o error. Última respuesta: {(last_text[:300] + '...') if last_text else '<sin cuerpo>'}")

    df = pd.DataFrame(data)
    needed = ["preferredName_A", "preferredName_B", "score"]
    for col in needed:
        if col not in df.columns:
            raise RuntimeError(f"Falta la columna '{col}' en la red descargada.")

    if df["score"].max() > 1.0:
        df["score"] = df["score"] / 1000.0

    df = df[df["preferredName_A"] != df["preferredName_B"]]
    df = df.drop_duplicates(subset=["preferredName_A", "preferredName_B"])
    return df


def build_graph_from_edges(df_edges: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for _, row in df_edges.iterrows():
        u = str(row["preferredName_A"])
        v = str(row["preferredName_B"])
        w = float(row["score"])
        if u and v:
            G.add_edge(u, v, weight=w)
    return G


# ----------------------------- Propagación (RWR) --------------------------------

def random_walk_with_restart(
    G: nx.Graph,
    seeds: Set[str],
    alpha: float = 0.5,
    tol: float = 1e-8,
    max_iter: int = 10_000
) -> pd.DataFrame:
    """
    p_{t+1} = (1 - alpha) * W * p_t + alpha * p0, con W columna-estocástica.
    Devuelve DataFrame con columnas: gene, score, is_seed.
    """
    if len(G) == 0:
        raise ValueError("El grafo está vacío.")

    nodes = sorted(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    deg = np.array([G.degree(n, weight="weight") for n in nodes], dtype=float)

    seeds = {s for s in seeds if s in idx}
    if not seeds:
        raise ValueError("Ninguna semilla está presente en la red.")
    p0 = np.zeros(n, dtype=float)
    for s in seeds:
        p0[idx[s]] = 1.0 / len(seeds)

    p = p0.copy()
    A = {u: {v: d.get("weight", 1.0) for v, d in G[u].items()} for u in nodes}

    for it in range(max_iter):
        p_next = np.zeros_like(p)
        for j, nj in enumerate(nodes):
            if deg[j] == 0 or p[j] == 0:
                continue
            for i, w in A[nj].items():
                ii = idx[i]
                p_next[ii] += (w / deg[j]) * p[j]
        p_next = (1 - alpha) * p_next + alpha * p0
        if np.linalg.norm(p_next - p, 1) < tol:
            p = p_next
            logging.info("RWR convergió en %d iteraciones.", it + 1)
            break
        p = p_next
    else:
        logging.warning("RWR alcanzó el máximo de iteraciones (%d) sin converger.", max_iter)

    out = pd.DataFrame({
        "gene": nodes,
        "score": p,
        "is_seed": [n in seeds for n in nodes]
    }).sort_values("score", ascending=False).reset_index(drop=True)
    return out


# ----------------------------- Enrichr helpers ----------------------------------

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"


def enrichr_submit_list(genes: List[str], description: str = "funcnet_pipeline") -> str:
    url = f"{ENRICHR_BASE}/addList"
    payload = {"list": "\n".join(genes), "description": description}
    r = requests.post(url, files=payload, timeout=60)
    r.raise_for_status()
    j = r.json()
    if "userListId" not in j:
        raise RuntimeError("No se obtuvo userListId de Enrichr.")
    return str(j["userListId"])


def enrichr_get_results(list_id: str, library: str) -> pd.DataFrame:
    url = f"{ENRICHR_BASE}/enrich"
    params = {"userListId": list_id, "backgroundType": library}
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    j = r.json()
    if library not in j or not j[library]:
        return pd.DataFrame()
    cols = [
        "rank", "term_name", "pvalue", "zscore", "combined_score",
        "overlapping_genes", "adjusted_pvalue", "old_pvalue", "old_adjusted_pvalue"
    ]
    rows = []
    for row in j[library]:
        if not isinstance(row, list) or len(row) < 9:
            continue
        rows.append({
            "rank": row[0],
            "term_name": row[1],
            "pvalue": row[2],
            "zscore": row[3],
            "combined_score": row[4],
            "overlapping_genes": row[5],
            "adjusted_pvalue": row[6],
            "old_pvalue": row[7],
            "old_adjusted_pvalue": row[8],
        })
    df = pd.DataFrame(rows, columns=cols)
    df = df.sort_values(["adjusted_pvalue", "pvalue", "combined_score"], ascending=[True, True, False])
    return df


def run_enrichment(genes: List[str], out_csv_dir: str, prefix: str) -> Dict[str, str]:
    """
    Ejecuta enriquecimiento en GO BP y KEGG.
    Devuelve dict {librería: ruta_csv}
    """
    os.makedirs(out_csv_dir, exist_ok=True)
    list_id = enrichr_submit_list(genes, description=f"{prefix}_enrichr")
    libraries = {
        "GO_Biological_Process_2021": f"{prefix}_enrichr_GO_BP.csv",
        "KEGG_2019_Human": f"{prefix}_enrichr_KEGG.csv",
    }
    saved: Dict[str, str] = {}
    for lib, fname in libraries.items():
        try:
            df = enrichr_get_results(list_id, lib)
            path = os.path.join(out_csv_dir, fname)
            df.to_csv(path, index=False)
            saved[lib] = path
            logging.info("Guardado enriquecimiento %s en %s", lib, path)
        except Exception as e:
            logging.exception("Fallo en enriquecimiento %s: %s", lib, e)
    return saved

# ----------------------------- Visualizaciones ----------------------------------

def plot_go_bar(
    enrich_csv: str,
    out_png: str,
    out_svg: str,
    top_n: int = 15,
    wrap_width: int = 35,
    fig_width: float = 16.0,
    font_size: int = 9,
    dpi: int = 300
) -> None:
    """
    Barplot legible para GO BP con etiquetas envueltas y altura dinámica. Guarda PNG+SVG.
    """
    import re, textwrap

    if not os.path.exists(enrich_csv):
        logging.warning("No existe %s para graficar GO.", enrich_csv)
        return
    df = pd.read_csv(enrich_csv)
    if df.empty:
        logging.warning("Enriquecimiento GO vacío, no se genera gráfico.")
        return

    df = df.head(top_n).copy()
    df["neglog10_adjP"] = -np.log10(np.clip(df["adjusted_pvalue"].replace(0, 1e-300), 1e-300, 1.0))
    df = df.iloc[::-1]

    def _wrap_label(s: str, width: int = 35) -> str:
        s = str(s)
        m = re.search(r"(.*)\s+\((GO:\d+)\)$", s)
        if m:
            name, go = m.groups()
            name_wrapped = textwrap.fill(name, width=width)
            return f"{name_wrapped}\n{go}"
        return textwrap.fill(s, width=width)

    labels = [_wrap_label(t, width=wrap_width) for t in df["term_name"]]
    n = len(labels)
    fig_height = max(3.5, min(0.8 * n + 1.5, 24.0))

    plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
    plt.barh(range(n), df["neglog10_adjP"])
    plt.yticks(range(n), labels, fontsize=font_size)
    plt.xlabel("-log10(FDR)", fontsize=font_size)
    plt.title("Top términos GO: Biological Process (Enrichr)", fontsize=font_size + 2, pad=12)
    plt.gcf().subplots_adjust(left=0.45, right=0.98, top=0.92, bottom=0.08)
    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_svg, bbox_inches="tight")
    plt.close()
    logging.info("Gráficos GO guardados en %s y %s", out_png, out_svg)


def plot_network_scores(G: nx.Graph, scores_df: pd.DataFrame, out_png: str) -> None:
    if len(G) == 0:
        logging.warning("Grafo vacío, no se grafica red.")
        return

    scores = dict(zip(scores_df["gene"], scores_df["score"]))
    nodes = list(G.nodes())
    vals = np.array([scores.get(n, 0.0) for n in nodes], dtype=float)
    sizes = 200 * (vals / vals.max()) + 50 if vals.max() > 0 else np.full_like(vals, 50) + 0.0

    pos = nx.spring_layout(G, seed=42, weight="weight", k=None)
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color=vals, cmap=plt.cm.viridis)
    top_nodes = list(scores_df.head(30)["gene"])
    labels = {n: n for n in top_nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    plt.title("Red STRING con puntuación de propagación (RWR)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=250)
    plt.close()
    logging.info("Gráfico de red guardado en %s", out_png)

# ----------------------------- DIAMOnD / GUILD (opcionales) ---------------------

def diamond_algorithm(G: nx.Graph, seed_genes: Set[str], num_nodes: int = 10):
    added_nodes, seed_set = [], set(seed_genes)
    candidate_nodes = set(G.nodes()) - seed_set
    for _ in range(num_nodes):
        best_node, best_p = None, 1
        for node in candidate_nodes:
            neighbors = set(G.neighbors(node))
            k_s = len(neighbors & seed_set)
            k = len(neighbors)
            K = len(seed_set)
            N = len(G)
            p = hypergeom.sf(k_s - 1, N, K, k)
            if p < best_p:
                best_node, best_p = node, p
        if best_node is None:
            break
        added_nodes.append((best_node, best_p))
        seed_set.add(best_node)
        candidate_nodes.remove(best_node)
    return added_nodes


def guild_netscore(G: nx.Graph, seed_genes: Set[str], max_iter: int = 5):
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


# ----------------------------- CLI principal ------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Análisis funcional con propagación en redes (STRING + RWR + Enrichr)."
    )
    parser.add_argument("--genes-file", type=str, default="data/genes_input.txt",
                        help="Ruta al archivo de genes (una ID por línea, p. ej. símbolos HGNC).")
    parser.add_argument("--outdir", type=str, default="results",
                        help="Directorio de salida para resultados.")
    parser.add_argument("--species", type=int, default=9606,
                        help="TaxID NCBI para STRING (9606=Human).")
    parser.add_argument("--string-score", type=int, default=700,
                        help="Puntuación mínima de STRING (400/700/900 ~ med/alta/máxima).")
    parser.add_argument("--neighbors", type=int, default=50,
                        help="Nº de vecinos a añadir desde STRING (first shell).")
    parser.add_argument("--alpha", type=float, default=0.5,
                        help="Parámetro de reinicio (restart) para RWR (0<alpha<1).")
    parser.add_argument("--top-k", type=int, default=200,
                        help="Top genes por score para enriquecimiento/tabla.")
    parser.add_argument("--no-enrich", action="store_true",
                        help="Si se indica, omite enriquecimiento funcional.")
    parser.add_argument("--strict-symbols", action="store_true",
                        help="Mantiene solo matches cuyo preferredName coincide con el símbolo de entrada (case-insensitive).")

    # Opciones de visualización (mini-mejora)
    parser.add_argument("--plot-top", type=int, default=15,
                        help="Nº de términos a mostrar en el barplot GO.")
    parser.add_argument("--plot-width", type=float, default=16.0,
                        help="Ancho de figura para barplots.")
    parser.add_argument("--plot-font", type=int, default=9,
                        help="Tamaño de fuente en barplots.")

    # Opcionales extra (bonus)
    parser.add_argument("--with-diamond", action="store_true",
                        help="Ejecuta DIAMOnD y guarda diamond_results.csv")
    parser.add_argument("--with-guild", action="store_true",
                        help="Ejecuta NetScore (GUILD) y guarda guild_results.csv")
    parser.add_argument("--combo-top", type=int, default=50,
                        help="Top para combinar DIAMOnD+GUILD con RWR (para enriquecimiento combinado)")

    args = parser.parse_args()

    # --- preparar carpetas de salida organizadas ---
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    DIR_CSV     = os.path.join(outdir, "csv")
    DIR_FIG     = os.path.join(outdir, "figures")
    DIR_RWR     = os.path.join(outdir, "rwr")
    DIR_GRAPHS  = os.path.join(outdir, "graphs")
    for d in (DIR_CSV, DIR_FIG, DIR_RWR, DIR_GRAPHS):
        os.makedirs(d, exist_ok=True)

    setup_logging(outdir)
    save_params(outdir, args)

    # 1) Cargar genes
    genes_in = read_gene_list(args.genes_file)
    logging.info("Genes de entrada (%d): %s", len(genes_in), ", ".join(genes_in))

    # 2) Mapear a IDs de STRING
    logging.info("Mapeando a IDs de STRING...")
    map_df = string_map_genes_to_ids(genes_in, species=args.species)

    # Guardar mapeos
    map_raw_csv  = os.path.join(DIR_CSV, "mapping_string_ids_raw.csv")
    map_clean_csv = os.path.join(DIR_CSV, "mapping_string_ids.csv")
    map_df.to_csv(map_raw_csv, index=False)

    # Avisar de desajustes símbolo→preferredName
    try:
        cmp_df = map_df.dropna(subset=["stringId"]).copy()
        _q = cmp_df["queryItem"].astype(str).str.upper().str.strip()
        _p = cmp_df["preferredName"].astype(str).str.upper().str.strip()
        mismatch = cmp_df[_q != _p]
        if not mismatch.empty:
            logging.warning("Se detectaron mapeos no idénticos símbolo→preferredName. Revisa %s", map_clean_csv)
    except Exception:
        pass

    mapped = map_df.dropna(subset=["stringId"]).copy()
    if args.strict_symbols:
        mapped["__q"] = mapped["queryItem"].astype(str).str.upper().str.strip()
        mapped["__p"] = mapped["preferredName"].astype(str).str.upper().str.strip()
        before = len(mapped)
        mapped = mapped[mapped["_q"] == mapped["p"]].drop(columns=["q", "_p"])
        after = len(mapped)
        if after < before:
            logging.warning("STRICT: %d mapeos descartados por no coincidir símbolo-preferredName. Conservados: %d",
                            before - after, after)

    if mapped.empty:
        raise RuntimeError("No se obtuvieron IDs de STRING para las semillas tras el mapeo/filtrado.")

    mapped.to_csv(map_clean_csv, index=False)
    logging.info("Guardado mapeo filtrado en %s", map_clean_csv)

    seed_names = set(mapped["preferredName"].astype(str).tolist())
    string_ids = mapped["stringId"].astype(str).tolist()
    logging.info("Semillas mapeadas en STRING (%d): %s", len(seed_names), ", ".join(sorted(seed_names)))

    # 3) Descargar red de STRING
    logging.info("Descargando red de STRING con score>=%d y %d vecinos...", args.string_score, args.neighbors)
    edges_df = string_fetch_network(
        string_ids=string_ids,
        species=args.species,
        required_score=args.string_score,
        add_neighbors=args.neighbors
    )
    edges_csv = os.path.join(DIR_CSV, "string_network_edges.csv")
    edges_tsv = os.path.join(DIR_CSV, "string_network_edges.tsv")
    edges_df.to_csv(edges_csv, index=False)
    edges_df.to_csv(edges_tsv, sep="\t", index=False)  # copia TSV para traza
    logging.info("Red: %d aristas entre %d genes. Guardada en %s",
                 len(edges_df), len(pd.unique(edges_df[["preferredName_A", "preferredName_B"]].values.ravel())), edges_csv)

    # 3b) Guardar edgelist con pesos
    G = build_graph_from_edges(edges_df)
    edgelist_path = os.path.join(DIR_GRAPHS, "network_weighted.edgelist")
    nx.write_weighted_edgelist(G, edgelist_path)
    logging.info("Edgelist (ponderada) guardada en %s", edgelist_path)

    # 4) RWR
    logging.info("Grafo con %d nodos y %d aristas.", G.number_of_nodes(), G.number_of_edges())
    scores_df = random_walk_with_restart(G, seeds=seed_names, alpha=args.alpha)
    scores_csv = os.path.join(DIR_RWR, "propagation_scores.csv")
    scores_df.to_csv(scores_csv, index=False)
    logging.info("Puntuaciones de propagación guardadas en %s", scores_csv)

    # 5) Top genes para enriquecimiento
    top_k = min(args.top_k, len(scores_df))
    top_df = scores_df.head(top_k).copy()
    top_csv = os.path.join(DIR_RWR, f"top_genes_{top_k}.csv")
    top_df.to_csv(top_csv, index=False)
    logging.info("Top %d genes guardados en %s", top_k, top_csv)

    # 6) Enriquecimiento base (RWR)
    saved_enrich: Dict[str, str] = {}
    if not args.no_enrich:
        try:
            saved_enrich = run_enrichment(
                genes=top_df["gene"].tolist(),
                out_csv_dir=DIR_CSV,
                prefix=f"top{top_k}"
            )
        except Exception as e:
            logging.exception("Error en enriquecimiento funcional: %s", e)

    # 7) Visualizaciones
    # 7a) GO BP (RWR)
    if "GO_Biological_Process_2021" in saved_enrich:
        go_png = os.path.join(DIR_FIG, "go_bp_barplot.png")
        go_svg = os.path.join(DIR_FIG, "go_bp_barplot.svg")
        plot_go_bar(
            saved_enrich["GO_Biological_Process_2021"],
            out_png=go_png,
            out_svg=go_svg,
            top_n=args.plot_top,
            fig_width=args.plot_width,
            font_size=args.plot_font
        )

    # 7b) Grafo con scores
    net_png = os.path.join(DIR_FIG, "network_scores.png")
    plot_network_scores(G, scores_df, net_png)

    # 8) Opcional: DIAMOnD/GUILD + enriquecimiento combinado
    if args.with_diamond or args.with_guild:
        combo_genes = set(top_df["gene"].tolist())

        if args.with_diamond:
            dres = diamond_algorithm(G, seed_names, num_nodes=args.combo_top)
            d_csv = os.path.join(DIR_CSV, "diamond_results.csv")
            pd.DataFrame(dres, columns=["node", "p_value"]).to_csv(d_csv, index=False)
            logging.info("DIAMOnD guardado en %s", d_csv)
            try:
                d = pd.read_csv(d_csv)
                combo_genes |= set(d.nsmallest(args.combo_top, "p_value")["node"].astype(str))
            except Exception:
                pass

        if args.with_guild:
            gscores = guild_netscore(G, seed_names)
            gdf = pd.DataFrame(sorted(gscores.items(), key=lambda x: x[1], reverse=True),
                               columns=["node", "score"])
            g_csv = os.path.join(DIR_CSV, "guild_results.csv")
            gdf.head(args.combo_top).to_csv(g_csv, index=False)
            logging.info("GUILD guardado en %s", g_csv)
            try:
                g = pd.read_csv(g_csv)
                combo_genes |= set(g.nlargest(args.combo_top, "score")["node"].astype(str))
            except Exception:
                pass

        if not args.no_enrich and len(combo_genes) > 0:
            saved_combo = run_enrichment(sorted(combo_genes), DIR_CSV, prefix=f"combo_top{len(combo_genes)}")
            if "GO_Biological_Process_2021" in saved_combo:
                go_combo_png = os.path.join(DIR_FIG, "go_bp_barplot_combo.png")
                go_combo_svg = os.path.join(DIR_FIG, "go_bp_barplot_combo.svg")
                plot_go_bar(
                    saved_combo["GO_Biological_Process_2021"],
                    out_png=go_combo_png,
                    out_svg=go_combo_svg,
                    top_n=args.plot_top,
                    fig_width=args.plot_width,
                    font_size=args.plot_font
                )

    # 9) README de resultados (referenciando subcarpetas)
    results_md = os.path.join(outdir, "README_results.md")
    with open(results_md, "w", encoding="utf-8") as f:
        f.write(f"# Resultados – FuncNet Pipeline\n\n")
        f.write(f"*Fecha:* {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Estructura de salidas\n")
        f.write("- csv/: mapeos, red, enriquecimientos, etc.\n")
        f.write("- figures/: gráficos en PNG y SVG (GO, red).\n")
        f.write("- rwr/: puntuaciones de propagación y top genes.\n")
        f.write("- graphs/: representaciones de la red (edgelist).\n\n")
        f.write("## Archivos clave\n")
        f.write(f"- csv/mapping_string_ids_raw.csv, csv/mapping_string_ids.csv\n")
        f.write(f"- csv/string_network_edges.csv (+ string_network_edges.tsv)\n")
        f.write(f"- graphs/network_weighted.edgelist\n")
        f.write(f"- rwr/propagation_scores.csv, rwr/top_genes_{top_k}.csv\n")
        f.write(f"- figures/go_bp_barplot.png|.svg, figures/network_scores.png\n")
        f.write(f"- pipeline.log, run_metadata.json\n")
        f.write("(si activa bonus) csv/diamond_results.csv, csv/guild_results.csv, "
                "figures/go_bp_barplot_combo.png|.svg\n\n")
        f.write("## Parámetros\njson\n")
        json.dump(vars(args), f, indent=2)
        f.write("\n\n")

    logging.info("Pipeline completado. Revisa %s", outdir)


if __name__ == "__main__":
    main()


# Get the project root directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, 'scripts')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')
NETWORK_SCRIPT = os.path.join(SCRIPTS_DIR, 'network_propagation_and_functional_analysis.py')

# Get virtual environment python
VENV_PYTHON = os.path.join(PROJECT_ROOT, 'venv', 'bin', 'python')
if not os.path.exists(VENV_PYTHON):
    VENV_PYTHON = sys.executable

def get_demo_args():
    """Returns demo arguments for testing"""
    # Look for STRING network file in data directory
    string_files = [f for f in os.listdir(DATA_DIR) if f.endswith('.tsv') and 'contrast' in f.lower()]
    
    if not string_files:
        print(f"Error: No network file found in {DATA_DIR}")
        print("Expected a file with pattern 'contrast.tsv' or similar")
        return None
    
    string_input = os.path.join(DATA_DIR, string_files[0])
    seed_genes = os.path.join(DATA_DIR, 'genes_input.txt')
    
    if not os.path.exists(seed_genes):
        print(f"Error: Seed genes file not found at {seed_genes}")
        return None
    
    return {
        'string_input': string_input,
        'seed_genes': seed_genes,
        'output_dir': RESULTS_DIR,
        'top_n': 10
    }

def run_analysis(args):
    """Run the network analysis script with given arguments"""
    cmd = [
        VENV_PYTHON,
        NETWORK_SCRIPT,
        '--string_input', args['string_input'],
        '--seed_genes', args['seed_genes'],
        '--output_dir', args['output_dir'],
        '--top_n', str(args['top_n'])
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    print("-" * 80)
    
    try:
        result = subprocess.run(cmd, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error running analysis: {e}")
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}")
        return 1

def main():
    parser = argparse.ArgumentParser(
        description="HAB Project CLI - Network Propagation and Functional Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run demo with default files
  python hab_cli.py --demo --input-genes data/genes_input.txt
  
  # Run with custom inputs
  python hab_cli.py --string-input <path> --seed-genes <path> --output-dir <path> --top-n 15
  
  # Show help
  python hab_cli.py --help
        """
    )
    
    # Demo mode
    parser.add_argument(
        '--demo',
        action='store_true',
        help='Run in demo mode using default data files'
    )
    
    # Custom arguments
    parser.add_argument(
        '--string-input',
        '--string_input',
        type=str,
        help='Path to STRING network file (protein links)'
    )
    parser.add_argument(
        '--seed-genes',
        '--seed_genes',
        type=str,
        help='Path to seed genes file (HUGO symbols)'
    )
    parser.add_argument(
        '--input-genes',
        '--input_genes',
        type=str,
        help='Alias for --seed-genes'
    )
    parser.add_argument(
        '--output-dir',
        '--output_dir',
        type=str,
        default=RESULTS_DIR,
        help=f'Output directory for results (default: {RESULTS_DIR})'
    )
    parser.add_argument(
        '--top-n',
        '--top_n',
        type=int,
        default=10,
        help='Number of top nodes for functional analysis (default: 10)'
    )
    
    args = parser.parse_args()
    
    # Handle demo mode
    if args.demo:
        print("Running in DEMO mode...")
        demo_args = get_demo_args()
        if demo_args is None:
            return 1
        
        # Override with command-line input-genes if provided
        if args.input_genes:
            demo_args['seed_genes'] = args.input_genes
        if args.seed_genes:
            demo_args['seed_genes'] = args.seed_genes
        if args.string_input:
            demo_args['string_input'] = args.string_input
        if args.top_n != 10:
            demo_args['top_n'] = args.top_n
        
        print(f"Using STRING network: {demo_args['string_input']}")
        print(f"Using seed genes: {demo_args['seed_genes']}")
        print(f"Output directory: {demo_args['output_dir']}")
        print(f"Top N genes: {demo_args['top_n']}")
        print()
        
        return run_analysis(demo_args)
    
    # Handle custom mode
    # Resolve aliases
    seed_genes = args.seed_genes or args.input_genes
    string_input = args.string_input
    
    if not string_input or not seed_genes:
        parser.print_help()
        print("\nError: --string-input and --seed-genes (or --input-genes) are required unless using --demo")
        return 1
    
    custom_args = {
        'string_input': string_input,
        'seed_genes': seed_genes,
        'output_dir': args.output_dir,
        'top_n': args.top_n
    }
    
    return run_analysis(custom_args)

if _name_ == '_main_':
    sys.exit(main())