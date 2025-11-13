# Resultados – FuncNet Pipeline

**Fecha:** 2025-11-11 20:22:04

## Estructura de salidas
- `csv/`: mapeos, red, enriquecimientos, etc.
- `figures/`: gráficos en PNG y SVG (GO, red).
- `rwr/`: puntuaciones de propagación y top genes.
- `graphs/`: representaciones de la red (edgelist).

## Archivos clave
- `csv/mapping_string_ids_raw.csv`, `csv/mapping_string_ids.csv`
- `csv/string_network_edges.csv` (+ `string_network_edges.tsv`)
- `graphs/network_weighted.edgelist`
- `rwr/propagation_scores.csv`, `rwr/top_genes_50.csv`
- `figures/go_bp_barplot.png|.svg`, `figures/network_scores.png`
- `pipeline.log`, `run_metadata.json`
(si activa bonus) `csv/diamond_results.csv`, `csv/guild_results.csv`, `figures/go_bp_barplot_combo.png|.svg`

## Parámetros
```json
{
  "genes_file": "data/genes_input.txt",
  "outdir": "results",
  "species": 9606,
  "string_score": 700,
  "neighbors": 50,
  "alpha": 0.5,
  "top_k": 50,
  "no_enrich": false,
  "strict_symbols": false,
  "plot_top": 15,
  "plot_width": 16.0,
  "plot_font": 9,
  "with_diamond": true,
  "with_guild": true,
  "combo_top": 50
}
```
