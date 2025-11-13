#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HAB CLI - Command Line Interface for the FuncNet pipeline.

Uso rápido:
  # Demo (usa data/genes_input.txt y results/)
  python scripts/hab_cli.py --demo

  # Personalizado
  python scripts/hab_cli.py \
    --genes-file data/genes_input.txt \
    --outdir results \
    --top-k 200 \
    --with-diamond --with-guild --combo-top 50
"""

import os
import sys
import argparse
import subprocess

# Rutas del proyecto
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR     = os.path.join(PROJECT_ROOT, "data")
SCRIPTS_DIR  = os.path.join(PROJECT_ROOT, "scripts")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results")

# Script principal (pipeline)
PIPELINE_SCRIPT = os.path.join(SCRIPTS_DIR, "funcnet_pipeline.py")

# Python del venv si existe
VENV_PYTHON = os.path.join(PROJECT_ROOT, "venv", "bin", "python")
if not os.path.exists(VENV_PYTHON):
    VENV_PYTHON = sys.executable  # usa el python actual

def get_demo_args():
    """Argumentos por defecto para correr un demo mínimo."""
    genes_file = os.path.join(DATA_DIR, "genes_input.txt")
    if not os.path.exists(genes_file):
        print(f"[ERROR] No se encuentra {genes_file}")
        return None
    return {
        "genes_file": genes_file,
        "outdir": RESULTS_DIR,
        "top_k": 50,
        # Puedes preactivar bonus aquí si quieres:
        # "with_diamond": True,
        # "with_guild": True,
        # "combo_top": 50,
    }

def run_pipeline(args_dict) -> int:
    """Construye el comando y llama al pipeline."""
    cmd = [
        VENV_PYTHON, PIPELINE_SCRIPT,
        "--genes-file", args_dict["genes_file"],
        "--outdir", args_dict.get("outdir", RESULTS_DIR),
        "--top-k", str(args_dict.get("top_k", 50)),
    ]

    # Parámetros opcionales que el CLI puede pasar al pipeline
    if args_dict.get("species"):       cmd += ["--species", str(args_dict["species"])]
    if args_dict.get("string_score"):  cmd += ["--string-score", str(args_dict["string_score"])]
    if args_dict.get("neighbors"):     cmd += ["--neighbors", str(args_dict["neighbors"])]
    if args_dict.get("alpha"):         cmd += ["--alpha", str(args_dict["alpha"])]
    if args_dict.get("no_enrich"):     cmd += ["--no-enrich"]
    if args_dict.get("strict_symbols"):cmd += ["--strict-symbols"]

    # Visualización
    if args_dict.get("plot_top"):   cmd += ["--plot-top", str(args_dict["plot_top"])]
    if args_dict.get("plot_width"): cmd += ["--plot-width", str(args_dict["plot_width"])]
    if args_dict.get("plot_font"):  cmd += ["--plot-font", str(args_dict["plot_font"])]

    # Bonus
    if args_dict.get("with_diamond"): cmd += ["--with-diamond"]
    if args_dict.get("with_guild"):   cmd += ["--with-guild"]
    if args_dict.get("combo_top"):    cmd += ["--combo-top", str(args_dict["combo_top"])]

    print(">> Ejecutando:", " ".join(cmd))
    try:
        return subprocess.run(cmd, check=True).returncode
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Pipeline devolvió código {e.returncode}")
        return e.returncode

def main() -> int:
    p = argparse.ArgumentParser(
        description="HAB Project CLI - Lanzador del FuncNet pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--demo", action="store_true", help="Ejecuta con datos por defecto")
    p.add_argument("--genes-file", "--genes_file", dest="genes_file", type=str,
                   help="Ruta a archivo con genes (uno por línea)")
    p.add_argument("--outdir", type=str, default=RESULTS_DIR, help="Directorio de salida")
    p.add_argument("--top-k", "--top_k", dest="top_k", type=int, default=50,
                   help="Top de genes para enriquecimiento")

    # Parámetros adicionales que propagamos al pipeline
    p.add_argument("--species", type=int)
    p.add_argument("--string-score", "--string_score", dest="string_score", type=int)
    p.add_argument("--neighbors", type=int)
    p.add_argument("--alpha", type=float)
    p.add_argument("--no-enrich", "--no_enrich", dest="no_enrich", action="store_true")
    p.add_argument("--strict-symbols", "--strict_symbols", dest="strict_symbols", action="store_true")
    p.add_argument("--plot-top", "--plot_top", dest="plot_top", type=int)
    p.add_argument("--plot-width", "--plot_width", dest="plot_width", type=float)
    p.add_argument("--plot-font", "--plot_font", dest="plot_font", type=int)

    # Bonus
    p.add_argument("--with-diamond", "--with_diamond", dest="with_diamond", action="store_true")
    p.add_argument("--with-guild", "--with_guild", dest="with_guild", action="store_true")
    p.add_argument("--combo-top", "--combo_top", dest="combo_top", type=int)

    args = p.parse_args()

    # DEMO
    if args.demo:
        demo = get_demo_args()
        if not demo:
            return 1
        # Permite sobreescribir algunos parámetros del demo
        if args.genes_file: demo["genes_file"] = args.genes_file
        demo["outdir"] = args.outdir
        if args.top_k:     demo["top_k"] = args.top_k
        # Pasa cualquier otro parámetro si quieres
        for k in ("species","string_score","neighbors","alpha","no_enrich","strict_symbols",
                  "plot_top","plot_width","plot_font","with_diamond","with_guild","combo_top"):
            v = getattr(args, k, None)
            if v is not None:
                demo[k] = v
        return run_pipeline(demo)

    # CUSTOM: requiere --genes-file
    if not args.genes_file:
        p.print_help()
        print("\n[ERROR] Debes indicar --genes-file o usar --demo")
        return 1

    cfg = vars(args)
    return run_pipeline(cfg)

if __name__ == "__main__":
    sys.exit(main())
