#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HAB CLI - Command Line Interface for HAB Project (Network Analysis)

Envuelve el script scripts/funcnet_pipeline.py
"""

import os
import sys
import argparse
import subprocess

# Directorios del proyecto
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
SCRIPTS_DIR = os.path.join(PROJECT_ROOT, "scripts")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results")

# Script de análisis de red (nuestro pipeline)
NETWORK_SCRIPT = os.path.join(SCRIPTS_DIR, "funcnet_pipeline.py")

# Python del entorno virtual (si existe)
VENV_PYTHON = os.path.join(PROJECT_ROOT, "venv", "bin", "python")
if not os.path.exists(VENV_PYTHON):
    VENV_PYTHON = sys.executable


def get_demo_args():
    """Devuelve argumentos por defecto para modo demo."""
    seed_genes = os.path.join(DATA_DIR, "genes_input.txt")
    if not os.path.exists(seed_genes):
        print(f"Error: Seed genes file not found at {seed_genes}")
        return None

    return {
        "seed_genes": seed_genes,
        "output_dir": RESULTS_DIR,
        "top_n": 10,
    }


def run_analysis(args):
    """Lanza el pipeline con los argumentos dados."""
    cmd = [
        VENV_PYTHON,
        NETWORK_SCRIPT,
        "--genes-file", args["seed_genes"],
        "--outdir", args["output_dir"],
        "--top-k", str(args["top_n"]),
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
  python hab_cli.py --demo
  
  # Run with custom input gene list
  python hab_cli.py --input-genes data/genes_input.txt --output-dir results --top-n 50
  
  # Show help
  python hab_cli.py --help
        """,
    )

    # Demo mode
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run in demo mode using default data files",
    )

    # Custom arguments
    parser.add_argument(
        "--string-input",
        "--string_input",
        type=str,
        help="(No usado por el pipeline actual, solo por compatibilidad)",
    )
    parser.add_argument(
        "--seed-genes",
        "--seed_genes",
        type=str,
        help="Ruta al archivo de genes semilla (HUGO symbols)",
    )
    parser.add_argument(
        "--input-genes",
        "--input_genes",
        type=str,
        help="Alias para --seed-genes",
    )
    parser.add_argument(
        "--output-dir",
        "--output_dir",
        type=str,
        default=RESULTS_DIR,
        help=f"Directorio de salida para resultados (default: {RESULTS_DIR})",
    )
    parser.add_argument(
        "--top-n",
        "--top_n",
        type=int,
        default=10,
        help="Número de top genes a considerar en el pipeline (default: 10)",
    )

    args = parser.parse_args()

    # Modo demo
    if args.demo:
        print("Running in DEMO mode...")
        demo_args = get_demo_args()
        if demo_args is None:
            return 1

        # Override con argumentos de línea de comandos si se dan
        if args.input_genes:
            demo_args["seed_genes"] = args.input_genes
        if args.seed_genes:
            demo_args["seed_genes"] = args.seed_genes
        if args.output_dir:
            demo_args["output_dir"] = args.output_dir
        if args.top_n != 10:
            demo_args["top_n"] = args.top_n

        print(f"Using seed genes: {demo_args['seed_genes']}")
        print(f"Output directory: {demo_args['output_dir']}")
        print(f"Top N genes: {demo_args['top_n']}")
        print()

        return run_analysis(demo_args)

    # Modo custom
    seed_genes = args.seed_genes or args.input_genes

    if not seed_genes:
        parser.print_help()
        print("\nError: necesitas --seed-genes o --input-genes (o usa --demo)")
        return 1

    custom_args = {
        "seed_genes": seed_genes,
        "output_dir": args.output_dir,
        "top_n": args.top_n,
    }

    return run_analysis(custom_args)


if __name__ == "__main__":
    sys.exit(main())
