# Proyecto de Análisis Funcional y Propagación en Redes (HAB / FuncNet)

En este repositorio hemos montado un **pipeline reproducible** para analizar listas de genes mediante:

- **Propagación en redes** (Random Walk with Restart, RWR) sobre **STRING**
- **Análisis funcional** con **Enrichr** (GO Biological Process y KEGG)
- Opcionalmente, **DIAMOnD** y **GUILD/NetScore** como métodos adicionales de priorización

El objetivo es que, dado un fichero con genes de interés, puedas ejecutar el script y obtener
tablas, gráficos y un resumen de parámetros listo para incluir en la memoria.

---

## Estructura del repositorio

```text
├── data/                        # Entradas (lista de genes, etc.)
│   └── genes_input.txt          # Lista de genes semilla (símbolos HUGO/HGNC)
├── results/                     # Salidas generadas al ejecutar
│   ├── csv/                     # Mapeos, red STRING, enriquecimientos…
│   ├── figures/                 # Gráficos (PNG, SVG)
│   ├── graphs/                  # Red (weighted edgelist)
│   └── rwr/                     # Scores de propagación y top genes
├── scripts/
│   ├── funcnet_pipeline.py      # Pipeline principal (RWR + Enrichr)
│   └── hab_cli.py               # CLI auxiliar
├── docs/                        # Documentación adicional
├── README.md                    # Este documento
└── requirements.txt             # Dependencias del proyecto
```

---

## Objetivos del proyecto

El proyecto implementa un **pipeline automatizado** que:

1. Mapea símbolos HUGO/HGNC → IDs de STRING.
2. Descarga la red PPI de STRING con filtros (score y vecinos).
3. Ejecuta **Random Walk with Restart (RWR)** desde genes semilla.
4. Ordena nodos por score de propagación.
5. Ejecuta **enriquecimiento funcional** (GO BP y KEGG, vía Enrichr).
6. Genera figuras y tablas listas para informe (PNG, SVG, CSV).
7. Documenta todo (logs, metadatos, README de resultados).
8. Opcionalmente ejecuta **DIAMOnD** y **GUILD/NetScore** con enriquecimiento combinado.

---

## Qué hace el pipeline (resumen)

1. **Lee** `data/genes_input.txt`.
2. **Mapea** genes → STRING IDs.
3. **Descarga** red PPI filtrada.
4. **Construye** el grafo y corre **RWR**.
5. **Selecciona** el *top-k* de genes más relevantes.
6. **Lanza Enrichr** (GO BP y KEGG).
7. **Genera gráficos**:

   * Barplot de GO BP
   * Red coloreada por score
8. *(Opcional)* Ejecuta DIAMOnD / GUILD.
9. Guarda logs, metadatos y un `README_results.md`.

---

## Metodología (breve)

* **STRING**: red PPI con evidencia integrada.
* **RWR**: difusión con reinicio periódico hacia las semillas.
* **Enrichr**: enriquecimiento funcional; aquí GO BP + KEGG.
* **DIAMOnD / GUILD (extra)**: estrategias alternativas de priorización.

---

## Requisitos e instalación

Dependencias principales:

```text
requests==2.31.0
pandas==2.2.2
numpy==1.26.4
networkx==3.2.1
matplotlib==3.8.4
scipy==1.11.4
```

# Instalación rápida:

## Crear entorno virtual
```bash
python -m venv venv
```
## Activarlo:
### Linux/macOS:
```bash
source venv/bin/activate
```
### Windows (PowerShell):
```bash
 .\venv\Scripts\activate
````

## Instalar dependencias
```bash
pip install -r requirements.txt
```

> Nota: se requiere conexión a Internet (STRING + Enrichr).

---

## Ejecución rápida

### Linux / macOS

```bash
python3 scripts/funcnet_pipeline.py \
  --genes-file data/genes_input.txt \
  --outdir results \
  --top-k 50 \
  --with-diamond --with-guild \
  --strict-symbols
```

### Windows

```bash
python scripts\funcnet_pipeline.py ^
  --genes-file data\genes_input.txt ^
  --outdir results ^
  --top-k 50 ^
  --with-diamond --with-guild ^
  --strict-symbols
---
```
## Uso detallado (`funcnet_pipeline.py`)

Ejemplo general:

```bash
python scripts/funcnet_pipeline.py [opciones]
```

## Parámetros principales:

* `--genes-file PATH` – Archivo con genes de entrada.
* `--outdir PATH` – Carpeta de salida.
* `--species INT` – TaxID para STRING (9606 = humano).
* `--string-score INT` – Score mínimo (400/700/900).
* `--neighbors INT` – Vecinos a añadir desde STRING.
* `--alpha FLOAT` – Parámetro RWR.
* `--top-k INT` – Nº de genes para Enrichr.
* `--no-enrich` – Omite análisis funcional.
* `--strict-symbols` – Solo mapeos exactos símbolo→preferredName.

### Visualización:

* `--plot-top INT`
* `--plot-width FLOAT`
* `--plot-font INT`

### Ejemplos útiles:

#### Ejecución básica
```bash
python scripts/funcnet_pipeline.py \
  --genes-file data/genes_input.txt \
  --outdir results
```
## Parámetros más estrictos
```bash
python scripts/funcnet_pipeline.py --string-score 900 --neighbors 100
```
---

## CLI auxiliar: `hab_cli.py`

### Modo demo

```bash
python scripts/hab_cli.py --demo
```

##### Con override de genes:

```bash
python scripts/hab_cli.py --demo --input-genes data/mis_genes.txt
````

### Modo personalizado

```bash
python scripts/hab_cli.py \
  --seed-genes data/genes_input.txt \
  --output-dir results \
  --top-n 15
```

---

## Entradas y salidas

### Entrada mínima

`data/genes_input.txt`:

```text
TP53
BRCA1
EGFR
```

### Salidas (en `results/`)

* **csv/** – Mapeos, red STRING, enriquecimientos, bonus.
* **rwr/** – Puntuaciones RWR y top genes.
* **figures/** – GO BP (PNG/SVG), red.
* **graphs/** – Weighted edgelist.
* `pipeline.log`, `run_metadata.json`, `README_results.md`.

---

## Reproducibilidad

* Parámetros guardados en `run_metadata.json`.
* Logging completo en `pipeline.log`.
* Dependencias fijadas en `requirements.txt`.
* Red y mapeos exportados en CSV para trazabilidad.
---

## Resumen de uso rápido

```bash
# 1) Activar entorno
source venv/bin/activate

# 2) Instalar dependencias
pip install -r requirements.txt

# 3) Ejecutar pipeline
python scripts/funcnet_pipeline.py \
  --genes-file data/genes_input.txt \
  --outdir results \
  --string-score 700 \
  --neighbors 50 \
  --alpha 0.5 \
  --top-k 200
``` 
---
