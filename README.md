<details>
<summary>Enunciado del Projecto</summary>

# Proyecto de Análisis Funcional y Propagación en Redes

Este repositorio contiene la plantilla base para el proyecto de análisis funcional.

## Objetivo
Desarrollar un script en Python que realice análisis funcional con propagación en redes a partir de una lista de genes diferencialmente expresados.

## Estructura del repositorio

```
├── data/                # Archivos de entrada (lista de genes, redes, etc.)
├── results/             # Resultados generados por el script
├── scripts/             # Código fuente del proyecto
├── docs/                # Documentación adicional (opcional)
├── README.md            # Este archivo
└── requirements.txt     # Dependencias del proyecto
```

## Instrucciones
1. Haz un fork de este repositorio.
2. Trabaja en tu fork con tu grupo.
3. Implementa el análisis funcional con propagación en redes.
4. Documenta tu código y resultados.
5. Sube tu proyecto a GitHub y asegúrate de que sea reproducible.

## Rúbrica de Evaluación

| Criterio | Descripción | Puntos |
|---------|-------------|--------|
| **1. Funcionalidad del script** | Correcta ejecución del análisis funcional y/o propagación. | 25 |
| **2. Elección y justificación de técnicas** | Adecuación y justificación de los métodos usados. | 15 |
| **3. Automatización y flujo de trabajo** | Entrada por CLI, conversión de IDs, descarga de datos, etc. | 15 |
| **4. Documentación y reproducibilidad** | README claro, dependencias, ejemplos. | 15 |
| **5. Calidad del código** | Estilo, modularidad, comentarios. | 10 |
| **6. Análisis y visualización de resultados** | Salidas interpretables, gráficos, tablas. | 10 |
| **Bonus** | Originalidad, integración de datos, visualizaciones interactivas. | +10 |

## Evaluación individual
Se tendrá en cuenta la participación de cada miembro del grupo a través del historial de commits en GitHub.
<\details>

# Proyecto de Análisis Funcional y Propagación en Redes (FuncNet)

Este repositorio implementa un flujo reproducible para análisis funcional (Enrichr) con propagación en redes (Random Walk with Restart, RWR) sobre genes de entrada. Opcionalmente integra DIAMOnD y GUILD/NetScore para bonus.

<pre> ```python 
├── data/                        # Entradas (lista de genes, etc.)
├── results/                     # Salidas (se generan al ejecutar)
│   ├── csv/                     # mapeos, red STRING, enriquecimientos…
│   ├── figures/                 # gráficos (PNG, SVG)
│   ├── graphs/                  # red exportada (edgelist)
│   └── rwr/                     # scores de propagación y top genes
├── scripts/
│   ├── funcnet_pipeline.py      # pipeline principal (RWR + Enrichr)
│   └── hab_cli.py               # CLI (lanzador práctico del pipeline)
├── README.md                    # (este documento)
└── requirements.txt             # dependencias
 ``` </pre>
</details>

# Proyecto de Análisis Funcional y Propagación en Redes (HAB)

Pipeline reproducible para analizar listas de genes mediante propagación en redes (Random Walk with Restart, RWR) sobre STRING y posterior enriquecimiento funcional con Enrichr (GO BP y KEGG). Incluye opciones bonus (DIAMOnD y GUILD/NetScore) y genera tablas y figuras listas para informe.

---

## Objetivos del proyecto

* Desarrollar un script en Python que:
  1. mapee símbolos HUGO/HGNC a IDs de STRING, 2) descargue la red PPI filtrada, 3) ejecute RWR desde genes semilla, 4) clasifique genes por proximidad funcional, 5) realice enriquecimiento funcional (GO BP, KEGG) y 6) produzca salidas estandarizadas (CSV, PNG, SVG, logs y metadatos).
* Proveer una CLI clara, reproducible y multiplataforma.
* (Opcional) Integrar DIAMOnD y/o GUILD (NetScore) para priorización adicional y enriquecimiento combinado.
* Cumplir los criterios de la rúbrica: funcionalidad, justificación técnica, automatización, documentación, calidad de código y visualización.

---

## Qué hace el pipeline

1. Lee una lista de genes (símbolos HUGO/HGNC).
2. Mapea símbolos → IDs de STRING (API).
3. Descarga la red de interacciones PPI desde STRING con filtros de score y vecinos.
4. Propaga la señal desde las semillas con RWR y rankea por puntuación.
5. Enriquece funcionalmente (GO: Biological Process y KEGG) con Enrichr.
6. Dibuja:
   * Barplot de GO BP (PNG + SVG) con etiquetas legibles.
   * Red coloreada por puntuación de propagación.
7. Opcional (bonus): ejecuta DIAMOnD y/o GUILD y permite enriquecimiento combinado.
8. Documenta la ejecución: logs, metadatos y README de resultados.

---

## Metodología (resumen)

* RWR: proceso de difusión con reinicio hacia las semillas; asigna scores globales de relevancia.
* STRING: red PPI curada que integra evidencia experimental, coexpresión y texto.
* Enrichr: enriquecimiento en GO BP y KEGG sobre el top de genes por RWR.
* DIAMOnD/GUILD (opcional): priorización basada en conectividad y propagación iterativa.

---

## Estructura del repositorio

data/ – Entradas (genes_input.txt)  
scripts/ – Código (funcnet_pipeline.py, hab_cli.py)  
results/ – Salidas generadas  
docs/ – Documentación adicional  
requirements.txt – Dependencias  

---

## Requisitos

requests==2.31.0  
pandas==2.2.2  
numpy==1.26.4  
networkx==3.2.1  
matplotlib==3.8.4  
scipy==1.11.4

Nota: se requiere conexión a Internet para STRING y Enrichr.

Instalación rápida:
python -m venv venv
source venv/bin/activate        # Linux/macOS
# .\venv\Scripts\activate   # Windows PowerShell
pip install -r requirements.txt
