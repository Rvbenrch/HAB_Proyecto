<details>
<summary>Enunciado del Projecto<summary>

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

 