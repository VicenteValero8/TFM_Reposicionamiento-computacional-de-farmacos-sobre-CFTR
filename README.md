# TFM_Reposicionamiento-computacional-de-farmacos-sobre-CFTR

Pipeline integrado de **PPI + docking molecular + QSAR** para la priorización de fármacos aprobados como candidatos a reposicionamiento sobre la proteína CFTR (fibrosis quística), con foco en la mutación ΔF508.

**Autor:** Vicente Valero Just
**Trabajo Final de Máster — Máster Universitario en Bioinformática y Bioestadística, 2025–2026**

---

## 1. Resumen del proyecto

El proyecto identifica y ordena fármacos ya aprobados que podrían unirse a CFTR, combinando tres capas de evidencia independientes y dos capas de caracterización farmacológica:

1. **Red PPI (STRING):** se construye la red de interacción proteína–proteína de CFTR y se identifican proteínas *hub*. Sobre las dianas con fármacos aprobados (DrugBank) se obtiene la lista de candidatos.
2. **Docking molecular (AutoDock Vina):** se calcula la afinidad de cada candidato por el bolsillo de unión, sobre la estructura experimental **6MSM** (CFTR WT) y sobre el modelo **AlphaFold2-ΔF508** (variante mutada).
3. **QSAR (Random Forest):** modelo de predicción de pIC50 entrenado sobre un dataset combinado de **ChEMBL + Papyrus + BindingDB** (1.082 compuestos únicos).
4. **Integración:** docking y pIC50 se normalizan y se combinan en un **score integrado** (60 % docking / 40 % QSAR), con análisis del **dominio de aplicabilidad** (Tanimoto).
5. **Caracterización complementaria:** **ADME** (SwissADME) y **perfil de interacciones** (PLIP) de los candidatos prioritarios.

**Resultado principal:** el bloque prioritario lo forman **Cromoglicic acid** (score 0,680) y **Glyburide** (0,655), separados del resto y por debajo de la referencia Ivacaftor (0,944).

---

## 2. Contenido del repositorio

```
.
├── CodigoTFM_VVJ.ipynb          # Notebook principal (dataset, QSAR, integración, ranking, figuras)
├── README.md
├── Data/
│   └── bindingdb_cftr.tsv       # Descarga manual de BindingDB (target CFTR)
├── Estructuras/
│   ├── Experimental/            # Estructura experimental CFTR-WT (PDB 6MSM)
│   └── AlphaFold/               # Modelo AlphaFold2-ΔF508
├── PLIP/                        # Resultados PLIP por candidato (exhaustividad 16)
│   ├── Glyburide_Results/
│   ├── Cromoglicic_acid_Results/
│   ├── Fostamatinib_Results/
│   ├── Desatinib_Results/
│   └── R406_Results/
├── Imagenes/                    # Figuras del trabajo
│   ├── fig_ranking_final
│   ├── fig_importancia_features
│   ├── fig_validacion_controles
│   ├── figura_plip_comparativa
│   └── STRING_CFTR
└── Outputs/                     # Datos y resultados en CSV
    ├── chembl_cftr_bruto.csv    # Descarga bruta de ChEMBL
    ├── dataset_completo.csv     # Dataset fusionado (ChEMBL + Papyrus + BindingDB)
    ├── comparacion_modelos.csv  # Modelo 1 (solo ChEMBL) vs Modelo 2 (ampliado)
    ├── ranking_final.csv        # Ranking integrado de candidatos
    └── validacion_controles.csv # Predicciones de candidatos y controles ±
```

> El notebook cubre los bloques de **dataset, QSAR, integración y visualización**. El **docking (AutoDock Vina)**, la **detección de bolsillos (PrankWeb)**, el **modelado AlphaFold2-ΔF508**, el **ADME (SwissADME)** y el **PLIP** se ejecutan con herramientas externas y sus resultados se incorporan al notebook como entradas.

---

## 3. Requisitos

- **Python 3.10+**
- Librerías:

```bash
pip install pandas numpy matplotlib scikit-learn rdkit chembl_webresource_client papyrus-scripts
```

| Librería | Uso |
|---|---|
| `pandas`, `numpy` | Manejo de datos |
| `rdkit` | Descriptores, Morgan fingerprints, Lipinski, Tanimoto |
| `scikit-learn` | Random Forest, validación cruzada, MinMaxScaler, métricas |
| `chembl_webresource_client` | Descarga de bioactividades de ChEMBL |
| `papyrus-scripts` | Descarga del dataset Papyrus |
| `matplotlib` | Figuras |

**Herramientas externas** (no Python): AutoDock Vina, AutoDockTools, PrankWeb, AlphaFold2, SwissADME, PLIP, STRING, DrugBank.

---

## 4. Datos de entrada

- **ChEMBL:** se descarga automáticamente desde el notebook (target `CHEMBL4051`, *Homo sapiens*).
- **Papyrus:** se descarga automáticamente con `papyrus-scripts` (v05.7, UniProt `P13569`).
- **BindingDB:** descarga **manual** del TSV de CFTR desde [bindingdb.org](https://www.bindingdb.org); colocar en `data/`.
- **Estructuras:** PDB **6MSM** (CFTR WT) y modelo **AlphaFold2-ΔF508** (generado a partir de la secuencia de la cadena A de 6MSM eliminando la Phe508).
- **Energías de docking:** calculadas externamente con AutoDock Vina (exhaustividad 8 para el cribado; 16 para el análisis PLIP de los candidatos prioritarios) y usadas como entrada del notebook.

---

## 5. Cómo ejecutarlo

1. Instalar las dependencias (sección 3).
2. Descargar el TSV de BindingDB y colocarlo en `data/`.
3. Abrir `CodigoTFM_VVJ.ipynb` y ejecutar las celdas en orden. El notebook está organizado en 16 secciones:

   1. Configuración y librerías (semilla `SEED = 42`, fingerprint 1024 bits, radio 2)
   2. Funciones reutilizables (Lipinski, fingerprints, etc.)
   3–5. Carga y limpieza de ChEMBL, Papyrus y BindingDB
   6. Fusión de las tres fuentes (deduplicación por SMILES canónico; pIC50 consolidado como mediana)
   7. Modelo 1 — baseline solo ChEMBL
   8. Modelo 2 — QSAR con dataset fusionado + comparación de modelos
   9. Candidatos y controles de validación (positivos y negativos)
   10. Dominio de aplicabilidad (Tanimoto)
   11. Predicción de pIC50 e integración con el docking (score integrado)
   12. Validación del modelo con los controles
   13. Ranking final
   14. Visualizaciones
   15. Exportación de resultados
   16. Figura comparativa de los perfiles PLIP

---

## 6. Salidas

- `ranking_final.csv` — ranking de los 7 candidatos (+ Ivacaftor de referencia) con docking, pIC50 predicho, Tanimoto, dominio de aplicabilidad y score final.
- `validacion_controles.csv` — predicciones para candidatos y controles positivos/negativos.
- Figuras: ranking coloreado por dominio de aplicabilidad, importancia de descriptores, espacio docking vs pIC50, y perfil comparativo de interacciones PLIP.

---

## 7. Reproducibilidad

- Semilla fija (`SEED = 42`) en NumPy, en el Random Forest y en las particiones.
- Random Forest: `n_estimators = 200`, `min_samples_leaf = 2`.
- Validación: k-fold (k = 5, *shuffle*) + *hold-out* independiente (20 %).
- Software exclusivamente de código abierto.

**Rendimiento del modelo QSAR ampliado:** R² = 0,649 ± 0,040 (CV k = 5); R² = 0,686 (*hold-out* 20 %); RMSE = 0,514.

---

## 8. Aviso

Todos los resultados son **computacionales** y constituyen hipótesis de priorización. Requieren validación experimental (electrofisiología / ensayos de función CFTR) antes de cualquier interpretación clínica.

---

## 9. Cómo citar

> Valero Just, V. (2026). *Reposicionamiento computacional de fármacos sobre CFTR: pipeline integrado de PPI, docking molecular y QSAR aplicado a la fibrosis quística.* Trabajo Final de Máster, Máster Universitario en Bioinformática y Bioestadística, Universitat Oberta de Catalunya (UOC) – Universitat de Barcelona (UB).
