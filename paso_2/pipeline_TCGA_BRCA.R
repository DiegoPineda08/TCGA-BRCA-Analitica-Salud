# ============================================================
# PIPELINE DE CIENCIA DE DATOS: TCGA-BRCA
# Oncología Traslacional — Universidad Externado de Colombia
# ============================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)

# 1. Consulta de 20 archivos para un análisis reproducible y manejable
resultados <- getResults(
  GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
)

casos <- resultados$cases[1:20]

query_small <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = casos
)

resultados_small <- getResults(query_small)
print(nrow(resultados_small))

# 2. Descarga
GDCdownload(query_small, method = "api", files.per.chunk = 5)

# 3. Preparación de la matriz
datos <- GDCprepare(query_small)
matriz_genes <- assay(datos)

# Evidencia de dimensiones
print(dim(matriz_genes))

# 4. Reducción por varianza
matriz_log <- log2(matriz_genes + 1)

varianzas_genes <- apply(
  matriz_log,
  1,
  var,
  na.rm = TRUE
)

top_10_indices <- order(
  varianzas_genes,
  decreasing = TRUE
)[1:10]

matriz_reducida <- matriz_genes[top_10_indices, ]

# 5. Identificación de genes
genes_ensembl <- sub("\\..*$", "", rownames(matriz_reducida))

anotacion <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes_ensembl,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)

print(anotacion)

# 6. Tabla final
tabla_top10 <- data.frame(
  ENSEMBL = genes_ensembl,
  SYMBOL = anotacion$SYMBOL,
  GENENAME = anotacion$GENENAME,
  VARIANZA = as.numeric(varianzas_genes[top_10_indices])
)

print(tabla_top10)

# Dimensiones finales
print(dim(matriz_reducida))
