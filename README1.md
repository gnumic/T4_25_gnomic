#
## ПЕРЕД ЗАПУСКОМ СКРИПТІВ ЗМІНІТЬ ШЛЯХИ В ФАЙЛІ PATHs.txt!!!!
#

## Огляд проекту

Цей проект відтворює аналіз даних РНК-секвенування мікроглії з зразків мозку пацієнтів із хворобою Альцгеймера (AD) та контрольних донорів з бета-амілоїдними бляшками (CTR+), опублікований у статті "Profiling Microglia From Alzheimer's Disease Donors and Non-demented Elderly in Acute Human Postmortem Cortical Tissue" (https://doi.org/10.3389/fnmol.2020.00134).

Проект використовує набір даних "Single-cell RNA Sequencing of human microglia from post mortem Alzheimers Disease CNS tissue" і фокусується на кластеризації клітин за експресією генів за допомогою Seurat та візуалізації результатів через UMAP.

## Структура проекту

```
project/
├── data/                          # Вихідні дані
│   ├── SingleCellsAlzheimers_*.tsv   # Метадані проекту
│   └── GSE146639_RAW/             # Вихідні файли експресії
├── export/                        # Результати аналізу
│   ├── analyze_*                  # Результати аналітичних скриптів
│   ├── build_*_objects_*          # Побудовані об'єкти Seurat
│   └── seurat_objects/            # Збережені об'єкти Seurat
├── modules/                       # Модулі обробки
│   ├── parser/                    # Скрипти парсингу вихідних даних
│   └── seurat_object/             # Скрипти створення та аналізу об'єктів Seurat
└── Supplementary Material/        # Додаткові матеріали
```

## Ключові результати

Проект виконує:
1. Обробку та нормалізацію даних scRNA-seq з використанням SCTransform
2. Кластеризацію клітин на основі експресії генів
3. Аналіз диференційно експресованих генів між групами AD та CTR+
4. UMAP візуалізацію кластерів клітин

# Інструкція для користувача

## Вимоги

- R версії 4.0 або вище
- Пакети:
  - Seurat (для аналізу одноклітинних даних)
  - glmGamPoi (для прискорення SCTransform)
  - presto (для швидких операцій з матрицями)
  - clustree (для визначення оптимальних параметрів кластеризації)
  - dplyr, ggplot2 (для маніпуляцій з даними та візуалізації)

## Встановлення залежностей

```r
# Встановлення пакетів CRAN
if (!require("Seurat")) install.packages("Seurat")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("clustree")) install.packages("clustree")

# Встановлення BioConductor пакетів
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("glmGamPoi")) BiocManager::install("glmGamPoi")

# Встановлення пакетів з GitHub
if (!require("devtools")) install.packages("devtools")
if (!require("presto")) devtools::install_github("immunogenomics/presto")
```

## Порядок використання

1. **Створення об'єктів Seurat**
   - Запустіть скрипт build_CTRplus_AD_objects.R
   - Скрипт обробить файли .txt з директорії даних і створить об'єкти Seurat
   - Результати збережуться в `export/build_CTRplus_AD_objects_[дата_час]/seurat_objects/`

2. **Аналіз об'єктів Seurat**
   - Запустіть скрипт analyze_CTRplus_AD.R
   - Скрипт виконає:
     * Завантаження створених об'єктів
     * Нормалізацію за допомогою SCTransform
     * Визначення варіабельних ознак
     * PCA та кластеризацію
     * UMAP візуалізацію
     * Пошук маркерів кластерів
     * Пошук диференційно експресованих генів між AD та CTR+
   - Результати аналізу збережуться в `export/analyze_CTRplus_AD_[дата_час]/`

3. **Перегляд результатів**
   - Візуалізації зберігаються в директорії з результатами аналізу
   - Об'єкти Seurat можна завантажити для подальшого аналізу:
   ```r
   library(Seurat)
   combined_seurat <- readRDS("path/to/combined_seurat.rds")
   ```

## Примітки до використання

1. Всі скрипти створюють окремі директорії для виводу з датою та часом, щоб уникнути перезапису попередніх результатів.

2. Логи виконання зберігаються у файли з назвою скрипта та датою/часом.
