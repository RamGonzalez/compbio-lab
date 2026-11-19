# Genómica de poblaciones: de la terminal al PCA y ADAMIXTURE

**Formato:** 7 sesiones de ~90 min, prácticas, cada estudiante en su propia laptop (Windows + WSL2/Ubuntu; una estudiante o estudiante con Mac usa la Terminal de macOS con los mismos comandos).
**Idioma del curso:** español.
**Público:** estudiantes de antropología física, biología y biotecnología. No se asume experiencia previa en programación ni en la línea de comandos.
**Presupuesto de hardware:** suponer **8 GB de RAM y 4 núcleos** por laptop. Todo comando del curso debe correr con holgura dentro de eso (usar `--threads 4`, `plink --memory 3000`, `adamixture -t 4`).
**Meta final:** cada estudiante puede correr PCA y ADAMIXTURE sobre un conjunto de datos pequeño y explicar qué dicen y qué no dicen los resultados.

## Principios de diseño

- **Primero la terminal, pero siempre al servicio de los datos.** Desde la sesión 2 se tocan datos genéticos reales, no ejercicios de juguete con `ls`.
- **Concepto, ejecución, interpretación.** Cada sesión de análisis tiene: intuición (sin matemáticas al inicio), una corrida y un ejercicio de interpretación.
- **Un solo conjunto de datos pequeño durante todo el curso.** Así, en la sesión 7 pueden predecir lo que mostrarán PCA y ADAMIXTURE y sorprenderse solo con lo interesante.
- **Todo es un script.** En la sesión 4 guardan sus comandos en un archivo `.sh`; en la sesión 7 pueden repetir todo el análisis con un solo comando.
- **El instructor prepara los datos.** Los estudiantes nunca descargan ni recortan datos de referencia.

## Software

Un entorno de conda en WSL (o macOS), creado por los estudiantes en la sesión 1 (instalado con Miniforge, entorno llamado `compbio-lab`):

- `plink` (1.9) y `plink2` (conda-forge/bioconda)
- `bcftools` (y `htslib` para bgzip/tabix) (bioconda)
- **ADAMIXTURE** (comando `adamixture`), instalado con **pip dentro del entorno** (está en PyPI, no en bioconda; requiere Python ≥ 3.10). Es una reimplementación rápida de ADMIXTURE que usa un optimizador tipo Adam dentro del algoritmo EM; ajusta el mismo modelo y escribe archivos `.P` y `.Q` compatibles con ADMIXTURE, con multihilo (`-t`), barrido de K integrado (`--min_k/--max_k`), validación cruzada integrada (`--cv`) y gráfica automática de Q. El modo CPU es el predeterminado y es el que usamos. Su licencia es gratuita solo para uso académico y de investigación.
- `r-base` con `ggplot2` y `tidyverse` (las gráficas se hacen en **R** en todo el curso)

**Estudiante con Mac:** usar el instalador de Miniforge para macOS (versión Apple Silicon o Intel). Antes de `pip install adamixture`, correr `brew install libomp` (se necesita OpenMP). Si algún paquete de bioconda no tiene versión nativa para Apple Silicon, crear el entorno con Rosetta (`CONDA_SUBDIR=osx-64`). Probar esto en una Mac antes de la sesión 1. Quedarse en CPU (no usar `--device mps`) para que los resultados coincidan con los del resto.

Mantener todo dentro de 8 GB / 4 núcleos: conjuntos pequeños de SNPs, pocas corridas en paralelo, ninguna herramienta que consuma mucha memoria. Fijar la versión de ADAMIXTURE en el archivo del entorno (el paquete evoluciona rápido y las opciones han cambiado entre versiones).

## Datos (se preparan más adelante)

Todavía no hacen falta. Plan para cuando llegue el momento: un subconjunto pequeño de 1000 Genomes, con **MXL** como grupo mezclado más anclas continentales (por ejemplo YRI, CEU, CHB), unos 40 individuos por población y solo unos pocos cromosomas, de modo que cada paso corra en segundos o minutos dentro del presupuesto de la laptop. Se entrega como un VCF, un conjunto PLINK ya filtrado y un `populations.tsv`. Opcionalmente una o dos variantes "trampa" para la sesión 7 (parientes que se dejan en la muestra, muestra desbalanceada).

---

## Sesión 1: Bienvenida, configuración y panorama del curso

**Objetivos de aprendizaje**
- Todos se conocen y saben qué espera aprender cada estudiante del curso.
- Todos pueden abrir WSL (o la Terminal en Mac) y tienen un entorno de conda funcionando con las herramientas del curso.
- Los estudiantes saben qué cubre el curso, qué podrán hacer al final y cómo se llevan las sesiones.

**Plan (90 min)**
- 0–20: Presentaciones (nombre, formación, qué saben de genética y cómputo, algo que quieran aprender). Presentación breve del curso: metas, las 7 sesiones, expectativas, cómo pedir ayuda.
- 20–35: Abrir WSL, primer vistazo al prompt, `pwd`, `ls`, `cd`, trabajar en el home de Linux (`~/`) y no en `/mnt/c`.
- 35–50: Instalar Miniforge/conda. Explicar qué es un gestor de paquetes y un entorno, y por qué los usamos.
- 50–80: Crear el entorno `compbio-lab` (plink, plink2, bcftools, R) y luego `pip install adamixture` dentro de él. Resolver y descargar el entorno puede tardar más de 10 minutos con conexiones lentas: empezarlo de inmediato y usar la espera para explicar lo que ocurre y mostrar el plan del curso.
- 80–90: Correr `check_setup.sh` (cada herramienta con `--version`/`--help`). Anotar quién sigue con problemas; resolverlos antes de la sesión 2 (horas de consulta o chat).

**Entregable/tarea:** una captura o la salida de `check_setup.sh` mostrando que todo está instalado. Dejar mucha holgura: los problemas de instalación siempre se alargan y las descargas con conexión lenta toman tiempo.

**Material:** [slides/sesion1_bienvenida_configuracion.pptx](slides/sesion1_bienvenida_configuracion.pptx)

## Sesión 2: La terminal y tu primer archivo genómico

**Objetivos de aprendizaje**
- Navegar el sistema de archivos y encontrar archivos de Windows desde Linux (`/mnt/c/...`).
- Usar `pwd ls cd mkdir cp mv rm less head tail wc cat`, autocompletado con tabulador e historial.
- Usar `grep`, `cut`, `sort`, `uniq`, pipes y redirección.
- Entender rutas, extensiones de archivo y archivos comprimidos.

**Plan**
- 0–10: Verificar rápidamente que el entorno se activa para todos.
- 10–40: Navegación guiada y comandos de archivos en una carpeta de práctica.
- 40–70: Pipes sobre datos reales: contar muestras por población en la tabla de metadatos con `cut | sort | uniq -c`; echar un vistazo a un VCF con `zcat | head`.
- 70–85: Primer script de shell que imprime el conteo de muestras por población.
- 85–90: Entrega de la guía rápida de comandos.

**Entregable/tarea:** `sesion2.sh` que produce la tabla de conteos.

## Sesión 3: Genética humana y variación genética, conceptos

**Objetivos de aprendizaje**
- Genoma, cromosoma, locus, alelo, genotipo, SNP; homocigoto vs heterocigoto.
- Frecuencia alélica; por qué la mayoría de las variantes son raras y unas pocas son comunes.
- Expectativa de Hardy-Weinberg (la idea y p², 2pq, q²).
- La variación humana es sobre todo *dentro* de las poblaciones; la estructura surge de la deriva, la migración y la historia demográfica (salida de África, efectos fundadores seriados, mezcla).
- La ancestría es una descripción estadística de similitud genética, no una categoría biológica de "raza". Dejar tiempo para la discusión; los estudiantes de antropología tendrán mucho que aportar.

**Plan**
- 0–30: Clase interactiva con dibujos (demostración de deriva en papel o con un simulador sencillo).
- 30–60: Práctica: extraer los genotipos de ~5 SNPs en unos pocos individuos con `bcftools query`; calcular a mano una frecuencia alélica y compararla con `plink --freq`.
- 60–80: Frecuencias por población de SNPs muy diferenciados (por ejemplo *SLC24A5*, *LCT*, Duffy) y por qué difieren.
- 80–90: Los estudiantes escriben predicciones: cuántos grupos habrá y dónde caerá MXL (se retoman en la sesión 7).

**Entregable/tarea:** hoja de una página con sus predicciones.

## Sesión 4: Datos genotípicos en la práctica (VCF, PLINK, control de calidad, poda por LD)

**Objetivos de aprendizaje**
- Leer un VCF: encabezado, columnas fijas, campos de genotipo (`0|1`, `./.`).
- Conocer los archivos PLINK `.bed/.bim/.fam` y convertir VCF → PLINK.
- Filtrar: SNPs bialélicos, datos faltantes, MAF, y por qué importa cada uno.
- Entender el desequilibrio de ligamiento (LD) y por qué PCA/ADAMIXTURE necesitan SNPs podados.
- Parentesco a nivel conceptual (y `plink2 --king-cutoff`).

**Plan**
- 0–20: Anatomía de un VCF (impreso anotado).
- 20–50: `plink --vcf ... --make-bed`; `--missing`, `--freq`, `--hardy`; revisar las salidas; un histograma en R.
- 50–70: Concepto de LD, `--indep-pairwise 50 5 0.2`, comparar el número de SNPs antes y después.
- 70–90: Reunir los comandos de control de calidad en `01_qc.sh`, que produce el conjunto limpio y podado que se usará de aquí en adelante.

**Entregable/tarea:** `01_qc.sh` corre desde cero; tabla con el número de SNPs tras cada paso.

## Sesión 5: PCA, conceptos y ejecución

**Objetivos de aprendizaje**
- Intuición: el PCA rota los datos hacia los ejes de mayor varianza; cada individuo recibe coordenadas (PCs).
- Valores propios / % de varianza explicada, cargas (brevemente).
- Correr `plink2 --pca` y graficar PC1 vs PC2 coloreado por población en R.
- Interpretar grupos, gradientes (individuos mezclados entre grupos) y valores atípicos.

**Plan**
- 0–25: Intuición sin álgebra matricial, con una analogía natural para biólogos; luego la matriz de genotipos (individuos × SNPs) → PCs.
- 25–45: Correr el PCA, revisar `.eigenvec` / `.eigenval`.
- 45–75: Gráficas en R: PC1–PC2, PC1–PC3, gráfica de sedimentación (scree), coloreadas por población.
- 75–90: Interpretación en parejas.

**Preguntas de interpretación**
- ¿Qué separa PC1? ¿Y PC2? ¿Qué poblaciones quedan entre otras?
- ¿Dónde caen los individuos MXL y por qué la dispersión?
- ¿Qué pasa si quitas una población y vuelves a correr? (Los PCs cambian de significado.)

**Entregable/tarea:** `02_pca.sh` + `plot_pca.R`; un párrafo interpretando la gráfica.

## Sesión 6: ADAMIXTURE, conceptos y ejecución

**Objetivos de aprendizaje**
- Intuición del modelo: cada genoma es una mezcla de K poblaciones ancestrales; las salidas son la matriz Q (proporciones) y la matriz P (frecuencias alélicas).
- K es una decisión de modelado; las "poblaciones ancestrales" son construcciones estadísticas, no grupos históricos reales.
- Correr `adamixture --cv` para K = 2…5; elegir K con el error de validación cruzada *y* con criterio biológico.
- ADAMIXTURE ajusta el mismo modelo que el ADMIXTURE clásico, solo que más rápido (optimizador Adam dentro del EM), así que lo aprendido se transfiere a los resultados de ADMIXTURE que aparecen en artículos.
- Gráficas de barras apiladas; las etiquetas y colores de los grupos son arbitrarios entre corridas.

**Plan**
- 0–25: Concepto (modelo de mezcla, K; relación con el PCA: la misma estructura, otra representación).
- 25–55: Correr K=2..5, con `--min_k 2 --max_k 5` o con un `for` en bash, usando `-t 4` como máximo; los datos son lo bastante pequeños para que cada corrida tarde de segundos a pocos minutos. Mirar primero el PNG automático de Q y luego reconstruirlo en R.
- 55–80: Extraer el error de validación cruzada por K del log (`tee`/`grep`), graficar CV vs K y hacer la gráfica de barras de los archivos `.Q` ordenada por población en R. (Revisar el formato exacto de la salida de CV en la versión fijada antes de la clase.)
- 80–90: Discusión: ¿qué muestra K=3 para MXL? ¿Qué aportan K=2 y K=4?

**Entregable/tarea:** `03_adamixture.sh` + `plot_adamixture.R`; gráficas de barras para K=2–5.

## Sesión 7: Integración, trampas comunes y mini proyecto

**Objetivos de aprendizaje**
- Combinar PCA y ADAMIXTURE: coincidencias, diferencias y límites.
- Reconocer trampas: muestreo desbalanceado, parientes, muy pocos SNPs, no podar por LD, sobreinterpretar K, leer la mezcla como historia o fechas.
- Reproducir todo el flujo con un solo script.

**Plan**
- 0–10: Revisar las predicciones de la sesión 3.
- 10–50: Mini proyecto en parejas con una variante "trampa" de los datos: volver a correr PCA y ADAMIXTURE y documentar qué cambió y por qué.
- 50–80: Presentaciones de 5 minutos por pareja: una figura, tres frases de interpretación y una salvedad.
- 80–90: Cierre: qué se puede y qué no se puede concluir; siguientes pasos (ancestría local, estadísticos f, paneles de referencia, ética de la inferencia de ancestría).

**Entregable:** informe final de 1–2 páginas con figuras de PCA y ADAMIXTURE y una interpretación.

---

## Notas transversales

- **Riesgos de ritmo:** las instalaciones de la sesión 1 y el tiempo de cómputo de la sesión 6 son lo que más se alarga. Tener una carpeta con resultados precalculados para que quien se atore pueda continuar.
- **Consejos de WSL:** trabajar en el home de Linux (`~/`) y no en `/mnt/c`, por velocidad; copiar las figuras a Windows para verlas (`explorer.exe .` abre la carpeta actual).
- **Evaluación (ligera):** un pequeño entregable por sesión + informe final; participación en las discusiones de interpretación.
- **Hilo de ética (2–3 momentos breves, no una sesión aparte):** ancestría vs raza, muestreo y consentimiento en conjuntos de datos genómicos, límites de los paneles de referencia.

## Materiales por construir (orden sugerido)

1. Diapositivas de la sesión 1 (listas) y `setup/` (pasos de instalación y `check_setup.sh`).
2. `data/`: conjunto de datos pequeño y README con su procedencia (más adelante).
3. `sessions/01…07/`: notas, ejercicios, scripts con solución, guías rápidas.
4. `results_reference/`: resultados precalculados de PCA/ADAMIXTURE del instructor.
