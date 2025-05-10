# 🧪 Essential Oils GC-MS Kovats Analysis

Este repositório contém um pipeline de análise de dados GC-MS de óleos essenciais, incluindo:

- Cálculo de Índice de Kovats
- Similaridade espectral por cosseno
- Análise de cromatogramas
- Merge de resultados espectrais e Kovats
- Visualizações interativas (.html)

## 🚩 Workflow configurado no MZmine
1. Converter os dados brutos à .CDF para importação no MZMine
2. Plotar o Chromatogram (como "Base peak intensity") e exportar em formato .XLSX
3. Iniciar o processamento seguindo as etapas abaixo.

Abaixo os parâmetros utilizados no arquivo `Essential_Oils_Batch.xml`: ESTES PARÂMETROS NÃO DEVEM SER SEGUIDOS, MAS SIM OTIMIZADOS DE ACORDO COM OS DADOS BRUTOS.

| **Etapa**                   | **Parâmetro**                           | **Valor**                      | **Descrição**                                                   |
|-----------------------------|----------------------------------------|--------------------------------|-----------------------------------------------------------------|
| **Mass Detection**           | Detector                                | Centroid                       | Detector de picos centroidizados                                |
|                             | Noise level                             | 10000.0                        | Limite mínimo de intensidade para detectar um pico              |
|                             | Isótopos abaixo do noise level          | Ativado                        | Permite isótopos mesmo abaixo do limite de ruído               |
|                             | Elementos para isótopos                 | H, C, N, O, S                  | Elementos considerados para detecção isotópica                  |
|                             | Tolerância m/z isótopos                 | 0.5 Da ou 0.0 ppm            | Tolerância para agrupamento de isótopos                         |
| **Chromatogram Builder**     | Mínimo de scans consecutivos            | 12                             | Número mínimo de pontos contínuos no cromatograma               |
|                             | Intensidade mínima por scan              | 1000                           | Intensidade mínima para contar um scan válido                   |
|                             | Altura mínima absoluta                  | 10000                          | Altura mínima do pico cromatográfico                            |
|                             | Tolerância m/z entre scans              | 0.5 Da                         | Tolerância de variação de m/z entre pontos                      |
| **ADAP Peak Picking**:      | S/N threshold                           | 5.0                            | Relação sinal/ruído mínima para detecção                        |
|                             | Mínima altura do pico                   | 1000                           | Altura mínima de um pico para ser considerado                   |
|                             | Coeficiente/área mínimo                 | 50                             | Mínimo coeficiente/área                                        |
|                             | Intervalo de tempo do pico              | 0 - 10 min                     | Intervalo aceito para duração do pico                           |
|                             | Pareamento com espectros MS/MS          | Ativado                        | Permite pareamento com espectros MS/MS                          |
|                             | Tolerância precursor m/z                | 0.01 Da / 10 ppm               | Tolerância para pareamento de precursor                         |
|                             | Filtro por tempo de retenção            | 0.2 min                        | Filtro de tempo de retenção                                     |
| **ADAP Hierarchical Clustering**     | Distância mínima de cluster             | 0.01                           | Distância mínima para formar um cluster                         |
|                             | Tamanho mínimo de cluster               | 2                              | Número mínimo de picos para formar um cluster                   |
|                             | Intensidade mínima do cluster           | 500                            | Intensidade mínima de um cluster                                |
|                             | Tolerância de similaridade de forma     | 60.0                           | Similaridade mínima para agrupar picos pela forma               |
## 📂 Estrutura do repositório

...
 

Sugestão de estrutura

/
├── essential_oils_kovats_processing.py   # Módulo Python com as funções auxiliares

├── README.md                             # Arquivo de descrição do projeto

└── Resources/

    ├── Hidrocarbon_RTs.csv/              # Dados obtidos da análise experimental da amostra padrão de hidrocarbonetos: "Name of Hydrocarbon"	"Number of Carbon Atoms"	"Seconds"	"Current RT in minutes"	"Typical RT in sec."	"Typical RT in minutes"

    ├── RI_Library.csv/                   # Banco de dados de Índice de Retenção para KOVATS
    ├── Essential_Oils_Filtered.msp/      # Banco de dados de Espectros de MS (filtrado)

└── Projeto_X/

    ├── dados/                            # Pasta com arquivos de dados brutos
    
    ├── GCMS_essentialOils_Kovats-Batch.ipynb  # Notebook principal de análise
    
    ├── results/                          # Resultados exportados (.xlsx, .html, etc.)
    
    └── images/                           # Gráficos gerados e salvos


