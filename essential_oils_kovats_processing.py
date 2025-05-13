"""
essential_oils_kovats_processing.py

Módulo com funções auxiliares para análise de óleos essenciais por GC-MS:
- Cálculo do Índice de Kovats
- Ajuste de curva de retenção
- Processamento de cromatogramas
- Parsing de arquivos MGF e MSP
- Cálculo de similaridade espectral
- Geração e salvamento de gráficos

Autor: [Seu nome ou laboratório]
Data: [opcional]
"""

# =========================
# 📐 FITTING / CURVAS
# =========================

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import os
import plotly.express as px

def fit_polynomial_trendline(df, x_col, y_col, degree=4, num_points=100):
    """
    Ajusta uma curva polinomial de grau especificado aos dados e calcula R².
    """
    z = np.polyfit(df[x_col], df[y_col], deg=degree)
    p = np.poly1d(z)
    x_fit = np.linspace(df[x_col].min(), df[x_col].max(), num_points)
    y_fit = p(x_fit)

    y_true = df[y_col]
    y_pred = p(df[x_col])
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    return x_fit, y_fit, p, r_squared

def save_plotly_figure(fig, filename="curva_retencao.html", output_dir="images"):
    """
    Salva uma figura Plotly como arquivo HTML interativo.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, filename)
    fig.write_html(output_path)
    print(f"✅ Figura salva em: {output_path}")

import os
import pandas as pd

def carregar_dados(input_dir, file_hidrocarbon_rts, file_ri_library):
    """
    Carrega os dados de RTs de alcanos, referência RI e os arquivos de cada amostra no diretório.

    Parâmetros:
    - input_dir: caminho da pasta com os arquivos de amostras
    - file_hidrocarbon_rts: caminho do arquivo Hidrocarbon_RTs.csv
    - file_ri_library: caminho do arquivo RI_Library.csv

    Retorna:
    - df_hidrocarbon: DataFrame da tabela de RTs dos alcanos
    - df_ri: DataFrame da tabela de referência RI
    - sample_prefixes: set com os prefixos das amostras
    - dfs_quant: dicionário prefix → DataFrame de quantificação
    - dfs_chrom: dicionário prefix → DataFrame do cromatograma
    """
    # ✅ Leitura da tabela de RTs dos alcanos
    df_hidrocarbon = pd.read_csv(file_hidrocarbon_rts, sep=',|\t', engine='python')
    print("✅ Tabela de RTs dos alcanos importada.")

    # ✅ Leitura da tabela de referência RI
    df_ri = pd.read_csv(file_ri_library, sep='\t')
    print("✅ Tabela de referência RI importada.")

    # 🔍 Procurar todos os arquivos na pasta
    all_files = os.listdir(input_dir)

    # 🎯 Identificar os prefixos das amostras
    sample_prefixes = set()
    for f in all_files:
        if f.endswith("_quant.csv"):
            prefix = f.replace("_quant.csv", "")
            sample_prefixes.add(prefix)

    # 📦 Criar dicionários para armazenar os DataFrames de cada amostra
    dfs_quant = {}
    dfs_chrom = {}

    for prefix in sample_prefixes:
        file_quant = os.path.join(input_dir, f"{prefix}_quant.csv")
        file_chrom = os.path.join(input_dir, f"{prefix}_Chromatogram.xlsx")
        
        if os.path.exists(file_quant) and os.path.exists(file_chrom):
            df_quant = pd.read_csv(file_quant)
            df_chrom = pd.read_excel(file_chrom)
            dfs_quant[prefix] = df_quant
            dfs_chrom[prefix] = df_chrom
            print(f"✅ Amostra {prefix} importada.")
        else:
            print(f"⚠️ Arquivos faltando para amostra {prefix}. Verifique!")

    return df_hidrocarbon, df_ri, sample_prefixes, dfs_quant, dfs_chrom



# =========================
# 🧮 ÍNDICE DE KOVATS
# =========================

def calculate_kovats(rt, alkane_table, rt_column):
    """
    Calcula o Índice de Kovats para um tempo de retenção dado.
    """
    for i in range(len(alkane_table) - 1):
        n = alkane_table.iloc[i]['Number of Carbon Atoms']
        rt_n = alkane_table.iloc[i][rt_column]
        rt_n1 = alkane_table.iloc[i + 1][rt_column]
        if rt_n <= rt <= rt_n1:
            return 100 * n + 100 * (rt - rt_n) / (rt_n1 - rt_n)
    return None

def identify_compound_multiple(kovats, ref_table, tolerance=10):
    """
    Identifica compostos por Índice de Kovats com tolerância.
    """
    match = ref_table.loc[np.abs(ref_table['RI'] - kovats) <= tolerance]
    if not match.empty:
        return "; ".join(match['Common Name'].unique())
    else:
        return 'N/A'

def process_sample_kovats(df_quant, alkane_table, rt_column, ref_table=None, tolerance=10):
    """
    Processa uma tabela de picos adicionando o Índice de Kovats e identificação.
    """
    kovats_list = [calculate_kovats(rt, alkane_table, rt_column) for rt in df_quant['row retention time']]
    df_quant = df_quant.copy()
    df_quant['Kovats Index'] = kovats_list

    if ref_table is not None:
        compounds = [identify_compound_multiple(k, ref_table, tolerance) for k in df_quant['Kovats Index']]
        df_quant['Kovats Identified Compound'] = compounds
        if "Unnamed: 4" in df_quant.columns:
            df_quant = df_quant.drop(columns="Unnamed: 4")
    
    return df_quant

# =========================
# 📊 CROMATOGRAMAS
# =========================

def plot_chromatogram_auto_multitrace(df_chrom, sample_prefix, output_dir="images"):
    """
    Plota e salva o cromatograma com múltiplos traços (um por par de colunas),
    salvando em HTML interativo e PNG de alta resolução.
    """
    if isinstance(df_chrom.columns, pd.MultiIndex):
        df_chrom.columns = [' | '.join([str(level).strip() for level in col if str(level) != 'nan']) for col in df_chrom.columns]
    else:
        df_chrom.columns = df_chrom.columns.str.strip()

    print(f"\n🔍 Colunas disponíveis para {sample_prefix}:")
    for idx, col in enumerate(df_chrom.columns):
        print(f"{idx}: {col}")

    fig = go.Figure()

    for i in range(0, len(df_chrom.columns)-1, 2):
        x_col = df_chrom.columns[i]
        y_col = df_chrom.columns[i+1]
        if df_chrom[x_col].isnull().all() or df_chrom[y_col].isnull().all():
            continue
        fig.add_trace(go.Scatter(
            x=df_chrom[x_col],
            y=df_chrom[y_col],
            mode='lines',
            name=x_col.split(':')[0]
        ))

    fig.update_layout(
        title=f"Cromatograma: {sample_prefix}",
        xaxis_title="Retention Time (min)",
        yaxis_title="Base Peak Intensity"
    )

    os.makedirs(output_dir, exist_ok=True)

    # ✅ Salva HTML
    output_html = os.path.join(output_dir, f"{sample_prefix}_chromatogram.html")
    fig.write_html(output_html)
    print(f"✅ Gráfico salvo como HTML em: {output_html}")

    # ✅ Salva PNG de alta resolução
    #output_png = os.path.join(output_dir, f"{sample_prefix}_chromatogram.png")
    #fig.write_image(output_png.replace(".png", ".svg"))#, scale=1.5, width=600, height=400)  # ajuste a resolução aqui
    #print(f"✅ Gráfico salvo como PNG em alta resolução em: {output_png}")

    return fig

# =========================
# 🗂️ PARSING DE ARQUIVOS
# =========================

def parse_mgf(file_path):
    """
    Parseia um arquivo .mgf e retorna um dicionário {FEATURE_ID: [(mz, intensidade), ...]}.
    """
    spectra = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        current_id = None
        peaks = []
        for line in f:
            line = line.strip()
            if line.startswith('FEATURE_ID='):
                current_id = line.split('=')[1]
            elif line.startswith('END IONS'):
                if current_id and peaks:
                    spectra[current_id] = peaks
                peaks, current_id = [], None
            elif line and not line.startswith(('BEGIN', 'PEPMASS', 'RTINSECONDS', 'SCANS', 'MSLEVEL', 'CHARGE')):
                parts = line.split()
                if len(parts) == 2:
                    mz, intensity = map(float, parts)
                    peaks.append((mz, intensity))
    return spectra

def parse_msp(file_path):
    """
    Parser tolerante para arquivos .msp no padrão NIST.
    Suporta blocos com 'Name:', 'Num Peaks:', e múltiplos formatos de
espaçamento.
    """
    spectra = {}
    current_name = None
    peaks = []
    in_peak_section = False
    total_lines = 0
    total_blocks = 0

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            total_lines += 1
            line = line.strip()
            if not line:
                continue

            if line.lower().startswith("name:"):
                if current_name and peaks:
                    spectra[current_name] = peaks
                    total_blocks += 1
                current_name = line.split(":", 1)[1].strip()
                peaks = []
                in_peak_section = False

            elif line.lower().startswith("num peaks"):
                in_peak_section = True

            elif in_peak_section:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                    except ValueError:
                        continue

        if current_name and peaks:
            spectra[current_name] = peaks
            total_blocks += 1

    print(f"📄 Linhas processadas: {total_lines}")
    print(f"✅ Blocos espectrais válidos: {total_blocks}")
    return spectra

import os
import pandas as pd

def processar_mgf_amostras(sample_prefixes, input_dir, spectra_lib, results_by_sample, tolerance=0.1, sim_threshold=0.7, output_dir="results"):
    """
    Processa os arquivos .mgf de cada amostra, busca matches na biblioteca e salva resultados.
    
    Parâmetros:
    - sample_prefixes: lista ou set de prefixos de amostra
    - input_dir: caminho da pasta com os arquivos .mgf
    - spectra_lib: biblioteca de espectros (já parseada)
    - results_by_sample: dicionário onde armazenar DataFrames de resultados por amostra
    - tolerance: tolerância de m/z para matching (default=0.1)
    - sim_threshold: similaridade mínima para aceitar match (default=0.7)
    - output_dir: diretório para salvar resultados .csv
    """
    os.makedirs(output_dir, exist_ok=True)

    for prefix in sample_prefixes:
        file_mgf = os.path.join(input_dir, f"{prefix}.mgf")
        
        if os.path.exists(file_mgf):
            print(f"\n🔍 Processando {file_mgf}...")

            # Parseia o arquivo MGF da amostra
            spectra_exp = parse_mgf(file_mgf)
            print(f"✅ {len(spectra_exp)} espectros lidos de {file_mgf}")

            # Aplica a busca de melhor match
            results = find_all_matches(spectra_exp, spectra_lib, tolerance=tolerance, sim_threshold=sim_threshold)

            # Converte para DataFrame
            df_results = pd.DataFrame(results)
            results_by_sample[prefix] = df_results

            # Exporta para CSV
            output_csv = os.path.join(output_dir, f"{prefix}_matches.csv")
            df_results.to_csv(output_csv, index=False)

            print(f"✅ Resultados salvos em {output_csv}")
        
        else:
            print(f"⚠️ Arquivo .mgf não encontrado para {prefix}")


# =========================
# 🔬 SIMILARIDADE ESPECTRAL
# =========================

def cosine_similarity(spec1, spec2, tolerance=0.1):
    """
    Calcula similaridade de cosseno entre dois espectros considerando tolerância de m/z.
    """
    matched_ints1, matched_ints2 = [], []
    mzs1, ints1 = zip(*spec1) if spec1 else ([], [])
    mzs2, ints2 = zip(*spec2) if spec2 else ([], [])
    
    for m1, i1 in zip(mzs1, ints1):
        best_match, best_delta = None, tolerance
        for m2, i2 in zip(mzs2, ints2):
            delta = abs(m1 - m2)
            if delta <= best_delta:
                best_delta = delta
                best_match = i2
        if best_match is not None:
            matched_ints1.append(i1)
            matched_ints2.append(best_match)

    if not matched_ints1 or not matched_ints2:
        return 0.0
    numerator = sum(i*j for i,j in zip(matched_ints1, matched_ints2))
    denom = np.sqrt(sum(i**2 for i in matched_ints1)) * np.sqrt(sum(j**2 for j in matched_ints2))
    return numerator / denom if denom else 0.0

def find_best_matches(spectra_exp, spectra_lib, tolerance=0.1):
    """
    Encontra o melhor match espectral para cada FEATURE_ID no espectro experimental.
    """
    results = []
    for feat_id, exp_spec in spectra_exp.items():
        best_score, best_name = 0, 'No match'
        for lib_name, lib_spec in spectra_lib.items():
            score = cosine_similarity(exp_spec, lib_spec, tolerance)
            if score > best_score:
                best_score, best_name = score, lib_name
        results.append({'FEATURE_ID': feat_id, 'Best Match Name': best_name, 'Similarity': best_score})
    return results

def find_all_matches(spectra_exp, spectra_lib, tolerance=0.1, sim_threshold=0.7):
    """
    Retorna todas as correspondências com similaridade >= limiar.
    """
    results = []
    for feat_id, exp_spec in spectra_exp.items():
        for lib_name, lib_spec in spectra_lib.items():
            score = cosine_similarity(exp_spec, lib_spec, tolerance)
            if score >= sim_threshold:
                results.append({'FEATURE_ID': feat_id, 'Best Match Name': lib_name, 'Similarity': round(score, 2)})
    return results

# =========================
# UTILITÁRIOS
# =========================

def read_first_lines(file_path, num_lines=20):
    """
    Lê as primeiras linhas de um arquivo texto.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        return [next(f) for _ in range(num_lines)]


from docxtpl import DocxTemplate
from docxtpl import InlineImage
from docx.shared import Cm
import pandas as pd
import os

import os
import pandas as pd

def processar_fusao_amostras(sample_prefixes, dfs_quant_processed, results_by_sample, prefix_ignorado="GT_PadraoHidroc_13-06-2024", output_dir="results"):
    """
    Processa a fusão dos DataFrames de quantificação e identificação para cada amostra,
    ignorando um prefixo especificado, e exporta o resultado em Excel.
    
    Parâmetros:
    - sample_prefixes: lista de nomes das amostras (prefixos)
    - dfs_quant_processed: dicionário de DataFrames quantificados por amostra
    - results_by_sample: dicionário de DataFrames de identificações por amostra
    - prefix_ignorado: nome da amostra padrão a ser ignorada
    - output_dir: diretório onde salvar os arquivos .xlsx
    """
    os.makedirs(output_dir, exist_ok=True)

    for prefix in sample_prefixes:
        if prefix == prefix_ignorado:
            print(f"⚠️ Ignorando amostra padrão {prefix}")
            continue

        print(f"\n🔍 Processando fusão para {prefix}...")

        df_quant = dfs_quant_processed[prefix].copy()
        df_matches = results_by_sample[prefix].copy()

        if 'FEATURE_ID' not in df_quant.columns:
            if 'row ID' in df_quant.columns:
                df_quant['FEATURE_ID'] = df_quant['row ID'].astype(str)
            else:
                raise ValueError(f"❌ Nenhuma coluna 'row ID' ou 'FEATURE_ID' em df_quant para {prefix}")
        else:
            df_quant['FEATURE_ID'] = df_quant['FEATURE_ID'].astype(str).str.replace('FEATURE_ID=', '', regex=False)

        df_quant['row ID'] = df_quant['FEATURE_ID'].astype(int)
        df_matches['row ID'] = df_matches['FEATURE_ID'].astype(int)

        for col in ['Kovats Index', 'Kovats Identified Compound']:
            if col not in df_quant.columns:
                df_quant[col] = None

        df_matches_renamed = df_matches.rename(columns={
            'Best Match Name': 'Spectral Match Name',
            'Similarity': 'Spectral Similarity'
        })

        df_combined = pd.merge(
            df_quant[['row ID', 'row m/z', 'row retention time', 'Kovats Index', 'Kovats Identified Compound']],
            df_matches_renamed[['row ID', 'Spectral Match Name', 'Spectral Similarity']],
            on='row ID',
            how='outer'
        )

        print(f"✅ Fusão completa para {prefix}")

        output_file = os.path.join(output_dir, f"{prefix}_Identificacao_Combinada.xlsx")
        df_combined.to_excel(output_file, index=False)
        print(f"✅ Arquivo salvo: {output_file}")

    print(f"\n🎉 Processamento de todas as amostras concluído.")



def gerar_relatorio_docx(tabela_resultados, amostra_nome, output_file, template_path="Relatorio_Analitico_Template.docx"):
    """
    Gera um relatório .docx baseado em um template Word com os dados da amostra.
    
    Parâmetros:
    - tabela_resultados: DataFrame com os resultados da amostra
    - amostra_nome: string com o nome da amostra
    - output_file: caminho completo do arquivo .docx de saída
    - template_path: caminho do arquivo template .docx (padrão: 'Relatorio_Analitico_Template.docx' na mesma pasta)
    """
    # Verifica se o template existe
    if not os.path.exists(template_path):
        raise FileNotFoundError(f"❌ Template não encontrado em: {template_path}")

    # Carrega o template
    doc = DocxTemplate(template_path)

    # Converte DataFrame para lista de listas (para usar no template)
    tabela_lista = [tabela_resultados.columns.tolist()] + tabela_resultados.values.tolist()

    # Define contexto de substituição
    context = {
        'amostra': amostra_nome,
        'tabela': tabela_lista,
        'data': pd.Timestamp.today().strftime('%d/%m/%Y')}

    # Renderiza o template com o contexto
    doc.render(context)

    # Salva o arquivo de saída
    doc.save(output_file)
    print(f"✅ Relatório gerado e salvo: {output_file}")
    
from docxtpl import DocxTemplate, InlineImage
from docx.shared import Cm
import os
import pandas as pd

def processar_amostra(prefix, dfs_quant_processed, results_by_sample, template_path="../Relatorio_Analitico_Template.docx"):
    if prefix == "GT_PadraoHidroc_13-06-2024":
        print(f"⚠️ Ignorando amostra padrão {prefix}")
        return

    print(f"\n🔍 Processando fusão e relatório para {prefix}...")
    
    df_quant = dfs_quant_processed[prefix].copy()
    df_matches = results_by_sample[prefix].copy()

    if 'FEATURE_ID' not in df_quant.columns:
        if 'row ID' in df_quant.columns:
            df_quant['FEATURE_ID'] = df_quant['row ID'].astype(str)
        else:
            raise ValueError(f"❌ Nenhuma coluna 'row ID' ou 'FEATURE_ID' em df_quant para {prefix}")
    else:
        df_quant['FEATURE_ID'] = df_quant['FEATURE_ID'].astype(str).str.replace('FEATURE_ID=', '', regex=False)

    df_quant['row ID'] = df_quant['FEATURE_ID'].astype(int)
    df_matches['row ID'] = df_matches['FEATURE_ID'].astype(int)

    # ✅ Garante que as colunas existam
    for col in ['Kovats Index', 'Kovats Identified Compound']:
        if col not in df_quant.columns:
            df_quant[col] = None

    df_matches_renamed = df_matches.rename(columns={
        'Best Match Name': 'Spectral Match Name',
        'Similarity': 'Spectral Similarity'
    })

    df_combined = pd.merge(
        df_quant[['row ID', 'row m/z', 'row retention time', 'Kovats Index', 'Kovats Identified Compound']],
        df_matches_renamed[['row ID', 'Spectral Match Name', 'Spectral Similarity']],
        on='row ID',
        how='outer'
    )

    print(f"✅ Fusão completa para {prefix}")
    #display(df_combined.head(10))

    # ✅ Salvar Excel
    os.makedirs("results", exist_ok=True)
    output_file = os.path.join("results", f"{prefix}_Identificacao_Combinada.xlsx")
    df_combined.to_excel(output_file, index=False)
    print(f"✅ Arquivo salvo: {output_file}")

    # ✅ Criar lista para relatório
    tabela_lista = []
    for _, row in df_combined.iterrows():
        tabela_lista.append({
            'row_mz': f"{row['row m/z']:.2f}" if pd.notnull(row['row m/z']) else "",
            'row_retention_time': f"{row['row retention time']:.2f}" if pd.notnull(row['row retention time']) else "",
            'kovats_index': f"{row['Kovats Index']:.2f}" if pd.notnull(row['Kovats Index']) else "",
            'kovats_identified_compound': row['Kovats Identified Compound'] if pd.notnull(row['Kovats Identified Compound']) else "",
            'spectral_match_name': row['Spectral Match Name'] if pd.notnull(row['Spectral Match Name']) else "",
            'spectral_similarity': f"{row['Spectral Similarity']:.2f}" if pd.notnull(row['Spectral Similarity']) else ""
        })

    # ✅ Gerar relatório Word
    doc = DocxTemplate(template_path)
    image_path = os.path.join("images", f"{prefix}_chromatogram.png")
    
    if os.path.exists(image_path):
        cromatograma_img = InlineImage(doc, image_path, width=Cm(12))
    else:
        cromatograma_img = None  # ou texto/placeholder

    context = {
        'amostra': prefix,
        'data': pd.Timestamp.today().strftime('%d/%m/%Y'),
        'tabela': tabela_lista,
        'cromatograma_img': cromatograma_img
    }

    try:
        doc.render(context)
        output_docx = os.path.join("results", f"{prefix}_relatorio.docx")
        doc.save(output_docx)
        print(f"✅ Relatório Word salvo: {output_docx}")
    except Exception as e:
        print(f"❌ Erro ao gerar relatório Word para {prefix}: {e}")

    return df_combined  # opcional, se quiser o dataframe resultante


# 📈 PLOTAR Número de Carbonos vs Retention Time
'''fig = px.scatter(df_alkanes, 
                 x='Number of Carbon Atoms', 
                 y=rt_column,
                 title="Curva de Retenção dos Alcanos",
                 labels={'Number of Carbon Atoms': 'Número de Carbonos', rt_column: 'Retention Time (min)'},
                # trendline='ols'
                )  # adiciona linha de tendência opcional

fig.update_traces(mode='markers+lines')
fig.show()'''
