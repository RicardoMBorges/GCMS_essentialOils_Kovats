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
    Plota e salva o cromatograma com múltiplos traços (um por par de colunas).
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
        fig.add_trace(go.Scatter(x=df_chrom[x_col], y=df_chrom[y_col], mode='lines', name=x_col.split(':')[0]))

    fig.update_layout(title=f"Cromatograma: {sample_prefix}", xaxis_title="Retention Time (min)", yaxis_title="Base Peak Intensity")

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{sample_prefix}_chromatogram.html")
    fig.write_html(output_file)
    print(f"✅ Gráfico salvo em: {output_file}")
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
    Parseia um arquivo .msp e retorna um dicionário {Name: [(mz, intensidade), ...]}.
    """
    spectra = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        current_name = None
        peaks = []
        for line in f:
            line = line.strip()
            if line.startswith('Name:'):
                if current_name and peaks:
                    spectra[current_name] = peaks
                current_name = line.split(':',1)[1].strip()
                peaks = []
            elif line and any(c.isdigit() for c in line.split()[0]):
                parts = line.split()
                if len(parts) >= 2:
                    mz, intensity = map(float, parts[:2])
                    peaks.append((mz, intensity))
        if current_name and peaks:
            spectra[current_name] = peaks
    return spectra

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
