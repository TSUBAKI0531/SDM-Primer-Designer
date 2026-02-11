import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

# --- ファイルアップロード ---
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リストをアップロード", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)

if not fasta_file or not mutations_file:
    st.info("サイドバーからファイルをアップロードしてください。")
    st.stop()

if st.button("プライマー設計を開始"):
    with st.spinner("解析中..."):
        try:
            fasta_content = fasta_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            
            # --- 追加: プラスミドパーツを自動検出 ---
            detected_features = designer.detect_features(str(record.seq))
            st.session_state['detected_features'] = detected_features
            
            df = pd.read_csv(mutations_file) if mutations_file.name.endswith('.csv') else pd.read_excel(mutations_file)
            
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows()]
            results = [r for r in results if r]
            
            if results:
                st.session_state['results'] = results
                result_df = pd.DataFrame(results).drop(['full_seq', 'mut_start', 'mut_end'], axis=1)
                st.subheader("✅ 設計結果")
                st.dataframe(result_df)
            else:
                st.warning("条件に合うプライマーが見つかりませんでした。")
        except Exception as e:
            st.error(f"エラーが発生しました: {e}")

# --- 視覚化セクション ---
if 'results' in st.session_state:
    st.divider()
    st.subheader("🖼️ ベクターマップ表示設定")
    
    col1, col2 = st.columns(2)
    with col1:
        selected_name = st.selectbox("詳細を表示する変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    with col2:
        view_mode = st.radio("表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

    res = next(r for r in st.session_state['results'] if r['mutation_name'] == selected_name)
    
    # 1. 特徴量（Feature）の集約
    features = []
    
    # 自動検出された主要領域を追加（青系）
    if 'detected_features' in st.session_state:
        for f in st.session_state['detected_features']:
            features.append(GraphicFeature(
                start=f['start'], end=f['end'], strand=f['strand'],
                color="#b3d9ff", label=f['name']
            ))

    # 変異箇所を追加（黄色）
    features.append(GraphicFeature(
        start=res['mut_start'], end=res['mut_end'], 
        color="#ffd700", label=f"Mutation: {res['mutation_name']}"
    ))
    
    # 新しい制限酵素サイトを追加（赤）
    if res['New_Sites'] != "None":
        for site in res['New_Sites'].split(", "):
            features.append(GraphicFeature(
                start=res['mut_start'], end=res['mut_start']+1, 
                color="#ff4b4b", label=site
            ))

    # 2. 描画
    if "Circular" in view_mode:
        record = CircularGraphicRecord(sequence_length=len(res['full_seq']), features=features)
    else:
        record = GraphicRecord(sequence_length=len(res['full_seq']), features=features)

    fig, ax = plt.subplots(figsize=(10, 8) if "Circular" in view_mode else (10, 3))
    record.plot(ax=ax)
    st.pyplot(fig)
    
    # ダウンロード
    fasta_out = f">{res['mutation_name']}_modified\n{res['full_seq']}"
    st.download_button(f"FASTA ({selected_name}) をダウンロード", data=fasta_out, file_name=f"{selected_name}_mod.fasta")