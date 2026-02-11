import streamlit as st
import pandas as pd
from io import StringIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

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
            
            # --- パーツ検出とデバッグ表示 ---
            detected = designer.detect_features(str(record.seq))
            st.session_state['detected_features'] = detected
            if detected:
                st.success(f"🔍 検出されたパーツ: {', '.join([f['name'] for f in detected])}")
            else:
                st.warning("⚠️ 主要パーツは見つかりませんでした（配列がライブラリと完全一致しません）")
            
            df = pd.read_csv(mutations_file) if mutations_file.name.endswith('.csv') else pd.read_excel(mutations_file)
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows() if designer.design(row, target_tm=target_tm)]
            
            if results:
                st.session_state['results'] = results
                st.dataframe(pd.DataFrame(results).drop(['full_seq', 'mut_start', 'mut_end'], axis=1))
        except Exception as e:
            st.error(f"エラー: {e}")

if 'results' in st.session_state:
    st.divider()
    col1, col2 = st.columns(2)
    with col1:
        sel = st.selectbox("変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    with col2:
        view = st.radio("表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

    res = next(r for r in st.session_state['results'] if r['mutation_name'] == sel)
    features = []
    
    # 検出パーツの描画（薄青）
    for f in st.session_state.get('detected_features', []):
        features.append(GraphicFeature(start=f['start'], end=f['end'], strand=f['strand'], color="#b3d9ff", label=f['name']))
    
    # 変異箇所（黄）と制限酵素（赤）
    features.append(GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=res['mutation_name']))
    if res['New_Sites'] != "None":
        for s in res['New_Sites'].split(", "):
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=s))

    record_cls = CircularGraphicRecord if "Circular" in view else GraphicRecord
    record = record_cls(sequence_length=len(res['full_seq']), features=features)
    fig, ax = plt.subplots(figsize=(10, 8) if "Circular" in view else (10, 3))
    record.plot(ax=ax)
    st.pyplot(fig)