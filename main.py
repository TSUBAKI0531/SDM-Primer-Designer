import streamlit as st
import pandas as pd
from io import StringIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

# カスタムパーツの一時保存用
if 'custom_features' not in st.session_state:
    st.session_state['custom_features'] = {}

# --- サイドバー ---
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リストをアップロード", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)

st.sidebar.divider()
with st.sidebar.expander("✨ カスタムパーツの登録"):
    new_f_name = st.text_input("パーツ名 (例: GFP)")
    new_f_seq = st.text_input("特徴的な配列 (20bp~)")
    if st.button("登録"):
        if new_f_name and new_f_seq:
            st.session_state['custom_features'][new_f_name] = new_f_seq.strip().upper()
            st.success(f"{new_f_name} を追加しました")
        else: st.error("入力を確認してください")
    
    # 登録済みリストの表示と削除
    if st.session_state['custom_features']:
        st.write("---")
        for n in list(st.session_state['custom_features'].keys()):
            c1, c2 = st.columns([4, 1])
            c1.caption(n)
            if c2.button("🗑️", key=f"del_{n}"):
                del st.session_state['custom_features'][n]
                st.rerun()

if not fasta_file or not mutations_file:
    st.info("サイドバーからファイルをアップロードしてください。")
    st.stop()

# --- 解析実行 ---
if st.button("プライマー設計を開始"):
    with st.spinner("解析中..."):
        try:
            fasta_content = fasta_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            
            # カスタムパーツを含めて検出
            detected = designer.detect_features(str(record.seq), custom_library=st.session_state['custom_features'])
            st.session_state['detected_features'] = detected
            
            df = pd.read_csv(mutations_file) if mutations_file.name.endswith('.csv') else pd.read_excel(mutations_file)
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows() if designer.design(row, target_tm=target_tm)]
            
            if results:
                st.session_state['results'] = results
                st.dataframe(pd.DataFrame(results).drop(['full_seq', 'mut_start', 'mut_end'], axis=1))
        except Exception as e:
            st.error(f"エラー: {e}")

# --- 視覚化 ---
if 'results' in st.session_state:
    st.divider()
    col1, col2 = st.columns(2)
    with col1:
        sel = st.selectbox("変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    with col2:
        view = st.radio("表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

    res = next(r for r in st.session_state['results'] if r['mutation_name'] == sel)
    features = []
    # 検出された全てのパーツを描画
    for f in st.session_state.get('detected_features', []):
        features.append(GraphicFeature(start=f['start'], end=f['end'], strand=f['strand'], color="#b3d9ff", label=f['name']))
    
    features.append(GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=res['mutation_name']))
    if res['New_Sites'] != "None":
        for s in res['New_Sites'].split(", "):
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=s))

    record_cls = CircularGraphicRecord if "Circular" in view else GraphicRecord
    record = record_cls(sequence_length=len(res['full_seq']), features=features)
    fig, ax = plt.subplots(figsize=(10, 8) if "Circular" in view else (10, 3))
    record.plot(ax=ax)
    st.pyplot(fig)