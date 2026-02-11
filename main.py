import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord

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
            
            df = pd.read_csv(mutations_file) if mutations_file.name.endswith('.csv') else pd.read_excel(mutations_file)
            
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows()]
            results = [r for r in results if r]
            
            if results:
                st.session_state['results'] = results # 視覚化のために保存
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
    st.subheader("🖼️ 変異後配列の視覚化とダウンロード")
    
    selected_name = st.selectbox("詳細を表示する変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    res = next(r for r in st.session_state['results'] if r['mutation_name'] == selected_name)
    
    # 1. マップの描画
    features = [GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=f"Mutation: {res['mutation_name']}")]
    if res['New_Sites'] != "None":
        for site in res['New_Sites'].split(", "):
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=site))

    record = GraphicRecord(sequence_length=len(res['full_seq']), features=features)
    fig, ax = plt.subplots(figsize=(10, 3))
    record.plot(ax=ax, with_ruler=True)
    st.pyplot(fig)
    
    # 2. 配列表示とダウンロード
    fasta_out = f">{res['mutation_name']}_modified\n{res['full_seq']}"
    st.text_area("変異後の全配列", fasta_out, height=150)
    st.download_button("FASTAをダウンロード", data=fasta_out, file_name=f"{res['mutation_name']}_mod.fasta")