import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO

st.set_page_config(page_title="SDM Primer Designer", layout="wide")

st.title("🧬 SDM Primer Designer")
st.markdown("農学・創薬研究のための部位特異的変異導入プライマー設計ツール")

# サイドバー設定
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAファイルをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リスト(CSV/Excel)をアップロード", type=["csv", "xlsx"])

target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68) # Tm値を低めにデフォルト設定
method = st.sidebar.selectbox("設計手法", ["overlapping"])

if not fasta_file or not mutations_file:
    st.info("👈 左側のサイドバーからファイルをアップロードしてください。")
    st.stop()

# 解析処理
if st.button("プライマー設計を開始"):
    with st.spinner("設計中..."):
        try:
            # FASTAの読み込み
            fasta_content = fasta_file.getvalue().decode("utf-8")
            if not fasta_content.strip():
                st.error("FASTAファイルが空です。")
                st.stop()
            
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            
            # 変異リストの読み込み
            if mutations_file.name.endswith('.csv'):
                df = pd.read_csv(mutations_file)
            else:
                df = pd.read_excel(mutations_file)
            
            # 設計実行
            results = []
            for _, row in df.iterrows():
                res = designer.design(row, method=method, target_tm=target_tm)
                if res:
                    results.append(res)
            
            if results:
                result_df = pd.DataFrame(results)
                st.subheader("✅ 設計結果")
                
                # 新しい制限酵素サイトがある行をハイライト
                def highlight_sites(val):
                    return 'background-color: #e6fffa' if val != "None" else ''
                
                st.dataframe(result_df.style.applymap(highlight_sites, subset=['New_Sites']))
                
                # Excelダウンロード
                output = BytesIO()
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    result_df.to_excel(writer, index=False, sheet_name='Primers')
                
                st.download_button(
                    label="結果をExcelでダウンロード",
                    data=output.getvalue(),
                    file_name="sdm_primer_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
            else:
                st.warning("条件に合うプライマーが見つかりませんでした。Tm値を下げるか、配列を確認してください。")
                
        except Exception as e:
            st.error(f"解析中にエラーが発生しました: {e}")