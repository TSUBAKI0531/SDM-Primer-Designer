import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO

# ページ設定は必ず一番最初に実行
st.set_page_config(page_title="SDM Primer Designer", layout="wide")

# タイトルは「if」の外に配置し、常に表示させる
st.title("🧬 SDM Primer Designer")
st.write("FASTAファイルと変異リストをアップロードしてプライマーを設計します。")

# サイドバーの設定
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAファイルをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リスト(CSV/Excel)をアップロード", type=["csv", "xlsx"])

# ファイルがまだアップロードされていない時の案内
if not fasta_file or not mutations_file:
    st.info("👈 左側のサイドバーからファイルをアップロードしてください。")
    # ここで処理を止める（以降の設計処理は行わない）
    st.stop()

# --- ここから下はファイルがアップロードされた時のみ実行される ---
st.success("ファイルの読み込みに成功しました！")

# 目標Tm値などのパラメータ
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 60)
method = st.sidebar.selectbox("設計手法", ["overlapping", "back-to-back"])

# --- main.py の設計開始ボタン内の処理 ---
if st.button("プライマー設計を開始"):
    with st.spinner("設計中..."):
        try:
            # FASTAの内容を読み込み
            fasta_content = fasta_file.getvalue().decode("utf-8")
            fasta_io = StringIO(fasta_content)
            record = SeqIO.read(fasta_io, "fasta")
            
            # クラスの初期化（配列文字列を渡す）
            designer = SDMPrimerDesigner(str(record.seq))
            
            # 変異リストの読み込み
            if mutations_file.name.endswith('.csv'):
                df = pd.read_csv(mutations_file)
            else:
                df = pd.read_excel(mutations_file)
            
            # 設計実行（結果の取得）
            results = []
            for _, row in df.iterrows():
                # run_design または design メソッド（作成したクラスに合わせて変更してください）
                res = designer.design(row, method=method, target_tm=target_tm)
                if res:
                    results.append(res)
            
            if results:
                result_df = pd.DataFrame(results)
                st.subheader("設計結果")
                st.dataframe(result_df)
                
                # Excelダウンロード機能
                output = BytesIO()
                result_df.to_excel(output, index=False)
                st.download_button(
                    label="結果をExcelでダウンロード",
                    data=output.getvalue(),
                    file_name="primer_results.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
            else:
                st.warning("条件に合うプライマーが見つかりませんでした。")
                
        except Exception as e:
            st.error(f"解析中にエラーが発生しました: {e}")