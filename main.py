import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
# CircularGraphicRecord を追加インポート
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

# --- サイドバー ---
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リストをアップロード", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)

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

# --- 視覚化セクション (強化版) ---
if 'results' in st.session_state:
    st.divider()
    st.subheader("🖼️ ベクターマップ表示設定")
    
    col1, col2 = st.columns(2)
    with col1:
        selected_name = st.selectbox("詳細を表示する変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    with col2:
        # マップ形式の切り替えスイッチ
        view_mode = st.radio("表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

    res = next(r for r in st.session_state['results'] if r['mutation_name'] == selected_name)
    
    # 特徴量（Feature）の準備
    features = [
        GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=f"Mutation: {res['mutation_name']}")
    ]
    
    # 新しい制限酵素サイトを特徴量として追加
    if res['New_Sites'] != "None":
        for site in res['New_Sites'].split(", "):
            # 重なり防止のため少しずつ位置をずらす等の処理はライブラリが自動で行うが、
            # ラベルを見やすくするために色を変える
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=site))

    # モードに応じたレコードの作成
    if "Circular" in view_mode:
        record = CircularGraphicRecord(sequence_length=len(res['full_seq']), features=features)
    else:
        record = GraphicRecord(sequence_length=len(res['full_seq']), features=features)

    # 図の描画
    # label_itinerary はラベルの重なりを調整するパラメータ
    fig, ax = plt.subplots(figsize=(10, 8) if "Circular" in view_mode else (10, 3))
    
    if "Circular" in view_mode:
        # 円形プロットの実行
        record.plot(ax=ax)
    else:
        # 直線プロットの実行（ruler付き）
        record.plot(ax=ax, with_ruler=True)
    
    st.pyplot(fig)
    
    # FASTAダウンロード
    fasta_out = f">{res['mutation_name']}_modified\n{res['full_seq']}"
    st.download_button(f"FASTA ({selected_name}) をダウンロード", data=fasta_out, file_name=f"{selected_name}_mod.fasta")