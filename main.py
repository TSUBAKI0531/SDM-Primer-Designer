import streamlit as st
import pandas as pd
import json
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
# xlsxwriterが明示的に必要になります
import xlsxwriter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

# --- ヘルパー関数: マップ画像の生成 ---
def create_map_image(res, detected_features, view_mode="Linear"):
    """指定された変異結果からマップ画像（PNGデータ）を生成してメモリ上に返す"""
    features = []
    # 検出パーツ（薄青）
    for f in detected_features:
        features.append(GraphicFeature(start=f['start'], end=f['end'], strand=f['strand'], color="#b3d9ff", label=f['name']))
    # 変異箇所（黄）
    features.append(GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=res['mutation_name']))
    # 新規制限酵素サイト（赤）
    if res['New_Sites'] != "None":
        for s in res['New_Sites'].split(", "):
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=s))

    # 描画モード設定（Excel貼付用は視認性が良いLinearを推奨ですが、選択も可能）
    is_circular = "Circular" in view_mode
    record_cls = CircularGraphicRecord if is_circular else GraphicRecord
    record = record_cls(sequence_length=len(res['full_seq']), features=features)

    # 図の生成
    fig, ax = plt.subplots(figsize=(8, 6) if is_circular else (10, 3))
    record.plot(ax=ax, with_ruler=not is_circular)
    
    # メモリバッファに画像を保存
    img_buf = BytesIO()
    fig.savefig(img_buf, format='png', bbox_inches='tight', dpi=100)
    img_buf.seek(0)
    plt.close(fig) # メモリ解放
    return img_buf

# --- セッション初期化 ---
if 'custom_features' not in st.session_state: st.session_state['custom_features'] = {}
if 'detected_features' not in st.session_state: st.session_state['detected_features'] = []

# --- サイドバー ---
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リストをアップロード", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)

st.sidebar.divider()
with st.sidebar.expander("✨ カスタムパーツの管理"):
    new_f_name = st.text_input("パーツ名 (例: GFP)")
    new_f_seq = st.text_input("特徴的な配列 (20bp~)")
    if st.button("登録"):
        if new_f_name and new_f_seq:
            st.session_state['custom_features'][new_f_name] = new_f_seq.strip().upper()
            st.success(f"{new_f_name} を追加しました")
    
    st.write("---")
    if st.session_state['custom_features']:
        st.download_button("JSON書き出し", json.dumps(st.session_state['custom_features'], indent=4), "custom_features.json", "application/json")
    uploaded_json = st.file_uploader("JSON読み込み", type=["json"])
    if uploaded_json:
        try:
            st.session_state['custom_features'].update(json.load(uploaded_json))
            st.success("読み込みました")
        except: st.error("読み込みエラー")

if not fasta_file or not mutations_file:
    st.info("サイドバーからファイルをアップロードしてください。")
    st.stop()

# --- 解析とExcel生成 ---
if st.button("プライマー設計を開始"):
    with st.spinner("解析とレポート生成中..."):
        try:
            # 1. 解析実行
            fasta_content = fasta_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            
            detected = designer.detect_features(str(record.seq), custom_library=st.session_state['custom_features'])
            st.session_state['detected_features'] = detected
            
            df = pd.read_csv(mutations_file) if mutations_file.name.endswith('.csv') else pd.read_excel(mutations_file)
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows() if designer.design(row, target_tm=target_tm)]
            
            if results:
                st.session_state['results'] = results
                result_df_clean = pd.DataFrame(results).drop(['full_seq', 'mut_start', 'mut_end'], axis=1)
                st.subheader("✅ 設計結果")
                st.dataframe(result_df_clean)

                # 2. 画像付きExcelの生成 (xlsxwriterエンジンを使用)
                output = BytesIO()
                # engine='xlsxwriter' を明示的に指定
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    # シート1: プライマーリスト
                    result_df_clean.to_excel(writer, index=False, sheet_name='Primer List')
                    
                    # シート2: ベクターマップ画像
                    workbook = writer.book
                    worksheet_maps = workbook.add_worksheet('Vector Maps')
                    worksheet_maps.set_column('A:A', 30) # A列を広げる
                    
                    current_row = 0
                    # 各変異についてループ
                    for i, res in enumerate(results):
                        # 見出し書き込み
                        worksheet_maps.write(current_row, 0, f"Mutation: {res['mutation_name']}", workbook.add_format({'bold': True, 'font_size': 12}))
                        
                        # 画像生成関数を呼び出し（Excel用はLinearの方が見やすいので固定）
                        img_buf = create_map_image(res, detected, view_mode="Linear")
                        
                        # 画像を挿入 (行間を空けて配置)
                        worksheet_maps.insert_image(current_row + 1, 0, f'map_{i}.png', {'image_data': img_buf, 'x_scale': 0.8, 'y_scale': 0.8})
                        current_row += 22 # 次の画像の配置位置へ移動（画像の高さに応じて調整）

                # ダウンロードボタン
                st.download_button(
                    label="結果Excelをダウンロード（マップ画像付き）",
                    data=output.getvalue(),
                    file_name="sdm_primers_with_maps.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

            else:
                st.warning("条件に合うプライマーが見つかりませんでした。")
                
        except Exception as e:
            st.error(f"エラーが発生しました: {e}")
            import traceback
            st.text(traceback.format_exc()) # 詳細なエラーログを表示

# --- 画面上での視覚化確認 (既存機能) ---
if 'results' in st.session_state and 'detected_features' in st.session_state:
    st.divider()
    col1, col2 = st.columns(2)
    with col1: sel = st.selectbox("詳細を表示する変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    with col2: view = st.radio("表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)
    res = next(r for r in st.session_state['results'] if r['mutation_name'] == sel)
    # 画面表示用に再度画像生成（非効率だがコードは単純化）
    img_buf = create_map_image(res, st.session_state['detected_features'], view_mode=view)
    st.image(img_buf, caption=f"{res['mutation_name']} Map", use_column_width=True)