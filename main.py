import streamlit as st
import pandas as pd
import json
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
import xlsxwriter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer", layout="wide")
st.title("🧬 SDM Primer Designer")

# --- マップ画像生成関数 ---
def create_map_image(res, detected_features, view_mode="Linear"):
    features = []
    for f in detected_features:
        features.append(GraphicFeature(start=f['start'], end=f['end'], strand=f['strand'], color="#b3d9ff", label=f['name']))
    features.append(GraphicFeature(start=res['mut_start'], end=res['mut_end'], color="#ffd700", label=res['mutation_name']))
    if res['New_Sites'] != "None":
        for s in res['New_Sites'].split(", "):
            features.append(GraphicFeature(start=res['mut_start'], end=res['mut_start']+1, color="#ff4b4b", label=s))

    is_circular = "Circular" in view_mode
    record_cls = CircularGraphicRecord if is_circular else GraphicRecord
    record = record_cls(sequence_length=len(res['full_seq']), features=features)

    # Excel用にサイズを調整
    fig, ax = plt.subplots(figsize=(6, 5) if is_circular else (8, 2))
    record.plot(ax=ax, with_ruler=not is_circular)
    
    img_buf = BytesIO()
    fig.savefig(img_buf, format='png', bbox_inches='tight', dpi=90)
    img_buf.seek(0)
    plt.close(fig)
    return img_buf

if 'custom_features' not in st.session_state: st.session_state['custom_features'] = {}

# --- サイドバー設定 ---
st.sidebar.header("1. 入力ファイルのアップロード")
fasta_file = st.sidebar.file_uploader("FASTAをアップロード", type=["fasta", "fa"])
mutations_file = st.sidebar.file_uploader("変異リストをアップロード", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)

# 表示モードをサイドバーへ移動（Excel出力にも反映させるため）
view_mode = st.sidebar.radio("ベクターマップ表示モード", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

st.sidebar.divider()
with st.sidebar.expander("✨ カスタムパーツの管理"):
    new_f_name = st.text_input("パーツ名")
    new_f_seq = st.text_input("配列")
    if st.button("登録"):
        if new_f_name and new_f_seq:
            st.session_state['custom_features'][new_f_name] = new_f_seq.strip().upper()

# --- 解析とExcel生成 ---
if st.button("プライマー設計を開始"):
    with st.spinner("解析とレポート生成中..."):
        try:
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
                
                # Excel生成処理
                output = BytesIO()
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    result_df_clean.to_excel(writer, index=False, sheet_name='SDM Report')
                    
                    workbook = writer.book
                    worksheet = writer.sheets['SDM Report']
                    
                    # 画像用の列(H列)と行の高さを調整
                    worksheet.set_column('H:H', 60) # 画像を表示する列の幅
                    header_format = workbook.add_format({'bold': True, 'bg_color': '#D9EAD3', 'border': 1})
                    worksheet.write(0, 7, 'Vector Map View', header_format)
                    
                    # 各行に画像を配置
                    row_height = 180 if "Circular" in view_mode else 80
                    for i, res in enumerate(results):
                        worksheet.set_row(i + 1, row_height)
                        img_buf = create_map_image(res, detected, view_mode=view_mode)
                        
                        # セル内に収まるようにスケール調整
                        y_offset = 5
                        worksheet.insert_image(i + 1, 7, f'map_{i}.png', 
                                               {'image_data': img_buf, 'x_scale': 0.5, 'y_scale': 0.5, 'y_offset': y_offset})

                st.success("解析完了！")
                st.dataframe(result_df_clean)
                st.download_button("Excelレポートをダウンロード", output.getvalue(), "sdm_analysis_report.xlsx")
        except Exception as e:
            st.error(f"エラー: {e}")

# --- 画面上での確認 ---
if 'results' in st.session_state:
    st.divider()
    sel = st.selectbox("詳細を表示する変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    res = next(r for r in st.session_state['results'] if r['mutation_name'] == sel)
    img_buf = create_map_image(res, st.session_state['detected_features'], view_mode=view_mode)
    st.image(img_buf, caption=f"Map: {sel} ({view_mode})", use_column_width=True)