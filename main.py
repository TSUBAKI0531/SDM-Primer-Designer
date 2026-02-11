import streamlit as st
import pandas as pd
import json
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner
from Bio import SeqIO
import matplotlib.pyplot as plt
import xlsxwriter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer Pro", layout="wide")
st.title("🧬 SDM Primer Designer Pro")

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
    fig, ax = plt.subplots(figsize=(6, 5) if is_circular else (8, 2))
    record.plot(ax=ax, with_ruler=not is_circular)
    img_buf = BytesIO(); fig.savefig(img_buf, format='png', bbox_inches='tight', dpi=90); img_buf.seek(0); plt.close(fig)
    return img_buf

if 'custom_features' not in st.session_state: st.session_state['custom_features'] = {}

st.sidebar.header("1. 解析の設定とデータ入力")
f_file = st.sidebar.file_uploader("FASTAファイル", type=["fasta", "fa"])
m_file = st.sidebar.file_uploader("変異リスト (CSV/Excel)", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)
view_mode = st.sidebar.radio("ベクターマップ表示形式", ["Linear (直線状)", "Circular (円形)"], horizontal=True)

st.sidebar.divider()
with st.sidebar.expander("✨ カスタムパーツの管理"):
    n_name, n_seq = st.text_input("パーツ名"), st.text_input("配列")
    if st.button("登録") and n_name and n_seq: st.session_state['custom_features'][n_name] = n_seq.strip().upper()
    if st.session_state['custom_features']:
        st.download_button("JSON保存", json.dumps(st.session_state['custom_features'], indent=4), "features.json", "application/json")
    up_json = st.sidebar.file_uploader("JSON読込", type=["json"], key="json_up")
    if up_json: 
        try: st.session_state['custom_features'].update(json.load(up_json))
        except: pass

if not f_file or not m_file:
    st.info("サイドバーからファイルをアップロードして解析を開始してください。")
    st.stop()

if st.button("🚀 プライマー設計と全解析を実行"):
    with st.spinner("DNA解析、物性計算、マップ描画を統合処理中..."):
        try:
            fasta_content = f_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            detected = designer.detect_features(str(record.seq), st.session_state['custom_features'])
            st.session_state['detected_features'] = detected
            
            df = pd.read_csv(m_file) if m_file.name.endswith('.csv') else pd.read_excel(m_file)
            results = [designer.design(row, target_tm=target_tm) for _, row in df.iterrows() if designer.design(row, target_tm=target_tm)]
            
            if results:
                st.session_state['results'] = results
                res_df = pd.DataFrame(results).drop(['full_seq', 'mut_start', 'mut_end'], axis=1)
                
                def style_dimer(val):
                    if val >= 50: return 'background-color: #ffcccc; color: red'
                    if val >= 40: return 'background-color: #fff4e6'
                    return ''
                
                st.subheader("✅ 総合解析レポート")
                st.dataframe(res_df.style.map(style_dimer, subset=['Dimer_Tm']))
                
                # Excel生成
                output = BytesIO()
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    res_df.to_excel(writer, index=False, sheet_name='SDM Report')
                    ws = writer.sheets['SDM Report']
                    ws.set_column('I:I', 60) # マップ用列
                    for i, res in enumerate(results):
                        ws.set_row(i + 1, 180 if "Circular" in view_mode else 80)
                        img = create_map_image(res, detected, view_mode=view_mode)
                        ws.insert_image(i + 1, 8, f'map_{i}.png', {'image_data': img, 'x_scale': 0.5, 'y_scale': 0.5})
                st.download_button("Excelレポート(画像・物性込)を保存", output.getvalue(), "sdm_full_report.xlsx")

                st.divider()
                st.subheader("📦 プライマー発注用リスト (Tab区切り)")
                order_lines = []
                for res in results:
                    order_lines.append(f"{res['mutation_name']}_F\t{res['fwd_primer']}")
                    order_lines.append(f"{res['mutation_name']}_R\t{res['rev_primer']}")
                st.text_area("そのままコピーして注文サイトへ", "\n".join(order_lines), height=150)
            else: st.warning("設計可能な条件が見つかりませんでした。")
        except Exception as e: st.error(f"解析エラー: {e}")

if 'results' in st.session_state:
    st.divider()
    st.subheader("🧪 実験ベンチ用溶解ガイド")
    prep_data = []
    for res in st.session_state['results']:
        prep_data.append({"Primer Name": f"{res['mutation_name']}_F", "nmol": 25.0})
        prep_data.append({"Primer Name": f"{res['mutation_name']}_R", "nmol": 25.0})
    edited_df = st.data_editor(pd.DataFrame(prep_data), use_container_width=True)
    edited_df['TE (μL) for 100 μM'] = edited_df['nmol'] * 10
    st.dataframe(edited_df)

    st.divider()
    sel = st.selectbox("詳細プレビューする変異を選択", [r['mutation_name'] for r in st.session_state['results']])
    res = next(r for r in st.session_state['results'] if r['mutation_name'] == sel)
    st.image(create_map_image(res, st.session_state['detected_features'], view_mode=view_mode), use_column_width=True)