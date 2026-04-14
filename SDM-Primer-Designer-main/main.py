import streamlit as st
import pandas as pd
import json
from io import StringIO, BytesIO
from sdm_designer import SDMPrimerDesigner, PCRProtocol
from Bio import SeqIO
import matplotlib.pyplot as plt
import xlsxwriter
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord

st.set_page_config(page_title="SDM Primer Designer Pro", layout="wide")
st.title("🧬 SDM Primer Designer Pro")


# =====================================================================
#  ユーティリティ関数
# =====================================================================
def create_map_image(res, detected_features, view_mode="Linear"):
    """ベクターマップ画像を生成して BytesIO で返す"""
    features = []
    for f in detected_features:
        features.append(
            GraphicFeature(start=f['start'], end=f['end'], strand=f['strand'],
                           color="#b3d9ff", label=f['name'])
        )
    features.append(
        GraphicFeature(start=res['mut_start'], end=res['mut_end'],
                       color="#ffd700", label=res['mutation_name'])
    )
    if res['New_Sites'] != "None":
        for s in res['New_Sites'].split(", "):
            features.append(
                GraphicFeature(start=res['mut_start'], end=res['mut_start'] + 1,
                               color="#ff4b4b", label=s)
            )
    is_circular = "Circular" in view_mode
    record_cls = CircularGraphicRecord if is_circular else GraphicRecord
    record = record_cls(sequence_length=len(res['full_seq']), features=features)
    fig, ax = plt.subplots(figsize=(6, 5) if is_circular else (8, 2))
    record.plot(ax=ax, with_ruler=not is_circular)
    img_buf = BytesIO()
    fig.savefig(img_buf, format='png', bbox_inches='tight', dpi=90)
    img_buf.seek(0)
    plt.close(fig)
    return img_buf


def render_protocol_section(protocol: PCRProtocol):
    """単一変異の PCR プロトコルを Streamlit UI に描画する"""
    size_kb = round(protocol.product_size_bp / 1000, 1)

    st.markdown(f"#### 🔬 {protocol.mutation_name} — PCR プロトコル")
    st.caption(f"産物サイズ: {protocol.product_size_bp} bp ({size_kb} kb) / "
               f"Primer Tm (NN法): {protocol.primer_tm}℃ / "
               f"Primer長: {protocol.primer_length} mer")

    # --- 反応液組成 & サイクリング条件 (2列) ---
    col_mix, col_cycle = st.columns(2)
    with col_mix:
        st.markdown("**PCR 反応液組成 (50 μL 系)**")
        mix_df = pd.DataFrame(protocol.to_reaction_table())
        st.dataframe(mix_df, hide_index=True, use_container_width=True)
        st.caption("※ PCR 反応液は室温調製可。酵素は氷上保管。")

    with col_cycle:
        st.markdown("**サイクリング条件**")
        st.code(protocol.to_cycling_text(), language=None)

    # --- 鋳型量ガイド / トラブルシューティング ---
    with st.expander("📋 鋳型 DNA 量ガイド / トラブルシューティング", expanded=False):
        guide_col, trouble_col = st.columns(2)
        with guide_col:
            st.markdown("**鋳型 DNA 量 (50 μL 反応系)**")
            template_data = [
                {"鋳型種類": "プラスミド DNA", "至適量": "10 pg 〜 1 ng"},
                {"鋳型種類": "λDNA", "至適量": "10 pg 〜 10 ng"},
                {"鋳型種類": "大腸菌ゲノム DNA", "至適量": "100 pg 〜 200 ng"},
                {"鋳型種類": "ヒトゲノム DNA", "至適量": "5 〜 200 ng"},
            ]
            st.dataframe(pd.DataFrame(template_data),
                         hide_index=True, use_container_width=True)
            st.caption("200 ng 超の場合は伸長時間を 30〜60 sec/kb に延長")

        with trouble_col:
            st.markdown("**トラブルシューティング**")
            trouble_data = [
                {"現象": "増幅しない",
                 "対策": "伸長 10〜60 sec/kb, サイクル 35〜40, アニーリング 15 sec"},
                {"現象": "増幅効率が悪い",
                 "対策": "アニーリング温度 2℃ずつ低下, 反応液 25 μL"},
                {"現象": "スメア",
                 "対策": "アニーリング 5 sec, 温度 58〜63℃, 2-step PCR"},
                {"現象": "エキストラバンド",
                 "対策": "鋳型DNA減量, サイクル 25〜30"},
            ]
            st.dataframe(pd.DataFrame(trouble_data),
                         hide_index=True, use_container_width=True)

    # --- 注意事項 ---
    if protocol.notes:
        with st.expander("⚠ この変異に対する注意事項", expanded=False):
            for note in protocol.notes:
                st.markdown(f"- {note}")


# =====================================================================
#  Session State 初期化
# =====================================================================
if 'custom_features' not in st.session_state:
    st.session_state['custom_features'] = {}

# =====================================================================
#  サイドバー
# =====================================================================
st.sidebar.header("1. 解析の設定とデータ入力")
f_file = st.sidebar.file_uploader("FASTAファイル", type=["fasta", "fa"])
m_file = st.sidebar.file_uploader("変異リスト", type=["csv", "xlsx"])
target_tm = st.sidebar.slider("目標 Tm値 (°C)", 50, 85, 68)
view_mode = st.sidebar.radio("ベクターマップ表示形式",
                              ["Linear (直線状)", "Circular (円形)"], horizontal=True)

st.sidebar.divider()
with st.sidebar.expander("✨ カスタムパーツの管理"):
    n_name = st.text_input("パーツ名")
    n_seq = st.text_input("配列")
    if st.button("登録") and n_name and n_seq:
        st.session_state['custom_features'][n_name] = n_seq.strip().upper()
    if st.session_state['custom_features']:
        st.download_button("JSON保存",
                           json.dumps(st.session_state['custom_features'], indent=4),
                           "features.json", "application/json")
    up_json = st.sidebar.file_uploader("JSON読込", type=["json"], key="json_up")
    if up_json:
        try:
            st.session_state['custom_features'].update(json.load(up_json))
        except Exception:
            pass

if not f_file or not m_file:
    st.info("ファイルをアップロードしてください。")
    st.stop()

# =====================================================================
#  メイン解析実行
# =====================================================================
if st.button("🚀 プライマー設計と全解析を実行"):
    with st.spinner("DNA解析、物性計算、マップ描画を統合処理中..."):
        try:
            fasta_content = f_file.getvalue().decode("utf-8")
            record = SeqIO.read(StringIO(fasta_content), "fasta")
            designer = SDMPrimerDesigner(str(record.seq))
            detected = designer.detect_features(
                str(record.seq), st.session_state['custom_features']
            )
            st.session_state['detected_features'] = detected
            st.session_state['designer'] = designer

            df = (pd.read_csv(m_file)
                  if m_file.name.endswith('.csv')
                  else pd.read_excel(m_file))

            # ★ 修正: design() の二重呼び出しバグを解消
            results = []
            for _, row in df.iterrows():
                res = designer.design(row, target_tm=target_tm)
                if res:
                    results.append(res)

            if results:
                st.session_state['results'] = results
                res_df = pd.DataFrame(results).drop(
                    ['full_seq', 'mut_start', 'mut_end'], axis=1
                )

                def style_dimer(val):
                    if val >= 50:
                        return 'background-color: #ffcccc; color: red'
                    if val >= 40:
                        return 'background-color: #fff4e6'
                    return ''

                st.subheader("✅ 総合解析レポート")
                st.dataframe(
                    res_df.style.map(style_dimer, subset=['Dimer_Tm'])
                )

                # Excel 生成
                output = BytesIO()
                with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                    res_df.to_excel(writer, index=False, sheet_name='SDM Report')
                    ws = writer.sheets['SDM Report']
                    ws.set_column('I:I', 60)
                    for i, res in enumerate(results):
                        ws.set_row(
                            i + 1, 180 if "Circular" in view_mode else 80
                        )
                        img = create_map_image(res, detected, view_mode=view_mode)
                        ws.insert_image(
                            i + 1, 8, f'map_{i}.png',
                            {'image_data': img, 'x_scale': 0.5, 'y_scale': 0.5}
                        )
                st.download_button("Excelレポート(画像・物性込)を保存",
                                   output.getvalue(), "sdm_full_report.xlsx")

                # --- 発注用フォーマット ---
                st.divider()
                st.subheader("📦 プライマー発注用リスト (Tab区切り)")
                order_lines = []
                for res in results:
                    order_lines.append(
                        f"{res['mutation_name']}_F\t{res['fwd_primer']}"
                    )
                    order_lines.append(
                        f"{res['mutation_name']}_R\t{res['rev_primer']}"
                    )
                st.text_area("そのままコピーして注文サイトへ",
                             "\n".join(order_lines), height=150)
            else:
                st.warning("設計可能な条件が見つかりませんでした。")
        except Exception as e:
            st.error(f"解析エラー: {e}")

# =====================================================================
#  到着後の溶解ガイド
# =====================================================================
if 'results' in st.session_state:
    st.divider()
    st.subheader("🧪 プライマー調製ガイド (到着後用)")
    prep_data = []
    for res in st.session_state['results']:
        prep_data.append({
            "Primer Name": f"{res['mutation_name']}_F", "nmol": 25.0
        })
        prep_data.append({
            "Primer Name": f"{res['mutation_name']}_R", "nmol": 25.0
        })
    edited_df = st.data_editor(pd.DataFrame(prep_data), use_container_width=True)
    edited_df['TE (μL) for 100 μM'] = edited_df['nmol'] * 10
    st.dataframe(edited_df)

    # ==================================================================
    #  🆕 PCR プロトコル自動出力 (PrimeSTAR Max 準拠)
    # ==================================================================
    st.divider()
    st.subheader("📝 PCR プロトコル (PrimeSTAR Max DNA Polymerase 準拠)")
    st.caption(
        "PrimeSTAR Max DNA Polymerase (R045A) の説明書に基づき、"
        "各変異の設計結果から反応液組成・サイクリング条件を自動生成します。"
    )

    designer = st.session_state.get('designer')
    if designer is None:
        designer = SDMPrimerDesigner("A")  # フォールバック用ダミー

    protocols = [designer.generate_protocol(res)
                 for res in st.session_state['results']]
    st.session_state['protocols'] = protocols

    mutation_names = [p.mutation_name for p in protocols]

    tab_all, tab_single = st.tabs(
        ["📄 一括プロトコル出力", "🔍 個別プロトコル表示"]
    )

    with tab_all:
        # 全変異の一括テキスト
        all_protocol_text = "\n\n".join([p.to_full_text() for p in protocols])
        st.download_button(
            "📥 全変異プロトコル (.txt) をダウンロード",
            all_protocol_text,
            file_name="sdm_pcr_protocols.txt",
            mime="text/plain",
        )
        # サマリーテーブル
        summary_data = []
        for p in protocols:
            summary_data.append({
                "変異名": p.mutation_name,
                "産物サイズ": f"{p.product_size_bp} bp",
                "変性": f"{p.denature_temp}℃ {p.denature_sec}s",
                "アニーリング": f"{p.annealing_temp}℃ {p.annealing_sec}s",
                "伸長": f"{p.extension_temp}℃ {p.extension_sec}s",
                "サイクル": p.cycles,
                "Primer Tm": f"{p.primer_tm}℃",
            })
        st.dataframe(pd.DataFrame(summary_data),
                     hide_index=True, use_container_width=True)

    with tab_single:
        sel_name = st.selectbox("プロトコルを表示する変異を選択",
                                mutation_names, key="protocol_select")
        sel_protocol = next(
            p for p in protocols if p.mutation_name == sel_name
        )
        render_protocol_section(sel_protocol)

        st.download_button(
            f"📥 {sel_name} のプロトコル (.txt)",
            sel_protocol.to_full_text(),
            file_name=f"protocol_{sel_name}.txt",
            mime="text/plain",
            key=f"dl_{sel_name}",
        )

    # ==================================================================
    #  ベクターマップ詳細プレビュー
    # ==================================================================
    st.divider()
    st.subheader("🗺 ベクターマップ詳細プレビュー")
    sel_map = st.selectbox(
        "詳細プレビューする変異を選択",
        [r['mutation_name'] for r in st.session_state['results']],
        key="map_select",
    )
    res = next(
        r for r in st.session_state['results'] if r['mutation_name'] == sel_map
    )
    st.image(
        create_map_image(res, st.session_state['detected_features'],
                         view_mode=view_mode),
        use_container_width=True,
    )
