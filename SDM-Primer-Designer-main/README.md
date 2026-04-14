# 🧬 SDM Primer Designer Pro

部位特異的変異導入（Site-Directed Mutagenesis; SDM）におけるプライマー設計から、実験ベンチでの試薬調製までを一貫してサポートする、バイオ・農学研究者向けのWebアプリケーションです。

## 🌟 主な機能

本アプリは、SDMのワークフローにおける「6つのステップ」を自動化し、ヒューマンエラーを最小限に抑えます。

1.  **高度なプライマー設計**
    * 置換 (Substitution)、挿入 (Insertion)、欠失 (Deletion) の全モードに対応。
    * オーバーラップ法を用いた Tm値 ($T_m$) の自動最適化。
    * **自己アニーリング (Self-Dimer $T_m$)** の算出と、PCR阻害リスクの自動警告。

2.  **インタラクティブな可視化**
    * `dna_features_viewer` を用いたベクターマップ表示（直線状/円形切り替え）。
    * 標準パーツ（AmpR, ori等）に加え、**カスタムパーツ（GFP, Tag等）** の自由な登録・JSON管理。

3.  **プロ仕様のレポート出力**
    * 設計データ、物理特性、制限酵素情報を網羅した **画像埋め込み型Excelレポート** の生成。
    * 表示モード（Linear/Circular）に同期したマップ画像の自動配置。

4.  **発注プロセスの自動化**
    * 合成会社（IDT, ユーロフィン等）のフォームにそのまま貼り付け可能な **タブ区切りテキスト形式** での一括出力。

5.  **実験ベンチ支援**
    * **PCR産物サイズ (Product Size)** の自動予測（電気泳動での確認用）。
    * 到着したプライマーの nmol 数から **$100 \mu M$ ストック液** 作成に必要な溶媒量を算出する溶解ガイド。

6.  **🆕 PCR プロトコル自動出力 (PrimeSTAR Max 準拠)**
    * TaKaRa PrimeSTAR Max DNA Polymerase (R045A) の説明書に基づく **PCR 反応液組成・サイクリング条件の自動生成**。
    * 産物サイズから **伸長時間 (5 sec/kb)** を自動算出、プライマー Tm から **アニーリング時間 (5 or 15 sec)** を自動判定。
    * 6 kb 超の産物に対する **2-step PCR 代替条件** の自動提示。
    * **鋳型 DNA 量ガイド**（プラスミド/ゲノム/λDNA/cDNA 別の至適量）と **トラブルシューティング表** を内蔵。
    * 変異ごとの個別プロトコルおよび **全変異一括プロトコル (.txt)** のダウンロード対応。

## 🛠 テックスタック

* **Language:** Python 3.12+
* **Framework:** Streamlit
* **Bioinformatics:** Biopython (Seq, SeqUtils, Restriction, SeqIO)
* **Visualization:** Matplotlib, dna_features_viewer
* **Reporting:** Pandas, XlsxWriter

## 🚀 ローカルでのセットアップ

```bash
# リポジトリのクローン
git clone [https://github.com/TSUBAKI0531/Vector2Fold.git](https://github.com/TSUBAKI0531/Vector2Fold.git)
cd Vector2Fold

# 仮想環境の作成と起動
python -m venv .venv
source .venv/bin/activate  # Windowsの場合は .venv\Scripts\activate

# 依存ライブラリのインストール
pip install -r requirements.txt
pip install xlsxwriter dna_features_viewer matplotlib biopython streamlit

# アプリの起動
streamlit run main.py
📝 使用方法
FASTAファイルのアップロード: 鋳型となるプラスミド配列を読み込みます。

変異リストのアップロード: CSV/Excel形式で変異情報を入力します。

カスタムパーツの登録: 必要に応じて、JSONから検索したい重要パーツを読み込みます。

解析実行: プライマーが自動設計され、物理特性とマップが表示されます。

レポート保存: Excelレポートをダウンロードし、発注用テキストをコピーします。


---

## 2. GitHub へのプッシュ用コード

READMEの更新を反映させます。

```bash
# 変更を反映
git add README.md

# コミット（「Pro版の機能を反映」）
git commit -m "Update README: Reflect all Pro features (Dimer check, Excel reporting, Dissolution guide)"

# プッシュ
git push origin main