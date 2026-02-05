SDM-Primer-Designer 🧬

SDM-Primer-Designer は、部位特異的変異導入（Site-Directed Mutagenesis）試験のためのプライマー設計を自動化する Python ツールです。
アミノ酸レベルでの変異指定、コドン最適化、制限酵素サイトの自動解析機能を備え、大規模な変異体ライブラリ作製を強力にサポートします。


✨ 主な特徴 (Features)
一括設計 (Batch Processing): ExcelやCSVリストから数百件の変異プライマーを一度に設計。
アミノ酸指定 (AA-level Specification): 塩基ではなくアミノ酸での変異指定が可能（例：3番目のAlaをAspに置換）。
コドン最適化 (Codon Optimization): 使用宿主（デフォルト：大腸菌）に合わせて頻度の高いコドンを自動選択。
制限酵素チェック (Restriction Site Analysis): 変異導入によって「新しく出現したサイト」や「消失したサイト」を自動判定。スクリーニングを容易にします。
複数の設計モード: * Overlapping: QuikChange法などに適した完全に重なるプライマー。
Back-to-back: NEB Q5法などに適した、外側を向いたプライマーペア。
詳細なログ出力: 各設計の成功・失敗理由をログファイルに自動保存。


🚀 インストール (Installation)
git clone https://github.com/あなたのユーザー名/SDM-Primer-Designer.git
cd SDM-Primer-Designer
pip install -r requirements.txt


🛠 使い方 (Usage)
1. 入力ファイルの準備
template.fasta: 鋳型となるDNA配列ファイル。
mutations.csv: 設計したい変異のリスト。以下の列を含めてください：
mutation_name,mode,aa_pos,target_aa,insert_seq,del_len
M1_Y3D,sub,3,D,,0
M2_Ins,ins,6,,GGGGGG,0
M3_Del,del,10,,,2


2. 実行
from sdm_designer import batch_process

# 実行
batch_process(
    fasta_path="examples/template.fasta",
    csv_path="examples/mutations.csv",
    method='overlapping' # または 'back-to-back'
)


📊 出力 (Output)
実行後、primer_results.xlsx が生成されます。

fwd/rev_primer: 設計されたプライマー配列
tm: 計算されたTm値
new_sites: 変異によって新しく作成された制限酵素サイト
lost_sites: 変異によって消失した制限酵素サイト


🧪 依存関係 (Dependencies)
Biopython
pandas
openpyxl


📄 ライセンス (License)
このプロジェクトは MIT License の下で公開されています。


🤝 貢献 (Contributing)
バグ報告や機能改善の提案は、Issues または Pull Request までお気軽にどうぞ！


💡 工夫したポイント（これを載せると評価が上がります）
研究者の視点: 「制限酵素サイトの消失・出現」をリストアップすることで、高額なシーケンス解析に回す前の「一次スクリーニング（酵素消化チェック）」を圧倒的に楽にします。
堅牢な設計: 実験条件に合うプライマーが見つからない場合、どのパラメータ（Tmなど）が原因かをログに出力し、試行錯誤の時間を短縮します。