from dataclasses import dataclass, field
from typing import Optional
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import molecular_weight
from Bio import Restriction
import math


@dataclass
class PCRProtocol:
    """PrimeSTAR Max DNA Polymerase ベースの PCR プロトコル"""
    mutation_name: str
    product_size_bp: int
    primer_tm: float
    primer_length: int
    fwd_primer: str
    rev_primer: str

    # --- 反応液組成 (50 μL 系) ---
    premix_vol: float = 25.0          # μL
    primer_pmol: float = 10.0         # pmol per primer
    primer_conc_uM: float = 0.2       # μM final
    template_ng: float = 1.0          # ng (プラスミド推奨)
    water_vol: float = 0.0            # 計算で埋める

    # --- サイクリング条件 ---
    denature_temp: float = 98.0
    denature_sec: int = 10
    annealing_temp: float = 55.0
    annealing_sec: int = 5
    extension_temp: float = 72.0
    extension_sec: int = 5            # 計算で埋める
    cycles: int = 30
    use_2step: bool = False

    # --- 注意事項 ---
    notes: list = field(default_factory=list)

    def __post_init__(self):
        self._calc_extension_time()
        self._calc_annealing()
        self._calc_water()
        self._generate_notes()

    # ---- 内部計算 ----
    def _calc_extension_time(self):
        """産物サイズから伸長時間を算出 (5 sec/kb, プラスミド≤200 ng 想定)"""
        size_kb = self.product_size_bp / 1000
        raw_sec = math.ceil(size_kb * 5)
        self.extension_sec = max(5, raw_sec)

    def _calc_annealing(self):
        """簡易 Tm 式で判定: Tm ≥ 55 → 5 sec, Tm < 55 → 15 sec"""
        simple_tm = self._wallace_tm(self.fwd_primer)
        if simple_tm >= 55:
            self.annealing_sec = 5
        else:
            self.annealing_sec = 15

    def _calc_water(self):
        primer_vol_each = self.primer_pmol / (10)  # 10 μM stock → μL
        self.water_vol = round(50.0 - self.premix_vol - primer_vol_each * 2 - 1.0, 1)

    @staticmethod
    def _wallace_tm(seq: str) -> float:
        """Wallace rule: Tm = 2(A+T) + 4(G+C) - 5 (≤25 mer 向け)"""
        s = seq.upper()[:25]
        at = s.count('A') + s.count('T')
        gc = s.count('G') + s.count('C')
        return 2 * at + 4 * gc - 5

    def _generate_notes(self):
        self.notes = []
        if self.product_size_bp > 10000:
            self.notes.append(
                "⚠ 増幅鎖長が 10 kb を超えています。伸長時間を 15〜30 sec/kb に延長し、"
                "サイクル数を 35 に増やすことを検討してください。"
            )
        if self.primer_length > 25:
            self.notes.append(
                "プライマーが 25 mer を超えるため、アニーリング時間は 5 sec を推奨します。"
            )
        if self.primer_tm > 75:
            self.notes.append(
                "Tm が高めです。非特異増幅が見られる場合は 2-step PCR (98℃/68℃) を検討してください。"
            )
        self.notes.append(
            "電気泳動には TAE Buffer を推奨します（TBE では裾広がりになる場合あり）。"
        )
        self.notes.append(
            "PCR 産物は平滑末端です。平滑末端ベクターへのクローニングが可能です。"
        )
        self.notes.append(
            "制限酵素処理・シーケンス前には PCR Clean-up でタンパク質を除去してください "
            "（3'→5' exonuclease 活性の残存を防ぐため）。"
        )

    # ---- 出力フォーマット ----
    def to_reaction_table(self) -> list[dict]:
        """反応液組成テーブル用の辞書リスト"""
        primer_vol = round(self.primer_pmol / 10, 1)  # 10 μM stock
        return [
            {"試薬": "PrimeSTAR Max Premix (2×)", "使用量": f"{self.premix_vol} μL",
             "最終濃度": "1×"},
            {"試薬": f"Primer F ({self.mutation_name}_F)", "使用量": f"{primer_vol} μL (={self.primer_pmol} pmol)",
             "最終濃度": f"{self.primer_conc_uM} μM"},
            {"試薬": f"Primer R ({self.mutation_name}_R)", "使用量": f"{primer_vol} μL (={self.primer_pmol} pmol)",
             "最終濃度": f"{self.primer_conc_uM} μM"},
            {"試薬": "Template DNA (プラスミド)", "使用量": f"{self.template_ng} ng",
             "最終濃度": "—"},
            {"試薬": "滅菌精製水", "使用量": f"up to 50 μL",
             "最終濃度": "—"},
        ]

    def to_cycling_text(self) -> str:
        """サイクリング条件のテキスト表現"""
        size_kb = round(self.product_size_bp / 1000, 1)
        lines = [
            f"【 3-Step Protocol — {self.mutation_name} ({size_kb} kb) 】",
            "",
            f"  {self.denature_temp}℃   {self.denature_sec} sec        ┐",
            f"  {self.annealing_temp}℃    {self.annealing_sec} sec         ├ {self.cycles} cycles",
            f"  {self.extension_temp}℃   {self.extension_sec} sec ({size_kb} kb × 5 sec/kb) ┘",
        ]
        if self.product_size_bp > 6000:
            lines.append("")
            lines.append("【 代替: 2-Step Protocol 】")
            ext_2step = max(30, math.ceil(size_kb * 30))
            lines.append(f"  98℃   10 sec   ┐")
            lines.append(f"  68℃   {ext_2step} sec   ├ {self.cycles} cycles")
            lines.append(f"                  ┘")
        return "\n".join(lines)

    def to_full_text(self) -> str:
        """プロトコル全文テキスト（ダウンロード用）"""
        sep = "=" * 60
        size_kb = round(self.product_size_bp / 1000, 1)
        primer_vol = round(self.primer_pmol / 10, 1)

        sections = []
        sections.append(f"{sep}")
        sections.append(f"  SDM PCR Protocol — {self.mutation_name}")
        sections.append(f"  PrimeSTAR Max DNA Polymerase (R045A)")
        sections.append(f"{sep}")
        sections.append("")
        sections.append(f"■ 産物サイズ: {self.product_size_bp} bp ({size_kb} kb)")
        sections.append(f"■ プライマー Tm (NN法): {self.primer_tm}℃")
        sections.append("")

        # 反応液組成
        sections.append("━━━ 1. PCR 反応液組成 (Total 50 μL) ━━━")
        sections.append("")
        sections.append(f"  {'試薬':<40} {'使用量':<20} {'最終濃度'}")
        sections.append(f"  {'─'*40} {'─'*20} {'─'*10}")
        for row in self.to_reaction_table():
            sections.append(f"  {row['試薬']:<40} {row['使用量']:<20} {row['最終濃度']}")
        sections.append("")
        sections.append("  ※ PCR 反応液の調製は室温でも可能です。")
        sections.append("    ただし酵素などの各試薬は氷上に置いてご使用ください。")
        sections.append("")

        # サイクリング
        sections.append("━━━ 2. PCR 条件 ━━━")
        sections.append("")
        sections.append(self.to_cycling_text())
        sections.append("")
        sections.append(f"  ● 変性条件: 98℃ 5〜10 sec 推奨 (94℃の場合は 10〜15 sec)")
        sections.append(f"  ● アニーリング温度: まず 55℃で試行")
        sections.append(f"  ● アニーリング時間: Wallace Tm ≥ 55℃ → 5 sec / < 55℃ → 15 sec")
        sections.append(f"    (25 mer 超のプライマーは常に 5 sec)")
        sections.append("")

        # テンプレート量ガイド
        sections.append("━━━ 3. 鋳型 DNA 量ガイド (50 μL 反応系) ━━━")
        sections.append("")
        sections.append("  プラスミド DNA       10 pg 〜 1 ng")
        sections.append("  λDNA               10 pg 〜 10 ng")
        sections.append("  大腸菌ゲノム DNA    100 pg 〜 200 ng")
        sections.append("  ヒトゲノム DNA       5 ng 〜 200 ng")
        sections.append("")
        sections.append("  ※ 200 ng 超の場合は伸長時間を 30〜60 sec/kb に延長")
        sections.append("")

        # ポスト PCR
        sections.append("━━━ 4. PCR 後の処理 ━━━")
        sections.append("")
        sections.append("  ● 電気泳動: TAE Buffer 推奨")
        sections.append(f"  ● 期待バンドサイズ: {self.product_size_bp} bp")
        sections.append("  ● クローニング: 平滑末端ベクターに直接挿入可能")
        sections.append("    (T-vector の場合は dA 付加反応が必要)")
        sections.append("  ● 制限酵素処理/シーケンス前: PCR Clean-up 必須")
        sections.append("    (3'→5' exo 活性の残存を除去)")
        sections.append("")

        # トラブルシューティング
        sections.append("━━━ 5. トラブルシューティング ━━━")
        sections.append("")
        sections.append("  【増幅しない / 効率が悪い】")
        sections.append("    → 伸長時間を 10〜60 sec/kb に延長")
        sections.append("    → サイクル数を 35〜40 に増加")
        sections.append("    → アニーリング時間を 15 sec に設定")
        sections.append("    → アニーリング温度を 2℃ずつ低下 (50〜53℃)")
        sections.append("    → 反応液量を 25 μL に変更")
        sections.append("")
        sections.append("  【スメア / エキストラバンド】")
        sections.append("    → アニーリング時間を 5 sec に短縮")
        sections.append("    → アニーリング温度を 58〜63℃ に上昇")
        sections.append("    → 2-step PCR (98℃/68℃) を試行")
        sections.append("    → 鋳型 DNA 量を減量")
        sections.append("    → サイクル数を 25〜30 に減少")
        sections.append("")

        # 注意事項
        if self.notes:
            sections.append("━━━ 6. この変異に対する注意事項 ━━━")
            sections.append("")
            for note in self.notes:
                sections.append(f"  • {note}")
            sections.append("")

        sections.append(f"{sep}")
        sections.append("  Generated by SDM Primer Designer Pro")
        sections.append(f"{sep}")
        return "\n".join(sections)


class SDMPrimerDesigner:
    # 標準的なプラスミドパーツのシグネチャー配列
    FEATURE_LIBRARY = {
        "AmpR (bla gene)": "TTTCCGTGTCGCCCTTATTCCCTTTTTTGC", 
        "pUC ori": "GGTGAGCGTGGGTCTCGCGGTATCATTGC", 
        "T7 promoter": "TAATACGACTCACTATAGGG",
        "CMV promoter": "TGACATTGATTATTGACTAGTTATTAATAG",
    }

    ECOLI_CODONS = {
        'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
        'Q': 'CAG', 'E': 'GAA', 'G': 'GGT', 'H': 'CAC', 'I': 'ATT',
        'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
        'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
        '*': 'TAA'
    }

    def __init__(self, template_sequence, host_codons=None):
        self.template_dna = str(template_sequence).upper().strip()
        self.host_codons = host_codons if host_codons else self.ECOLI_CODONS
        self.enzymes = Restriction.CommOnly

    def _calculate_gc_content(self, seq_str):
        g, c = seq_str.count('G'), seq_str.count('C')
        return round((g + c) / len(seq_str) * 100, 1)

    def _calculate_self_dimer_tm(self, seq_str):
        seq = str(seq_str).upper()
        rc = str(Seq(seq).reverse_complement())
        n, max_dimer_tm = len(seq), 0
        for shift in range(-n + 4, n - 3):
            s1, s2 = (seq[:n+shift], rc[-shift:]) if shift < 0 else (seq[shift:], rc[:n-shift])
            if len(s1) >= 4:
                try:
                    res_tm = mt.Tm_NN(Seq(s1), cseq=Seq(s2))
                    if res_tm > max_dimer_tm: max_dimer_tm = res_tm
                except: continue
        return round(max_dimer_tm, 2)

    def detect_features(self, sequence, custom_library=None):
        found = []
        seq_str = str(sequence).upper()
        lib = self.FEATURE_LIBRARY.copy()
        if custom_library: lib.update(custom_library)
        for name, sig in lib.items():
            start = seq_str.find(sig.upper())
            if start != -1:
                found.append({"name": name, "start": start, "end": start + len(sig), "strand": 1})
                continue
            rc_sig = str(Seq(sig).reverse_complement()).upper()
            start_rc = seq_str.find(rc_sig)
            if start_rc != -1:
                found.append({"name": name, "start": start_rc, "end": start_rc + len(rc_sig), "strand": -1})
        return found

    def design(self, row, method='overlapping', target_tm=68):
        name = row['mutation_name']
        dna_idx = (int(row['aa_pos']) - 1) * 3
        mode = row['mode']
        if mode == 'sub':
            mut_dna, ref_len = self.host_codons.get(row['target_aa'], "NNN"), 3
        elif mode == 'ins':
            mut_dna, ref_len = str(row['insert_seq']).upper(), 0
        elif mode == 'del':
            mut_dna, ref_len = "", int(row.get('del_len', 1)) * 3
        else: return None
        
        mod_seq = self.template_dna[:dna_idx] + mut_dna + self.template_dna[dna_idx + ref_len:]
        win = 20
        c_s, c_e = max(0, dna_idx - win), min(len(self.template_dna), dna_idx + len(mut_dna) + win)
        a_orig = Restriction.Analysis(self.enzymes, Seq(self.template_dna[c_s:c_e]))
        a_mod = Restriction.Analysis(self.enzymes, Seq(mod_seq[c_s:c_e]))
        new_sites = set([str(e) for e, p in a_mod.full().items() if p]) - set([str(e) for e, p in a_orig.full().items() if p])

        if method == 'overlapping':
            for length in range(15, 45):
                start, end = dna_idx - length, dna_idx + len(mut_dna) + length
                if start < 0 or end > len(mod_seq): break
                fwd = str(mod_seq[start:end])
                tm = mt.Tm_NN(Seq(fwd))
                if tm >= target_tm:
                    f_obj = Seq(fwd)
                    return {
                        "mutation_name": name, "fwd_primer": fwd, "rev_primer": str(f_obj.reverse_complement()),
                        "Tm": round(tm, 2), "GC_%": self._calculate_gc_content(fwd),
                        "MW": round(molecular_weight(f_obj, seq_type='DNA'), 1),
                        "Dimer_Tm": self._calculate_self_dimer_tm(fwd),
                        "Product_Size(bp)": len(mod_seq),
                        "New_Sites": ", ".join(new_sites) if new_sites else "None",
                        "full_seq": mod_seq, "mut_start": dna_idx, "mut_end": dna_idx + len(mut_dna)
                    }
        return None

    def generate_protocol(self, result: dict) -> PCRProtocol:
        """設計結果辞書から PCRProtocol オブジェクトを生成"""
        return PCRProtocol(
            mutation_name=result['mutation_name'],
            product_size_bp=result['Product_Size(bp)'],
            primer_tm=result['Tm'],
            primer_length=len(result['fwd_primer']),
            fwd_primer=result['fwd_primer'],
            rev_primer=result['rev_primer'],
        )