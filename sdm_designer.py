import pandas as pd
import logging
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Data import CodonTable

# --- ログ設定 ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.FileHandler("design.log"), logging.StreamHandler()]
)

class SDMPrimerDesigner:
    def __init__(self, template_sequence, host_codons=None):
        # ファイルパスではなく、文字列を直接セットするように修正
        self.template_dna = str(template_sequence).upper().strip()
        self.host_codons = host_codons if host_codons else self.ECOLI_CODONS
        self.enzymes = Restriction.CommOnly

    def _get_codon(self, aa):
        return self.host_codons.get(aa, "NNN")

    def analyze_sites(self, modified_seq):
        orig_analysis = Restriction.Analysis(self.enzymes, self.template_dna)
        mod_analysis = Restriction.Analysis(self.enzymes, modified_seq)
        new = set(mod_analysis.with_sites().keys()) - set(orig_analysis.with_sites().keys())
        lost = set(orig_analysis.with_sites().keys()) - set(mod_analysis.with_sites().keys())
        return ", ".join([str(e) for e in new]), ", ".join([str(e) for e in lost])

    def design_primers(self, row, method='overlapping', target_tm=78):
        name = row['mutation_name']
        dna_idx = (int(row['aa_pos']) - 1) * 3
        mode = row['mode']

        # 変異DNAの構築
        if mode == 'sub':
            mut_dna = self._get_codon(row['target_aa'])
            ref_len = 3
        elif mode == 'ins':
            mut_dna = str(row['insert_seq']).upper()
            ref_len = 0
        elif mode == 'del':
            mut_dna = ""
            ref_len = int(row.get('del_len', 1)) * 3
        else:
            logging.error(f"Unknown mode in {name}")
            return None

        # 変異後配列の生成とサイトチェック
        mod_seq = self.template_dna[:dna_idx] + mut_dna + self.template_dna[dna_idx + ref_len:]
        new_s, lost_s = self.analyze_sites(mod_seq)

        # プライマー設計ロジック (Overlapping)
        if method == 'overlapping':
            for length in range(15, 35):
                start, end = dna_idx - length, dna_idx + len(mut_dna) + length
                if start < 0 or end > len(mod_seq): break
                fwd = str(mod_seq[start:end])
                tm = mt.Tm_NN(Seq(fwd))
                if tm >= target_tm:
                    return {"name":name, "fwd":fwd, "rev":str(Seq(fwd).reverse_complement()), "tm":round(tm,2), "new_sites":new_s, "lost_sites":lost_s}
        
        # Back-to-back等はここに同様に実装
        return None