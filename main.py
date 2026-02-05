import pandas as pd
import logging
from Bio import SeqIO, Restriction
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Data import CodonTable

# --- 簡易的な大腸菌コドン頻度テーブル (Most frequent codons) ---
# 実用時は外部のJSONやBio.SeqUtils.CodonUsageから読み込むとより正確です
E_COLI_BEST_CODONS = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGT', 'H': 'CAC', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
    '*': 'TAA'
}

class SDMPrimerDesigner:
    def __init__(self, template_dna, host_codons=E_COLI_BEST_CODONS):
        self.template_dna = Seq(template_dna.upper().strip().replace(" ", ""))
        self.host_codons = host_codons
        # 一般的な制限酵素セット (EcoRI, BamHI, HindIII, XhoIなど)
        self.enzyme_set = Restriction.CommOnly

    def get_optimized_codon(self, amino_acid):
        return self.host_codons.get(amino_acid, None)

    def analyze_restriction(self, original_seq, modified_seq):
        """変異前後の制限酵素サイトの変化を解析"""
        orig_analysis = Restriction.Analysis(self.enzyme_set, original_seq)
        mod_analysis = Restriction.Analysis(self.enzyme_set, modified_seq)
        
        orig_sites = set(orig_analysis.with_sites().keys())
        mod_sites = set(mod_analysis.with_sites().keys())
        
        new_sites = mod_sites - orig_sites
        lost_sites = orig_sites - mod_sites
        
        return list(new_sites), list(lost_sites)

    def design(self, row, method='overlapping', target_tm=78):
        name = row['mutation_name']
        dna_idx = (int(row['aa_pos']) - 1) * 3
        
        # 1. 変異DNAの決定（コドン最適化適用）
        if row['mode'] == 'sub':
            mut_dna = self.get_optimized_codon(row['target_aa'])
            target_len = 3
        elif row['mode'] == 'ins':
            mut_dna = str(row['insert_seq']).upper()
            target_len = 0
        elif row['mode'] == 'del':
            mut_dna = ""
            target_len = int(row.get('del_len', 1)) * 3
        else:
            return None

        # 2. 変異後の全配列を作成（制限酵素チェック用）
        modified_full_seq = self.template_dna[:dna_idx] + mut_dna + self.template_dna[dna_idx + target_len:]
        new_sites, lost_sites = self.analyze_restriction(self.template_dna, modified_full_seq)

        # 3. プライマー設計（前回のロジックを流用）
        # ※ここではoverlappingを例に簡略化
        length = 15
        while length < 35:
            start = dna_idx - length
            end = dna_idx + len(mut_dna) + length
            if start < 0 or end > len(modified_full_seq): break
            
            fwd = str(modified_full_seq[start:end])
            tm = mt.Tm_NN(Seq(fwd))
            if tm >= target_tm:
                res = {
                    "name": name,
                    "fwd": fwd,
                    "rev": str(Seq(fwd).reverse_complement()),
                    "tm": round(tm, 2),
                    "new_sites": ", ".join([str(e) for e in new_sites]),
                    "lost_sites": ", ".join([str(e) for e in lost_sites])
                }
                logging.info(f"[{name}] 設計成功 (New: {res['new_sites']})")
                return res
            length += 1
        return None