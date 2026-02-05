from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Restriction

class SDMPrimerDesigner:
    # 1. クラス変数としてコドン表を定義（これがないとエラーになります）
    ECOLI_CODONS = {
        'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
        'Q': 'CAG', 'E': 'GAA', 'G': 'GGT', 'H': 'CAC', 'I': 'ATT',
        'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
        'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
        '*': 'TAA'
    }

    def __init__(self, template_sequence, host_codons=None):
        # 文字列として配列を受け取り、大文字に統一
        self.template_dna = str(template_sequence).upper().strip()
        self.host_codons = host_codons if host_codons else self.ECOLI_CODONS
        self.enzymes = Restriction.CommOnly

    def _get_codon(self, aa):
        return self.host_codons.get(aa, "NNN")

    def design(self, row, method='overlapping', target_tm=78):
        """
        メインの設計ロジック
        """
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
            return None

        # 変異後配列の生成（制限酵素サイト解析用など）
        mod_seq = self.template_dna[:dna_idx] + mut_dna + self.template_dna[dna_idx + ref_len:]

        # プライマー設計（Overlapping形式）
        if method == 'overlapping':
            for length in range(15, 35):
                start = dna_idx - length
                end = dna_idx + len(mut_dna) + length
                if start < 0 or end > len(mod_seq): break
                
                fwd = str(mod_seq[start:end])
                tm = mt.Tm_NN(Seq(fwd))
                if tm >= target_tm:
                    return {
                        "mutation_name": name,
                        "fwd_primer": fwd,
                        "rev_primer": str(Seq(fwd).reverse_complement()),
                        "Tm": round(tm, 2),
                        "length": len(fwd)
                    }
        return None