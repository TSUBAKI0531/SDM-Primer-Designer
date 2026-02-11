from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Restriction

class SDMPrimerDesigner:
    # 標準パーツライブラリ
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

    def _calculate_self_dimer_tm(self, seq_str):
        """プライマーの自己結合における最大 Tm値を算出する"""
        seq = str(seq_str).upper()
        rc = str(Seq(seq).reverse_complement())
        n = len(seq)
        max_dimer_tm = 0
        
        # 配列をスライドさせて最も安定な結合(Tm)を探す
        for shift in range(-n + 4, n - 3):
            if shift < 0:
                s1, s2 = seq[:n+shift], rc[-shift:]
            else:
                s1, s2 = seq[shift:], rc[:n-shift]
            
            if len(s1) >= 4:
                try:
                    # 二本鎖としての安定性を計算
                    res_tm = mt.Tm_NN(Seq(s1), cseq=Seq(s2))
                    if res_tm > max_dimer_tm:
                        max_dimer_tm = res_tm
                except: continue
        return round(max_dimer_tm, 2)

    def detect_features(self, sequence, custom_library=None):
        found_features = []
        seq_str = str(sequence).upper()
        search_library = self.FEATURE_LIBRARY.copy()
        if custom_library: search_library.update(custom_library)
        
        for name, sig_seq in search_library.items():
            start = seq_str.find(sig_seq.upper())
            if start != -1:
                found_features.append({"name": name, "start": start, "end": start + len(sig_seq), "strand": 1})
                continue
            rc_sig = str(Seq(sig_seq).reverse_complement()).upper()
            start_rc = seq_str.find(rc_sig)
            if start_rc != -1:
                found_features.append({"name": name, "start": start_rc, "end": start_rc + len(rc_sig), "strand": -1})
        return found_features

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
        window = 20
        check_start, check_end = max(0, dna_idx - window), min(len(self.template_dna), dna_idx + len(mut_dna) + window)
        a_orig = Restriction.Analysis(self.enzymes, Seq(self.template_dna[check_start:check_end]))
        a_mod = Restriction.Analysis(self.enzymes, Seq(mod_seq[check_start:check_end]))
        new_sites = set([str(e) for e, p in a_mod.full().items() if p]) - set([str(e) for e, p in a_orig.full().items() if p])

        if method == 'overlapping':
            for length in range(15, 45):
                start, end = dna_idx - length, dna_idx + len(mut_dna) + length
                if start < 0 or end > len(mod_seq): break
                fwd = str(mod_seq[start:end])
                tm = mt.Tm_NN(Seq(fwd))
                if tm >= target_tm:
                    return {
                        "mutation_name": name, "fwd_primer": fwd, 
                        "rev_primer": str(Seq(fwd).reverse_complement()),
                        "Tm": round(tm, 2), 
                        "Dimer_Tm": self._calculate_self_dimer_tm(fwd), # ここで計算
                        "New_Sites": ", ".join(new_sites) if new_sites else "None",
                        "full_seq": mod_seq, "mut_start": dna_idx, "mut_end": dna_idx + len(mut_dna)
                    }
        return None