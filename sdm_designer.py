from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import molecular_weight
from Bio import Restriction

class SDMPrimerDesigner:
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
            if start != -1: found.append({"name": name, "start": start, "end": start + len(sig), "strand": 1})
            else:
                rc_sig = str(Seq(sig).reverse_complement()).upper()
                start_rc = seq_str.find(rc_sig)
                if start_rc != -1: found.append({"name": name, "start": start_rc, "end": start_rc + len(rc_sig), "strand": -1})
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
        a_orig, a_mod = Restriction.Analysis(self.enzymes, Seq(self.template_dna[c_s:c_e])), Restriction.Analysis(self.enzymes, Seq(mod_seq[c_s:c_e]))
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
                        "Product_Size(bp)": len(mod_seq), # 追加
                        "New_Sites": ", ".join(new_sites) if new_sites else "None",
                        "full_seq": mod_seq, "mut_start": dna_idx, "mut_end": dna_idx + len(mut_dna)
                    }
        return None