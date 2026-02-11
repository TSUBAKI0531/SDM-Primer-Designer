from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Restriction

class SDMPrimerDesigner:
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

    def _get_restriction_diff(self, original_segment, modified_segment):
        analysis_orig = Restriction.Analysis(self.enzymes, Seq(original_segment))
        analysis_mod = Restriction.Analysis(self.enzymes, Seq(modified_segment))
        orig_sites = set([str(e) for e, pos in analysis_orig.full().items() if pos])
        mod_sites = set([str(e) for e, pos in analysis_mod.full().items() if pos])
        return list(mod_sites - orig_sites), list(orig_sites - mod_sites)

    def design(self, row, method='overlapping', target_tm=68):
        name = row['mutation_name']
        dna_idx = (int(row['aa_pos']) - 1) * 3
        mode = row['mode']

        if mode == 'sub':
            mut_dna = self.host_codons.get(row['target_aa'], "NNN")
            ref_len = 3
        elif mode == 'ins':
            mut_dna = str(row['insert_seq']).upper()
            ref_len = 0
        elif mode == 'del':
            mut_dna = ""
            ref_len = int(row.get('del_len', 1)) * 3
        else: return None

        mod_seq = self.template_dna[:dna_idx] + mut_dna + self.template_dna[dna_idx + ref_len:]

        window = 20
        check_start = max(0, dna_idx - window)
        check_end = min(len(self.template_dna), dna_idx + len(mut_dna) + window)
        new_sites, lost_sites = self._get_restriction_diff(
            self.template_dna[check_start:check_end], mod_seq[check_start:check_end]
        )

        if method == 'overlapping':
            for length in range(15, 45):
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
                        "New_Sites": ", ".join(new_sites) if new_sites else "None",
                        "full_seq": mod_seq,
                        "mut_start": dna_idx,
                        "mut_end": dna_idx + len(mut_dna)
                    }
        return None