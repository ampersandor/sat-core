import os
from Bio import SeqIO
from collections import Counter
from dataclasses import dataclass

Nucleotide_common = ["A", "T", "G", "C"]
Nucleotide_IUPAC_code = {
    "A": "A",
    "G": "G",
    "C": "C",
    "T": "T",
    "CT": "Y",
    "AG": "R",
    "AT": "W",
    "GT": "K",
    "CG": "S",
    "AC": "M",
    "AGT": "D",
    "ACG": "V",
    "ACT": "H",
    "CGT": "B",
    "": "None",
    "ACGT": "X",
}
DEG_list = ["Y", "R", "W", "K", "S", "M", "N", "D", "V", "H", "B", "-", "."]


@dataclass
class Statistic:
    total_seq: int
    gap_seq_count: int
    gap_count: int
    gap_frequency: float
    gap_sum_length: int
    gap_length: int
    sum_of_blue_bases: int
    no_blue_bases: int
    no_miss_bases: int
    blue_base_ratio: float
    blue_base_count: dict[float, int]

class BlueBase:
    def __init__(self, input_file, output_dir):
        self.input_file = input_file
        self.output_dir = output_dir
        self.base_name = os.path.basename(self.input_file)


    def get_statistics(self, stat_file):
        (
            stat_values,
            gap_statheader,
            gap_stat_result,
            pct_id_cutoff,
            blue_base_count,
        ) = self.align_to_statistics()

        with open(stat_file, "w") as wf:
            wf.write("\n".join(stat_values) + "\n")
        
        blue_base_count = {
            pct_id_cutoff[i]: blue_base_count[pct_id_cutoff[i]] for i in range(len(pct_id_cutoff))
        }
        statistic = Statistic(
            total_seq=gap_stat_result[0],
            gap_seq_count=gap_stat_result[1],
            gap_count=gap_stat_result[2],
            gap_frequency=gap_stat_result[3],
            gap_sum_length=gap_stat_result[4],
            gap_length=gap_stat_result[5],
            sum_of_blue_bases=gap_stat_result[6],
            no_blue_bases=gap_stat_result[7],
            no_miss_bases=gap_stat_result[8],
            blue_base_ratio=gap_stat_result[9],
            blue_base_count=blue_base_count,
        )
        return statistic

    def align_to_statistics(self):
        bp_color = dict()
        fasta_list = list()
        miss_front_list = list()
        miss_end_list = list()
        gap_seq_count = 0
        gap_count = 0
        gap_sumlength = 0

        for record in SeqIO.parse(self.input_file, "fasta"):
            fasta_list.append(str(record.seq).replace(".", "-").upper())  # append fasta sequence
            seq = str(record.seq).replace(".", "-").upper()  # sequence for miss base count

            miss_front = 0
            miss_front_seq = ""
            miss_continue = False
            base_start = 0

            for i in range(len(seq)):  # front -> end read
                if seq[i] != "-":
                    if base_start == 0:
                        base_start = i
                if miss_continue:
                    miss_front_seq += "|"
                elif i == 0 and seq[i] == "-":
                    miss_front += 1
                    miss_front_seq += "*"
                elif i == 0 and seq[i] != "-":
                    miss_front_seq += "|"
                    miss_continue = True
                elif i != 0 and seq[i] != "-":
                    miss_continue = True
                    miss_front_seq += "|"
                elif miss_front != 0 and seq[i] == "-":
                    miss_front += 1
                    miss_front_seq += "*"
                    
            miss_front_list.append(miss_front_seq)

            miss_end = 0
            miss_end_seq = ""
            miss_continue = False
            base_end = 0

            for i in reversed(range(len(seq))):  # end -> front read
                if seq[i] != "-":
                    if base_end == 0:
                        base_end = i
                if miss_continue:
                    miss_end_seq = "|" + miss_end_seq
                elif i == len(seq) - 1 and seq[i] == "-":
                    miss_end += 1
                    miss_end_seq = "*" + miss_end_seq
                elif i == len(seq) - 1 and seq[i] != "-":
                    miss_end_seq = "|" + miss_end_seq
                    miss_continue = True
                elif i != len(seq) - 1 and seq[i] != "-":
                    miss_continue = True
                    miss_end_seq = "|" + miss_end_seq
                elif miss_end != 0 and seq[i] == "-":
                    miss_end += 1
                    miss_end_seq = "*" + miss_end_seq
            miss_end_list.append(miss_end_seq)

            nomissgap_seq = seq[base_start : base_end + 1]
            if "-" in nomissgap_seq:
                gap_seq_count += 1

                flag = 0
                gap_starts = list()
                gap_ends = list()
                for n in range(len(nomissgap_seq)):
                    if nomissgap_seq[n] == "-":
                        if flag == 0:
                            gap_starts.append(n)
                            flag = 1
                    else:
                        if flag == 1:
                            gap_ends.append(n)
                            flag = 0
                gap_count += len(gap_starts)

                for g in range(len(gap_starts)):
                    gap_start = gap_starts[g]
                    gap_end = gap_ends[g]
                    gap_length = gap_end - gap_start
                    gap_sumlength += gap_length

        # gap statistics
        if gap_seq_count == 0:
            gap_freq = 0
            gap_avg_length = 0
        else:
            gap_freq = gap_count / float(gap_seq_count)
            gap_avg_length = gap_sumlength / float(gap_seq_count)

        list_length = len(fasta_list)  # Extract total fasta count
        sequence_length = len(fasta_list[0])  # Extract sequence length
        a_data = []
        t_data = []
        g_data = []
        c_data = []
        etc_data = []
        miss_data = []
        real_gap_data = []
        miss_gap_data = []
        total_data = []
        header_data = ["Position"]
        apcr_data = []
        a_freq_data = []
        t_freq_data = []
        g_freq_data = []
        c_freq_data = []
        iupac_data = []
        max_seq_data = []
        max_seq_count_data = []

        cnt = 1
        for i in range(0, sequence_length):
            data = []
            mdata = []  # base miss data
            not_gap_data = []
            header_data.append(str(cnt))
            for j in range(0, list_length):
                try:
                    if fasta_list[j][i] in Nucleotide_common:
                        not_gap_data.append(fasta_list[j][i])
                    if miss_front_list[j][i] == "*" and miss_end_list[j][i] == "|":
                        mdata.append("*")
                    elif miss_front_list[j][i] == "|" and miss_end_list[j][i] == "|":
                        mdata.append("|")
                    elif miss_front_list[j][i] == "|" and miss_end_list[j][i] == "*":
                        mdata.append("*")
                    data.append(fasta_list[j][i])
                except Exception as e:
                    data.append("-")

            nuc_counter = Counter(data)


            a_count = nuc_counter["A"]
            t_count = nuc_counter["T"]
            g_count = nuc_counter["G"]
            c_count = nuc_counter["C"]
            y_count = nuc_counter["Y"]
            r_count = nuc_counter["R"]
            w_count = nuc_counter["W"]
            s_count = nuc_counter["S"]
            k_count = nuc_counter["K"]
            m_count = nuc_counter["M"]
            d_count = nuc_counter["D"]
            v_count = nuc_counter["V"]
            h_count = nuc_counter["H"]
            b_count = nuc_counter["B"]
            ambi_count = nuc_counter["N"]
            miss_count = mdata.count("*")

            if len(not_gap_data) != 0:
                max_nucleotide, num_max_nucleotide = Counter(not_gap_data).most_common(1)[0]
            else:
                max_nucleotide = "-"
                num_max_nucleotide = 0

            etc_count = (
                y_count
                + r_count
                + w_count
                + s_count
                + k_count
                + m_count
                + d_count
                + v_count
                + h_count
                + b_count
                + ambi_count
            )
            miss_gap_count = data.count("-") + data.count(".")
            real_gap_count = miss_gap_count - miss_count
            sequence_count = a_count + t_count + g_count + c_count  # + etc_count
            total_count = a_count + t_count + g_count + c_count + etc_count + miss_gap_count
            apcr_count = (sequence_count / (total_count * 1.0)) * 100

            delete_duplicate = list(set(data))
            delete_duplicate.sort()
            temp = "".join(delete_duplicate)

            for deg_item in DEG_list:
                temp = temp.replace(deg_item, "")

            a_freq = round((float(a_count) / total_count) * 100)
            t_freq = round((float(t_count) / total_count) * 100)
            g_freq = round((float(g_count) / total_count) * 100)
            c_freq = round((float(c_count) / total_count) * 100)

            iupac_bases = Nucleotide_IUPAC_code[temp]
            a_data.append(str(a_count))
            t_data.append(str(t_count))
            g_data.append(str(g_count))
            c_data.append(str(c_count))
            etc_data.append(str(etc_count))
            miss_data.append(str(miss_count))
            real_gap_data.append(str(real_gap_count))
            miss_gap_data.append(str(miss_gap_count))
            total_data.append(str(total_count))
            apcr_data.append(str(round(apcr_count)))
            a_freq_data.append(str(a_freq))
            t_freq_data.append(str(t_freq))
            g_freq_data.append(str(g_freq))
            c_freq_data.append(str(c_freq))
            iupac_data.append(iupac_bases)
            max_seq_data.append(max_nucleotide)
            max_seq_count_data.append(str(num_max_nucleotide))

            freq_max = max([a_freq, t_freq, g_freq, c_freq])
            bp_color[cnt] = [max_nucleotide, freq_max]
            cnt = cnt + 1

        # Blue base max.
        blue_base_count = dict()
        nomiss_base_count = 0
        noblue_base_count = 0
        pctid_cutoff = [90.0, 80.0, 70.0, 60.0, 50.0, 40.0]
        for seq in fasta_list:
            for j in range(len(seq)):
                if seq[j] != "-":
                    nomiss_base_count += 1
                if seq[j] == bp_color[j + 1][0]:
                    for c in range(len(pctid_cutoff)):
                        pctid = pctid_cutoff[c]
                        if c == 0:
                            if pctid + 10 >= bp_color[j + 1][1] >= pctid:
                                if pctid not in blue_base_count:
                                    blue_base_count[pctid] = 0
                                blue_base_count[pctid] += 1
                        else:
                            if pctid + 10 > bp_color[j + 1][1] >= pctid:
                                if pctid not in blue_base_count:
                                    blue_base_count[pctid] = 0
                                blue_base_count[pctid] += 1
                    if bp_color[j + 1][1] < 40.0:
                        pctid = 40.0
                        if pctid not in blue_base_count:
                            blue_base_count[pctid] = 0
                        blue_base_count[pctid] += 1
                else:
                    if seq[j] != "-":
                        noblue_base_count += 1

        blue_base_ratio = sum(blue_base_count.values()) / float(nomiss_base_count)

        gap_stat_header = [
            "Total seqs",
            "Gap seq. count",
            "Gap count",
            "Gap frequency",
            "Sum of gap length",
            "Gap length",
            "Sum of blue bases",
            "No blue bases",
            "No miss bases",
            "Blue base ratio",
        ]
        gap_stat_result = [
            len(fasta_list),
            gap_seq_count,
            gap_count,
            gap_freq,
            gap_sumlength,
            gap_avg_length,
            sum(blue_base_count.values()),
            noblue_base_count,
            nomiss_base_count,
            blue_base_ratio,
        ]

        # Stat
        header = "\t".join(header_data)

        value_names = (
            list(
                map(
                    lambda x: "{} count".format(x),
                    ["A", "T", "G", "C", "Miss", "Gap", "etc", "MissGap"],
                )
            )
            + list(map(lambda x: "{} freq".format(x), ["A", "T", "G", "C"]))
            + ["Total count", "Coverage", "IUPAC", "Major base", "Major base count"]
        )
        value_data_list = [
            a_data,
            t_data,
            g_data,
            c_data,
            miss_data,
            real_gap_data,
            etc_data,
            miss_gap_data,
            a_freq_data,
            t_freq_data,
            g_freq_data,
            c_freq_data,
            total_data,
            apcr_data,
            iupac_data,
            max_seq_data,
            max_seq_count_data,
        ]

        values = [header]
        for value_index in range(len(value_names)):
            value = "{}\t{}".format(value_names[value_index], "\t".join(value_data_list[value_index]))
            values.append(value)

        return values, gap_stat_header, gap_stat_result, pctid_cutoff, blue_base_count

    def main(self):
        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"File not found: {self.input_file}")
        stat_file = os.path.join(self.output_dir, self.base_name + ".txt")
        statistic = self.get_statistics(stat_file)

        return stat_file, statistic


cur_dir = os.path.dirname(os.path.abspath(__file__))
BlueBase(os.path.join(cur_dir, "aligned.fasta"), cur_dir).main()