import os
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import numpy as np

project_root = os.path.dirname(os.path.abspath(__file__))
fasta_directory = os.path.join(project_root, 'viruses', 'data')


def get_six_frames(seq):
    frames = []
    reverse_complement = seq.reverse_complement()
    for i in range(3):
        frames.append(seq[i:])
        frames.append(reverse_complement[i:])
    return frames


def find_orfs(reading_frame):
    start_codons = ['ATG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    seq_len = len(reading_frame)

    for i in range(0, seq_len - 2, 3):
        codon = reading_frame[i:i + 3]
        if codon in start_codons:
            for j in range(i + 3, seq_len - 2, 3):
                stop_codon = reading_frame[j:j + 3]
                if stop_codon in stop_codons:
                    orf = reading_frame[i:j + 3]
                    if len(orf) > 100:
                        orfs.append(orf)
                    break
    return orfs


def count_and_translate_codons_dicodons(orfs):
    codon_count = defaultdict(int)
    dicodon_count = defaultdict(int)

    total_codons = 0
    total_dicodons = 0

    for orf in orfs:
        for i in range(0, len(orf) - 2, 3):
            codon = orf[i:i + 3]
            translated_codon = str(Seq(codon).translate())
            if translated_codon != '*':
                codon_count[translated_codon] += 1
                total_codons += 1

            if i + 6 <= len(orf):
                dicodon = orf[i:i + 6]
                translated_dicodon = str(Seq(dicodon).translate())
                if '*' not in translated_dicodon:
                    dicodon_count[translated_dicodon] += 1
                    total_dicodons += 1

    codon_frequency = {codon: count / total_codons for codon, count in codon_count.items()} if total_codons > 0 else {}
    dicodon_frequency = {dicodon: count / total_dicodons for dicodon, count in
                         dicodon_count.items()} if total_dicodons > 0 else {}

    return codon_frequency, dicodon_frequency


def calculate_distance_matrix(codon_frequencies):
    all_files = list(codon_frequencies.keys())
    num_files = len(all_files)

    all_codons = set()
    for freqs in codon_frequencies.values():
        all_codons.update(freqs.keys())

    distance_matrix = np.zeros((num_files, num_files))

    for i in range(num_files):
        for j in range(num_files):
            if i != j:
                vec_i = np.array([codon_frequencies[all_files[i]].get(codon, 0) for codon in all_codons])
                vec_j = np.array([codon_frequencies[all_files[j]].get(codon, 0) for codon in all_codons])

                squared_diff = np.abs(vec_i ** 2 - vec_j ** 2)
                distance_matrix[i][j] = np.sum(squared_diff)

    return distance_matrix


def calculate_dicodon_distance_matrix(dicodon_frequencies):
    all_files = list(dicodon_frequencies.keys())
    num_files = len(all_files)

    all_dicodons = set()
    for freqs in dicodon_frequencies.values():
        all_dicodons.update(freqs.keys())

    distance_matrix = np.zeros((num_files, num_files))

    for i in range(num_files):
        for j in range(num_files):
            if i != j:
                vec_i = np.array([dicodon_frequencies[all_files[i]].get(dicodon, 0) for dicodon in all_dicodons])
                vec_j = np.array([dicodon_frequencies[all_files[j]].get(dicodon, 0) for dicodon in all_dicodons])
                squared_diff = np.abs(vec_i ** 2 - vec_j ** 2)
                distance_matrix[i][j] = np.sum(squared_diff)

    return distance_matrix


def write_phylip_format(filename, distance_matrix, file_names):
    with open(filename, 'w') as f:
        f.write(f"{len(file_names)}\n")
        for i in range(len(file_names)):
            distances = " ".join(f"{distance_matrix[i][j]:.3f}" for j in range(len(file_names)))
            name = file_names[i][:10].ljust(10)
            f.write(f"{name} {distances}\n")


def get_top_frequencies(frequency_dict, top_n=10):
    sorted_items = sorted(frequency_dict.items(), key=lambda item: item[1], reverse=True)
    return sorted_items[:top_n]


def process_fasta_files():
    codon_frequencies = defaultdict(lambda: defaultdict(int))
    dicodon_frequencies = defaultdict(lambda: defaultdict(int))
    mammalian_codon_frequencies = defaultdict(int)
    bacterial_codon_frequencies = defaultdict(int)
    mammalian_dicodon_frequencies = defaultdict(int)
    bacterial_dicodon_frequencies = defaultdict(int)

    for file_name in os.listdir(fasta_directory):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(fasta_directory, file_name)
            print(f"Processing {file_name}...")

            for record in SeqIO.parse(file_path, "fasta"):
                seq = record.seq
                six_frames = get_six_frames(seq)

                all_orfs = []
                for frame in six_frames:
                    orfs = find_orfs(frame)
                    all_orfs.extend(orfs)

                codon_frequency_total, dicodon_frequency_total = count_and_translate_codons_dicodons(all_orfs)

                if "mamalian" in file_name.lower():
                    for codon, freq in codon_frequency_total.items():
                        mammalian_codon_frequencies[codon] += freq
                    for dicodon, freq in dicodon_frequency_total.items():
                        mammalian_dicodon_frequencies[dicodon] += freq
                elif "bacterial" in file_name.lower():
                    for codon, freq in codon_frequency_total.items():
                        bacterial_codon_frequencies[codon] += freq
                    for dicodon, freq in dicodon_frequency_total.items():
                        bacterial_dicodon_frequencies[dicodon] += freq

                for codon, freq in codon_frequency_total.items():
                    codon_frequencies[file_name][codon] += freq
                for dicodon, freq in dicodon_frequency_total.items():
                    dicodon_frequencies[file_name][dicodon] += freq

    codon_distance_matrix = calculate_distance_matrix(codon_frequencies)
    dicodon_distance_matrix = calculate_dicodon_distance_matrix(dicodon_frequencies)

    write_phylip_format("codon_distance_matrix.phy", codon_distance_matrix, list(codon_frequencies.keys()))
    write_phylip_format("dicodon_distance_matrix.phy", dicodon_distance_matrix, list(dicodon_frequencies.keys()))

    print(f"\nTop 10 codons in mammalian viruses:")
    for codon, count in get_top_frequencies(mammalian_codon_frequencies):
        print(f"{codon}: {count:.4f}")

    print(f"\nTop 10 codons in bacterial viruses:")
    for codon, count in get_top_frequencies(bacterial_codon_frequencies):
        print(f"{codon}: {count:.4f}")

    print(f"\nTop 10 dicodons in mammalian viruses:")
    for dicodon, count in get_top_frequencies(mammalian_dicodon_frequencies):
        print(f"{dicodon}: {count:.4f}")

    print(f"\nTop 10 dicodons in bacterial viruses:")
    for dicodon, count in get_top_frequencies(bacterial_dicodon_frequencies):
        print(f"{dicodon}: {count:.4f}")





if __name__ == "__main__":
    process_fasta_files()
