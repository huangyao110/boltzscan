import os
import argparse
from Bio import SeqIO

def read_a3m(a3m_file):
    return list(SeqIO.parse(a3m_file, 'fasta'))[:-1]

def combine_unpaired_a3m_for_colabfold(chain_a_path, chain_b_path, output_path):
    records_a = read_a3m(chain_a_path)
    records_b = read_a3m(chain_b_path)

    L1 = len(str(records_a[0].seq))
    L2 = len(str(records_b[0].seq))

    header_line = f"#{L1},{L2}\t1,1"
    combined_lines = [header_line]

    # Add combined reference line
    combined_lines.append(">101\t102")
    combined_lines.append(str(records_a[0].seq) + str(records_b[0].seq))

    # Chain A block
    for record in records_a:
        combined_lines.append(f">{record.id}")
        combined_lines.append(str(record.seq) + '-' * L2)

    # Chain B block
    for i, record in enumerate(records_b):
        header = ">102" if i == 0 else f">{record.id.replace(chr(9), ' ').split()[0]}"
        combined_lines.append(header)
        combined_lines.append('-' * L1 + str(record.seq))

    with open(output_path, "w") as f:
        f.write("\n".join(combined_lines) + "\n")

    print(f" Combined and written to: {output_path}")

def crop_a3m_file(a3m_file, out_path, interval=None):
    records = read_a3m(a3m_file)
    # crop_res = []
    with open(out_path, 'w') as f:
        for idx, rec in enumerate(records):
            if interval:
                crop_seq = str(rec.seq[interval[0]: interval[1]]).upper()
            else:
                crop_seq= str(rec.seq).upper()
            crop_info = rec.description
            if all([i=='-' for i in crop_seq]):
                continue
            else:
                if idx == 0:
                    f.write('>' + crop_info + '\n')
                    f.write(crop_seq)
                else:
                    f.write('\n' + '>' + crop_info + '\n')
                    f.write(crop_seq)
    return out_path