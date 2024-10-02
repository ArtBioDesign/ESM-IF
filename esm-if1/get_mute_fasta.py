import pandas as pd
import os
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Data import CodonTable

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

def pdb_to_chain_fasta(pdb_file, chain_id):
    """
    从PDB文件中提取指定链的氨基酸序列.

    参数:
        pdb_file (str): PDB文件路径
        chain_id (str): 链ID

    返回:
        str: 氨基酸序列
    """
    parser = PDBParser(QUIET=True)
    sequence = ""

    try:
        structure = parser.get_structure('structure', pdb_file)
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.get_resname() != 'HOH':  # 忽略水分子
                            res_name = residue.get_resname()
                            try:
                                # 将三字母代码转换为一字母代码
                                amino_acid = seq1(res_name)
                                sequence += amino_acid
                            except ValueError:
                                print(f"警告: 残基 {res_name} 无法转换为氨基酸缩写")

                    if not sequence:
                        raise ValueError(f"未找到链 {chain_id} 的氨基酸序列")

                    return sequence

    except Exception as e:
        print(f"读取PDB文件时出错: {e}")
        return None

def apply_single_mutation(sequence, mutation):
    """
    根据单个突变修改序列.

    参数:
        sequence (str): 原始氨基酸序列
        mutation (str): 单个突变, 如 'A10G'

    返回:
        Seq: 突变后的序列，作为Biopython的Seq对象
    """
    sequence = list(sequence)
    original, pos, mutated = mutation[0], int(mutation[1:-1]) - 1, mutation[-1]
    if sequence[pos] == original:
        sequence[pos] = mutated
    else:
        print(f"突变错误: 位置 {pos+1} 的氨基酸应为 {original}，实际为 {sequence[pos]}")
    return Seq(''.join(sequence))

def mutate_sequence_from_csv(sequence, csv_file, output_fasta_file):
    """
    根据CSV中的突变信息分别应用到序列，并将每条突变后的序列写入FASTA文件。

    参数:
        sequence (str): 原始氨基酸序列
        csv_file (str): 包含突变信息的CSV文件路径
        output_fasta_file (str): 输出FASTA文件路径
    """
    try:
        df = pd.read_csv(csv_file)
        mutations = df['mutant'].tolist()
        seq_records = []

        # 对每个突变独立应用，不相互影响
        for i, mutation in enumerate(mutations):
            mutated_sequence = apply_single_mutation(sequence, mutation)
            seq_record = SeqRecord(mutated_sequence, id=f"mut_{i+1}", description=f"Mutation {mutation}")
            seq_records.append(seq_record)
            print(f"Mutation {mutation} 突变后的序列: {mutated_sequence}")

        # 将所有突变后的序列写入一个FASTA文件
        os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)
        with open(output_fasta_file, "w") as fasta_out:
            SeqIO.write(seq_records, fasta_out, "fasta")

        print(f"所有突变后的序列已写入: {output_fasta_file}")

    except Exception as e:
        print(f"处理突变时出错: {e}")


def main(pdb_file, chain_id, csv_file, output_fasta_file):

    chain_sequence = pdb_to_chain_fasta(pdb_file, chain_id)
    print(chain_sequence)
    mutate_sequence_from_csv(chain_sequence, csv_file, output_fasta_file)


if __name__ == '__main__':
    pdb_file = "data/5YH2.pdb"
    chain_id = "C"
    mute_csv_file = "data/input.csv"
    output_fasta_file = "data/mutated_sequence.fasta"
    main( pdb_file, chain_id, mute_csv_file, output_fasta_file )



# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="从PDB文件中提取链序列，应用突变并输出FASTA文件")
#     parser.add_argument("--pdb_file", help="PDB文件路径")
#     parser.add_argument("--chain_id", help="要提取的链ID")
#     parser.add_argument("--csv_file", help="包含突变信息的CSV文件路径")
#     parser.add_argument("--output_fasta_file", help="输出的FASTA文件路径")

#     args = parser.parse_args()
#     main(args.pdb_file, args.chain_id, args.csv_file, args.output_fasta_file)


    # 示例使用
    # pdb_file = "/hpcfs/fhome/yangchh/ai/esm-if1/data/5YH2.pdb"
    # chain_id = "C"
    # mute_csv_file = "/hpcfs/fhome/yangchh/ai/esm-if1/input.csv"
    # output_fasta_file = "/hpcfs/fhome/yangchh/ai/esm-if1/data/mutated_sequence.fasta"

    # chain_sequence = pdb_to_chain_fasta(pdb_file, chain_id)
    # if chain_sequence:
    #     print(f"链 {chain_id} 序列: {chain_sequence}")
    #     mutate_sequence_from_csv(chain_sequence, mute_csv_file, output_fasta_file)

