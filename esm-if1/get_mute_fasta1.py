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

def apply_mutations(sequence, mutations):
    """根据突变列表对序列进行修改"""
    sequence = list(sequence)
    for mutation in mutations:
        original, pos, mutated = mutation[0], int(mutation[1:-1]) - 1, mutation[-1]
        if sequence[pos] == original:
            sequence[pos] = mutated
        else:
            print(f"突变错误: 位置 {pos+1} 的氨基酸应为 {original}，实际为 {sequence[pos]}")
    return Seq(''.join(sequence))

def mutate_sequence_from_csv(sequence, csv_file, output_fasta_file):
    """根据CSV中的突变信息应用到序列并写入FASTA文件"""
    try:
        df = pd.read_csv(csv_file)
        mutations = df['mutant'].tolist()
        mutated_sequence = apply_mutations(sequence, mutations)
        os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)
        seq_record = SeqRecord(mutated_sequence, id="mutated_sequence", description="Mutated sequence")
        with open(output_fasta_file, "w") as fasta_out:
            SeqIO.write(seq_record, fasta_out, "fasta")
        print(f"突变后的序列已写入: {output_fasta_file}")
        return mutated_sequence
    except Exception as e:
        print(f"处理突变时出错: {e}")

def main(pdb_file, chain_id, csv_file, output_fasta_file):
    """主函数，执行PDB到FASTA转换及突变应用"""
    chain_sequence = pdb_to_chain_fasta(pdb_file, chain_id)
    if chain_sequence:
        print(f"链 {chain_id} 序列: {chain_sequence}")
        mutate_sequence_from_csv(chain_sequence, csv_file, output_fasta_file)


     


    # 示例使用
    # pdb_file = "/hpcfs/fhome/yangchh/ai/esm-if1/data/5YH2.pdb"
    # chain_id = "C"
    # mute_csv_file = "/hpcfs/fhome/yangchh/ai/esm-if1/input.csv"
    # output_fasta_file = "/hpcfs/fhome/yangchh/ai/esm-if1/data/mutated_sequence.fasta"

    # chain_sequence = pdb_to_chain_fasta(pdb_file, chain_id)
    # if chain_sequence:
    #     print(f"链 {chain_id} 序列: {chain_sequence}")
    #     mutate_sequence_from_csv(chain_sequence, mute_csv_file, output_fasta_file)


