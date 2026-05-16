from Bio import SeqIO
from Bio.Seq import Seq
from typing import Dict, List, Generator, Any
import sys
import re
import gzip
import argparse


def parse_gff3(file_path: str) -> Generator[List[Dict[str, Any]], None, None]:
    """
    解析GFF3文件，返回包含基因组特征的生成器
    每个特征以字典形式返回，包含所有标准字段和解析后的属性字段
    按基因分组处理，返回每个基因的所有特征 (assuming groups are separated by ### or implicitly by file end)
    """
    current_gene = []
    open_func = gzip.open if file_path.endswith('.gz') else open
    read_mode = 'rt' if file_path.endswith('.gz') else 'r'
    print(f"Parsing GFF using generator: {file_path}")

    try:
        with open_func(file_path, read_mode, encoding='utf-8') as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue

                # Check for standard GFF separator, handle end of file later
                if line == '###':
                    if current_gene:
                        yield current_gene
                        current_gene = []
                    continue

                if line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) != 9:
                    # print(f"Skipping malformed line: {line}", file=sys.stderr)
                    continue

                seqid, source, type_, start_str, end_str, score_str, strand, phase_str, attributes_str = fields

                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    # print(f"Skipping line with invalid coordinates: {line}", file=sys.stderr)
                    continue

                feature = {
                    'seqid': seqid,
                    'source': source,
                    'type': type_,
                    'start': start,
                    'end': end,
                    'score': None if score_str == '.' else float(score_str),
                    'strand': strand,
                    'phase': None if phase_str == '.' else int(phase_str),
                    'attributes': {}
                }

                for pair in attributes_str.split(';'):
                    if '=' in pair:
                        try:
                            key, value = pair.strip().split('=', 1)
                            feature['attributes'][key] = value
                        except ValueError:
                            # Handle cases like flag attributes without '='
                            # print(f"Skipping attribute without '=': {pair} in line: {line}", file=sys.stderr)
                            pass


                # --- Grouping Logic ---
                # Decide how to group. Simplest: group by features sharing the same *first* field (seqid)?
                # Or rely *only* on '###' separator? Let's rely on '###' for now.
                current_gene.append(feature)

            # Yield the last gene group if file doesn't end with ###
            if current_gene:
                yield current_gene
        print("Finished parsing GFF.")
    except FileNotFoundError:
        print(f"Error: GFF file not found: {file_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading GFF file {file_path}: {type(e).__name__} - {e}", file=sys.stderr)
        sys.exit(1)


def find_tss(gene_lst):
    # Extract gene basic information
    gene_features = [feature for feature in gene_lst if feature['type'] == 'gene']
    if not gene_features:
        return None

    first_gene = gene_features[0]
    chr_num = first_gene['seqid']
    strand = first_gene['strand']
    gene_name = first_gene['attributes'].get('ID', None)

    cds_starts = []
    cds_ends = []

    # Collect all CDS start and end positions
    for feature in gene_lst:
        if feature['type'] == 'CDS':
            cds_starts.append(feature['start'])
            cds_ends.append(feature['end'])

    # Determine TSS based on strand direction
    if strand == '+':
        tss = min(cds_starts) if cds_starts else None
    elif strand == '-':
        tss = max(cds_ends) if cds_ends else None
    else:
        tss = None  # Handle non-standard strand values

    if tss is not None and gene_name is not None:
        return [gene_name.split(':')[-1], chr_num, tss, strand]
    return None


def extract_bed_regions(gff_file, output_file, feature_type="gene", name_attribute="ID"):
    """
    从GFF文件中提取指定类型的特征并生成BED格式文件
    
    参数:
        gff_file: GFF文件路径
        output_file: 输出BED文件路径
        feature_type: 要提取的特征类型，默认为"gene"
        name_attribute: 用作BED文件第4列名称的属性，默认为"ID"
    """
    print(f"从GFF文件 {gff_file} 中提取 {feature_type} 特征到BED文件 {output_file}")
    
    try:
        with open(output_file, 'w') as bed_file:
            for gene_features in parse_gff3(gff_file):
                for feature in gene_features:
                    if feature['type'] == feature_type:
                        # BED格式: chrom, chromStart, chromEnd, name, score, strand
                        chrom = feature['seqid']
                        start = feature['start'] - 1  # BED使用0-based坐标系统
                        end = feature['end']
                        name = feature['attributes'].get(name_attribute, f"{feature_type}_{start}_{end}")
                        score = feature['score'] if feature['score'] is not None else 0
                        strand = feature['strand']
                        
                        # 写入BED格式行
                        bed_file.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        
        print(f"成功生成BED文件: {output_file}")
    except Exception as e:
        print(f"生成BED文件时出错: {type(e).__name__} - {e}", file=sys.stderr)
        sys.exit(1)


def extract_promoters_both(args):
    """
    从GFF文件中提取启动子区域并保存为BED和FASTA两种格式
    
    参数:
        args: 包含必要参数的命名空间对象
            - gff: GFF文件路径
            - genome: 基因组FASTA文件路径
            - output: 输出文件前缀，将生成.bed和.fasta两个文件
            - upstream: 上游区域长度
            - downstream: 下游区域长度
    """
    # --- Load Genome into Memory ---
    print(f"\n加载基因组文件 '{args.genome}' 到内存...")
    try:
        genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
        print(f"从基因组文件加载了 {len(genome_dict)} 个序列。")
        if not genome_dict:
             print("错误: 从基因组文件未加载到序列。", file=sys.stderr)
             sys.exit(1)
    except FileNotFoundError:
        print(f"错误: 未找到基因组文件 '{args.genome}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"读取基因组文件 '{args.genome}' 时出错: {type(e).__name__} - {e}", file=sys.stderr)
        sys.exit(1)

    # --- 准备输出文件 ---
    bed_output = args.output + ".bed"
    fasta_output = args.output + ".fasta"
    
    print(f"\n处理GFF并提取启动子区域到: {bed_output} 和 {fasta_output}...")
    success_count = 0
    error_count = 0
    gff_gene_groups_processed = 0

    try:
        with open(bed_output, 'w') as bed_file, open(fasta_output, 'w') as fasta_file:
            # 遍历由解析器生成的基因特征列表
            for gene_feature_list in parse_gff3(args.gff):
                gff_gene_groups_processed += 1
                if not gene_feature_list: continue # 跳过空组

                # 使用函数找到锚点(TSS代理)
                tss_info = find_tss(gene_feature_list)

                if tss_info is None:
                    error_count += 1
                    continue # 如果无法确定TSS锚点则跳过

                gene_id, chrom, tss_anchor, strand = tss_info

                # --- 核心提取逻辑 ---
                if chrom not in genome_dict:
                    print(f"警告: 染色体 '{chrom}' 对基因 '{gene_id}' 在基因组中未找到。跳过。", file=sys.stderr)
                    error_count += 1
                    continue
                if tss_anchor is None:
                    print(f"警告: 基因 '{gene_id}' 未找到TSS。跳过。", file=sys.stderr)
                    error_count += 1
                    continue

                genome_seq_record = genome_dict[chrom]
                chrom_len = len(genome_seq_record.seq)
                promoter_seq_final_case = None
                promoter_start, promoter_end = 0, 0
                split_index_in_extracted = -1

                # --- 使用TSS锚点计算坐标 (1-based) ---
                if strand == '+':
                    # tss_anchor是min(cds_starts)
                    promoter_start = max(1, tss_anchor - args.upstream)
                    promoter_end = min(chrom_len, tss_anchor + args.downstream - 1)
                    if promoter_start > promoter_end: promoter_end = promoter_start
                    split_index_in_extracted = max(0, tss_anchor - promoter_start)

                elif strand == '-':
                    # tss_anchor是max(cds_ends)
                    promoter_start = max(1, tss_anchor - args.downstream + 1) # 下游区域具有较低的坐标
                    promoter_end = min(chrom_len, tss_anchor + args.upstream)    # 上游区域具有较高的坐标
                    if promoter_start > promoter_end: promoter_start = promoter_end
                    # 反向互补后，分割索引计算在概念上保持不变
                    split_index_in_extracted = max(0, promoter_end - tss_anchor)
                else:
                     # 应该已经被find_tss捕获，但再次检查
                    error_count += 1
                    print('链信息错误')
                    continue

                # --- 序列提取 (0-based切片) ---
                try:
                    seq_0_start = promoter_start - 1
                    seq_0_end = promoter_end

                    if seq_0_start < 0 or seq_0_end > chrom_len or seq_0_start >= seq_0_end:
                        raise IndexError(f"无效的切片坐标 [{seq_0_start}:{seq_0_end}] 染色体长度 {chrom_len}")

                    promoter_seq_slice = genome_seq_record.seq[seq_0_start:seq_0_end]
                    extracted_seq = str(promoter_seq_slice)

                    if not extracted_seq:
                        error_count += 1
                        continue

                    # --- 处理链和大小写 ---
                    if strand == '-':
                        extracted_seq = str(Seq(extracted_seq).reverse_complement())

                    len_extracted = len(extracted_seq)
                    if split_index_in_extracted >= len_extracted:
                         promoter_seq_final_case = extracted_seq.lower()
                    elif split_index_in_extracted <= 0:
                         promoter_seq_final_case = extracted_seq.upper()
                    else:
                         upstream_part = extracted_seq[:split_index_in_extracted].upper()
                         downstream_part = extracted_seq[split_index_in_extracted:].lower()
                         promoter_seq_final_case = upstream_part + downstream_part

                    # --- 写入BED文件 ---
                    # BED格式使用0-based坐标系统，所以需要将起始位置减1
                    bed_start = promoter_start - 1
                    name = f"{gene_id}_promoter"
                    score = 1000  # 默认分数
                    bed_file.write(f"{chrom}\t{bed_start}\t{promoter_end}\t{name}\t{score}\t{strand}\n")
                    
                    # --- 写入FASTA文件 ---
                    # 清理gene_id用于FASTA标题（替换空格等）- 可选
                    safe_gene_id = re.sub(r'\s+', '_', str(gene_id)) # 基本清理
                    fasta_file.write(f">{safe_gene_id} promoter_cds_anchor {chrom}:{promoter_start}-{promoter_end}({strand}) anchor:{tss_anchor}\n")
                    for i in range(0, len(promoter_seq_final_case), 60):
                        fasta_file.write(promoter_seq_final_case[i:i+60] + "\n")
                    
                    success_count += 1

                except IndexError as e:
                    print(f"基因 '{gene_id}' ({chrom}:{promoter_start}-{promoter_end}) 切片序列时出错: {e}。跳过。", file=sys.stderr)
                    error_count += 1
                except Exception as e:
                    print(f"处理 '{gene_id}' 启动子时意外错误: {type(e).__name__} - {e}", file=sys.stderr)
                    error_count += 1
                # --- 核心提取逻辑结束 ---

    except IOError as e:
        print(f"错误: 无法写入输出文件: {e}", file=sys.stderr)
        # 如果文件写入中途失败，无法准确报告计数

    print(f"\n处理了 {gff_gene_groups_processed} 个GFF基因组。")
    print(f"成功将 {success_count} 个启动子区域写入BED文件: {bed_output}")
    print(f"成功将 {success_count} 个启动子序列写入FASTA文件: {fasta_output}")
    if error_count > 0:
         print(f"遇到 {error_count} 个错误/警告 (跳过的基因)。")


def extract_promoters(args):
    # --- Load Genome into Memory ---
    print(f"\nLoading genome file '{args.genome}' into memory...")
    try:
        genome_dict = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
        print(f"Loaded {len(genome_dict)} sequences from genome file.")
        if not genome_dict:
             print("Error: No sequences loaded from the genome file.", file=sys.stderr)
             sys.exit(1)
    except FileNotFoundError:
        print(f"Error: Genome file not found at '{args.genome}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading genome file '{args.genome}': {type(e).__name__} - {e}", file=sys.stderr)
        sys.exit(1)

    # --- Process GFF and Extract Promoters Gene by Gene ---
    print(f"\nProcessing GFF and extracting promoters to: {args.output}...")
    success_count = 0
    error_count = 0
    gff_gene_groups_processed = 0

    try:
        with open(args.output, 'w') as outfile:
            # Iterate through gene feature lists yielded by the parser
            for gene_feature_list in parse_gff3(args.gff):
                gff_gene_groups_processed += 1
                if not gene_feature_list: continue # Skip empty groups

                # Find the anchor point (TSS proxy) using your function
                tss_info = find_tss(gene_feature_list)

                if tss_info is None:
                    error_count += 1
                    continue # Skip if TSS anchor couldn't be determined

                gene_id, chrom, tss_anchor, strand = tss_info

                # --- Core Extraction Logic (adapted from previous function) ---
                if chrom not in genome_dict:
                    print(f"Warning: Chromosome '{chrom}' for gene '{gene_id}' not found in genome. Skipping.", file=sys.stderr)
                    error_count += 1
                    continue
                if tss_anchor is None:
                    print(f"Warning: gene '{gene_id}' not found tts. Skipping.", file=sys.stderr)
                    error_count += 1
                    continue

                genome_seq_record = genome_dict[chrom]
                chrom_len = len(genome_seq_record.seq)
                promoter_seq_final_case = None
                promoter_start, promoter_end = 0, 0
                split_index_in_extracted = -1

                # --- Coordinate Calculation (1-based) using tss_anchor ---
                if strand == '+':
                    # tss_anchor is min(cds_starts)
                    promoter_start = max(1, tss_anchor - args.upstream)
                    promoter_end = min(chrom_len, tss_anchor + args.downstream - 1)
                    if promoter_start > promoter_end: promoter_end = promoter_start
                    split_index_in_extracted = max(0, tss_anchor - promoter_start)

                elif strand == '-':
                    # tss_anchor is max(cds_ends)
                    promoter_start = max(1, tss_anchor - args.downstream + 1) # Downstream biological region has lower coordinates
                    promoter_end = min(chrom_len, tss_anchor + args.upstream)    # Upstream biological region has higher coordinates
                    if promoter_start > promoter_end: promoter_start = promoter_end
                    # Split index calculation remains the same conceptually after reverse complement
                    split_index_in_extracted = max(0, promoter_end - tss_anchor)
                else:
                     # Should have been caught by find_tss, but double check
                    error_count += 1
                    print('stand errotr')
                    continue

                # --- Sequence Extraction (0-based slicing) ---
                try:
                    seq_0_start = promoter_start - 1
                    seq_0_end = promoter_end

                    if seq_0_start < 0 or seq_0_end > chrom_len or seq_0_start >= seq_0_end:
                        raise IndexError(f"Invalid slice coords [{seq_0_start}:{seq_0_end}] for chrom len {chrom_len}")

                    promoter_seq_slice = genome_seq_record.seq[seq_0_start:seq_0_end]
                    extracted_seq = str(promoter_seq_slice)

                    if not extracted_seq:
                        error_count += 1
                        continue

                    # --- Handle Strand and Case ---
                    if strand == '-':
                        extracted_seq = str(Seq(extracted_seq).reverse_complement())

                    len_extracted = len(extracted_seq)
                    if split_index_in_extracted >= len_extracted:
                         promoter_seq_final_case = extracted_seq.lower()
                    elif split_index_in_extracted <= 0:
                         promoter_seq_final_case = extracted_seq.upper()
                    else:
                         upstream_part = extracted_seq[:split_index_in_extracted].upper()
                         downstream_part = extracted_seq[split_index_in_extracted:].lower()
                         promoter_seq_final_case = upstream_part + downstream_part

                    # --- Write to single file ---
                    # Clean up gene_id for FASTA header (replace spaces, etc.) - Optional
                    safe_gene_id = re.sub(r'\s+', '_', str(gene_id)) # Basic cleaning
                    outfile.write(f">{safe_gene_id} promoter_cds_anchor {chrom}:{promoter_start}-{promoter_end}({strand}) anchor:{tss_anchor}\n")
                    for i in range(0, len(promoter_seq_final_case), 60):
                        outfile.write(promoter_seq_final_case[i:i+60] + "\n")
                    success_count += 1

                except IndexError as e:
                    print(f"Error slicing sequence for gene '{gene_id}' ({chrom}:{promoter_start}-{promoter_end}): {e}. Skipping.", file=sys.stderr)
                    error_count += 1
                except Exception as e:
                    print(f"Unexpected error processing promoter for '{gene_id}': {type(e).__name__} - {e}", file=sys.stderr)
                    error_count += 1
                # --- End Core Extraction Logic ---

    except IOError as e:
        print(f"Error: Could not write to output file '{args.output}': {e}", file=sys.stderr)
        # No good way to report counts accurately if file writing fails mid-way

    print(f"\nProcessed {gff_gene_groups_processed} gene groups from GFF.")
    print(f"Successfully wrote {success_count} promoter sequences to: {args.output}")
    if error_count > 0:
         print(f"Encountered {error_count} errors/warnings (genes skipped).")


def main():
    """
    主函数，解析命令行参数并执行相应操作
    """
    parser = argparse.ArgumentParser(description='基因组区域处理工具')
    subparsers = parser.add_subparsers(dest='command', help='可用命令')
    
    # 添加BED提取子命令
    bed_parser = subparsers.add_parser('extract-bed', help='从GFF文件提取区域并生成BED文件')
    bed_parser.add_argument('-i', '--input', required=True, help='输入GFF文件路径')
    bed_parser.add_argument('-o', '--output', required=True, help='输出BED文件路径')
    bed_parser.add_argument('-t', '--type', default='gene', help='要提取的特征类型 (默认: gene)')
    bed_parser.add_argument('-n', '--name', default='ID', help='用作名称的属性字段 (默认: ID)')
    
    # 添加启动子提取子命令
    promoter_parser = subparsers.add_parser('extract-promoters', help='提取启动子序列')
    promoter_parser.add_argument('-g', '--gff', required=True, help='GFF文件路径')
    promoter_parser.add_argument('--genome', required=True, help='基因组FASTA文件路径')
    promoter_parser.add_argument('-o', '--output', required=True, help='输出FASTA文件路径')
    promoter_parser.add_argument('-u', '--upstream', type=int, default=2000, help='上游区域长度 (默认: 2000)')
    promoter_parser.add_argument('-d', '--downstream', type=int, default=200, help='下游区域长度 (默认: 200)')
    
    # 添加启动子BED和FASTA同时提取子命令
    both_parser = subparsers.add_parser('extract-promoters-both', help='同时提取启动子序列的BED文件和FASTA文件')
    both_parser.add_argument('-g', '--gff', required=True, help='GFF文件路径')
    both_parser.add_argument('--genome', required=True, help='基因组FASTA文件路径')
    both_parser.add_argument('-o', '--output', required=True, help='输出文件前缀，将生成.bed和.fasta两个文件')
    both_parser.add_argument('-u', '--upstream', type=int, default=2000, help='上游区域长度 (默认: 2000)')
    both_parser.add_argument('-d', '--downstream', type=int, default=200, help='下游区域长度 (默认: 200)')
    
    args = parser.parse_args()
    
    if args.command == 'extract-bed':
        extract_bed_regions(args.input, args.output, args.type, args.name)
    elif args.command == 'extract-promoters':
        extract_promoters(args)
    elif args.command == 'extract-promoters-both':
        extract_promoters_both(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()