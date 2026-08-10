from Bio import SeqIO
from Bio.Seq import Seq
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, List, Generator, Any
import sys
import re
import gzip


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
                # Prefer explicit '###' separators (handled above). When a GFF
                # has none (e.g. coreset.gff), group implicitly per gene: a new
                # top-level 'gene' feature closes the previous group. Without
                # this, the whole file collapses into one group and find_tss()
                # emits at most a single (wrong) promoter.
                if type_ == 'gene' and current_gene:
                    yield current_gene
                    current_gene = []

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


def extract_promoter_regions(args, output_format=None):
    """
    从GFF文件中提取启动子区域，输出 BED、FASTA 或两者。
    
    参数:
        args: 包含必要参数的命名空间对象
            - gff: GFF文件路径
            - genome: 基因组FASTA文件路径
            - output: 输出文件前缀
            - format: bed、fasta 或 both
            - upstream: 上游区域长度
            - downstream: 下游区域长度
    """
    output_format = output_format or getattr(args, 'format', 'both')
    if output_format not in {'bed', 'fasta', 'both'}:
        raise ValueError(f"Unsupported promoter output format: {output_format}")

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
    # Accept either a clean prefix or a familiar .bed/.fa/.fasta filename.
    output_path = Path(args.output)
    output_prefix = (
        output_path.with_suffix('')
        if output_path.suffix.lower() in {'.bed', '.fa', '.fasta'}
        else output_path
    )
    bed_output = str(output_prefix) + ".bed" if output_format in {'bed', 'both'} else None
    fasta_output = str(output_prefix) + ".fasta" if output_format in {'fasta', 'both'} else None
    
    outputs = [path for path in (bed_output, fasta_output) if path is not None]
    print(f"\n处理GFF并提取启动子区域到: {', '.join(outputs)}...")
    success_count = 0
    error_count = 0
    gff_gene_groups_processed = 0

    try:
        with ExitStack() as stack:
            bed_file = stack.enter_context(open(bed_output, 'w')) if bed_output else None
            fasta_file = stack.enter_context(open(fasta_output, 'w')) if fasta_output else None
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
                    if bed_file is not None:
                        bed_file.write(f"{chrom}\t{bed_start}\t{promoter_end}\t{name}\t{score}\t{strand}\n")
                    
                    # --- 写入FASTA文件 ---
                    # 清理gene_id用于FASTA标题（替换空格等）- 可选
                    if fasta_file is not None:
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
    if bed_output:
        print(f"成功将 {success_count} 个启动子区域写入BED文件: {bed_output}")
    if fasta_output:
        print(f"成功将 {success_count} 个启动子序列写入FASTA文件: {fasta_output}")
    if error_count > 0:
         print(f"遇到 {error_count} 个错误/警告 (跳过的基因)。")
