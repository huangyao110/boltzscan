import sys
import re

def meme_to_transfac(meme_file, output_file):
    print(f"Converting {meme_file} to {output_file} ...")
    
    with open(meme_file, 'r') as f_in, open(output_file, 'w') as f_out:
        in_matrix = False
        motif_id = ""
        motif_alt = ""
        matrix_rows = []
        counter = 1
        
        # 正则表达式匹配 MEME 的 Motif 定义行
        # 格式通常是: MOTIF <id> <alt_name>
        motif_header_re = re.compile(r"^MOTIF\s+(\S+)(?:\s+(.*))?")
        
        for line in f_in:
            line = line.strip()
            if not line: continue
            
            # 1. 识别 Motif 开始
            match = motif_header_re.match(line)
            if match:
                # 如果上一个 Motif 还没写完（理论上不会，除非没有矩阵数据），先处理
                if in_matrix and matrix_rows:
                    write_transfac_entry(f_out, motif_id, motif_alt, matrix_rows)
                    matrix_rows = []
                
                motif_id = match.group(1)
                motif_alt = match.group(2) if match.group(2) else ""
                # 如果 ID 重复或者为空，用计数器兜底
                if not motif_id:
                    motif_id = f"M{counter:04d}"
                counter += 1
                in_matrix = False # 等待 letter-probability matrix 标记
                continue
            
            # 2. 识别矩阵数据块开始
            if line.startswith("letter-probability"):
                in_matrix = True
                matrix_rows = []
                continue
            
            # 3. 读取矩阵行
            if in_matrix:
                # MEME 矩阵行通常是以数字或小数点开头的
                # 例如: 0.123  0.456  0.000  0.421
                parts = line.split()
                # 简单的合法性检查：必须是4个数值
                if len(parts) == 4 and is_number(parts[0]):
                    matrix_rows.append(parts)
                else:
                    # 遇到非数字行（如 URL 或下一个 MOTIF），说明当前矩阵结束
                    if matrix_rows:
                        write_transfac_entry(f_out, motif_id, motif_alt, matrix_rows)
                        matrix_rows = []
                    in_matrix = False
                    
                    # 重新检查这行是不是新 Motif
                    # (虽然上面的 loop 逻辑会处理下一行，但这里防止丢失当前行是 MOTIF 的情况)
                    match_new = motif_header_re.match(line)
                    if match_new:
                        motif_id = match_new.group(1)
                        motif_alt = match_new.group(2) if match_new.group(2) else ""
                        in_matrix = False

        # 4. 处理文件末尾最后一个 Motif
        if in_matrix and matrix_rows:
            write_transfac_entry(f_out, motif_id, motif_alt, matrix_rows)

    print("Conversion Done!")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def write_transfac_entry(f, m_id, m_alt, rows):
    """写入单个 TRANSFAC 块"""
    # AC: Accession number
    f.write(f"AC  {m_id}\n")
    f.write("XX\n")
    
    # ID: Identifier (可以是 alt name，也可以是 id)
    # TRANSFAC 习惯 ID 放名字，NA 放描述
    display_id = m_id
    if m_alt:
        display_id = f"{m_id}_{m_alt}"
        # 去掉可能的非法字符
        display_id = display_id.replace(" ", "_")
    
    f.write(f"ID  {display_id}\n")
    f.write("XX\n")
    
    # P0: Position Matrix (A C G T)
    # TRANSFAC 的标准列顺序是 A, C, G, T
    f.write("PO      A      C      G      T\n")
    
    for i, row in enumerate(rows):
        # 格式化：位置编号(从1开始)  A  C  G  T
        # MEME 是概率(0-1)，TRANSFAC 也可以接受概率，或者 Count
        # 我们直接保留 MEME 的原值
        pos = i + 1
        # 使用制表符分隔
        f.write(f"{pos:02d}\t{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\n")
        
    f.write("XX\n")
    f.write("//\n") # 结束标记

if __name__ == "__main__":
    # 用法: python meme2transfac.py <input.meme> <output.transfac>
    if len(sys.argv) < 3:
        print("Usage: python meme2transfac.py input.meme output.transfac")
        # 默认测试
        # meme_to_transfac("all_meme.meme", "all_motifs.transfac")
    else:
        meme_to_transfac(sys.argv[1], sys.argv[2])