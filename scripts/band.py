import sys, re

hapblock = sys.argv[1]

def hapcut2_to_ideogram(input_file, output_file):
    # 初始化变量
    current_block = {"chr": None, "positions": [], "value": 1}
    block_counter = 0  # 用于交替颜色
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        # 写入表头
        f_out.write("Chr\tStart\tEnd\tValue\n")
        
        for line in f_in:
            line = line.strip()
            
            # 检测BLOCK起始行
            if line.startswith("BLOCK:"):
                # 如果上一个区块未处理则先处理（理论上不会出现）
                if current_block["positions"]:
                    _write_block(current_block, f_out)
                # 初始化新block
                current_block = {"chr": None, "positions": [], "value": block_counter % 2 + 1}
                block_counter += 1  # 颜色交替
                
            # 检测区块结束分隔符
            elif line.startswith("********"):
                if current_block["positions"]:
                    _write_block(current_block, f_out)
                current_block = {"chr": None, "positions": [], "value": None}
                
            # 处理数据行
            elif line and not line.startswith("#"):  # 忽略注释行
                parts = re.split(r'\s+', line)  # 兼容制表符和空格分隔
                if len(parts) >= 5:
                    chr = parts[3].replace("chr", "")  # 移除chr前缀（根据karyotype文件调整）
                    pos = int(parts[4])
                    if not current_block["chr"]:
                        current_block["chr"] = chr
                    # 确保同一区块内染色体一致
                    if current_block["chr"] == chr:
                        current_block["positions"].append(pos)
                    else:
                        print(f"Warning: Chromosome mismatch in block {block_counter}")

def _write_block(block, f_out):
    if block["positions"]:
        start = min(block["positions"])
        end = max(block["positions"])
        f_out.write(f"{block['chr']}\t{start}\t{end}\t{block['value']}\n")

# 调用示例
hapcut2_to_ideogram(hapblock, "contigs.tsv")
