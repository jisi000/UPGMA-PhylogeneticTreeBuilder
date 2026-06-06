"""
序列解析与距离矩阵计算

  - FASTA 格式解析
  - Hamming 距离（直接碱基差异比例）
  - Jukes-Cantor 距离（进化模型修正）
"""

import re
import math
from typing import Optional

#FASTA 解析
def parse_fasta(text: str) -> dict[str, str]:
    """
    解析 FASTA 格式文本，返回 {名称: 序列} 字典。

    示例输入：
        >Human
        ATGCATGC
        >Chimp
        ATGCTTGC
    """
    sequences: dict[str, str] = {}
    current_name: Optional[str] = None
    current_seq: list[str] = []

    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_name is not None:
                sequences[current_name] = "".join(current_seq).upper()
            current_name = line[1:].split()[1]   # 取第二个 token 作为名称，因为Genbank下载的fasta文件通常是accession-name
            current_seq = []
        else:
            # 只保留合法碱基字符
            current_seq.append(re.sub(r"[^ATGCatgcNn-]", "", line))

    if current_name is not None:
        sequences[current_name] = "".join(current_seq).upper()

    return sequences


#距离计算

#预处理
def trim_to_shortest(sequences: dict[str, str]) -> dict[str, str]:
    """截取所有序列到最短长度（仅建议用于长度差异 < 5% 的情况）"""
    min_len = min(len(s) for s in sequences.values())
    trimmed = {name: seq[:min_len] for name, seq in sequences.items()}

    max_len = max(len(s) for s in sequences.values())
    if (max_len - min_len) / max_len > 0.05:
        import warnings
        warnings.warn(f"序列长度差异达 {max_len - min_len}bp "f"({(max_len - min_len) / max_len:.1%})，截断结果可能有偏差，"f"建议使用其他方法比对。")
    return trimmed
#计算汉明距离（只能处理长度一致的序列）
def hamming_distance(seq1: str, seq2: str) -> float:
    """
    Hamming 距离：两序列中碱基不同的位点比例。
    仅统计双方均非 gap('-') 且非 N 的位点。

    时间复杂度：O(L)，L 为序列长度
    """
    if len(seq1) != len(seq2):
        raise ValueError(f"序列长度不一致：{len(seq1)} vs {len(seq2)}")

    mismatches = 0
    valid_sites = 0
    for a, b in zip(seq1, seq2):
        if a in "-N" or b in "-N":
            continue
        valid_sites += 1
        if a != b:
            mismatches += 1

    if valid_sites == 0:
        return 0.0
    return mismatches / valid_sites

def jukes_cantor_distance(seq1: str, seq2: str) -> float:
    """
    Jukes-Cantor 距离：对 Hamming 距离做进化模型修正。
    公式：d = -3/4 * ln(1 - 4/3 * p)
    其中 p = Hamming 距离

    当 p >= 0.75 时模型不适用，返回 Hamming 距离作为近似。
    """
    p = hamming_distance(seq1, seq2)
    if p == 0.0:
        return 0.0
    if p >= 0.75:
        # 超出模型适用范围，退回 Hamming
        return p
    return -0.75 * math.log(1 - (4.0 / 3.0) * p)


#距离矩阵

def build_distance_matrix(sequences: dict[str, str],method: str) -> tuple[list[str], list[list[float]]]:
    """
    计算所有序列对之间的距离，返回 (名称列表, 距离矩阵)。

    method: "hamming" | "jukes_cantor"

    时间复杂度：O(n² * L)，n = 序列数，L = 序列长度
    空间复杂度：O(n²)
    """
    names = list(sequences.keys())
    n = len(names)
    dist_fn = jukes_cantor_distance if method == "jukes_cantor" else hamming_distance
    legal_sequences=trim_to_shortest(sequences)
    matrix: list[list[float]] = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            d = dist_fn(legal_sequences[names[i]], legal_sequences[names[j]])
            matrix[i][j] = d
            matrix[j][i] = d   # 对称矩阵

    return names, matrix


def print_distance_matrix(names: list[str], matrix: list[list[float]]) -> None:
    """打印距离矩阵（调试用）"""
    width = max(len(n) for n in names) + 2
    header = " " * width + "  ".join(f"{n:>{width}}" for n in names)
    print(header)
    for i, name in enumerate(names):
        row = "  ".join(f"{matrix[i][j]:>{width}.4f}" for j in range(len(names)))
        print(f"{name:<{width}}{row}")