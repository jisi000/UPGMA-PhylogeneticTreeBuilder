"""
UPGMA（非加权配对组平均法）建树算法
核心假设：假设进化速率恒定（分子钟假说）
算法：小顶堆优化距离查找
时间复杂度：朴素实现 O(n³)，堆优化后 O(n² log n)

UPGMA 步骤
  每一步找到距离最近的两个簇，将它们合并为一个新簇，
  新簇与其他簇的距离为两子簇距离的加权平均（按大小加权）。
  重复 n-1 次，直到只剩一个根节点。
"""
import heapq
from tree import TreeNode


def upgma(names: list[str], dist_matrix: list[list[float]]) -> TreeNode:
    n = len(names)
    # 深拷贝距离矩阵
    d = [row[:] for row in dist_matrix]
    # 当前簇高度
    height = [0.0] * n
    # 当前簇大小
    size = [1] * n
    # 树节点
    nodes = [TreeNode(name=name)for name in names]
    # 是否仍为活跃簇
    active = [True] * n
    # 小顶堆
    heap = []
    for i in range(n):
        for j in range(i + 1, n):
            heapq.heappush(heap,(d[i][j], i, j))
    new_idx = n
    for _ in range(n - 1):
        # 找最近的两个活跃簇
        while True:
            dist, i, j = heapq.heappop(heap)
            if active[i] and active[j]:
                break
        # 枝长计算
        branch_i = dist / 2.0 - height[i]
        branch_j = dist / 2.0 - height[j]
        nodes[i].branch_length = max(branch_i, 0.0)
        nodes[j].branch_length = max(branch_j, 0.0)
        # 新内部节点
        new_node = TreeNode(name=f"internal_{new_idx}")
        new_node.left = nodes[i]
        new_node.right = nodes[j]
        new_height = dist / 2.0
        new_size = size[i] + size[j]
        # 计算新簇距离
        new_distances = [0.0] * new_idx
        for k in range(new_idx):
            if active[k] and k != i and k != j:
                new_distances[k] = (size[i] * d[i][k]+ size[j] * d[j][k]) / new_size
        # 扩展矩阵
        # 保持对称
        for row_idx in range(len(d)):
            d[row_idx].append(new_distances[row_idx])
        d.append(new_distances + [0.0])
        # 更新辅助数组
        height.append(new_height)
        size.append(new_size)
        nodes.append(new_node)
        active.append(True)
        active[i] = False
        active[j] = False
        # 新簇入堆
        for k in range(new_idx):
            if active[k]:
                heapq.heappush(heap,(new_distances[k],min(k, new_idx),max(k, new_idx)))
        new_idx += 1
    # 找最后活跃簇
    root_idx = -1
    for i in range(len(active)):
        if active[i]:
            root_idx = i
            break
    return nodes[root_idx]