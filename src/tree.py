"""
树节点定义与 Newick 格式序列化

数据结构：二叉树（TreeNode）
核心操作：后序递归遍历生成 Newick 字符串
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class TreeNode:
    """
    系统发育树节点。
    - 叶节点：name 为物种名，left/right 为 None
    - 内部节点：name 为空，left/right 为子节点
    """
    name: str = ""
    branch_length: float = 0.0
    left: Optional["TreeNode"] = field(default=None, repr=False)
    right: Optional["TreeNode"] = field(default=None, repr=False)

    def is_leaf(self) -> bool:
        return self.left is None and self.right is None

    def to_newick(self) -> str:
        """
        后序递归遍历，生成标准 Newick 格式字符串。
        例：((A:0.1,B:0.2):0.3,C:0.4)

        时间复杂度：O(n)，n 为节点总数
        """
        if self.is_leaf():
            return f"{self.name}:{self.branch_length:.4f}"
        left_str  = self.left.to_newick()  if self.left  else ""
        right_str = self.right.to_newick() if self.right else ""
        inner = f"({left_str},{right_str})"
        if self.branch_length > 0:
            return f"{inner}:{self.branch_length:.4f}"
        return inner

    def get_leaves(self) -> list[str]:
        """前序遍历收集所有叶节点名称"""
        if self.is_leaf():
            return [self.name]
        result = []
        if self.left:
            result += self.left.get_leaves()
        if self.right:
            result += self.right.get_leaves()
        return result

    def depth(self) -> int:
        """递归计算树的最大深度"""
        if self.is_leaf():
            return 0
        left_d  = self.left.depth()  if self.left  else 0
        right_d = self.right.depth() if self.right else 0
        return 1 + max(left_d, right_d)
