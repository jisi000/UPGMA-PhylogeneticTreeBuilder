# UPGMA-PhylogeneticTreeBuilder 🌿

> 基于UPGMA算法的系统发育树构建与可视化工具

生命科学学院，孙婧斯

------

## 项目背景

**系统发育树（Phylogenetic Tree）**，也称进化树，是生命科学中描述物种或基因进化关系的核心工具。

进化论是生物学的核心理论。从20世纪60年代，分子进化的研究开始。分子进化上的一个讨论热点是分子钟理论，该理论认为在特定的系统发育世系中核苷酸替换速率是一个恒定的常数。因此DNA分子之间的核苷酸差异数即可代表该同源基因分开的时间。根据这一理论，分子钟被用来测定物种是从何时从某一个祖先开始分开的。分子钟可根据分子的遗传信息决定生物间的分类，推算系统发育事件的发生时间，是分子系统发育分析中的重要技术。

目前用的最多的系统发育树构建方法有三种：距离法、最大简约法、最大似然法。其中距离法速度最快，在序列由较高相似性时适用。距离法是首先算出序列间的遗传距离，再根据这些距离将序列分别依次合并的聚类分析方法。

距离的计算有多种方法，本项目使用最基础的hamming距离与常用的Jukes-Can-tor距离K

**非加权组平均法**（unweighted pair group method with arithmetic mean，简称**UPGMA**）是最早的距离法。

本项目运用生物学原理和”数据结构与算法“中**二叉树**、**括号嵌套树**、**堆**、**贪心算法**等知识，实现根据输入的一组fasta格式的基因序列构建系统发生树，并提供可视化和交互界面。展现了对计算机基础知识的跨学科应用和对生物信息学的初步探索。

------

## 核心算法

### 1. 数据结构：二叉树

- 使用 `TreeNode` 类存储进化树。每个叶子节点包含物种名、枝长。内部节点代表系统发育树上的分支点，没有物种名，有左右子节点
- **后序递归遍历**生成标准 Newick 格式字符串：`((A:0.1,B:0.2):0.3,C:0.4)`，类似**括号嵌套树**，数字代表该节点到父节点的枝长

### 2. 距离计算

DNA序列有ATGC四种可能，对于genbank中的基因序列资源，可能含有N（未定）或者 -（空位）等无效位点，计算距离时将两个**长度相同**的序列共同遍历，统计不同点位和有效点位数量。

| 方法              | 公式                          | 适用场景                         |
| ----------------- | ----------------------------- | -------------------------------- |
| Hamming 距离      | `d = 不同位点数 / 有效位点数` | 序列差异较小，认为位点仅单次突变 |
| Jukes-Cantor 距离 | `K = -3/4 · ln(1 - 4d/3)`     | 考虑多次突变的进化模型修正       |

时间复杂度：O(n² · L)，n 为序列数，L 为序列长度。

### 3. UPGMA 算法

思路：类似构建**哈夫曼编码树**的思路（贪心 + 小顶堆）

- 核心思想：每步合并距离最近的两个簇，新簇与其他簇的距离取加权平均
- 时间复杂度：朴素实现为 O(n³)，用 `heapq` 维护所有簇对距离，避免每步 O(n²) 全局扫描，堆优化后 O(n² log n)
- 假设：分子钟（进化速率恒定），适合同源性高的序列

```bash
每轮迭代：
  1. heappop 取最小距离对 (i, j)
  2. 计算枝长：branch_i = d(i,j)/2 - height[i] #根据系统发生树本身的结构特点
  3. 创建内部节点，更新距离矩阵（加权均值公式）#①
  new_distances[k] = (size[i] * d[i][k]+ size[j] * d[j][k]) / new_size
  4. heappush 新簇与所有活跃簇的距离
重复 n-1 次 → 得到根节点
```

> 对①加权均值公式的解释：
>
> **UPGMA 叫 “非加权”，但计算时反而要加权**：
>
> 它的 “非加权” 是指**不对合并后的簇赋予额外权重**，而是**按每个原始物种的数量来加权**。
>
> 比如合并了含 2 个物种的簇 A 和含 3 个物种的簇 B，新簇与簇 C 的距离是：
>
> \(d(AB,C) = \frac{2 \times d(A,C) + 3 \times d(B,C)}{2+3}\)
>
> 这样每个原始物种的贡献是均等的，没有偏向大簇或小簇。

### 4.可视化

- 使用biopython和matplotlib，将Newick格式的字符串转化为系统发生树的图
- 使用Tkinter制作用户界面

## 运行指南

### 环境要求

**Python 版本：** 3.9 及以上

**推荐 IDE：** PyCharm（Community 或 Professional 版均可）

**安装依赖：**
在 PyCharm 中打开终端，执行：
    pip install -r requirements.txt

（文件应在根目录下）

**PyCharm 插件：**
首次用 PyCharm 打开 `.fasta` 格式文件时，编辑器会自动弹出插件安装提示，按提示安装即可获得语法高亮支持。

### 操作步骤

1.在本地准备.fasta格式文件

格式展示

![image-20260607004603530](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260607004603530.png)

2.将 `tree.py`、`upgma.py`、`sequence.py`、`gui.py`放在根目录下，确保文件之间的import可执行

![image-20260606232743492](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260606232743492.png)

3.运行`gui.py`

![image-20260606232852273](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260606232852273.png)

4.点击**Browse FASTA**打开目标文件（会显示文件名），选择一个**Distance Model**（model介绍见上文）

![image-20260606233438991](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260606233438991.png)

5.点击**Build Tree**建树

Phylogenetic Tree, Distance Matrix, Newick都会出现内容，同时底部状态变为‘Tree construction completed.'

 ![image-20260606233947070](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260606233947070.png)

其中Distance Matrix可以滚动查看，Phologenetic Tree可以点击**Vieew Full Size**查看大图和保存图片

![image-20260606234538384](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260606234538384.png)

6.说明

由于UPGMA算法适合对序列差异小的基因序列进行分析和建树，并且UPGMA算法本身不支持长度不同的序列比较，本项目对输入序列进行了截取到最短序列的长度的暴力处理，当序列长度差异不超过5%时可认为是可信的。但是序列长度差异更大的基因不建议用此方法建树。本项目遇到这类序列组会发出警告。

![image-20260607002615263](C:\Users\嵇兕\AppData\Roaming\Typora\typora-user-images\image-20260607002615263.png)

用户可以使用提供的数据进行体验，也可以到Genbank下载fasta格式的基因序列，这里提供NCBI的整合链接[Home - Nucleotide - NCBI](https://www.ncbi.nlm.nih.gov/nucleotide/)

构造数据时，建议选取一个比较保守的同源基因，选取亲缘关系较近的物种，也建议选取一个外类群作为参照点

例如上图所示的例子，使用 GAPDH（甘油醛 - 3 - 磷酸脱氢酶）基因，选取了*Neptis*（环蛱蝶属）的20个物种，还有外类群*Pantoporia*（蟠蛱蝶属）

**提供的四个测试数据的文件名即为基因名称。**

------

## AI 工具与参考资源声明

本项目使用**ChatGPT** 生成 GUI 的分模块代码，由本人进行整合和调整。

核心算法逻辑（`upgma.py`、`tree.py`、`sequence.py`）由本人在参考以下教材后编写，但使用了**Claud**润色（添加注释、规范函数格式）和debug。

项目背景与原理解释都参考了以下教材。

参考教材：《生物信息学》许忠能 P213-240

------

