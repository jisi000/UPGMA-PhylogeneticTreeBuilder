# =========================
# 构建主窗口
# =========================
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

root = tk.Tk()
root.title("UPGMA Phylogenetic Tree Builder")
root.geometry("1100x800")
root.minsize(900, 600)
root.rowconfigure(0,weight=0)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=1)
root.rowconfigure(3, weight=0)
root.rowconfigure(4,weight=0)
root.columnconfigure(0, weight=1)
# =========================
# ttk样式
# =========================
style = ttk.Style()
style.theme_use("clam")
style.configure("Treeview",rowheight=25,font=("Segoe UI", 9))
style.configure("Treeview.Heading",font=("Segoe UI", 10, "bold"))

# =========================
# 标题
# =========================
title_label = tk.Label(root,text="UPGMA Phylogenetic Tree Builder",font=("Segoe UI", 18, "bold"),fg="#1E3A5F")
title_label.grid(row=0,column=0,pady=10,sticky="ew")

# =========================
# 页面左右布局
# =========================
top_frame = tk.Frame(root)
top_frame.grid(row=1,column=0,sticky='nsew',padx=10)
root.rowconfigure(0,weight=1)
top_frame.columnconfigure(0, weight=1, minsize=260)
top_frame.columnconfigure(1, weight=1)

left_frame = tk.LabelFrame(top_frame,text="Controls",padx=10,pady=10)
left_frame.grid(row=0,column=0,sticky='nw',padx=(0,10))
left_frame.columnconfigure(0, weight=1)
tree_frame = tk.LabelFrame(top_frame,text="Phylogenetic Tree",padx=5,pady=5)
tree_frame.grid(row=0,column=1,sticky='nsew')
view_btn = tk.Button(tree_frame, text="View Full Size")
view_btn.pack(side="top", pady=3)
matrix_frame = tk.LabelFrame(root,text="Distance Matrix")
matrix_frame.grid(row=2,column=0,sticky='nsew',padx=10,pady=5)

# =========================
# 文件选择
# =========================
selected_file = ""
def select_file():
    global selected_file
    selected_file = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta*")])
    if selected_file:
        file_label.config(text=selected_file)

# 浏览按钮
browse_button = tk.Button(left_frame,text="Browse FASTA",command=select_file,bg="#2196F3",fg="white",font=("Segoe UI", 10),width=20)
browse_button.pack(pady=5)

# 文件状态标签
file_label = tk.Label(left_frame,text="No file selected",wraplength=220,justify="left")
file_label.pack(pady=5)


# =========================
# 距离模型选择
# =========================
distance_method = tk.StringVar(value="jukes_cantor")
model_label = tk.Label(left_frame,text="Distance Model")
model_label.pack(pady=(15, 2))
model_box = ttk.Combobox(left_frame,textvariable=distance_method,state="readonly",width=18)
model_box["values"] = ("hamming","jukes_cantor")
model_box.pack(pady=(0, 10))

# =========================
# 导入建树模块
# =========================
from sequence import parse_fasta
from sequence import build_distance_matrix
from upgma import upgma

from Bio import Phylo
from io import StringIO

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

tree_canvas = None
tree_figure=None
# =========================
# 距离矩阵显示
# =========================
def show_distance_matrix(names, matrix):
    table.delete(*table.get_children())
    table["columns"] = ["Species"] + names
    table["show"] = "headings"
    table.heading("Species",text="Species")
    table.column("Species",width=120,anchor="center",stretch=False)
    for name in names:
        table.heading(name,text=name)
        table.column(name,width=90,anchor="center",stretch=False)
    for i in range(len(names)):
        row = []
        for value in matrix[i]:
            row.append(f"{value:.4f}")
        table.insert("","end",values=[names[i]] + row)


# =========================
# 插入树图
# =========================
def display_tree(newick):
    global tree_canvas
    global tree_figure
    if tree_canvas:
        tree_canvas.get_tk_widget().destroy()
    if tree_figure:
        plt.close(tree_figure)
    tree = Phylo.read(StringIO(newick),"newick")
    tree_figure = plt.Figure(figsize=(7, 5),dpi=100)
    ax = tree_figure.add_subplot(111)
    Phylo.draw(tree,axes=ax,do_show=False)
    tree_canvas = FigureCanvasTkAgg(tree_figure,master=tree_frame)
    tree_canvas.draw()
    tree_canvas.get_tk_widget().pack(fill="both",expand=True,after=view_btn)
    def open_full_view():
        fig, ax = plt.subplots(figsize=(12, 8))
        Phylo.draw(tree, axes=ax, do_show=False)
        plt.tight_layout()
        plt.show()

    view_btn.config(command=open_full_view)
# =========================
# 建树主函数
# =========================
def build_tree():
    if not selected_file:
        status_label.config(text="Please select a FASTA file first.")
        return
    status_label.config(text="Building tree...")
    root.update_idletasks()
    try:
        # 读取文件
        with open(selected_file) as f:
            fasta_text = f.read()
        # 计算距离矩阵
        sequences = parse_fasta(fasta_text)
        names, matrix = build_distance_matrix(sequences,method=distance_method.get())
        # 显示距离矩阵
        show_distance_matrix(names,matrix)
        # UPGMA建树
        root_node = upgma(names,matrix)
        # Newick
        newick = root_node.to_newick() + ";"
        result_text.delete("1.0",tk.END)
        result_text.insert(tk.END,newick)
        # 状态更新
        status_label.config(text="Tree construction completed.")
        # 可视化树
        display_tree(newick)
    except Exception as e:
        status_label.config(text=f"Error: {str(e)}")

# 建树按钮
build_button = tk.Button(left_frame,text="Build Tree",command=build_tree,bg="#4CAF50",fg="white",font=("Segoe UI", 10, "bold"),width=20)
build_button.pack(pady=15)

# =========================
# 距离矩阵区域
# =========================
table_frame = tk.Frame(matrix_frame)
table_frame.pack(fill="both",expand=True)
table_frame.grid_rowconfigure(0,weight=1)
table_frame.grid_columnconfigure(0,weight=1)

# Treeview
table = ttk.Treeview(table_frame)

# 滚动条
scroll_y = ttk.Scrollbar(table_frame,orient="vertical",command=table.yview)
scroll_x = ttk.Scrollbar(table_frame,orient="horizontal",command=table.xview)
table.configure(yscrollcommand=scroll_y.set,xscrollcommand=scroll_x.set)
# 布局
table.grid(row=0,column=0,sticky="nsew")
scroll_y.grid(row=0,column=1,sticky="ns")
scroll_x.grid(row=1,column=0,sticky='ew')

# =========================
# 鼠标滚轮支持
# =========================
def on_mousewheel(event):
    table.yview_scroll(
        int(-1 * (event.delta / 120)),"units")
table.bind_all("<MouseWheel>",on_mousewheel)

# =========================
# Newick输出区域
# =========================
result_frame = tk.LabelFrame(root,text="Newick")
result_frame.grid(row=3,column=0,sticky='ew',padx=10,pady=5)

result_text = tk.Text(result_frame,height=8,font=("Consolas", 10),wrap="word")
result_text.pack(fill="both",expand=True,padx=5,pady=5)

# =========================
# 状态栏
# =========================
status_label = tk.Label(root,text="Ready",anchor="w",relief="sunken")
status_label.grid(row=4,column=0,sticky='ew')

# =========================
# 启动
# =========================
root.mainloop()