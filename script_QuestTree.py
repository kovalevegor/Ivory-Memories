# -*- coding: utf-8 -*-
"""
Quest Tree Visualization Script
"""

import matplotlib.pyplot as plt
import networkx as nx
from textwrap import fill

# Описываем дерево квеста вручную на основе твоего текста:
quest_tree = {
    0: {"desc": "Fetch item 0 and bring to knot 128",    "children": [1]},
    1: {"desc": "Obtain item 0",                          "children": [2, 4, 6]},
    2: {"desc": "Defeat enemy 6 to get item 0",           "children": [3]},
    3: {"desc": "Go to knot 128 (after defeating)",       "children": []},
    4: {"desc": "Pick up item 0 from knot 73",            "children": [5]},
    5: {"desc": "Go to knot 128 (after pickup)",          "children": []},
    6: {"desc": "Talk to character 0 to get item 0",     "children": [7]},
    7: {"desc": "Go to knot 128 (after talking)",         "children": []},
}

# Строим ориентированный граф
G = nx.DiGraph()
for nid, data in quest_tree.items():
    G.add_node(nid, label=data["desc"])
    for child in data["children"]:
        G.add_edge(nid, child)

# Расположение: пробуем graphviz, иначе fall back на spring_layout
try:
    pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
except Exception:
    pos = nx.spring_layout(G, seed=42)

# Подготовим обёрнутые подписи (максимум 20 символов в строке)
raw_labels = nx.get_node_attributes(G, "label")
labels = {nid: fill(text, width=20) for nid, text in raw_labels.items()}

# Рисуем
plt.figure(figsize=(12, 8))
# Узлы
nx.draw_networkx_nodes(
    G, pos,
    node_size=2000,
    node_color="#FFDD99",
    edgecolors="#333333"
)
# Стрелки
nx.draw_networkx_edges(
    G, pos,
    arrowstyle="->",
    arrowsize=20,
    edge_color="#555555",
    width=2
)
# Подписи (без wrap, используем уже обёрнутые)
nx.draw_networkx_labels(
    G, pos,
    labels,
    font_size=9,
    font_color="#000000"
)

plt.title("Quest Tree")
plt.axis("off")
plt.tight_layout()
plt.savefig("quest_tree.png", dpi=300)
plt.show()
