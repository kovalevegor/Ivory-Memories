# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 18:28:21 2025

@author: gogak
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

def read_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split())
            points.append((x, y))
    return np.array(points)

def read_triangles(filename):
    triangles = []
    with open(filename, 'r') as file:
        for line in file:
            coords = list(map(float, line.strip().split()))
            # Каждая строка должна содержать координаты трех вершин: (x1,y1), (x2,y2), (x3,y3)
            triangle = [(coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])]
            triangles.append(triangle)
    return triangles

def plot_triangles(points, triangles, point_size=30, grid_width=0.8, font_size=14):
    # Настройка цветовой схемы
    background_color = "#131416"
    grid_color = "#1f2022"
    point_color = "#86a67d"
    text_color = "#ffffff"
    frame_color = "#1f2022"
    triangle_edge_color = "#e07a5f"  # цвет границы треугольника
    triangle_face_color = "none"      # прозрачная заливка

    # Создание фигуры и осей
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')

    # Установка пределов осей по данным из knots.txt
    x_min, x_max = np.min(points[:, 0]), np.max(points[:, 0])
    y_min, y_max = np.min(points[:, 1]), np.max(points[:, 1])
    ax.set_xlim(x_min - 0.1, x_max + 0.1)
    ax.set_ylim(y_min - 0.1, y_max + 0.1)

    # Настройка рамки и сетки
    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Отображение точек (фиксированное изображение)
    ax.scatter(points[:, 0], points[:, 1], s=point_size, color=point_color)

    # Список для хранения добавленных патчей (треугольников)
    drawn_triangles = []

    # Функция анимации: на каждом кадре добавляется один треугольник
    def animate(i):
        if i < len(triangles):
            coords = triangles[i]
            triangle = patches.Polygon(coords, closed=True,
                                       edgecolor=triangle_edge_color,
                                       facecolor=triangle_face_color,
                                       linewidth=2)
            ax.add_patch(triangle)
            drawn_triangles.append(triangle)
        return drawn_triangles

    ani = FuncAnimation(
        fig=fig,
        func=animate,
        frames=len(triangles),
        interval=0.001,  # время между появлением треугольников в мс (можно подстроить)
        blit=False,
        repeat=False
    )

    ani.save('animation_T.gif', writer='pillow', dpi=100)
    plt.close()

def main():
    points = read_points('knots.txt')
    triangles = read_triangles('triangles.txt')
    plot_triangles(points, triangles)

if __name__ == '__main__':
    main()
