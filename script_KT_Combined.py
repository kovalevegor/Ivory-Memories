# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 19:50:31 2025

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
            triangle = [(coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])]
            triangles.append(triangle)
    return triangles

def plot_animation(points, triangles, point_size=30, grid_width=0.8, font_size=14):
    # Настройки цветовой схемы
    background_color = "#131416"
    grid_color = "#1f2022"
    point_color = "#86a67d"
    text_color = "#ffffff"
    frame_color = "#1f2022"
    triangle_edge_color = "#e07a5f"
    triangle_face_color = "none"

    # Создание фигуры и осей
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')

    # Установка пределов осей
    x_min, x_max = np.min(points[:, 0]), np.max(points[:, 0])
    y_min, y_max = np.min(points[:, 1]), np.max(points[:, 1])
    ax.set_xlim(x_min - 0.1, x_max + 0.1)
    ax.set_ylim(y_min - 0.1, y_max + 0.1)

    # Настройка сетки и осей
    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Инициализация графических объектов
    scat = ax.scatter([], [], s=point_size, color=point_color)  # Анимация точек
    drawn_triangles = []  # Список для хранения треугольников

    total_frames = len(points) + len(triangles)  # Всего кадров

    def animate(i):
        if i < len(points):  # Анимация точек
            scat.set_offsets(points[:i+1])
        else:  # Анимация треугольников
            t_index = i - len(points)
            if t_index < len(triangles):
                triangle = patches.Polygon(triangles[t_index], closed=True,
                                           edgecolor=triangle_edge_color,
                                           facecolor=triangle_face_color,
                                           linewidth=2)
                ax.add_patch(triangle)
                drawn_triangles.append(triangle)

        return [scat] + drawn_triangles

    ani = FuncAnimation(
        fig=fig,
        func=animate,
        frames=total_frames,
        interval=100,  # Время между кадрами (можно подстроить)
        blit=False,
        repeat=False
    )

    ani.save('combined_animation.gif', writer='pillow', dpi=100)
    plt.close()

def main():
    points = read_points('knots.txt')
    triangles = read_triangles('triangles.txt')
    plot_animation(points, triangles)

if __name__ == '__main__':
    main()
