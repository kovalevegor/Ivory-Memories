# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 21:15:19 2025

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

def read_segments(filename):
    segments = []
    with open(filename, 'r') as file:
        for line in file:
            coords = list(map(float, line.strip().split()))
            # Каждая строка должна содержать координаты двух точек: (x1,y1), (x2,y2)
            segment = [(coords[0], coords[1]), (coords[2], coords[3])]
            segments.append(segment)
    return segments

def plot_segments(points, segments, point_size=100, grid_width=0.5, font_size=20):
    # Настройка цветовой схемы
    background_color = "#ffffff"
    grid_color = "#000000"
    point_color = "#000000"
    text_color = "#000000"
    frame_color = "#000000"
    segment_color = "#000000"
    

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

    # Список для хранения добавленных патчей (отрезков)
    drawn_segments = []

    # def animate(i):
    #     if i < len(segments):
    #         start, end = segments[i]
    #         segment = plt.Line2D([start[0], end[0]], [start[1], end[1]], color=segment_color, linewidth=2)
    #         ax.add_line(segment)
    #         drawn_segments.append(segment)
    #     return drawn_segments

    # # Сохраняем анимацию
    # ani = FuncAnimation(
    #     fig=fig,
    #     func=animate,
    #     frames=len(segments),
    #     interval=0.001,
    #     blit=False,
    #     repeat=False
    # )

    # ani.save('polygon.gif', writer='pillow', dpi=100)

    # 🔽 Создаём финальное изображение
    # Заново рисуем всё с нуля, чтобы результат был чистым
    fig2, ax2 = plt.subplots(figsize=(10, 10))
    fig2.patch.set_facecolor(background_color)
    ax2.set_facecolor(background_color)
    ax2.set_aspect('equal')
    ax2.set_xlim(x_min - 0.1, x_max + 0.1)
    ax2.set_ylim(y_min - 0.1, y_max + 0.1)

    for spine in ax2.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax2.grid(True, color=grid_color, linewidth=grid_width)
    ax2.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax2.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Рисуем финальные точки и все сегменты
    ax2.scatter(points[:, 0], points[:, 1], s=point_size, color=point_color)
    for start, end in segments:
        ax2.plot([start[0], end[0]], [start[1], end[1]], color=segment_color, linewidth=2)

    # Сохраняем результат
    plt.savefig("polygon.png", dpi=300)
    plt.close(fig2)


def main():
    points = read_points('knots.txt')
    segments = read_segments('maze.txt')
    plot_segments(points, segments)

if __name__ == '__main__':
    main()
