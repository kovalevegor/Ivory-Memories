# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 21:32:00 2025

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

def plot_animation(points, segments, point_size=30, grid_width=0.8, font_size=14):
    # Настройки цветовой схемы
    background_color = "#131416"
    grid_color = "#1f2022"
    point_color = "#86a67d"
    text_color = "#ffffff"
    frame_color = "#1f2022"
    segment_color = "#e07a5f"  # Цвет отрезков

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
    drawn_segments = []  # Список для хранения отрезков

    total_frames = len(points) + len(segments)  # Всего кадров

    def animate(i):
        if i < len(points):  # Анимация точек
            scat.set_offsets(points[:i+1])
        else:  # Анимация отрезков
            s_index = i - len(points)
            if s_index < len(segments):
                start, end = segments[s_index]
                segment = plt.Line2D([start[0], end[0]], [start[1], end[1]], color=segment_color, linewidth=2)
                ax.add_line(segment)
                drawn_segments.append(segment)

        return [scat] + drawn_segments

    ani = FuncAnimation(
        fig=fig,
        func=animate,
        frames=total_frames,
        interval=100,  # Время между кадрами (можно подстроить)
        blit=False,
        repeat=False
    )

    ani.save('knots_polygon.gif', writer='pillow', dpi=100)
    plt.close()

def main():
    points = read_points('knots.txt')
    segments = read_segments('polygon.txt')  # Используем новый файл для отрезков
    plot_animation(points, segments)

if __name__ == '__main__':
    main()
