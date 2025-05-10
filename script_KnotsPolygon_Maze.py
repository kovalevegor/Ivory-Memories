# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 21:32:00 2025

@author: gogak
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

def read_knots(filename):
    knots = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split())
            knots.append((x, y))
    return np.array(knots)

def read_segments(filename):
    segments = []
    with open(filename, 'r') as file:
        for line in file:
            coords = list(map(float, line.strip().split()))
            # Каждая строка должна содержать координаты двух точек: (x1,y1), (x2,y2)
            segment = [(coords[0], coords[1]), (coords[2], coords[3])]
            segments.append(segment)
    return segments

def plot_animation(knots, segments, spawn_exit_knots=None, asylum_knot=None, knot_size=100, grid_width=0.5, font_size=20):
    # Цвета
    background_color = "#0e141a"
    grid_color = "#cceb5b"
    knot_color = "#ffffff"
    spawn_exit_knots_color = "#f8c818"
    asylum_knot_color = "#cceb5b"
    segment_color = "#d6d0d0"
    text_color = "#cceb5b"
    frame_color = "#0e141a"

    # Границы графика
    all_x = knots[:, 0]
    all_y = knots[:, 1]
    if spawn_exit_knots is not None:
        all_x = np.concatenate([all_x, spawn_exit_knots[:, 0]])
        all_y = np.concatenate([all_y, spawn_exit_knots[:, 1]])
    if asylum_knot is not None and len(asylum_knot) > 0:
        all_x = np.concatenate([all_x, asylum_knot[:, 0]])
        all_y = np.concatenate([all_y, asylum_knot[:, 1]])

    x_min, x_max = np.min(all_x), np.max(all_x)
    y_min, y_max = np.min(all_y), np.max(all_y)
    width, height = x_max - x_min, y_max - y_min

    fig_width = 10
    fig_height = 10 * height / width if width > 0 else 10  # сохранение пропорций

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')
    # ax.set_xlim(x_min - 0.1, x_max + 0.1)
    # ax.set_ylim(y_min - 0.1, y_max + 0.1)
    x_padding = (x_max - x_min) * 0.1  # 10% от ширины
    y_padding = (y_max - y_min) * 0.1  # 10% от высоты
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding, y_max + y_padding)

    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Отображаем специальные точки заранее
    if spawn_exit_knots is not None:
        ax.scatter(spawn_exit_knots[:, 0], spawn_exit_knots[:, 1], s=knot_size * 1.5,
                   color=spawn_exit_knots_color, zorder=23)
    if asylum_knot is not None and len(asylum_knot) > 0:
        ax.scatter(asylum_knot[:, 0], asylum_knot[:, 1],
                   s=knot_size * 2, color=asylum_knot_color, marker='D', zorder=21)

    scat = ax.scatter([], [], s=knot_size, color=knot_color, zorder=20)
    drawn_segments = []
    total_frames = len(knots) + len(segments)

    # def animate(i):
    #     if i < len(knots):
    #         scat.set_offsets(knots[:i+1])
    #     else:
    #         s_index = i - len(knots)
    #         if s_index < len(segments):
    #             start, end = segments[s_index]
    #             segment = plt.Line2D([start[0], end[0]], [start[1], end[1]],
    #                                   color=segment_color, linewidth=2, zorder=10)
    #             ax.add_line(segment)
    #             drawn_segments.append(segment)
    #     return [scat] + drawn_segments

    # ani = FuncAnimation(
    #     fig=fig,
    #     func=animate,
    #     frames=total_frames,
    #     interval=100,
    #     blit=False,
    #     repeat=False
    # )
    # ani.save('maze.gif', writer='pillow', dpi=100)

    # Финальный рендер
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')
    # ax.set_xlim(x_min - 0.1, x_max + 0.1)
    # ax.set_ylim(y_min - 0.1, y_max + 0.1)
    x_padding = (x_max - x_min) * 0.1  # 10% от ширины
    y_padding = (y_max - y_min) * 0.1  # 10% от высоты
    
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding, y_max + y_padding)

    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    for start, end in segments:
        ax.plot([start[0], end[0]], [start[1], end[1]], color=segment_color, linewidth=2, zorder=10)

    ax.scatter(knots[:, 0], knots[:, 1], s=knot_size, color=knot_color, zorder=20)
    if spawn_exit_knots is not None:
        ax.scatter(spawn_exit_knots[:, 0], spawn_exit_knots[:, 1], s=knot_size * 1.5,
                   color=spawn_exit_knots_color, zorder=22)
    if asylum_knot is not None and len(asylum_knot) > 0:
        ax.scatter(asylum_knot[:, 0], asylum_knot[:, 1],
                   s=knot_size * 2, color=asylum_knot_color, marker='D', zorder=21)

    plt.savefig("maze.png", dpi=300, bbox_inches='tight', pad_inches=0.5)
    plt.close()



def main():
    spawn_exit_knots = read_knots('spawn_exit_knots.txt')
    asylum_knot = read_knots('asylum.txt')
    knots = read_knots('knots.txt')
    segments = read_segments('maze.txt') 
    plot_animation(knots, segments, spawn_exit_knots, asylum_knot)

if __name__ == '__main__':
    main()
