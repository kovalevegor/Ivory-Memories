# -*- coding: utf-8 -*-
"""
Created on Sat Apr  5 21:32:00 2025

@author: gogak
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from scipy.spatial import ConvexHull
import re

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
            segment = [(coords[0], coords[1]), (coords[2], coords[3])]
            segments.append(segment)
    return segments

def read_clusters(filename):
    clusters = {}
    with open(filename, 'r') as file:
        for line in file:
            # Обновленное регулярное выражение для дробных чисел
            match = re.match(r"Knot \(([0-9]*\.?[0-9]+), ([0-9]*\.?[0-9]+)\) -> Cluster (\d+)", line.strip())
            if match:
                x, y, cluster = float(match.group(1)), float(match.group(2)), int(match.group(3))
                clusters[(x, y)] = cluster
                print(f"Read: Knot ({x}, {y}) -> Cluster {cluster}")
    print(f"Total clusters read: {len(clusters)}")
    return clusters

def plot_animation(knots, segments, spawn_exit_knots=None, asylum_knot=None, 
                   clusters=None, knot_size=60, grid_width=0.1, font_size=10):
    # Цвета
    background_color = "#0e141a"
    grid_color = "#cceb5b"
    knot_color = "#ffffff"
    spawn_exit_knots_color = "#f8c818"
    asylum_knot_color = "#cceb5b"
    segment_color = "#d6d0d0"
    text_color = "#cceb5b"
    frame_color = "#0e141a"
    
    # background_color = "#ffffff"
    # grid_color = "#000000"
    # knot_color = "#000000"
    # text_color = "#000000"
    # frame_color = "#000000"
    # segment_color = "#000000"
    # spawn_exit_knots_color = "#ffffff"
    # asylum_knot_color = "#ffffff"

    # Цвета и альфа для оболочек сообществ
    cluster_colors = [
        '#f75049',  # Зарезервированный цвет 1
        '#1ded83',  # Зарезервированный цвет 2
        '#2570d4',  # Зарезервированный цвет 3
        '#fb932e',  # Зарезервированный цвет 4
        '#9d2bf5',  # Зарезервированный цвет 5
        '#ff4d4d',  # Запасной цвет 1
        '#4dff4d',  # Запасной цвет 2
        '#4d4dff',  # Запасной цвет 3
        '#ff4dff',  # Запасной цвет 4
        '#4dffff',  # Запасной цвет 5
        '#ffff4d'   # Запасной цвет 6
    ]
    cluster_alpha = 0.25  # Прозрачность оболочек

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
    fig_height = 10 * height / width if width > 0 else 10

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.1
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding, y_max + y_padding)

    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Отображаем специальные точки
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
    #                                  color=segment_color, linewidth=2, zorder=10)
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
    x_padding = (x_max - x_min) * 0.1
    y_padding = (y_max - y_min) * 0.1
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

    # Отрисовка выпуклых оболочек для сообществ
    if clusters:
        # Группируем точки по сообществам
        cluster_points = {}
        for knot in knots:
            matched = False
            for (x, y), cluster in clusters.items():
                if np.isclose(knot[0], x, atol=1e-4) and np.isclose(knot[1], y, atol=1e-4):
                    if cluster not in cluster_points:
                        cluster_points[cluster] = []
                    cluster_points[cluster].append(knot)
                    matched = True
                    break  # нашли совпадение — дальше не ищем
            if not matched:
                # Если точка ни в какой кластер не попала — можно залогировать
                pass  # или print(f"Point {knot} not matched to any cluster")
        
        print(f"Clusters found: {len(cluster_points)}")
        for cluster, points in cluster_points.items():
            print(f"Cluster {cluster}: {len(points)} points")

        # Рисуем выпуклые оболочки
        for cluster, points in cluster_points.items():
            points = np.array(points)
            if len(points) >= 3:  # Для построения оболочки нужно минимум 3 точки
                hull = ConvexHull(points)
                color = cluster_colors[cluster % len(cluster_colors)]
                # Заполненная область
                plt.fill(points[hull.vertices, 0], points[hull.vertices, 1],
                         color=color, alpha=cluster_alpha, zorder=31)
                # Контур
                hull_points = np.vstack([points[hull.vertices], points[hull.vertices[0]]])
                plt.plot(hull_points[:, 0], hull_points[:, 1], color=color, linewidth=1, zorder=30)

    # Отрисовка всех точек
    ax.scatter(knots[:, 0], knots[:, 1], s=knot_size, color=knot_color, zorder=20)
    if spawn_exit_knots is not None:
        ax.scatter(spawn_exit_knots[:, 0], spawn_exit_knots[:, 1], s=knot_size * 1.5,
                   color=spawn_exit_knots_color, zorder=21)
    if asylum_knot is not None and len(asylum_knot) > 0:
        ax.scatter(asylum_knot[:, 0], asylum_knot[:, 1],
                   s=knot_size * 2, color=asylum_knot_color, marker='D', zorder=22)

    plt.savefig("clusters.png", dpi=300, bbox_inches='tight', pad_inches=0.5)
    plt.close()

def main():
    spawn_exit_knots = read_knots('spawn_exit_knots.txt')
    asylum_knot = read_knots('asylum.txt')
    knots = read_knots('knots.txt')
    segments = read_segments('maze.txt')
    clusters = read_clusters('clustering.txt')
    plot_animation(knots, segments, spawn_exit_knots, asylum_knot, clusters)

if __name__ == '__main__':
    main()