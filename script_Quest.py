# -*- coding: utf-8 -*-
"""
Created on Tue May 13 20:53:20 2025

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
            match = re.match(r"Knot \(([0-9]*\.?[0-9]+), ([0-9]*\.?[0-9]+)\) -> Cluster (\d+)", line.strip())
            if match:
                x, y, cluster = float(match.group(1)), float(match.group(2)), int(match.group(3))
                clusters[(x, y)] = cluster
    return clusters

# Новая функция для чтения заданий квестов
def read_quest_assignments(filename):
    quests = []
    with open(filename, 'r') as file:
        for line in file:
            match = re.match(r"Tag: (.+), Node: (\d+) \((.+), (.+)\)", line.strip())
            if match:
                tag = match.group(1)
                node = int(match.group(2))
                x = float(match.group(3))
                y = float(match.group(4))
                quests.append({'coords': (x, y), 'tag': tag, 'node': node})
    return quests

def plot_animation(knots, segments, spawn_exit_knots=None, asylum_knot=None, 
                   clusters=None, quest_assignments=None, knot_size=60, grid_width=0.1, font_size=10):
    background_color = "#0e141a"
    grid_color = "#cceb5b"
    knot_color = "#ffffff"
    spawn_exit_knots_color = "#f8c818"
    asylum_knot_color = "#cceb5b"
    segment_color = "#d6d0d0"
    text_color = "#cceb5b"
    frame_color = "#0e141a"
    quest_color = "#f8c818"
    
    cluster_colors = [
        '#f75049', '#1ded83', '#2570d4', 
        '#fb932e', '#9d2bf5', '#ff4d4d',
        '#4dff4d', '#4d4dff', '#ff4dff',
        '#4dffff', '#ffff4d'
    ]
    cluster_alpha = 0.25

    # Обновляем границы с учетом точек квестов
    all_x = knots[:, 0]
    all_y = knots[:, 1]
    
    additional_points = []
    if spawn_exit_knots is not None:
        additional_points.append(spawn_exit_knots)
    if asylum_knot is not None and len(asylum_knot) > 0:
        additional_points.append(asylum_knot)
    if quest_assignments is not None:
        quest_coords = np.array([q['coords'] for q in quest_assignments])
        additional_points.append(quest_coords)
    
    for points in additional_points:
        all_x = np.concatenate([all_x, points[:, 0]])
        all_y = np.concatenate([all_y, points[:, 1]])

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

    # Отрисовка сегментов
    for start, end in segments:
        ax.plot([start[0], end[0]], [start[1], end[1]], color=segment_color, linewidth=2, zorder=10)

    # Отрисовка кластеров (остается без изменений)
    if clusters:
        cluster_points = {}
        for knot in knots:
            matched = False
            for (x, y), cluster in clusters.items():
                if np.isclose(knot[0], x, atol=1e-4) and np.isclose(knot[1], y, atol=1e-4):
                    if cluster not in cluster_points:
                        cluster_points[cluster] = []
                    cluster_points[cluster].append(knot)
                    matched = True
                    break
            if not matched:
                pass
        
        for cluster, points in cluster_points.items():
            points = np.array(points)
            if len(points) >= 3:
                hull = ConvexHull(points)
                color = cluster_colors[cluster % len(cluster_colors)]
                plt.fill(points[hull.vertices, 0], points[hull.vertices, 1],
                          color=color, alpha=cluster_alpha, zorder=31)
                hull_points = np.vstack([points[hull.vertices], points[hull.vertices[0]]])
                plt.plot(hull_points[:, 0], hull_points[:, 1], color=color, linewidth=1, zorder=25)

    # Отрисовка специальных точек
    ax.scatter(knots[:, 0], knots[:, 1], s=knot_size, color=knot_color, zorder=20)
    if spawn_exit_knots is not None:
        ax.scatter(spawn_exit_knots[:, 0], spawn_exit_knots[:, 1], s=knot_size * 1.5,
                   color=spawn_exit_knots_color, zorder=21)
    if asylum_knot is not None and len(asylum_knot) > 0:
        ax.scatter(asylum_knot[:, 0], asylum_knot[:, 1],
                   s=knot_size * 2, color=asylum_knot_color, marker='D', zorder=22)

    # Новая секция: отрисовка точек квестов и подписей
    if quest_assignments is not None:
        quest_coords = np.array([q['coords'] for q in quest_assignments])
        ax.scatter(quest_coords[:, 0], quest_coords[:, 1], 
                   s=knot_size*1.2, color=quest_color, zorder=60)
        
        for quest in quest_assignments:
            x, y = quest['coords']
            label = f"{quest['tag']}"
            ax.annotate(
                label,
                (x, y),
                textcoords="offset points",
                xytext=(10, 5),
                ha='left',
                fontsize=font_size,
                color=quest_color,
                weight='bold', 
                zorder=100
            )

    plt.savefig("quests.png", dpi=300, bbox_inches='tight', pad_inches=0.5)
    plt.close()

def main():
    spawn_exit_knots = read_knots('spawn_exit_knots.txt')
    asylum_knot = read_knots('asylum.txt')
    knots = read_knots('knots.txt')
    segments = read_segments('maze.txt')
    clusters = read_clusters('clustering.txt')
    quest_assignments = read_quest_assignments('quest_assignments.txt')  # Читаем новые данные
    
    plot_animation(
        knots=knots,
        segments=segments,
        spawn_exit_knots=spawn_exit_knots,
        asylum_knot=asylum_knot,
        clusters=clusters,
        quest_assignments=quest_assignments  # Передаем в функцию визуализации
    )

if __name__ == '__main__':
    main()