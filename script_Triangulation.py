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
            triangle = [(coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])]
            triangles.append(triangle)
    return triangles

def plot_triangles(points, triangles, point_size=100, grid_width=0.0, font_size=20):
    # Цветовая схема
    background_color = "#ffffff"
    grid_color = "#000000"
    point_color = "#000000"
    text_color = "#000000"
    frame_color = "#000000"
    triangle_edge_color = "#000000"
    triangle_face_color = "none"
    
    # background_color = "#ffffff"
    # grid_color = "#1f2022"
    # point_color = "#100f14"
    # text_color = "#100f14"
    # frame_color = "#100f14"
    # triangle_edge_color = "#100f14"
    # triangle_face_color = "none"

    # Создание фигуры
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')

    x_min, x_max = np.min(points[:, 0]), np.max(points[:, 0])
    y_min, y_max = np.min(points[:, 1]), np.max(points[:, 1])
    ax.set_xlim(x_min - 0.1, x_max + 0.1)
    ax.set_ylim(y_min - 0.1, y_max + 0.1)

    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    ax.scatter(points[:, 0], points[:, 1], s=point_size, color=point_color)

    drawn_triangles = []

    # def animate(i):
    #     if i < len(triangles):
    #         coords = triangles[i]
    #         triangle = patches.Polygon(coords, closed=True,
    #                                    edgecolor=triangle_edge_color,
    #                                    facecolor=triangle_face_color,
    #                                    linewidth=2)
    #         ax.add_patch(triangle)
    #         drawn_triangles.append(triangle)
    #     return drawn_triangles

    # ani = FuncAnimation(
    #     fig=fig,
    #     func=animate,
    #     frames=len(triangles),
    #     interval=0.001,
    #     blit=False,
    #     repeat=False
    # )

    # ani.save('triangulation.gif', writer='pillow', dpi=100)
    # plt.close()

    # 🔽 Финальный кадр в PNG
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    ax.set_aspect('equal')
    ax.set_xlim(x_min - 0.1, x_max + 0.1)
    ax.set_ylim(y_min - 0.1, y_max + 0.1)

    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    ax.scatter(points[:, 0], points[:, 1], s=point_size, color=point_color)

    for triangle in triangles:
        poly = patches.Polygon(triangle, closed=True,
                               edgecolor=triangle_edge_color,
                               facecolor=triangle_face_color,
                               linewidth=2)
        ax.add_patch(poly)

    plt.savefig("triangulation.png", dpi=300)
    plt.close()

def main():
    points = read_points('knots.txt')
    triangles = read_triangles('triangles.txt')
    plot_triangles(points, triangles)

if __name__ == '__main__':
    main()
