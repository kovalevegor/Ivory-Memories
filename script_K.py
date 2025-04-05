import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

def read_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split())
            points.append((x, y))
    return np.array(points)

def plot_knots(points, point_size=30, grid_width=0.8, font_size=14):
    # Настройка цветовой схемы
    background_color = "#131416"
    grid_color = "#1f2022"
    point_color = "#86a67d"
    text_color = "#ffffff"
    frame_color = "#1f2022"

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

    # Настройка рамки и сетки
    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    ax.grid(True, color=grid_color, linewidth=grid_width)
    
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)

    # Инициализация пустого графика
    scat = ax.scatter([], [], s=point_size, color=point_color)
    
    # Функция анимации
    def animate(i):
        scat.set_offsets(points[:i+1])  # Добавляем точки по одной
        return (scat,)

    # Создание анимации
    ani = FuncAnimation(
        fig=fig,
        func=animate,
        frames=len(points),
        interval=0.0001,  # мс между точками
        blit=True,
        repeat=False
    )

    # Сохранение в GIF
    ani.save('knots.gif', writer='pillow', dpi=100)
    plt.close()

def main():
    points = read_points('knots.txt')
    # plot_triangulation(points)
    plot_knots(points)

if __name__ == '__main__':
    main()