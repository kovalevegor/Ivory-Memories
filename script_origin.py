import matplotlib.pyplot as plt
import numpy as np

def read_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split())
            points.append((x, y))
    return np.array(points)

def plot_triangulation(points, point_size=20, line_width=1.5, grid_width=0.5, font_size=12):
    # Цветовая схема
    background_color = "#131416"
    grid_color = "#1f2022"
    point_color = "#86a67d"
    line_color = "#f65447"
    text_color = "#ffffff"
    frame_color = "#1f2022"  # Цвет рамки графика
    
    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor(background_color)  # Устанавливаем цвет всего фона
    ax.set_facecolor(background_color)  # Цвет области графика
    ax.set_aspect('equal')
    
    # Устанавливаем цвет рамки графика
    for spine in ax.spines.values():
        spine.set_edgecolor(frame_color)
        spine.set_linewidth(1.5)
    
    # Отображение точек
    ax.scatter(points[:, 0], points[:, 1], color=point_color, s=point_size, label='Точки')
    
    # Настройки сетки
    ax.grid(True, linewidth=grid_width, color=grid_color)
    
    # Настройки подписей
    # ax.set_title('Триангуляция Делоне', fontsize=font_size, color=text_color)
    # ax.set_xlabel('X', fontsize=font_size, color=text_color)
    # ax.set_ylabel('Y', fontsize=font_size, color=text_color)
    ax.tick_params(axis='x', colors=text_color, labelsize=font_size)
    ax.tick_params(axis='y', colors=text_color, labelsize=font_size)
    # ax.legend()
    
    plt.show()

def main():
    points = read_points('knots.txt')
    plot_triangulation(points, point_size=30, line_width=2, grid_width=0.8, font_size=14)
    
if __name__ == '__main__':
    main()
