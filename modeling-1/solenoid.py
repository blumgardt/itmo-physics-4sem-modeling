"""
Моделирование магнитного поля цилиндрического соленоида.
Принцип суперпозиции: поле = сумма полей N круговых витков с током.
Поле каждого витка — закон Био-Савара через эллиптические интегралы K(k), E(k).

Интерактивное окно: слайдеры D, L, N, I — пересчёт и перерисовка в реальном времени.
"""

import numpy as np
from scipy.special import ellipk, ellipe
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

MU0 = 4 * np.pi * 1e-7  # Тл·м/А


def single_loop_field(R, I, r, z):
    """
    Поле одного кругового витка радиуса R с током I
    в точках (r, z) (цилиндрические координаты относительно центра витка).
    Возвращает (Br, Bz).
    """
    Br = np.zeros_like(r)
    Bz = np.zeros_like(r)

    on_axis = np.abs(r) < 1e-10
    off_axis = ~on_axis

    if np.any(on_axis):
        z_ax = z[on_axis]
        Bz[on_axis] = MU0 * I * R**2 / (2 * (R**2 + z_ax**2)**1.5)

    if np.any(off_axis):
        r_o = r[off_axis]
        z_o = z[off_axis]
        alpha2 = (R + r_o)**2 + z_o**2
        beta2 = (R - r_o)**2 + z_o**2
        k2 = np.clip(4 * R * r_o / alpha2, 0, 1 - 1e-12)
        K = ellipk(k2)
        E = ellipe(k2)
        alpha = np.sqrt(alpha2)
        c = MU0 * I / (2 * np.pi)
        Bz[off_axis] = c / alpha * (K + (R**2 - r_o**2 - z_o**2) / beta2 * E)
        Br[off_axis] = c * z_o / (r_o * alpha) * (-K + (R**2 + r_o**2 + z_o**2) / beta2 * E)

    return Br, Bz


def solenoid_field(D, L, N, I, r_grid, z_grid):
    """Суперпозиция полей N витков, расположенных от -L/2 до +L/2."""
    R = D / 2
    positions = np.linspace(-L / 2, L / 2, N)
    Br_total = np.zeros_like(r_grid)
    Bz_total = np.zeros_like(z_grid)
    for z_coil in positions:
        dBr, dBz = single_loop_field(R, I, r_grid, z_grid - z_coil)
        Br_total += dBr
        Bz_total += dBz
    return Br_total, Bz_total


def compute(D, L, N, I):
    """Вычислить поле на фиксированной сетке, вернуть всё нужное для отрисовки."""
    N = max(int(N), 2)
    R = D / 2
    nz, ny = 200, 140
    z_lin = np.linspace(-1.5 * L, 1.5 * L, nz)
    y_lin = np.linspace(-1.5 * D, 1.5 * D, ny)
    Z, Y = np.meshgrid(z_lin, y_lin)
    r_abs = np.abs(Y)

    Br, Bz = solenoid_field(D, L, N, I, r_abs, Z)
    By = Br * np.sign(Y)
    B_mag = np.sqrt(Bz**2 + By**2)
    B_theory = MU0 * (N / L) * I

    return z_lin, y_lin, Z, Y, Bz, By, B_mag, B_theory, R


def main():
    # Начальные параметры
    D0, L0, N0, I0 = 0.1, 0.3, 100, 1.0

    # Первый расчёт
    z_lin, y_lin, Z, Y, Bz, By, B_mag, B_theory, R = compute(D0, L0, N0, I0)

    # --- Окно ---
    fig, ax = plt.subplots(figsize=(13, 8))
    plt.subplots_adjust(bottom=0.32)

    # Цветовая карта
    levels = np.linspace(0, 1.5 * B_theory, 50)
    cf = ax.contourf(Z * 100, Y * 100, B_mag, levels=levels, cmap='inferno', extend='max')
    cbar = fig.colorbar(cf, ax=ax, label='|B|, Тл')

    # Линии поля
    sp = ax.streamplot(
        z_lin * 100, y_lin * 100, Bz, By,
        color='white', linewidth=0.6, density=1.5, arrowsize=0.7,
        broken_streamlines=False,
    )

    # Контур соленоида
    sol_top, = ax.plot(
        [-L0/2*100, L0/2*100], [R*100, R*100], 'c-', linewidth=2.5
    )
    sol_bot, = ax.plot(
        [-L0/2*100, L0/2*100], [-R*100, -R*100], 'c-', linewidth=2.5
    )

    ax.set_xlabel('z, см', fontsize=12)
    ax.set_ylabel('y, см', fontsize=12)
    title = ax.set_title(
        f'Магнитное поле соленоида  |  '
        f'D={D0*100:.1f} см, L={L0*100:.1f} см, N={int(N0)}, I={I0:.1f} А\n'
        f'B(центр) = {B_mag[B_mag.shape[0]//2, B_mag.shape[1]//2]*1e6:.1f} мкТл   '
        f'B_теор(∞) = {B_theory*1e6:.1f} мкТл',
        fontsize=11,
    )
    ax.set_aspect('equal')

    # --- Слайдеры ---
    slider_color = '#e0e0e0'
    ax_D = plt.axes([0.15, 0.20, 0.65, 0.03], facecolor=slider_color)
    ax_L = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor=slider_color)
    ax_N = plt.axes([0.15, 0.10, 0.65, 0.03], facecolor=slider_color)
    ax_I = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=slider_color)

    s_D = Slider(ax_D, 'D, см', 2, 30, valinit=D0 * 100, valstep=0.5)
    s_L = Slider(ax_L, 'L, см', 5, 100, valinit=L0 * 100, valstep=1)
    s_N = Slider(ax_N, 'N витков', 5, 500, valinit=N0, valstep=5)
    s_I = Slider(ax_I, 'I, А', 0.1, 10, valinit=I0, valstep=0.1)

    # Кнопка «Сброс»
    ax_reset = plt.axes([0.85, 0.05, 0.1, 0.04])
    btn_reset = Button(ax_reset, 'Сброс', color=slider_color, hovercolor='#c0c0c0')

    def update(val=None):
        D = s_D.val / 100
        L = s_L.val / 100
        N = int(s_N.val)
        I = s_I.val

        z_lin, y_lin, Z, Y, Bz, By, B_mag, B_theory, R = compute(D, L, N, I)

        ax.clear()

        levels = np.linspace(0, max(1.5 * B_theory, 1e-10), 50)
        ax.contourf(Z * 100, Y * 100, B_mag, levels=levels, cmap='inferno', extend='max')

        ax.streamplot(
            z_lin * 100, y_lin * 100, Bz, By,
            color='white', linewidth=0.6, density=1.5, arrowsize=0.7,
            broken_streamlines=False,
        )

        ax.plot([-L/2*100, L/2*100], [R*100, R*100], 'c-', linewidth=2.5)
        ax.plot([-L/2*100, L/2*100], [-R*100, -R*100], 'c-', linewidth=2.5)

        B_center = B_mag[B_mag.shape[0]//2, B_mag.shape[1]//2]

        ax.set_xlabel('z, см', fontsize=12)
        ax.set_ylabel('y, см', fontsize=12)
        ax.set_title(
            f'Магнитное поле соленоида  |  '
            f'D={D*100:.1f} см, L={L*100:.1f} см, N={N}, I={I:.1f} А\n'
            f'B(центр) = {B_center*1e6:.1f} мкТл   '
            f'B_теор(∞) = {B_theory*1e6:.1f} мкТл',
            fontsize=11,
        )
        ax.set_aspect('equal')
        fig.canvas.draw_idle()

    s_D.on_changed(update)
    s_L.on_changed(update)
    s_N.on_changed(update)
    s_I.on_changed(update)

    def reset(event):
        s_D.reset()
        s_L.reset()
        s_N.reset()
        s_I.reset()

    btn_reset.on_clicked(reset)

    plt.show()


if __name__ == '__main__':
    main()
