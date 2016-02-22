# -*- coding: utf-8 -*-
"""This module contains some setting functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['setAxesEqual']


def setAxesEqual(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    from numpy import mean

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]
    x_mean = mean(x_limits)
    y_range = y_limits[1] - y_limits[0]
    y_mean = mean(y_limits)
    z_range = z_limits[1] - z_limits[0]
    z_mean = mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])

    return None
