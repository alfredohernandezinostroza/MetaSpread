import numpy as np
import matplotlib.pyplot as plt
#import random

# Returns an array with 1 in the region that contains the N centre-most grid points and 0 on the others
def find_quasi_circle(N, grid_width, grid_height):
    '''
    Returns the coordinate of the centre-most points and also an
    array with 1 in that region and 0 on the rest
    '''

    current_grid = np.zeros((grid_height, grid_width))
    less_than_ideal_radius = np.sqrt(N/np.pi) / 1.2   
    n_cells_current = 0

    # Center of the grid in 0 based indexing
    a = grid_width/2 - 0.5
    b = grid_height/2 - 0.5

    # Increases the radius of the circle until if finds the correct size
    for radius in np.arange(less_than_ideal_radius, (min(grid_width, grid_height)/2), 0.1):
        last_grid = current_grid.copy()
        n_cells_last = n_cells_current

        # fill the circle
        for x in range(grid_width):
            for y in range(grid_height):
                if (x-a)**2 + (y-b)**2 < radius**2:
                    current_grid[x,y] = 1

        n_cells_current = np.count_nonzero(current_grid == 1)

        if n_cells_current >= N:

            if n_cells_current == N:
                last_grid = current_grid

            if n_cells_current > N:
                n_cells_to_add = N - n_cells_last
                circle_border = np.ma.masked_where(last_grid == 1, current_grid)

                # Fill the outer cells in a diametrically opposit counter-clocksise way
                outer_p = np.where(circle_border == 1)
                coords_outer_p = list(zip(outer_p[0], outer_p[1]))
                sign = 1
                for i in range(n_cells_to_add):
                    coord_to_add = coords_outer_p[i*sign]
                    last_grid[coord_to_add] = 1
                    sign *= -1

            # Create a list with possible places
            possible_places = np.where(last_grid == 1)
            coords = [list(tup) for tup in zip(possible_places[0], possible_places[1], (0 for i in possible_places[0]))]
            return last_grid, coords


def find_quasi_sphere(N, grid_width, grid_height, grid_depth):
    '''
    3D analogue of find_quasi_circle. Returns a boolean-ish mask with 1 in the
    region containing the N centre-most grid points, and the list of those
    points as [x, y, z, 0] (the trailing 0 is an occupancy counter used by the
    seeding loop, mirroring the 2D helper).
    '''
    a = grid_width / 2 - 0.5
    b = grid_height / 2 - 0.5
    c = grid_depth / 2 - 0.5
    xs, ys, zs = np.indices((grid_width, grid_height, grid_depth))
    dist2 = (xs - a) ** 2 + (ys - b) ** 2 + (zs - c) ** 2

    ideal_radius = (N * 3 / (4 * np.pi)) ** (1 / 3) / 1.2
    grid = np.zeros((grid_width, grid_height, grid_depth))
    n_current = 0
    for radius in np.arange(ideal_radius, min(grid_width, grid_height, grid_depth) / 2, 0.1):
        last_grid = grid.copy()
        n_last = n_current
        grid = (dist2 < radius ** 2).astype(float)
        n_current = int(grid.sum())
        if n_current >= N:
            if n_current == N:
                last_grid = grid.copy()
            elif n_current > N:
                n_to_add = N - n_last
                border = np.argwhere((grid == 1) & (last_grid == 0))
                sign = 1
                for i in range(n_to_add):
                    coord = border[(i * sign) % len(border)]
                    last_grid[tuple(coord)] = 1
                    sign *= -1
            places = np.argwhere(last_grid == 1)
            coords = [[int(p[0]), int(p[1]), int(p[2]), 0] for p in places]
            return last_grid, coords

    places = np.argwhere(grid == 1)
    coords = [[int(p[0]), int(p[1]), int(p[2]), 0] for p in places]
    return grid, coords