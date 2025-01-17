"""Module for plasma (3d solver) initialization routines."""
import numpy as np

from ..config.config import Config
from .data import Fields, Particles, Currents, Const_Arrays
from .weights import initial_deposition

ELECTRON_CHARGE = -1
ELECTRON_MASS = 1


# Solving Laplace equation with Dirichlet boundary conditions (Ez) #

def dirichlet_matrix(grid_steps, grid_step_size):
    """
    Calculate a magical matrix that solves the Laplace equation
    if you elementwise-multiply the RHS by it "in DST-space".
    See Samarskiy-Nikolaev, p. 187.
    """
    # mul[i, j] = 1 / (lam[i] + lam[j])
    # lam[k] = 4 / h**2 * sin(k * pi * h / (2 * L))**2, where L = h * (N - 1)
    k = np.arange(1, grid_steps - 1)
    lam = 4 / grid_step_size**2 * np.sin(k * np.pi / (2 * (grid_steps - 1)))**2
    lambda_i, lambda_j = lam[:, None], lam[None, :]
    mul = 1 / (lambda_i + lambda_j)
    return mul / (2 * (grid_steps - 1))**2  # additional 2xDST normalization


# Solving Laplace or Helmholtz equation with mixed boundary conditions #

def mixed_matrix(grid_steps, grid_step_size, subtraction_trick):
    """
    Calculate a magical matrix that solves the Helmholtz or Laplace equation
    (subtraction_trick=True and subtraction_trick=False correspondingly)
    if you elementwise-multiply the RHS by it "in DST-DCT-transformed-space".
    See Samarskiy-Nikolaev, p. 189 and around.
    """
    # mul[i, j] = 1 / (lam[i] + lam[j])
    # lam[k] = 4 / h**2 * sin(k * pi * h / (2 * L))**2, where L = h * (N - 1)
    # but k for lam_i spans from 1..N-2, while k for lam_j covers 0..N-1
    ki, kj = np.arange(1, grid_steps - 1), np.arange(grid_steps)
    li = 4 / grid_step_size**2 * np.sin(ki * np.pi / (2 * (grid_steps - 1)))**2
    lj = 4 / grid_step_size**2 * np.sin(kj * np.pi / (2 * (grid_steps - 1)))**2
    lambda_i, lambda_j = li[:, None], lj[None, :]
    mul = 1 / (lambda_i + lambda_j + (1 if subtraction_trick else 0))
    return mul / (2 * (grid_steps - 1))**2  
    # return additional 2xDST normalization


# Solving Laplace equation with Neumann boundary conditions (Bz) #

def neumann_matrix(grid_steps, grid_step_size):
    """
    Calculate a magical matrix that solves the Laplace equation
    if you elementwise-multiply the RHS by it "in DST-space".
    See Samarskiy-Nikolaev, p. 187.
    """
    # mul[i, j] = 1 / (lam[i] + lam[j])
    # lam[k] = 4 / h**2 * sin(k * pi * h / (2 * L))**2, where L = h * (N - 1)

    k = np.arange(0, grid_steps)
    lam = 4 / grid_step_size**2 * np.sin(k * np.pi / (2 * (grid_steps - 1)))**2
    lambda_i, lambda_j = lam[:, None], lam[None, :]
    mul = 1 / (lambda_i + lambda_j)  # WARNING: zero division in mul[0, 0]!
    mul[0, 0] = 0  # doesn't matter anyway, just defines constant shift
    return mul / (2 * (grid_steps - 1))**2  # additional 2xDST normalization


# Plasma particles initialization #

def make_plasma_grid(steps, step_size, fineness):
    """
    Create initial plasma particles coordinates
    (a single 1D grid for both x and y).
    Avoids positioning particles at the cell edges and boundaries, example:
    `fineness=3`:
        +-----------+-----------+-----------+-----------+
        | .   .   . | .   .   . | .   .   . | .   .   . |
        |           |           |           |           |
        | .   .   . | .   .   . | .   .   . | .   .   . |
        |           |           |           |           |
        | .   .   . | .   .   . | .   .   . | .   .   . |
        +-----------+-----------+-----------+-----------+
    `fineness=2`:
        +-------+-------+-------+-------+-------+
        | .   . | .   . | .   . | .   . | .   . |
        |       |       |       |       |       |
        | .   . | .   . | .   . | .   . | .   . |
        +-------+-------+-------+-------+-------+
    """
    plasma_step = step_size / fineness
    if fineness % 2:  # some on zero axes, none on cell corners
        right_half = np.arange(steps // 2 * fineness) * plasma_step
        left_half = -right_half[:0:-1]  # invert, reverse, drop zero
    else:  # none on zero axes, none on cell corners
        right_half = (.5 + np.arange(steps // 2 * fineness)) * plasma_step
        left_half = -right_half[::-1]  # invert, reverse
    plasma_grid = np.concatenate([left_half, right_half])
    assert(np.array_equal(plasma_grid, -plasma_grid[::-1]))
    return plasma_grid


def make_plasma(steps, cell_size, fineness=2):
    """
    Initialize default plasma state, fineness**2 particles per cell.
    """
    pl_grid = make_plasma_grid(steps, cell_size, fineness)

    Np = len(pl_grid)

    y_init = np.broadcast_to(pl_grid, (Np, Np))
    x_init = y_init.T

    x_offt = np.zeros((Np, Np))
    y_offt = np.zeros((Np, Np))
    px = np.zeros((Np, Np))
    py = np.zeros((Np, Np))
    pz = np.zeros((Np, Np))
    q = np.ones((Np, Np)) * ELECTRON_CHARGE / fineness**2
    m = np.ones((Np, Np)) * ELECTRON_MASS / fineness**2

    return x_init, y_init, x_offt, y_offt, px, py, pz, q, m


def init_plasma(config: Config):
    """
    Initialize all the arrays needed (for what?).
    """
    grid_steps            = config.getint('window-width-steps')
    grid_step_size        = config.getfloat('window-width-step-size')
    reflect_padding_steps = config.getint('reflect-padding-steps')
    plasma_padding_steps  = config.getint('plasma-padding-steps')
    plasma_fineness       = config.getint('plasma-fineness')
    solver_trick          = config.getint('field-solver-subtraction-trick')

    # for convenient diagnostics, a cell should be in the center of the grid
    assert grid_steps % 2 == 1

    # particles should not reach the window pre-boundary cells
    if reflect_padding_steps <= 2:
        raise Exception("'reflect_padding_steps' parameter is too low.\n" +
                        "Details: 'reflect_padding_steps' must be bigger than" +
                        " 2. By default it is 5.")
    
    x_init, y_init, x_offt, y_offt, px, py, pz, q, m = \
        make_plasma(grid_steps - plasma_padding_steps * 2,
                    grid_step_size,
                    fineness=plasma_fineness)

    def zeros():
        return np.zeros((grid_steps, grid_steps), dtype=np.float64)

    ro_initial = initial_deposition(grid_steps, grid_step_size,
                                    x_init, y_init, x_offt, y_offt,
                                    px, py, pz, m, q)
    dir_matrix = dirichlet_matrix(grid_steps, grid_step_size)
    mix_matrix = mixed_matrix(grid_steps, grid_step_size, solver_trick)
    neu_matrix = neumann_matrix(grid_steps, grid_step_size)

    fields = Fields(grid_steps)
    particles = Particles(x_init, y_init, x_offt, y_offt, px, py, pz, q, m)
    currents = Currents(zeros(), zeros(), zeros(), zeros())
    const_arrays = Const_Arrays(ro_initial, dir_matrix, mix_matrix, neu_matrix)

    return fields, particles, currents, const_arrays


def load_plasma(config: Config, path_to_plasmastate: str):
    fields, particles, currents, const_arrays = init_plasma(config)

    with np.load(file=path_to_plasmastate) as state:
        fields.Ex, fields.Ey, fields.Ez =\
             state['Ex'], state['Ey'], state['Ez']
        fields.Bx, fields.By, fields.Bz =\
             state['Bx'], state['By'], state['Bz']
        fields.Phi = state['Phi']

        particles.x_offt, particles.x_offt =\
            state['x_offt'], state['y_offt']
        particles.px, particles.py, particles.pz =\
            state['px'], state['py'], state['pz']

        currents.ro, currents.jx, currents.jy, currents.jz = \
            state['ro'], state['jx'], state['jy'], state['jz']

    return fields, particles, currents, const_arrays
