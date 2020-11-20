"""
    This file is part of Avaframe.
"""

# Load modules
import logging
import numpy as np

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def getNormal(x, y, Nx, Ny, Nz, csz):
    # by default bilinear interpolation of the Nx, Ny, Nz of the grid
    nx = geoTrans.projectOnRasterRoot(x, y, Nx, csz=csz)
    ny = geoTrans.projectOnRasterRoot(x, y, Ny, csz=csz)
    nz = geoTrans.projectOnRasterRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalArray(x, y, Nx, Ny, Nz, csz):
    nrow, ncol = np.shape(Nx)
    # by default bilinear interpolation of the Nx, Ny, Nz of the grid
    nx, _ = geoTrans.projectOnRasterVectRoot(x, y, Nx, csz=csz)
    ny, _ = geoTrans.projectOnRasterVectRoot(x, y, Ny, csz=csz)
    nz, _ = geoTrans.projectOnRasterVectRoot(x, y, Nz, csz=csz)
    nx, ny, nz = normalize(nx, ny, nz)
    return nx, ny, nz


def getNormalMesh(z, csz, num=4):
    n, m = np.shape(z)
    Nx = np.ones((n, m))
    Ny = np.ones((n, m))
    Nz = np.ones((n, m))
    # first and last row, first and last column are inacurate
    if num == 4:
        # filling the inside of the matrix
        # normal calculation with 4 triangles
        # (Zl - Zr) / csz
        Nx[1:n-1, 1:m-1] = (z[1:n-1, 0:m-2] - z[1:n-1, 2:m]) / csz
        # (Zd - Zu) * csz
        Ny[1:n-1, 1:m-1] = (z[0:n-2, 1:m-1] - z[2:n, 1:m-1]) / csz
        Nz = 2 * Nz
        # filling the first col of the matrix
        # -2*(Zr - Zp) / csz
        Nx[1:n-1, 0] = - 2*(z[1:n-1, 1] - z[1:n-1, 0]) / csz
        # (Zd - Zu) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0]) / csz
        # filling the last col of the matrix
        # 2*(Zl - Zp) / csz
        Nx[1:n-1, m-1] = 2*(z[1:n-1, m-2] - z[1:n-1, m-1]) / csz
        # (Zd - Zu) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1]) / csz
        # filling the first row of the matrix
        # (Zl - Zr) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m]) / csz
        # -2*(Zu - Zp) / csz
        Ny[0, 1:m-1] = - 2*(z[1, 1:m-1] - z[0, 1:m-1]) / csz
        # filling the last row of the matrix
        # (Zl - Zr) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m]) / csz
        # 2*(Zd - Zp) / csz
        Ny[n-1, 1:m-1] = 2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) / csz
        # filling the corners of the matrix
        Nx[0, 0] = -(z[0, 1] - z[0, 0]) / csz
        Ny[0, 0] = -(z[1, 0] - z[0, 0]) / csz
        Nz[0, 0] = 1
        Nx[n-1, 0] = -(z[n-1, 1] - z[n-1, 0]) / csz
        Ny[n-1, 0] = (z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 1
        Nx[0, m-1] = (z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = -(z[1, m-1] - z[0, m-1]) / csz
        Nz[0, m-1] = 1
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1]) / csz
        Nz[n-1, m-1] = 1

    if num == 6:
        # filling the inside of the matrix
        # normal calculation with 6 triangles
        # (2*(Zl - Zr) - Zur + Zdl + Zu - Zd) / csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            + z[2:n, 1:m-1] - z[0:n-2, 1:m-1]) / csz
        # (2*(Zd - Zu) - Zur + Zdl - Zl + Zr) / csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
                            - z[2:n, 2:m] + z[0:n-2, 0:m-2]
                            - z[1:n-1, 0:m-2] + z[1:n-1, 2:m]) / csz
        Nz = 6 * Nz
        # filling the first col of the matrix
        # (- 2*(Zr - Zp) + Zu - Zur ) / csz
        Nx[1:n-1, 0] = (- 2*(z[1:n-1, 1] - z[1:n-1, 0]) + z[2:n, 0] - z[2:n, 1]) / csz
        # (Zd - Zu + Zr - Zur) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0] + z[1:n-1, 1] - z[2:n, 1]) / csz
        Nz[1:n-1, 0] = 3
        # filling the last col of the matrix
        # (2*(Zl - Zp) + Zdl - Zd) / csz
        Nx[1:n-1, m-1] = (2*(z[1:n-1, m-2] - z[1:n-1, m-1]) + z[0:n-2, m-2] - z[0:n-2, m-1]) / csz
        # (Zd - Zu + Zdl - Zl) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1] + z[0:n-2, m-2] - z[1:n-1, m-2]) / csz
        Nz[1:n-1, m-1] = 3
        # filling the first row of the matrix
        # (Zl - Zr + Zu - Zur) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m] + z[1, 1:m-1] - z[1, 2:m]) / csz
        # (-2*(Zu - Zp) + Zr - Zur) / csz
        Ny[0, 1:m-1] = (- 2*(z[1, 1:m-1] - z[0, 1:m-1]) + z[0, 2:m] - z[1, 2:m]) / csz
        Nz[0, 1:m-1] = 3
        # filling the last row of the matrix
        # (Zl - Zr + Zdl - Zd) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m] + z[n-2, 0:m-2] - z[n-2, 1:m-1]) / csz
        # (2*(Zd - Zp) + Zdl - Zl) / csz
        Ny[n-1, 1:m-1] = (2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) + z[n-2, 0:m-2] - z[n-1, 0:m-2]) / csz
        Nz[n-1, 1:m-1] = 3
        # filling the corners of the matrix
        Nx[0, 0] = (z[1, 0] - z[1, 1] - (z[0, 1] - z[0, 0])) / csz
        Ny[0, 0] = (z[0, 1] - z[1, 1] - (z[1, 0] - z[0, 0])) / csz
        Nz[0, 0] = 2
        Nx[n-1, 0] = -(z[n-1, 1] - z[n-1, 0]) / csz
        Ny[n-1, 0] = (z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 1
        Nx[0, m-1] = (z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = -(z[1, m-1] - z[0, m-1]) / csz
        Nz[0, m-1] = 1
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1] + z[n-2, m-2] - z[n-2, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1] + z[n-2, m-2] - z[n-1, m-2]) / csz
        Nz[n-1, m-1] = 2

    if num == 8:
        # filling the inside of the matrix
        # normal calculation with 8 triangles
        # (2*(Zl - Zr) + Zul - Zur + Zdl - Zdr) / csz
        Nx[1:n-1, 1:m-1] = (2 * (z[1:n-1, 0:m-2] - z[1:n-1, 2:m])
                            + z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] - z[0:n-2, 2:m]) / csz
        # (2*(Zd - Zu) - Zul - Zur + Zdl + Zdr) / csz
        Ny[1:n-1, 1:m-1] = (2 * (z[0:n-2, 1:m-1] - z[2:n, 1:m-1])
                            - z[2:n, 0:m-2] - z[2:n, 2:m]
                            + z[0:n-2, 0:m-2] + z[0:n-2, 2:m]) / csz
        Nz = 8 * Nz
        # filling the first col of the matrix
        # (- 2*(Zr - Zp) + Zu - Zur + Zd - Zdr) / csz
        Nx[1:n-1, 0] = (- 2*(z[1:n-1, 1] - z[1:n-1, 0]) + z[2:n, 0] - z[2:n, 1] + z[0:n-2, 0] - z[0:n-2, 1]) / csz
        # (Zd - Zu + Zdr - Zur) / csz
        Ny[1:n-1, 0] = (z[0:n-2, 0] - z[2:n, 0] + z[0:n-2, 1] - z[2:n, 1]) / csz
        Nz[1:n-1, 0] = 4
        # filling the last col of the matrix
        # (2*(Zl - Zp) + Zdl - Zd + Zul - Zu) / csz
        Nx[1:n-1, m-1] = (2*(z[1:n-1, m-2] - z[1:n-1, m-1]) + z[0:n-2, m-2] - z[0:n-2, m-1] + z[2:n, m-2] - z[2:n, m-1]) / csz
        # (Zd - Zu + Zdl - Zul) / csz
        Ny[1:n-1, m-1] = (z[0:n-2, m-1] - z[2:n, m-1] + z[0:n-2, m-2] - z[2:n, m-2]) / csz
        Nz[1:n-1, m-1] = 4
        # filling the first row of the matrix
        # (Zl - Zr + Zul - Zur) / csz
        Nx[0, 1:m-1] = (z[0, 0:m-2] - z[0, 2:m] + z[1, 0:m-2] - z[1, 2:m]) / csz
        # (-2*(Zu - Zp) + Zr - Zur + Zl - Zul) / csz
        Ny[0, 1:m-1] = (- 2*(z[1, 1:m-1] - z[0, 1:m-1]) + z[0, 2:m] - z[1, 2:m] + z[0, 0:m-2] - z[1, 0:m-2]) / csz
        Nz[0, 1:m-1] = 4
        # filling the last row of the matrix
        # (Zl - Zr + Zdl - Zdr) / csz
        Nx[n-1, 1:m-1] = (z[n-1, 0:m-2] - z[n-1, 2:m] + z[n-2, 0:m-2] - z[n-2, 2:m]) / csz
        # (2*(Zd - Zp) + Zdl - Zl + Zdr - Zr) / csz
        Ny[n-1, 1:m-1] = (2*(z[n-2, 1:m-1] - z[n-1, 1:m-1]) + z[n-2, 0:m-2] - z[n-1, 0:m-2] + z[n-2, 2:m] - z[n-1, 2:m]) / csz
        Nz[n-1, 1:m-1] = 4
        # filling the corners of the matrix
        Nx[0, 0] = (z[1, 0] - z[1, 1] - (z[0, 1] - z[0, 0])) / csz
        Ny[0, 0] = (z[0, 1] - z[1, 1] - (z[1, 0] - z[0, 0])) / csz
        Nz[0, 0] = 2
        Nx[n-1, 0] = (-(z[n-1, 1] - z[n-1, 0]) + z[n-2, 0] - z[n-2, 1]) / csz
        Ny[n-1, 0] = (z[n-2, 1] - z[n-1, 1] + z[n-2, 0] - z[n-1, 0]) / csz
        Nz[n-1, 0] = 2
        Nx[0, m-1] = (z[1, m-2] - z[1, m-1] + z[0, m-2] - z[0, m-1]) / csz
        Ny[0, m-1] = (-(z[1, m-1] - z[0, m-1]) + z[0, m-2] - z[1, m-2]) / csz
        Nz[0, m-1] = 2
        Nx[n-1, m-1] = (z[n-1, m-2] - z[n-1, m-1] + z[n-2, m-2] - z[n-2, m-1]) / csz
        Ny[n-1, m-1] = (z[n-2, m-1] - z[n-1, m-1] + z[n-2, m-2] - z[n-1, m-2]) / csz
        Nz[n-1, m-1] = 2

    if num == 1:
        z1 = np.append(z, z[:, -2].reshape(n, 1), axis=1)
        n1, m1 = np.shape(z1)
        z2 = np.append(z1, z1[-2, :].reshape(1, m1), axis=0)
        n2, m2 = np.shape(z2)
        Nx = - (z2[0:n2-1, 1:m2] - z2[0:n2-1, 0:m2-1]) / csz
        Ny = - (z2[1:n2, 0:m2-1] - z2[0:n2-1, 0:m2-1]) / csz
        Ny[n-1, 0:m-1] = -Ny[n-1, 0:m-1]
        Nx[0:n-1, m-1] = -Nx[0:n-1, m-1]
        Ny[n-1, m-1] = -Ny[n-1, m-1]
        Nx[n-1, m-1] = -Nx[n-1, m-1]

        # Nx[0:n-1, 0:m-1] = - (z[0:n-1, 1:m] - z[0:n-1, 0:m-1]) / csz
        # Ny[0:n-1, 0:m-1] = - (z[1:n, 0:m-1] - z[0:n-1, 0:m-1]) / csz

    Nx, Ny, Nz = normalize(Nx, Ny, Nz)

    return Nx, Ny, Nz


def getAreaMesh(Nx, Ny, Nz, csz):
    A = csz * csz / Nz
    A = np.where(A > 3*csz*csz, 3*csz*csz, A)
    return A


def norm(x, y, z):
    norme = np.sqrt(x*x + y*y + z*z)
    return norme


def norm2(x, y, z):
    norme2 = (x*x + y*y + z*z)
    return norme2


def normalize(x, y, z):
    # TODO : avoid error message when input vector is zero and make sure
    # to return zero
    norme = norm(x, y, z)
    xn = x / norme
    xn = np.where(np.isnan(xn), 0, xn)
    yn = y / norme
    yn = np.where(np.isnan(yn), 0, yn)
    zn = z / norme
    zn = np.where(np.isnan(zn), 0, zn)

    return xn, yn, zn
