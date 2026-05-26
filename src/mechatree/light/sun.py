"""Sun direction grid for hemispherical light integration.

Ports the elevation / azimuth construction from ``legacy/fortran/tree.f90:139``::

    elev(k) = acos((i - 0.5) / Nelev)        i = 1..Nelev
    azim(k) = 2*pi * (j - 1) / Nazim          j = 1..Nazim

The ``acos((i-0.5)/N)`` schedule places elevations so the *cosines* are
uniformly spaced — a discrete sampling of the Lambert (cosine-weighted)
hemisphere, which is the standard isotropic-sky integration scheme.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class Sun:
    """A hemispherical grid of sun directions for ``intercept``.

    Parameters
    ----------
    n_elevations:
        Number of elevation bands (Fortran ``Nelev``; default 4).
    n_azimuths:
        Number of azimuth bins per band (Fortran ``Nazim``; default 8).
    size_leaf:
        Cell size of the 2D shadow grid in world units (Fortran
        ``SizeLeaf``; default 1.0). A leaf falls into one cell.
    """

    n_elevations: int = 4
    n_azimuths: int = 8
    size_leaf: float = 1.0

    # Computed grids (n_dir = n_elevations * n_azimuths).
    elev: np.ndarray = field(init=False, repr=False)
    azim: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if self.n_elevations <= 0 or self.n_azimuths <= 0:
            raise ValueError("n_elevations and n_azimuths must be positive")
        if self.size_leaf <= 0.0:
            raise ValueError("size_leaf must be positive")

        i = np.arange(1, self.n_elevations + 1, dtype=np.float64)
        j = np.arange(1, self.n_azimuths + 1, dtype=np.float64)
        elev_band = np.arccos((i - 0.5) / self.n_elevations)
        azim_band = 2.0 * np.pi * (j - 1.0) / self.n_azimuths
        # Match the Fortran walk order: outer loop i over elevations, inner
        # loop j over azimuths — k goes 1..Nelev*Nazim.
        elev, azim = np.meshgrid(elev_band, azim_band, indexing="ij")
        object.__setattr__(self, "elev", elev.ravel())
        object.__setattr__(self, "azim", azim.ravel())

    @property
    def n_directions(self) -> int:
        return len(self.elev)

    @classmethod
    def from_arrays(cls, elev: np.ndarray, azim: np.ndarray, size_leaf: float = 1.0) -> Sun:
        """Construct a Sun with arbitrary elev/azim arrays.

        Bypasses the Fortran ``acos((i-0.5)/N)`` schedule — useful for tests
        that need a sun straight overhead (``elev=[0.0]``) or for users
        targeting a non-default integration scheme.
        """
        elev = np.asarray(elev, dtype=np.float64)
        azim = np.asarray(azim, dtype=np.float64)
        if elev.shape != azim.shape or elev.ndim != 1:
            raise ValueError("elev and azim must be 1-D arrays of the same length")
        if size_leaf <= 0.0:
            raise ValueError("size_leaf must be positive")
        s = cls.__new__(cls)
        object.__setattr__(s, "n_elevations", len(elev))
        object.__setattr__(s, "n_azimuths", 1)
        object.__setattr__(s, "size_leaf", size_leaf)
        object.__setattr__(s, "elev", elev)
        object.__setattr__(s, "azim", azim)
        return s
