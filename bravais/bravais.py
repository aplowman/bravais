"""`bravais.bravais.py`

Module defining a BravaisLattice class.

"""

from pathlib import Path
import pkgutil
from enum import Enum

import numpy as np
import yaml

from bravais.validator import NumericValidator

# Lattice sites (as column vectors) in fractional coordinates, for each
# centring type:
CENTRING_LATTICE_SITES = {
    'primitive': np.array([
        [0.0],
        [0.0],
        [0.0],
    ]),
    'base': np.array([
        [0.0, 0.5],
        [0.0, 0.5],
        [0.0, 0.0],
    ]),
    'body': np.array([
        [0.0, 0.5],
        [0.0, 0.5],
        [0.0, 0.5],
    ]),
    'face': np.array([
        [0.0, 0.5, 0.5, 0.0],
        [0.0, 0.5, 0.0, 0.5],
        [0.0, 0.0, 0.5, 0.5],
    ]),
    'rhombohedral': np.array([
        [0.0, 1/3, 2/3],
        [0.0, 2/3, 1/3],
        [0.0, 1/3, 2/3],
    ]),
}


class CentringType(Enum):
    """Class to represent the centring type."""

    primitive = 'P'
    base = 'B'
    body = 'I'
    face = 'F'
    rhombohedral = 'R'


class BravaisLattice(object):
    """Class to represent a Bravais lattice unit cell.

    Attributes
    ----------
    lattice_system : str
        Lattice system. One of "cubic", "rhombohedral", "orthorhombic",
        "tetragonal", "monoclinic", "triclinic", "hexagonal".
    centring_type : str
        Centring type also known as lattice type.
    a : float
        Lattice parameter, magnitude of the first unit cell edge vector.
    b : float
        Lattice parameter, magnitude of the second unit cell edge vector.
    c : float
        Lattice parameter, magnitude of the third unit cell edge vector.
    alpha : float
        Lattice parameter, angle in degrees between the second and third unit
        cell edge vectors.
    beta : float
        Lattice parameter, angle in degrees between the first and third unit
        cell edge vectors.
    gamma : float
        Lattice parameter, angle in degrees between the first and second unit
        cell edge vectors.
    unit_cell : ndarray of shape (3, 3)
        Array of column or row (depending on `row_or_column`) vectors defining
        the lattice unit cell.
    lat_sites_std: ndarray of shape (3, N)
        Array of column or row (depending on `row_or_column`) vectors defining
        the Cartesian positions of lattice sites within the lattice unit cell.
    lat_sites_frac: ndarray of shape (3, N)
        Array of column or row (depending on `row_or_column`) vectors defining
        the fractional positions of lattice sites within the lattice unit cell.
    row_or_column : str ("row" or "column")
        Specified whether `unit_cell`, `lat_sites_std` and `lat_sites_frac` are
        arrays of row or column vectors.

    Notes
    -----
    Lattice vectors are formed by aligning the crystallographic x-axis (i.e
    with magnitude `a`) along the Cartesian x-axis and aligning the
    crystallographic xy-plane (i.e. the vectors with magnitudes `a` and `b`)
    parallel to the Cartesian xy-plane.

    Conventional unit cells are generated and expected in the lattice parameter
    parameters. For instance the rhombohedral lattice system is represented
    with a rhombohedrally-centred hexagonal unit cell (the hexagonal axes/
    "setting"), rather than a primitive rhombohedral (cube stretched along its
    diagonal) cell.

    References
    ----------
    TODO: <Insert reference here to the 14 Bravais lattices in 3D.>

    Examples
    --------
    TODO: <Insert example here>

    TODO:
    -   Add primitive centring type to allowed centring types for rhombohedral
        lattice system.
    -   Regarding restrictions/validations on lattice types, need to allow
        restrictions to be "any n parameters must be..." rather than "the first
        two parameters must be". E.g. monoclinic: "two of the angle parameters
        must be 90 deg" rather than "parameters 3 and 5 must be 90 deg".
    -   Add align option ('ax' or 'cz').

    """

    def __init__(self, lattice_system=None, centring_type=None, a=None, b=None,
                 c=None, α=None, β=None, γ=None, alpha=None, beta=None,
                 gamma=None, degrees=True, alignment='ax',
                 row_or_column='column'):
        """Constructor method for BravaisLattice object.

        Parameters
        ----------
        lattice_system : str
            Lattice system is one of: cubic, rhombohedral, orthorhombic,
            tetragonal, monoclinic, triclinic, hexagonal.
        centring_type : str, optional
            The centring type of the lattice, also known as the lattice type,
            is one of P (primitive), B (base-centred), I (body-centred), F
            (face-centred) or R (rhombohedrally-centred). Not all centring
            types are compatible with all lattice systems. Default is None, in
            which case the rhombohedrally-centred centring type (for the
            rhombohedral lattice system) or primitive centring type (for all
            other lattice systems) will be chosen.
        a : float, optional
            Lattice parameter, magnitude of the first unit cell edge vector.
        b : float, optional
            Lattice parameter, magnitude of the second unit cell edge vector.
        c : float, optional
            Lattice parameter, magnitude of the third unit cell edge vector.
        α : float, optional
            Lattice parameter, angle in degrees between the second and third
            unit cell edge vectors. Specify `alpha` or `α`, but not both.
        β : float, optional
            Lattice parameter, angle in degrees between the first and third
            unit cell edge vectors. Specify `beta` or `β`, but not both.
        γ : float, optional
            Lattice parameter, angle in degrees between the first and second
            unit cell edge vectors. Specify `gamma` or `γ`, but not both.
        alpha : float, optional
            Lattice parameter, angle in degrees between the second and third
            unit cell edge vectors. Specify `alpha` or `α`, but not both.
        beta : float, optional
            Lattice parameter, angle in degrees between the first and third
            unit cell edge vectors. Specify `beta` or `β`, but not both.
        gamma : float, optional
            Lattice parameter, angle in degrees between the first and second
            unit cell edge vectors. Specify `gamma` or `γ`, but not both.
        row_or_column : str ("row" or "column"), optional
            Defines vector direction in array attributes. For instance, if
            "row", unit cell edge vectors are represented as row vectors in
            `unit_cell`. If "column", unit cell edge vectors are represented as
            column vectors in `unit_cell`. Default is "column".

        """

        lat_cent = self._validate_lattice_system(lattice_system, centring_type)
        self.lattice_system = lat_cent[0]
        self.centring_type = lat_cent[1]

        alpha, beta, gamma = self._normalise_angle_spec(
            alpha, beta, gamma, α, β, γ)

        lengths, angles = self._validate_params(a, b, c, alpha, beta, gamma)
        self.a = lengths['a']
        self.b = lengths['b']
        self.c = lengths['c']
        self.alpha = angles['alpha']
        self.beta = angles['beta']
        self.gamma = angles['gamma']

        self.row_or_column = row_or_column
        self.alignment = alignment
        self._unit_cell = self._compute_unit_cell(alignment)
        self.lattice_sites_frac = CENTRING_LATTICE_SITES[self.centring_type.name]

        # TODO: add a `degrees` boolean attribute for the angles.
        # TODO: mustn't be able to assign centring_type and lattice_system or
        # lengths/angles after init

    def _normalise_angle_spec(self, alpha, beta, gamma, α, β, γ):
        """Check angles are not specified as both spelled-out and greek
        symbols."""

        for i, j in zip([alpha, beta, gamma], [α, β, γ]):
            if i is not None and j is not None:
                msg = ('For each angle, specify either the spelled-out version'
                       ' (e.g. "alpha"), or the Greek letter version (e.g. '
                       '"α"), but not both!.')
                raise ValueError(msg)

        alpha = alpha or α
        beta = beta or β
        gamma = gamma or γ

        return alpha, beta, gamma

    def _compute_unit_cell(self, alignment):
        """Use the lattice parameters to form the unit cell with a particular
        alignment with the Cartesian axes.

        Parameters
        ----------
        alignment : str
            Alignment of the unit cell with repsect to the Cartesian axes, must
            be one either "ax" (align the `a` lattice vector with the x-axis,
            and the `b` lattice vector in the xy-plane) or "cz" (align the `c`
            lattice vector with the z-axis and the `b` lattice vector in the
            yz-plane).

        Returns
        -------
        unit_cell : ndarray of shape (3, 3)
            Array of column vectors representing the unit cell edge vectors.

        TODO: add some references for these equations. Consider doing one
        original alignment, and then rotating from that? Might want to support
        arbitrary orientations, which would involve specifying the orientation
        (as a vector) of one of the lattice vectors, plus a plane (specified
        as a normal vector?) in which a second lattice vector should lie.

        """

        alignment_opts = ['ax', 'cz']
        if alignment not in alignment_opts:
            msg = ('"{}" is not a valid axes alignment option. `align` must be'
                   ' one of: {}.'.format(alignment, alignment_opts))
            raise ValueError(msg)

        alpha_r = np.deg2rad(self.alpha)
        beta_r = np.deg2rad(self.beta)
        gamma_r = np.deg2rad(self.gamma)

        if alignment == 'ax':
            # Align `a` lattice vector with the x-axis and `b` lattice vector
            # in the xy-plane.

            a_x = self.a
            b_x = self.b * np.cos(gamma_r)
            b_y = self.b * np.sin(gamma_r)
            c_x = self.c * np.cos(beta_r)
            c_y = (abs(self.c) * abs(self.b) *
                   np.cos(alpha_r) - b_x * c_x) / b_y
            c_z = np.sqrt(self.c**2 - c_x**2 - c_y**2)

            unit_cell = np.array([
                [a_x, b_x, c_x],
                [0,  b_y, c_y],
                [0,  0, c_z]
            ])

        elif alignment == 'cz':
            # Align `c` lattice vector with the z-axis and `b` lattice vector
            # in the yz-plane.

            f = (1 - (np.cos(alpha_r))**2 - (np.cos(beta_r))**2 - (np.cos(gamma_r))**2
                 + 2 * np.cos(alpha_r) * np.cos(beta_r) * np.cos(gamma_r))**0.5
            a_x = self.a * f / np.sin(alpha_r)
            a_y = self.a * (np.cos(gamma_r) - np.cos(alpha_r)
                            * np.cos(beta_r)) / np.sin(alpha_r)
            a_z = self.a * np.cos(beta_r)
            b_y = self.b * np.sin(alpha_r)
            b_z = self.b * np.cos(alpha_r)
            c_z = self.c

            unit_cell = np.array([
                [a_x, 0, 0],
                [a_y, b_y, 0],
                [a_z, b_z, c_z]
            ])

        return unit_cell

    def _validate_lattice_system(self, lattice_system, centring_type):
        """Validate the specified lattice system and centring type."""

        if centring_type is None:
            if lattice_system == 'rhombohedral':
                centring_type = 'R'
            else:
                centring_type = 'P'

        # List valid Bravais lattice systems and centring types:
        all_lat = [
            'triclinic',
            'monoclinic',
            'orthorhombic',
            'tetragonal',
            'rhombohedral',
            'hexagonal',
            'cubic'
        ]
        centring_compat = {
            'P': list(set(all_lat) - set(['rhombohedral'])),
            'B': ['monoclinic', 'orthorhombic'],
            'I': ['orthorhombic', 'tetragonal'],
            'F': ['orthorhombic', 'cubic'],
            'R': ['rhombohedral'],
        }
        all_cent = list(centring_compat.keys())

        lattice_system = lattice_system.lower()
        centring_type = centring_type.upper()

        # Check specified lattice system and centring type are compatible:
        if lattice_system not in all_lat:
            msg = ('"{}" is not a valid lattice system. `lattice_system` must '
                   'be one of: {}.')
            raise ValueError(msg.format(lattice_system, all_lat))

        if centring_type not in all_cent:
            msg = ('"{}" is not a valid centering type. `centring_type` must '
                   'be one of {}.')
            raise ValueError(msg.format(centring_type, all_cent))

        if lattice_system not in centring_compat[centring_type]:
            msg = 'Lattice system {} and centring type {} are not compatible.'
            raise ValueError(msg.format(lattice_system, centring_type))

        centring_type = CentringType(centring_type)

        return lattice_system, centring_type

    def _validate_params(self, a, b, c, alpha, beta, gamma):
        """Validate lattice parameters against the specified lattice system."""

        data = pkgutil.get_data('bravais', 'valid_parameters.yml')
        valid_parameters = yaml.load(data)

        valid_lat_params = valid_parameters[self.lattice_system]
        length_validator = NumericValidator(**valid_lat_params['lengths'])
        angles_validator = NumericValidator(**valid_lat_params['angles'])

        lengths_valid = length_validator.validate(a=a, b=b, c=c)
        angles_valid = angles_validator.validate(
            alpha=alpha, beta=beta, gamma=gamma)

        return lengths_valid, angles_valid

    @property
    def row_or_column(self):
        return self._row_or_column

    @row_or_column.setter
    def row_or_column(self, row_or_column):
        if row_or_column not in ['row', 'column']:
            msg = ('`row_or_column` must be specified as a string, either '
                   '"row" or "column".')
            raise ValueError(msg)
        self._row_or_column = row_or_column

    @property
    def unit_cell(self):
        if self.row_or_column == 'column':
            return self._unit_cell
        elif self.row_or_column == 'row':
            return self._unit_cell.T

    @property
    def lattice_sites(self):
        """Get the position of the lattice sites in Cartesian coordinates."""
        return np.dot(self._unit_cell, self.lattice_sites_frac)

    @property
    def α(self):
        return self.alpha

    @property
    def β(self):
        return self.beta

    @property
    def γ(self):
        return self.gamma

    def __repr__(self):
        return (
            '<{}('
            'lattice_system="{}", '
            'centring_type="{}", '
            'a={:.4f}, '
            'b={:.4f}, '
            'c={:.4f}, '
            'alpha={:.2f}, '
            'beta={:.2f}, '
            'gamma={:.2f}'
            ')>'.format(
                self.__class__.__name__,
                self.lattice_system,
                self.centring_type.value,
                self.a,
                self.b,
                self.c,
                self.alpha,
                self.beta,
                self.gamma,
            )
        )

    def __str__(self):
        return (
            '{}-centred {} lattice ('
            'a={:.4f}, b={:.4f}, c={:.4f}, '
            'alpha={:.2f}, beta={:.2f}, gamma={:.2f}'
            ')'.format(
                self.centring_type.value,
                self.lattice_system,
                self.a,
                self.b,
                self.c,
                self.alpha,
                self.beta,
                self.gamma,
            )
        )
