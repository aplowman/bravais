"""`bravais.bravais.py`

Module defining a BravaisLattice class.

"""

from pathlib import Path
import pkgutil

import numpy as np
import yaml

from bravais.validator import NumericValidator

# Lattice sites (as column vectors) in fractional coordinates, for each
# centring type:
centring_lattice_sites = {
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
centring_type_code = {
    'P': 'primitive',
    'B': 'base',
    'I': 'body',
    'F': 'face',
    'R': 'rhombohedral',
}


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
    with a rhombohedrally-centred hexagonal unit cell, rather than a primitive
    rhombohedral cell.

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
        self._unit_cell = self._get_unit_cell()

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

    def _get_unit_cell(self):
        """Use the lattice parameters to form the unit cell with a particular
        alignment with the Cartesian axes."""
        return np.array([])

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
            'P': all_lat,
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
    def α(self):
        return self.alpha

    @property
    def β(self):
        return self.beta

    @property
    def γ(self):
        return self.gamma
