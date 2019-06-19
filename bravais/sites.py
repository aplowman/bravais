"""`bravais.sites.py`

Module defining a Sites class that represents a set of points in space.

"""

import warnings
import re

import numpy as np

from bravais.utils import check_indices


class SitesLabel(object):
    """Class to represent the labelling of a set of points in space.

    Attributes
    ----------
    name : str
    unique_values : ndarray
    values_idx : ndarray of int
    values : ndarray

    """

    def __init__(self, name, values=None, unique_values=None, values_idx=None):

        args = [values, unique_values, values_idx]
        msg = ('Specify either `values` or both `unique_values` and '
               '`values_idx`')
        if all([i for i in args]) or all([not i for i in args]):
            raise ValueError(msg)

        if values and not isinstance(values, np.ndarray):
            values = np.array(values)

        if unique_values and not isinstance(unique_values, np.ndarray):
            unique_values = np.array(unique_values)

        if values_idx and not isinstance(values_idx, np.ndarray):
            values_idx = np.array(values_idx)

        if values is not None:
            # Get unique `values` and place indices in `values_idx`:
            unique_values, values_idx = np.unique(values, return_inverse=True)

        else:
            # Check unique values are all unique:
            if len(np.unique(unique_values)) != len(unique_values):
                msg = ('Not all of the values in `unique_values` are unique.')
                raise ValueError(msg)

            # Check all `values_idx` do index `unique_values`:
            check_indices(unique_values, values_idx)

        self._validate_name(name)
        self.name = name
        self.unique_values = unique_values
        self.values_idx = values_idx

    @property
    def values(self):
        return self.unique_values[self.values_idx]

    @property
    def dtype(self):
        return self.unique_values.dtype

    def __len__(self):
        return len(self.values)

    def _validate_name(self, name):
        """Ensure name is safe to use as an object attribute."""
        pattern = r'^(?![0-9])[a-zA-Z0-9_]+$'
        if not re.match(pattern, name):
            msg = ('SitesLabel name "{}" is not valid since it cannot be '
                   'used as an object attribute name. Names must match the '
                   'regular expression "{}".')
            raise ValueError(msg.format(name, pattern))


class Sites(object):
    """Class to represent a set of points in space.

    Attributes
    ----------
    sites : ndarray
    dimension : int
    vector_direction : str
    labels : dict

    """

    def __init__(self, sites, labels=None, vector_direction='column',
                 dimension=3):

        self.vector_direction = vector_direction
        self._sites = self._validate(sites, self.vector_direction, dimension)
        self._dimension = dimension
        self.labels = self._init_labels(labels)

    def _init_labels(self, labels):
        """Set labels as attributes for easy access."""

        lable_objs = {}
        for k, v in (labels or {}).items():

            msg = ('Specify site labels as either a single list/tuple of '
                   'values, or as a list/tuple of length two, whose first '
                   'element is a list/tuple of unique values, and whose '
                   'second element is a list/tuple of indices that index the '
                   'first element.')
            values = None
            unique_values = None
            values_idx = None

            if isinstance(v[0], (np.ndarray, list, tuple)):
                if len(v) == 2:
                    unique_values, values_idx = v
                else:
                    raise ValueError(msg)
            else:
                values = v

            sites_label = SitesLabel(
                k,
                values=values,
                unique_values=unique_values,
                values_idx=values_idx
            )

            msg = ('Length of site labels named "{}" ({}) does not match '
                   'the number of sites ({}).')
            vals = sites_label.values
            if len(vals) != len(self):
                raise ValueError(msg.format(k, len(vals), len(self)))

            setattr(self, k, vals)
            lable_objs.update({k: sites_label})

        return lable_objs

    def __len__(self):
        """Get how many sites there are in this Sites objects."""
        return self._sites.shape[1]

    @property
    def dimension(self):
        return self._dimension

    @property
    def vector_direction(self):
        return self._vector_direction

    @vector_direction.setter
    def vector_direction(self, vector_direction):

        if vector_direction not in ['row', 'column', 'col']:
            msg = ('`vector_direction` must be specified as a string, either '
                   '"row" or "column" (or "col").')
            raise ValueError(msg)

        if vector_direction == 'col':
            vector_direction = 'column'

        if getattr(self, '_vector_direction', None):
            if vector_direction == getattr(self, '_vector_direction'):
                msg = '`vector_direction` is already set to "{}"'
                warnings.warn(msg.format(vector_direction))

        self._vector_direction = vector_direction

    def _validate(self, sites, vector_direction, dimension):
        """Validate inputs."""

        if dimension not in [2, 3]:
            msg = '`dimension` must be an integer: 2 or 3.'
            raise ValueError(msg)

        if not isinstance(sites, np.ndarray):
            sites = np.array(sites)

        if sites.ndim != 2:
            raise ValueError('`sites` must be a 2D array.')

        vec_len_idx = 0 if vector_direction == 'column' else 1
        vec_len = sites.shape[vec_len_idx]

        if vec_len != dimension:
            msg = ('The length of {}s in `sites` ({}) must be equal to '
                   '`dimension` ({}). Change `vector_direction` to "{}" if '
                   'you would like an individual site to be represented as '
                   'a {}-vector')
            non_vec_dir = 'row' if vector_direction == 'column' else 'column'
            raise ValueError(
                msg.format(
                    vector_direction,
                    vec_len,
                    dimension,
                    non_vec_dir,
                    non_vec_dir,
                )
            )

        if self.vector_direction == 'row':
            return sites.T
        else:
            return sites

    @property
    def sites(self):
        if self.vector_direction == 'column':
            return self._sites
        else:
            return self._sites.T

    def _validate_label_filter(self, **kwargs):
        """Validation for the `index` method."""

        if not self.labels:
            raise ValueError(
                'No labels are associated with this Sites object.')

        if not kwargs:
            msg = ('Provide a label condition to filter the sites. Available '
                   'labels are: {}')
            raise ValueError(msg.format(list(self.labels.keys())))

        if len(kwargs) > 1:
            msg = 'Only one label condition is currently supported by `whose`.'
            raise NotImplementedError(msg)

        for match_label, match_val in kwargs.items():
            try:
                getattr(self, match_label)
            except AttributeError as err:
                msg = 'No Sites label called "{}" was found.'
                raise ValueError(msg.format(match_label))

            return match_label, match_val

    def index(self, **kwargs):
        """Filter site indices by a label with a particular value."""

        match_label, match_val = self._validate_label_filter(**kwargs)
        label_vals = getattr(self, match_label)
        match_idx = np.where(label_vals == match_val)[0]

        return match_idx

    def whose(self, **kwargs):
        """Filter sites by a label with a particular value."""

        match_idx = self.index(**kwargs)
        match_sites = self._sites[:, match_idx]

        if self.vector_direction == 'row':
            match_sites = match_sites.T

        return match_sites
