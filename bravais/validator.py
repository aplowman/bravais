"""`bravais.validator.py`"""

import copy
import random


class ValidationError(Exception):
    """Raised when NumericValidator.validate fails."""


class NumericValidator(object):
    """Class to represent a set of validation conditions to be applied to
    a set of numeric parameters."""

    def __init__(self, equal_to=None, not_equal_to=None, equal_groups=None,
                 not_equal_groups=None, default_min=0, default_max=1):
        """
        Initialise a NumericValidator with a set of conditions.

        Parameters
        ----------
        equal_to : dict, optional
        not_equal_to : dict, optional
        equal_groups : list of list of str, optional
        not_equal_groups : list of list of str, optional
        default_min : int or float, optional
        default_max : int or float, optional

        """

        if equal_to is not None:
            if not isinstance(equal_to, dict):
                msg = '`equal_to` must be a dict.'
                raise ValueError(msg)

        if not_equal_to is not None:
            if not isinstance(not_equal_to, dict):
                msg = '`not_equal_to` must be a dict.'

        if equal_groups is not None:
            msg = '`equal_groups` must be a list of list of str.'
            if isinstance(equal_groups, list):
                for i in equal_groups:
                    if isinstance(i, list):
                        for j in i:
                            if not isinstance(j, str):
                                raise ValueError(msg)
                    else:
                        raise ValueError(msg)
            else:
                raise ValueError(msg)

        if not_equal_groups is not None:
            msg = '`not_equal_groups` must be a list of list of str.'
            if isinstance(not_equal_groups, list):
                for i in not_equal_groups:
                    if isinstance(i, list):
                        for j in i:
                            if not isinstance(j, str):
                                raise ValueError(msg)
                    else:
                        raise ValueError(msg)
            else:
                raise ValueError(msg)

        self.equal_to = equal_to or {}
        self.not_equal_to = not_equal_to or {}
        self.equal_groups = equal_groups or [[]]
        self.not_equal_groups = not_equal_groups or [[]]
        self.default_min = default_min
        self.default_max = default_max

        self._validate_conditions()

    def _get_parameter_equal_group(self, param):
        """Get the `equal_group` in which a parameter exists, if at all."""
        for i in self.equal_groups:
            if param in i:
                return i

    def _get_parameter_not_equal_group(self, param):
        """Get the `not_equal_group` in which a parameter exists, if at all."""
        for i in self.not_equal_groups:
            if param in i:
                return i

    def _validate_conditions(self):
        """Check conditions make sense."""

        msg = ('Parameter names should not be repeated within equal_groups` '
               'elements, nor within `not_equal_groups` elements')
        for i in (self.equal_groups + self.not_equal_groups):
            if len(set(i)) != len(i):
                raise ValueError(msg)

        msg = ('The parameter "{}" appears in more than one sublist of '
               '`equal_groups`; these two sublists should be merged.')
        eq_groups_params = []
        for i in self.equal_groups:
            for j in i:
                if j in eq_groups_params:
                    raise ValueError(msg.format(j))
                eq_groups_params.append(j)

        msg = ('The parameter "{}" appears in more than one sublist of '
               '`not_equal_groups`; these two sublist should be merged.')
        neq_groups_params = []
        for i in self.not_equal_groups:
            for j in i:
                if j in neq_groups_params:
                    raise ValueError(msg.format(j))
                neq_groups_params.append(j)

        msg = ('The parameters "{}" are specified in both `equal_to` and '
               '`not_equal_to`. Specify the parameter in only one of these.')
        eq_to_params = self.equal_to.keys()
        neq_to_params = self.not_equal_to.keys()
        intersect = list(set(eq_to_params) & set(neq_to_params))
        if intersect:
            raise ValueError(msg.format(intersect))

        msg = ('Only one of the parameters specified within an `equal_group`'
               'element may be specified in `equal_to`.')
        for i in self.equal_groups:
            if len(list(set(i) & set(self.equal_to.keys()))) > 1:
                raise ValueError(msg)

        msg = ('Parameters specified within a `not_equal_group` element cannot'
               ' be assigned to the same value in `equal_to`.')
        for i in self.not_equal_groups:
            eq_to_vals = [self.equal_to.get(j) for j in i
                          if self.equal_to.get(j) is not None]
            if len(set(eq_to_vals)) < len(eq_to_vals):
                raise ValueError(msg)

    def _set_default_value(self, all_params, params_to_set):
        """Assign a default value to a parameter.

        Parameters
        ----------
        all_params : dict
            All parameters. The the value of the parameter specified will be
            changed within this dict.
        params_to_set : list of str
            Names of the parameters to be set to the same value.

        Returns
        -------
        None

        """

        assert all([i in all_params for i in params_to_set])

        # Collect values param must not be:
        bad_vals = []
        for i in params_to_set:
            if i in self.not_equal_to:
                bad_vals.append(self.not_equal_to[i])

            ne_group = self._get_parameter_not_equal_group(i) or []
            for j in ne_group:
                if j != i:
                    ne_val = all_params[j]
                    if ne_val:
                        bad_vals.append(ne_val)

        # Find a suitable default value:
        default_val = None
        count = 0
        while default_val is None or default_val in bad_vals:

            default_val = random.random()
            default_val *= (self.default_max - self.default_min)
            default_val += self.default_min

            count += 1
            if count > 10:
                msg = ('Failed to find suitable default value for '
                       'parameter {}.')
                raise RuntimeError(msg.format(params_to_set))

        for i in params_to_set:
            all_params[i] = default_val

    def validate(self, **kwargs):
        """Validate named parameters."""

        params = kwargs.copy()

        # Firstly, validate according to rules, but ignore parameters set to
        # `None`:

        msg = 'Parameter "{}" must be equal to "{}".'
        for param_name, val in self.equal_to.items():
            if kwargs[param_name] is not None and kwargs[param_name] != val:
                raise ValidationError(msg.format(param_name, val))

        msg = 'Parameter "{}" must not be equal to "{}".'
        for param_name, val in self.not_equal_to.items():
            if kwargs[param_name] is not None and kwargs[param_name] == val:
                raise ValidationError(msg.format(param_name, val))

        msg = 'Parameters "{}" must all have the same value.'
        for i in self.equal_groups:
            first_val = None
            for j in i:
                if first_val is None and kwargs[j] is not None:
                    first_val = kwargs[j]
                    continue
                if kwargs[j] is not None and kwargs[j] != first_val:
                    raise ValidationError(msg.format(i))

        msg = 'Parameters "{}" must all have distinct values.'
        for i in self.not_equal_groups:
            all_vals = []
            for j in i:
                if kwargs[j] is not None:
                    if kwargs[j] in all_vals:
                        raise ValidationError(msg.format(i))
                    all_vals.append(kwargs[j])

        none_params = [k for k, v in kwargs.items() if v is None]

        # Fill in `None`s.
        # First set any `equal_to`s:
        for i in copy.copy(none_params):
            if i in self.equal_to:
                params[i] = self.equal_to[i]
                none_params.remove(i)

        # Then try using `equal_groups`:
        for eq_group in self.equal_groups:

            num_none = sum([params.get(i) is None for i in eq_group])

            if len(eq_group) > num_none > 0:
                # At least one param in the eq_group has a value, use this one.

                val = None
                for i in eq_group:
                    if params[i] is not None:
                        val = params[i]
                        break

                for i in eq_group:
                    if params[i] is None:
                        params[i] = val
                        none_params.remove(i)

            elif len(eq_group) == num_none:
                # None of the equal group values are set. Use a default.
                self._set_default_value(params, eq_group)
                for i in eq_group:
                    none_params.remove(i)

        # Try to set remaining `None`s using defaults. Need a default value
        # that doesn't conflict with `not_equal_to`, nor with
        # `not_equal_groups`.
        if none_params:
            for i in copy.copy(none_params):
                self._set_default_value(params, [i])
                none_params.remove(i)

        return params

    def __repr__(self):
        return (
            '<{}('
            'equal_to={}, '
            'not_equal_to={}, '
            'equal_groups={}, '
            'not_equal_groups={}'
            ')>'.format(
                self.__class__.__name__,
                self.equal_to,
                self.not_equal_to,
                self.equal_groups,
                self.not_equal_groups,
            )
        )
