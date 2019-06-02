"""Tests on the `bravais.validator` module."""

import unittest

from bravais.validator import NumericValidator, ValidationError


class NumericValidatorTestCase(unittest.TestCase):
    """Tests on the NumericValidator class."""

    def test_raise_on_repeated_intra_equal_group(self):
        """Test ValueError is raised when a parameter name is repeated within
        an `equal_group` element."""

        with self.assertRaises(ValueError):
            NumericValidator(equal_groups=[['a', 'a', 'b']])

    def test_raise_on_repeated_intra_not_equal_group(self):
        """Test ValueError is raised when a parameter name is repeated within
        a `not_equal_group` element."""

        with self.assertRaises(ValueError):
            NumericValidator(not_equal_groups=[['a', 'a', 'b']])

    def test_raise_on_repeated_inter_equal_group(self):
        """Test ValueError is raised when a parameter name is repeated across
        multiple `equal_group` elements."""

        with self.assertRaises(ValueError):
            NumericValidator(equal_groups=[
                ['a', 'b'],
                ['a', 'c'],
            ])

    def test_raise_on_repeated_inter_not_equal_group(self):
        """Test ValueError is raised when a parameter name is repeated across
        multiple `not_equal_group` elements."""

        with self.assertRaises(ValueError):
            NumericValidator(not_equal_groups=[
                ['a', 'b'],
                ['a', 'c'],
            ])

    def test_raise_on_inconsistent_assignment(self):
        """Test ValueError is raised when a parameter name appears in both
        `equal_to` and `not_equal_to`."""

        with self.assertRaises(ValueError):
            NumericValidator(
                equal_to={'a': 1},
                not_equal_to={'a': 2},
            )

    def test_raise_on_multiple_equal_group_assignment(self):
        """Test ValueError is raised when distinct parameters that belong
        to an `equal_group` element are assigned in `equal_to`."""

        with self.assertRaises(ValueError):
            NumericValidator(
                equal_groups=[
                    ['a', 'b']
                ],
                equal_to={
                    'a': 1,
                    'b': 2,
                }
            )

    def test_raise_on_multiple_not_equal_group_equal_assignment(self):
        """Test ValueError is raised when distinct parameters that belong
        to a `not_equal_group` element are assigned the same value in
        `equal_to`."""

        with self.assertRaises(ValueError):
            NumericValidator(
                not_equal_groups=[
                    ['a', 'b']
                ],
                equal_to={
                    'a': 1,
                    'b': 1,
                }
            )

    def test_invalid_on_paramater_equal(self):
        """Test ValidationError is raised when a parameter's value is other
        than specified in `equal_to`."""

        validator = NumericValidator(equal_to={'a': 1})
        with self.assertRaises(ValidationError):
            validator.validate(a=2)

    def test_invalid_on_paramater_not_equal(self):
        """Test ValidationError is raised when a parameter's value is equal to
        that specified in `not_equal_to`."""

        validator = NumericValidator(not_equal_to={'a': 1})
        with self.assertRaises(ValidationError):
            validator.validate(a=1)

    def test_invalid_on_paramater_group_not_equal(self):
        """Test ValidationError is raised when a set of parameters that belong
        to an `equal_group` element do not have equal values."""

        validator = NumericValidator(equal_groups=[
            ['a', 'b', 'c']
        ])
        with self.assertRaises(ValidationError):
            validator.validate(a=1, b=1, c=2)

    def test_invalid_on_parameter_group_equal(self):
        """Test ValidationError is raised when a set of parameters that belong
        to a `not_equal_group` element have equal values."""

        validator = NumericValidator(not_equal_groups=[
            ['a', 'b', 'c']
        ])
        with self.assertRaises(ValidationError):
            validator.validate(a=1, b=2, c=2)
