"""bravais.utils.py`

Utility functions.

"""


def check_indices(seq, seq_idx):
    """
    Given a sequence (e.g. list, tuple, ndarray) which is indexed by another,
    check the indices are sensible.

    TODO: check integers as well?

    Parameters
    ----------
    seq : sequence
    seq_idx : sequence of int

    """

    # Check: minimum index is greater than zero
    if min(seq_idx) < 0:
        raise IndexError('Found index < 0.')

    # Check maximum index is equal to length of sequence - 1
    if max(seq_idx) > len(seq) - 1:
        raise IndexError('Found index larger than sequence length.')
