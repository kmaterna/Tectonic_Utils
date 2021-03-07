"""
Calculations that deal with seismic moment tensors
Notes from Lay and Wallace Chapter 8:
Decomposition 1: Mij = isotropic + deviatoric
Decomposition 2: Mij = isotropic + 3 vector dipoles
Decomposition 3: Mij = isotropic + 3 double couples
Decomposition 4: Mij = isotropic + 3 CLVDs
Decomposition 5: Mij = isotropic + major DC + minor DC
Decomposition 6: Mij = isotropic + DC + CLVD
The most useful in practical terms are Decomposition 1 and Decomposition 6.
"""

import numpy as np

def get_MT(mrr, mtt, mpp, mrt, mrp, mtp):
    """build a matrix from the six components of the moment tensor"""
    MT = np.array([[mrr, mrt, mrp], [mrt, mtt, mtp], [mrp, mtp, mpp]]);
    return MT;

def diagonalize_MT(MT):
    """Return a diagonal matrix"""
    eigvals, eigvecs = np.linalg.eig(MT);
    return np.multiply(eigvecs, eigvals);

def get_deviatoric_MT(MT):
    """Get deviatoric MT"""
    iso_MT = get_iso_matrix(MT);
    M_dev = np.subtract(MT, iso_MT);
    return M_dev;

def get_iso_matrix(MT):
    """Return the isotropic part of a moment tensor (returns a matrix)"""
    x = (1 / 3) * np.trace(MT);
    iso_MT = np.multiply(np.eye(3), x);
    return iso_MT

def decompose_iso_dc_clvd(MT):
    """A useful function to decompose a moment tensor into an isotropic part, a double-couple, and a clvd component"""
    diag_MT = diagonalize_MT(MT)
    iso_MT = get_iso_matrix(diag_MT);
    _ = get_deviatoric_MT(diag_MT);
    # Next: assert components of deviatoric moment tensor are in the right order
    # Next: calculate epsilon
    # Next: Construct DC and CLVD matrices, return them.
    return iso_MT;

def get_scalar_moment():
    """Not sure how to do this yet"""
    return;
