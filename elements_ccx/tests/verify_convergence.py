#!/usr/bin/env python3
"""Manufactured-solution convergence of F-barES-FEM-T4 against C3D4.

The exact solution is a smooth, divergence-free field on the unit cube,

    u* = curl A,   A = (sin 2pi y sin 2pi z,
                        sin 2pi z sin 2pi x,
                        sin 2pi x sin 2pi y),

so that div u* = 0 exactly and, for an isotropic material,

    -div sigma(u*) = 8 pi^2 G u*.

Divergence-free data is what makes this a volumetric-locking test: the
exact solution carries no volumetric energy, and an element whose discrete
incompressibility constraint is too strong cannot approximate it as the
bulk modulus grows.  The body force is G times the shape of u* (no K), so
the same exact field solves the single-phase problem for every K/G and the
two-phase problem when the shear modulus is shared and only K differs
between the phases.

Two layouts are run:

  single   one material, K/G in the sweep;
  two      a cube inclusion [0.25, 0.75]^3 (material 1) in a matrix
           (material 0), same G, K_m = 2.17 G against K_i = (K/G) G.
           Every interface plane is a mesh plane for n divisible by 4, so
           the geometry is exact on every mesh and only h changes.

Errors are the exact discrete norms of the P1 field e = u_h - I_h u*:
L2 through the consistent tetrahedral mass matrix, H1 through the piecewise
constant gradients.  Both the prototype (`assemble`) and the CalculiX MMS
driver (`mms_ccx.py`) feed the same `errors` routine, so the two levels are
directly comparable.

    python3 elements_ccx/tests/verify_convergence.py --ns 8 12 16
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import scipy.spatial as spat

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smoothing_proto as proto  # noqa: E402

TWOPI2 = (2.0 * np.pi) ** 2
LAP_COEFF = 8.0 * np.pi ** 2

E_MATRIX, NU_MATRIX = 9.37e9, 0.33
G_INCLUSION = 4.4e5

Geom = Callable[[float, float, float, float, float], int]


def u_star(x: np.ndarray) -> np.ndarray:
    """Exact manufactured displacement at points x of shape (..., 3).

    Curl of A = (sin 2pi y sin 2pi z, sin 2pi z sin 2pi x,
    sin 2pi x sin 2pi y), hence divergence-free and periodic.

    Args:
        x: Point coordinates, last axis the three components.

    Returns:
        The exact displacement at every point.
    """
    s, c = np.sin(2.0 * np.pi * x), np.cos(2.0 * np.pi * x)
    out = np.empty_like(x)
    out[..., 0] = 2.0 * np.pi * s[..., 0] * (c[..., 1] - c[..., 2])
    out[..., 1] = 2.0 * np.pi * s[..., 1] * (c[..., 2] - c[..., 0])
    out[..., 2] = 2.0 * np.pi * s[..., 2] * (c[..., 0] - c[..., 1])
    return out


def geom_single(x: float, y: float, z: float, lo: float, hi: float) -> int:
    """One material everywhere.

    Args:
        x: Cube-centre coordinate.
        y: Cube-centre coordinate.
        z: Cube-centre coordinate.
        lo: Unused.
        hi: Unused.

    Returns:
        Always 1.
    """
    return 1


def geom_cube(x: float, y: float, z: float, lo: float, hi: float) -> int:
    """Cube inclusion [0.25, 0.75]^3 in a matrix, on cube centres.

    Args:
        x: Cube-centre coordinate.
        y: Cube-centre coordinate.
        z: Cube-centre coordinate.
        lo: Unused.
        hi: Unused.

    Returns:
        1 inside the inclusion, 0 in the matrix.
    """
    if 0.25 <= x < 0.75 and 0.25 <= y < 0.75 and 0.25 <= z < 0.75:
        return 1
    return 0


def iso(kg: float, g: float = G_INCLUSION) -> Tuple[float, float]:
    """Bulk and shear modulus from a ratio K/G.

    Args:
        kg: The ratio K/G.
        g: The shear modulus.

    Returns:
        The pair (K, G).
    """
    return kg * g, g


def matrix_props() -> Tuple[float, float]:
    """Isotropic moduli of the compliant matrix.

    Returns:
        The pair (K, G) of the matrix.
    """
    g = E_MATRIX / (2.0 * (1.0 + NU_MATRIX))
    return E_MATRIX / (3.0 * (1.0 - 2.0 * NU_MATRIX)), g


def boundary_mask(nodes: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Boolean mask of the nodes on the six faces of the unit cube.

    Args:
        nodes: Node coordinates, shape (n, 3).
        tol: Face-coordinate tolerance.

    Returns:
        True for every boundary node.
    """
    return ((nodes <= tol) | (nodes >= 1.0 - tol)).any(axis=1)


def delaunay_mesh(
    n: int, radius: float, seed: int = 3
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Jittered-grid Delaunay tetrahedralisation with a spherical inclusion.

    Args:
        n: Cells per edge of the underlying grid.
        radius: Radius of the inclusion, centred at (0.5, 0.5, 0.5).
        seed: Random seed of the jitter.

    Returns:
        The tuple (points, tets, mat).
    """
    g = np.linspace(0.0, 1.0, n + 1)
    x, y, z = np.meshgrid(g, g, g, indexing='ij')
    pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    h = 1.0 / n
    rng = np.random.default_rng(seed)
    interior = ~((pts <= 0.0) | (pts >= 1.0)).any(axis=1)
    pts[interior] += (rng.random((int(interior.sum()), 3)) - 0.5) * 0.4 * h
    tri = spat.Delaunay(pts)
    tets = tri.simplices
    p = pts[tets]
    det = np.linalg.det(p[:, 1:, :] - p[:, :1, :])
    # Qhull mixes orientations; keep the non-degenerate simplices and let
    # `smoothing_proto.grads` flip the negatively oriented ones.
    tets = tets[np.abs(det) > 1e-14]
    centre = pts[tets].mean(axis=1)
    mat = (np.sum((centre - 0.5) ** 2, axis=1) < radius ** 2).astype(int)
    return pts, tets, mat


def build_mesh(
    mesh: str, n: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the single-phase mesh of one case.

    Args:
        mesh: 'box' for the structured six-tet split, 'delaunay' for the
            unstructured tetrahedralisation of the jittered grid.
        n: Cells per edge.

    Returns:
        The tuple (nodes, tets, mat) with one material.
    """
    if mesh == 'delaunay':
        nodes, tets, _ = delaunay_mesh(n, radius=0.0)
        return nodes, tets, np.ones(len(tets), dtype=int)
    nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0, geom=geom_single)
    return nodes, tets, mat


def interior_elements(
    nodes: np.ndarray, tets: np.ndarray, n: int, layers: int = 2
) -> np.ndarray:
    """Elements with every node at least `layers` cells from the boundary.

    Args:
        nodes: Node coordinates.
        tets: Element connectivity.
        n: Cells per edge of the underlying grid.
        layers: Distance from the boundary, in cells.

    Returns:
        Boolean mask over the elements.
    """
    lo = layers / n
    inside = ((nodes >= lo) & (nodes <= 1.0 - lo)).all(axis=1)
    return inside[tets].all(axis=1)


def body_loads(
    nodes: np.ndarray, tets: np.ndarray, mat: np.ndarray, vol: np.ndarray,
    props: Dict[int, Tuple[float, float]],
) -> np.ndarray:
    """Consistent nodal loads of f = 8 pi^2 G u*.

    The nodal (vertex) quadrature integrates the linear shape functions
    exactly against a constant; for the smooth f here the load error is
    O(h^2), two orders below the H1 error under test.

    Args:
        nodes: Node coordinates.
        tets: Element connectivity.
        mat: Material index per element.
        vol: Element volumes.
        props: Material moduli, keyed by material index.

    Returns:
        The global load vector of length 3 * n_nodes.
    """
    coeff = np.array([LAP_COEFF * props[int(m)][1] for m in mat])
    weight = vol * coeff / 4.0
    ux = u_star(nodes)
    idx = (3 * tets[..., None] + np.arange(3)).reshape(-1)
    vals = (weight[:, None, None] * ux[tets]).reshape(-1)
    f = np.zeros(3 * len(nodes))
    np.add.at(f, idx, vals)
    return f


def solve_dirichlet(
    K: sp.csr_matrix, f: np.ndarray, fixed: np.ndarray, values: np.ndarray
) -> np.ndarray:
    """Solve K u = f with the Dirichlet degrees of freedom eliminated.

    Args:
        K: Assembled stiffness, shape (ndof, ndof).
        f: Load vector.
        fixed: Boolean mask of the prescribed degrees of freedom.
        values: Prescribed values on the full degree-of-freedom vector.

    Returns:
        The full solution vector, with the prescribed entries set to `values`.
    """
    free = ~fixed
    u = np.zeros(K.shape[0])
    u[fixed] = values[fixed]
    Kff = K[free][:, free].tocsc()
    rhs = f[free] - K[free][:, fixed] @ u[fixed]
    u[free] = spl.spsolve(Kff, rhs)
    return u


def errors(
    nodes: np.ndarray, tets: np.ndarray, mat: np.ndarray, g: np.ndarray,
    vol: np.ndarray, u: np.ndarray, phase: Optional[int] = None,
    element_mask: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """Relative L2 and H1 errors of u against I_h u*.

    Args:
        nodes: Node coordinates.
        tets: Element connectivity.
        mat: Material index per element.
        g: Shape-function gradients, shape (n_elem, 4, 3).
        vol: Element volumes.
        u: The discrete solution, length 3 * n_nodes.
        phase: If given, restrict the norms to the elements of that material.
        element_mask: Additional element restriction, or None.

    Returns:
        The pair (relative L2 error, relative H1 error), both measured
        against the P1 interpolant of the exact solution.
    """
    ex = u_star(nodes)
    en = u.reshape(-1, 3) - ex
    mask = np.ones(len(tets), dtype=bool)
    if phase is not None:
        mask = mat == phase
    if element_mask is not None:
        mask = mask & element_mask
    vol = vol[mask]
    en4 = en[tets[mask]]
    ex4 = ex[tets[mask]]
    gg = g[mask]

    mass = (vol / 20.0)[:, None, None] * (np.ones((4, 4)) + np.eye(4))
    dot_e = np.einsum('eac,ebc->eab', en4, en4)
    dot_x = np.einsum('eac,ebc->eab', ex4, ex4)
    l2_e = float((mass * dot_e).sum())
    l2_x = float((mass * dot_x).sum())

    grad_e = np.einsum('eac,eaj->ecj', en4, gg)
    grad_x = np.einsum('eac,eaj->ecj', ex4, gg)
    h1_e = float((vol * np.einsum('ecj,ecj->e', grad_e, grad_e)).sum())
    h1_x = float((vol * np.einsum('ecj,ecj->e', grad_x, grad_x)).sum())
    return float(np.sqrt(l2_e / l2_x)), float(np.sqrt(h1_e / h1_x))


def study(
    name: str, geom: Geom, props: Dict[int, Tuple[float, float]],
    scheme: str, ns: Sequence[int], phase: Optional[int],
    mesh: str = 'box', interior: int = 0,
) -> List[Tuple[int, float, float, float]]:
    """Run one scheme on a mesh sequence.

    Args:
        name: Label printed in the table.
        geom: Material geometry callback for `mesh_box`.
        props: Material moduli.
        scheme: Scheme name understood by `assemble`.
        ns: Mesh sizes.
        phase: Material to restrict the norms to, or None for the whole cell.
        mesh: 'box' or 'delaunay'; the latter ignores `geom` and is
            single-phase only.
        interior: Element layers to strip from the boundary when reporting
            the interior errors, or 0 for no interior row.

    Returns:
        Rows (n, h, relative L2 error, relative H1 error).
    """
    rows: List[Tuple[int, float, float, float]] = []
    inner: List[Tuple[int, float, float]] = []
    for n in ns:
        if mesh == 'delaunay':
            nodes, tets, mat = build_mesh(mesh, n)
        else:
            nodes, tets, mat = proto.mesh_box(n, 0.0, 1.0, 0.0, geom=geom)
        g, vol = proto.grads(nodes, tets)
        faces, patch = proto.topology(tets, mat)
        K = proto.assemble(scheme, nodes, tets, mat, g, vol, faces, patch,
                           props)
        f = body_loads(nodes, tets, mat, vol, props)
        fixed = np.repeat(boundary_mask(nodes), 3)
        u = solve_dirichlet(K, f, fixed, u_star(nodes).reshape(-1))
        l2, h1 = errors(nodes, tets, mat, g, vol, u, phase)
        rows.append((n, 1.0 / n, l2, h1))
        line = '  %-14s n=%-3d  L2 %.4e  H1 %.4e' % (name, n, l2, h1)
        if interior:
            mask = interior_elements(nodes, tets, n, interior)
            l2i, h1i = errors(nodes, tets, mat, g, vol, u, phase,
                              element_mask=mask)
            inner.append((n, l2i, h1i))
            line += '   | interior L2 %.4e  H1 %.4e' % (l2i, h1i)
        print(line)
    if len(rows) >= 2:
        (n1, _, l2a, h1a), (n2, _, l2b, h1b) = rows[-2], rows[-1]
        rate_l2 = np.log(l2a / l2b) / np.log(n2 / n1)
        rate_h1 = np.log(h1a / h1b) / np.log(n2 / n1)
        print('  %-14s last-pair rates: L2 %.2f  H1 %.2f' %
              (name, rate_l2, rate_h1))
    if len(inner) >= 2:
        (n1, l2a, h1a), (n2, l2b, h1b) = inner[-2], inner[-1]
        print('  %-14s interior rates: L2 %.2f  H1 %.2f'
              % ('', np.log(l2a / l2b) / np.log(n2 / n1),
                 np.log(h1a / h1b) / np.log(n2 / n1)))
    return rows


def main() -> int:
    """Run the single- and two-phase studies.

    Returns:
        Process exit status.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument('--ns', type=int, nargs='+', default=[8, 12, 16],
                    help='elements per edge (multiples of 4)')
    ap.add_argument('--kg', type=float, nargs='+', default=[100.0, 5000.0],
                    help='inclusion K/G values')
    ap.add_argument('--phase', choices=['single', 'two', 'all'], default='all')
    ap.add_argument('--schemes', nargs='+',
                    default=['c3d4', 'fbar_0', 'fbar_1'])
    ap.add_argument('--csv', default=None,
                    help='append rows to this CSV instead of only printing')
    ap.add_argument('--mesh', choices=['box', 'delaunay'], default='box',
                    help='structured six-tet split or an unstructured '
                         'Delaunay tetrahedralisation (single phase only)')
    ap.add_argument('--interior', type=int, default=0,
                    help='boundary layers stripped for the interior rates')
    a = ap.parse_args()
    if a.mesh == 'delaunay' and a.phase == 'two':
        raise SystemExit('verify_convergence: the two-phase cube geometry '
                         'needs the structured mesh')
    if a.csv:
        with open(a.csv, 'w') as handle:
            handle.write('layout,kg,scheme,phase,n,h,l2,h1\n')

    def record(layout: str, kg: float, scheme: str, phase: Optional[int],
               rows: List[Tuple[int, float, float, float]]) -> None:
        """Append rows to the optional CSV.

        Args:
            layout: 'single' or 'two'.
            kg: The K/G ratio.
            scheme: Scheme label.
            phase: Material the norms were restricted to, or None.
            rows: The rows returned by `study`.
        """
        if not a.csv:
            return
        with open(a.csv, 'a') as handle:
            for n, h, l2, h1 in rows:
                handle.write('%s,%g,%s,%s,%d,%.8e,%.8e,%.8e\n'
                             % (layout, kg, scheme,
                                '' if phase is None else phase, n, h, l2, h1))

    if a.phase in ('single', 'all'):
        print('\nSINGLE PHASE  (one material, divergence-free u*)')
        for kg in a.kg:
            props = {1: iso(kg)}
            print(' K/G = %g' % kg)
            for scheme in a.schemes:
                record('single', kg, scheme, None,
                       study(scheme, geom_single, props, scheme, a.ns, None,
                             mesh=a.mesh, interior=a.interior))

    if a.phase in ('two', 'all'):
        print('\nTWO PHASE  (cube inclusion, same G, only K differs)')
        km_ratio = 2.1666667          # nu = 0.30 in the matrix
        for kg in a.kg:
            props = {0: (km_ratio * G_INCLUSION, G_INCLUSION),
                     1: iso(kg)}
            print(' inclusion K/G = %g' % kg)
            for scheme in a.schemes:
                record('two', kg, scheme, None,
                       study(scheme, geom_cube, props, scheme, a.ns, None))
                record('two', kg, scheme, 1,
                       study(scheme + ' (incl)', geom_cube, props, scheme,
                             a.ns, 1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
