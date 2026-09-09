"""
se3.py : the Special Euclidean group (Lie group, also algebraic).
    SE(3) = {[[R,t],[0,1]] in R^4 x R^4 | R in SO(3), p in R^3}.
    Topologically homeomorphic to SO(3)xR^3 which is homeomorphic to a 6d manifold. 
    As Lie groups, isomorphic to SO(3) semidirprod R^3 (the group of translations).
    Elements of SE(3) are sometimes called 'pose', describing a position and an orientation.
    It is often called the group of rigid-body motion.
"""

import numpy as np
import so3
_EPS = 1e-10  # used to avoid errors when computing with small angle thetas


def make(R : np.ndarray, p : np.ndarray) -> np.ndarray:
    # assemble a 4x4 matrix in SE(3) from a rotation R in SO(3) and a translation p in R^3
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = np.asarray(p, dtype=float).reshape(3)
    return T

def split(T : np.ndarray) -> np.ndarray:
    # returns (R,p) from an element of SE(3)
    T = np.asarray(T, dtype=float)
    return (T[:3,:3].copy(), T[:3,3].copy())

def inverse(T : np.ndarray) -> np.ndarray:
    # returns elements [[R^T, -R^Tp],[0, 0, 0, 1]]
    R, p = split(T)
    return make(R.T, -R.T @ p)


# ---------------------------------------------------------------
# mappings between Lie algebra se(3) and Lie group SE(3) (exp and log maps)
# ---------------------------------------------------------------
def left_jacobian(omega: np.ndarray) -> np.ndarray:
    # in closed-form expression I + (1-cos(theta))/theta^2 omega^ + (theta - sin(theta))/theta^3 (omega^)^2
    # turns translational twist omega into actual displacement on curved surface during exp mapping
    omega = np.asarray(omega, dtype=float).reshape(3)
    theta = float(np.linalg.norm(omega))
    K = so3.hat(omega)
    if theta < _EPS: 
        # for small theta, cos and sin are approximated
        return np.eye(3) + 0.5 * K + (1.0 / 6.0) * (K @ K)
    A = (1.0 - np.cos(theta)) / (theta * theta)
    B = (theta - np.sin(theta)) / (theta ** 3)
    return np.eye(3) + A * K + B * (K @ K)

def exp(v: np.ndarray) -> np.ndarray:
    # twist vector xi in se(3) ->  pose matrix T in SE(3)
    v = np.asarray(v, dtype=float).reshape(6)
    rho, phi = v[:3], v[3:]
    R = so3.exp(phi) # we can use methods from so(3)->SO(3)
    p = left_jacobian(phi) @ rho
    return make(R, p)
 
def log(T: np.ndarray) -> np.ndarray:
    # pose matrix T in SE(3)  ->  twist vector xi in se(3) (inverse of exp map)
    R, p = split(T)
    phi = so3.log(R) # w can use methods from SO(3)->so(3)
    V = left_jacobian(phi)
    rho = np.linalg.solve(V, p) # numpy implicitly computs rho = V^-1 p (equation is rho = Vp)
    return np.concatenate([rho, phi])

def adjoint(T: np.ndarray) -> np.ndarray:
    # adjoint matrices map twists vectors between coordinate frames (xi' = Ad_T * xi) 
    # a 6x6 matrix : Ad_T = [[ R,      [p]x R ],
    #                        [ 0(3x3), R      ]]
    R, p = split(T)
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[:3, 3:] = so3.hat(p) @ R
    Ad[3:, 3:] = R
    return Ad

# ---------------------------------------------------------------
# interpolation 
# ---------------------------------------------------------------
# imagine the vector orientation is like a screw
def screw_interpolate(T0: np.ndarray, T1: np.ndarray, t: float) -> np.ndarray:
    # SE(3) analogue of slerp, geodesic pose interpolation, with continuous evolution of rotation
    # log(T0^{-1} T1) the "constant" velo. (screw motion) to travel from T0 to T1 along shortest SE(3) geod.
    # exp(t log(T0^{-1} T1)) scalling and projects back onto curved manifold
    return T0 @ exp(t * log(inverse(T0) @ T1))

def decoupled_interpolate(T0: np.ndarray, T1: np.ndarray, t: float) -> np.ndarray:
    # not a SE(3) geod., actually a straight line (often cheaper)
    # but keeps track of rotating SO(3) using slerp, so that it matches with rotation of geod_interpolation
    R0, p0 = split(T0)
    R1, p1 = split(T1)
    q = so3.quat_slerp(so3.matrix_to_quat(R0), so3.matrix_to_quat(R1), t)
    R = so3.quat_to_matrix(q)
    p = (1.0 - t) * p0 + t * p1
    return make(R, p)

def geodesic_distance(T0: np.ndarray, T1: np.ndarray, w_rot: float = 1.0, w_trans: float = 1.0) -> float:
    # computes geodesic length from linear velocity and angular velocity components
    # log(T0^{-1} T1) the "constant" velo. to travel from T0 to T1 along shortest SE(3) geod.
    # from thm 3.14 of Lee (Existence of Bi-invariant Metrics), SE(3) admits no bi-invariant metric because Ad(SE(3)) is unbounded (so closure cannot be compact)
    # hence the need to specify an arbitrary unit scale factor (e.g. matching 1 radian or 1 meter)
    xi = log(inverse(T0) @ T1)
    rho, phi = xi[:3], xi[3:]
    return float(np.sqrt(w_trans * rho @ rho + w_rot * phi @ phi))
