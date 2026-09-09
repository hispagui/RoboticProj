"""
so3.py : the Special Orthogonal group (Lie group, also algebraic).
    SO(3) = {M in R^3 x R^3 | MM^t = id, det(M) = 1}
    Used for rotations 
    Has a double covering by S^3 = {q = (a,b,c,d) in R^4 | a^2+b^2+c^2+d^2 = 1} 
    exponential map goes from Lie algebra so(3) -> Lie group SO(3)
    so(n) is the tangent space at identity (nxn skew-matrices)
"""

import numpy as np
_EPS = 1e-10  # used to avoid errors when computing with small angle thetas



# ---------------------------------------------------------------
# Lie algebra : so(3) is isomorphic to R^3 
# ---------------------------------------------------------------
def hat(w : np.ndarray) -> np.ndarray:
    # 3-vector => 3x3 symmetric-skew matrix
    w = np.asarray(w, dtype = float).reshape(3)
    return np.array([[0.0, -w[2], w[1]],
                    [w[2], 0.0, -w[0]],
                    [-w[1], w[0], 0.0]])

def vee(W : np.ndarray) -> np.ndarray:
    # inverse of hat
    W = np.asarray(W, dtype=float)
    return np.array([W[2, 1], W[0, 2], W[1, 0]])

# ---------------------------------------------------------------
# mapping between Lie algebra so(3) and Lie group SO(3) (exp and log maps)
# ---------------------------------------------------------------
def exp(w : np.ndarray) -> np.ndarray:
    # rotation vector so(3) -> rotation matrix SO(3)
    # e^w = I + w + w^2/2! + w^3/3! + ...
    # recall w represents a rotation, w = theta * u for u a unit vector
    # => Rodrigues' rotation formula : I + sin(theta)u + (1-cos(theta))u^2    
    w = np.asarray(w, dtype=float).reshape(3) 
    theta = float(np.linalg.norm(w)) # L2 norm
    K = hat(w) # do not forget, it represents a rotation by theta.
    if theta < _EPS:
        # for small theta, cos theta ~ 1 - cos theta ^2 / 2 + Taylor expansion (removing thetas)
        A = 1.0 - theta * theta / 6.0
        B = 0.5 - theta * theta / 24.0
    else:
        A = np.sin(theta) / theta
        B = (1 - np.cos(theta)) /(theta * theta)
    return np.eye(3) + A*K + B*(K@K)

def log(R : np.ndarray) -> np.ndarray:
    # rotation matrix SO(3) -> rotation vector so(3) (inverse of exp map)
    # log R = (R-I) - (R-I)^2/2 + (R-I)^3/3 - ...
    # Rodrigues' rotation formula : log R = theta (R - R^T) / 2 sin(theta)
    # with theta = cos^-1( (tr(R)-1)/2 ) 
    R = np.asarray(R, dtype = float)
    cos_theta = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    theta = float(np.arccos(cos_theta))
    if theta < 1e-8: # theta is near 0, here (R - R^T)/2 ~ hat(w) to first order.
        return vee(0.5 * (R - R.T))  # vee(hat(w)) = w
    if theta > np.pi - 1e-6: # theta is near pi, we recover the axis from the eigenvector of R
                             # whose eigenvalue is +1 (the rotation axis is fixed by R)
        vals, vecs = np.linalg.eig(R)
        idx = int(np.argmin(np.abs(vals - 1.0)))
        axis = np.real(vecs[:, idx])
        axis = axis / np.linalg.norm(axis)
        # Fix the sign so it agrees with the (tiny) skew part when available.
        skew = vee(R - R.T)
        if np.dot(skew, axis) < 0:
            axis = -axis
        return theta * axis
    return (theta / (2.0 * np.sin(theta))) * vee(R - R.T)

# ---------------------------------------------------------------
# interpolation
# ---------------------------------------------------------------
def geodesic_distance(P: np.ndarray, Q: np.ndarray) -> float:
    # (Riemannian) metric in SO(3) : the rotation angle of P^T Q
    # length of the shortest great_circle path between P and Q (lies in [0,pi])
    # d_R(P,Q) = 1/sqrt(2)*|log(P^-1*Q)
    return float(np.linalg.norm(log(P.T @ Q)))

def interpolate(R0: np.ndarray, R1: np.ndarray, t: float) -> np.ndarray:
    # spherical linear interpolation in SO(3), slerp(R0,R1,t) = (1-t)R0 + tR1 with t in [0,1]
    # constant-speed geodesic between rotations(points)
    return R0 @ exp(t * log(R0.T @ R1)) # slerp(R0, R1, t)

# ---------------------------------------------------------------
# quaternions interpretation (the double cover S^3 = SU(2))
# ---------------------------------------------------------------
def quat_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    # rotation, axis -> unit quaternion
    axis = np.asarray(axis, dtype=float).reshape(3)
    n = np.linalg.norm(axis)
    if n < _EPS:
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = axis / n
    h = angle / 2.0
    return np.concatenate(([np.cos(h)], np.sin(h) * axis))

def quat_to_matrix(q: np.ndarray) -> np.ndarray:
    # unit quaternion q = (w, x, y, z) -> rotation matrix in SO(3)
    # note that R(q) == R(-q).
    w, x, y, z = np.asarray(q, dtype=float) / np.linalg.norm(q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y)],
        [2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y)],
    ])

def matrix_to_quat(R: np.ndarray) -> np.ndarray:
    # rotation matrix -> unit quaternion q (Shepperd's method; returns w >= 0)
    # since double covering, it only returns +q since -q is the other point above (in the double covering)
    R = np.asarray(R, dtype=float)
    tr = np.trace(R)
    if tr > 0:
        S = np.sqrt(tr + 1.0) * 2.0
        w = 0.25 * S
        x = (R[2, 1] - R[1, 2]) / S
        y = (R[0, 2] - R[2, 0]) / S
        z = (R[1, 0] - R[0, 1]) / S
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        S = np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2]) * 2.0
        w = (R[2, 1] - R[1, 2]) / S
        x = 0.25 * S
        y = (R[0, 1] + R[1, 0]) / S
        z = (R[0, 2] + R[2, 0]) / S
    elif R[1, 1] > R[2, 2]:
        S = np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2]) * 2.0
        w = (R[0, 2] - R[2, 0]) / S
        x = (R[0, 1] + R[1, 0]) / S
        y = 0.25 * S
        z = (R[1, 2] + R[2, 1]) / S
    else:
        S = np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1]) * 2.0
        w = (R[1, 0] - R[0, 1]) / S
        x = (R[0, 2] + R[2, 0]) / S
        y = (R[1, 2] + R[2, 1]) / S
        z = 0.25 * S
    q = np.array([w, x, y, z])
    if q[0] < 0:            # canonicalize to the w >= 0 hemisphere
        q = -q
    return q / np.linalg.norm(q)

def quat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    # Hamilton product of two quaternions (this is composition of rotations)
    w0, x0, y0, z0 = a
    w1, x1, y1, z1 = b
    return np.array([
        w0 * w1 - x0 * x1 - y0 * y1 - z0 * z1,
        w0 * x1 + x0 * w1 + y0 * z1 - z0 * y1,
        w0 * y1 - x0 * z1 + y0 * w1 + z0 * x1,
        w0 * z1 + x0 * y1 - y0 * x1 + z0 * w1,
    ])

def quat_conjugate(q: np.ndarray) -> np.ndarray:
    # conjugate of a quaternion q^* = a - bi -cj - dk
    w, x, y, z = q
    return np.array([w, -x, -y, -z])

def quat_slerp(q0: np.ndarray, q1: np.ndarray, t: float) -> np.ndarray:
    # spherical linear interpolation on S^3, along the "shorter" arc.

    q0 = q0 / np.linalg.norm(q0)
    q1 = q1 / np.linalg.norm(q1)
    dot = float(np.dot(q0, q1))
    if dot < 0.0: # manages the double cover, to interpolate the short way round
        q1 = -q1
        dot = -dot
    if dot > 0.9995:                      # nearly parallel: linearize
        q = q0 + t * (q1 - q0)
        return q / np.linalg.norm(q)
    theta_0 = np.arccos(dot)
    theta = theta_0 * t
    q_perp = q1 - q0 * dot # versor form
    q_perp /= np.linalg.norm(q_perp)
    return q0 * np.cos(theta) + q_perp * np.sin(theta) # slerp(q0, q1, t)

# --------------------------------------------------------------------------
# monodoromy handling, pi_1(SO(3)) = Z/2 (the universal cover of SO(3) is SU(2) and its a 2-fold cover)
# --------------------------------------------------------------------------
def loop_lift(axis: np.ndarray, turns: float, n: int = 400) -> np.ndarray:
    # lift a loop in SO(3) to a path of S^3 (lift a loop of rotations to a path of quaternions)
    # in SO(3) spinning about the axis through a total of "turns"*2*pi is a loop
    # lifted to S^3 the continuous path is q(theta) = (cos(theta/2), sin(theta/2) * axis)
    # return an (n, 4) array of quaternions along the path
    thetas = np.linspace(0.0, turns * 2.0 * np.pi, n) # [0.0, turn*2*pi/n, 2*turn*2*pi/n, ..., turn*2*pi]
    return np.array([quat_from_axis_angle(axis, th) for th in thetas])

def loop_lift_closes(axis: np.ndarray, turns: float) -> bool:
    # is the lifted quaternion path closed?
    # a single, 2*pi spin in SO(3) does NOT lift to a loop in S^3 
    # a double, 4*pi spin does lift to a loop in S^3
    path = loop_lift(axis, turns)
    return bool(np.allclose(path[0], path[-1], atol=1e-6)) # equality whithin a given tolerance
