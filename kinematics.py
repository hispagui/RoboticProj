"""
kinameatics.py : 6-DOF manipulator
    class Manipulator keeps track of screws, whereeach screw is a pose in SE(3)
    Forward Kinematics (from angle of each joint, determines the end effector pose)
    Inverse Kinematics (from end effector pose, determines angle of each joint)
    Computes a few Jacobians and checks for validity
"""


import numpy as np
from maths import so3
from maths import se3
import math

def screw_axis(omega : np.ndarray, q : np.ndarray) -> np.ndarray:
    # omega = unit rotation axis and q = any point on that axis 
    # basically an element in SE(3)
    omega = np.asarray(omega, float)
    q = np.asarray(q, float)
    return np.concatenate([-np.cross(omega, q), omega])

class Manipulator :
    def __init__(self, screws : np.ndarray, home : np.ndarray, name = "arm", joint_limits = None, draw_points = None):
        self.screws = np.asarray(screws, float).reshape(-1, 6) # revolute joints
        self.M = np.asarray(home,float) # pose of ee when every joint angle is 0
        if self.M.shape != (4,4):
            raise ValueError("home pose must be a 4x4 matrix (in se3)")
        self.dof = len(self.screws)
        self.name = name
        self.joint_limits = (None if joint_limits is None
                             else np.asarray(joint_limits, float).reshape(-1,2)) # possibility to define restriction on joint angles for more realistic results
        self.draw_points = draw_points

    @classmethod
    def from_revolute(cls, axes_and_points : np.ndarray, home, **kwargs) -> "Manipulator":
        # builds objects from a list of (omega, q), one per revolution joint
        screws = [screw_axis(w, q) for (w, q) in axes_and_points]
        return cls(screws, home, **kwargs)

    def _theta(self, thetas : np.ndarray) -> np.ndarray:
        # angles for each frame (also a sanity check)
        thetas = np.asarray(thetas, float).reshape(-1)
        if thetas.size != self.dof:
            raise ValueError(f"{self.name}: expected {self.dof} joint angles, " f"got {thetas.size}")
        return thetas

    def __repr__(self) -> str:
        return f"Manipulator(name={self.name!r}, dof={self.dof})"


    # ---------------------------------------------------------------
    # forward kinematics (thetas => EE pose)
    # ---------------------------------------------------------------
    def fk(self, thetas : np.ndarray) -> np.ndarray: # forward kinematics
        # PoE product for given angles and home pose (basically a chain of se3.exp calls)
        T = np.eye(4)
        for S, th in zip(self.screws, self._theta(thetas)):
            T = T @ se3.exp(th * S)  # forward kinematics equation
        final_eq = T @ self.M 
        return final_eq
 
    def ee_position(self, thetas : np.ndarray) -> np.ndarray:
        # end-effector position, p in (R,p) for element of se3
        ee_pos = se3.split(self.fk(thetas))[1]
        return ee_pos

    # ---------------------------------------------------------------
    # skeleton tracking and reachability 
    # ---------------------------------------------------------------
    def skeleton(self, thetas : np.ndarray) -> np.ndarray:
        # nb of screws = nb of vertices in skeleton
        """  check again  ! """
        thetas = self._theta(thetas)
        G = [np.eye(4)]                
        T = np.eye(4)
        # G[k] = exp(th_1 * S_1) * ... * exp(th_k * S_k)
        for S, th in zip(self.screws, thetas):
            T = T @ se3.exp(th * S)
            G.append(T.copy())
        dp = self.draw_points or [(0, (0, 0, 0)), (self.dof, se3.split(self.M)[1])] 
        pts = []
        # k is pt index + coord and p is deg of freedom + orientation
        for k, p in dp: 
            hom = np.array([p[0], p[1], p[2], 1.0])
            pts.append((G[k] @ hom)[:3])
        return np.array(pts)

    def max_reach(self) -> float:
        # a revolute serial chain can't extend past the sum of its link lengths
        pts = self.skeleton(np.zeros(self.dof))
        p_ee = self.ee_position(np.zeros(self.dof))
        pts = np.vstack([pts, p_ee])  # make sure tool tip is accounted for
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        return float(seg.sum())

    def _within_limits(self, theta : np.ndarray, tol=1e-9) -> bool:
        # takes joint_limits into account
        if self.joint_limits is None:
            return True
        lo, hi = self.joint_limits[:, 0], self.joint_limits[:, 1]
        return bool(np.all(theta >= lo - tol) and np.all(theta <= hi + tol))

    def is_reachable(self, T : np.ndarray, pos_tol=1e-4, restarts=8, return_solution=False):
        # decides reachability of pose, (not a proof but strong evidence)
        # first check radius of sum of lengths of manipulator 
        # then numerical IK check.
        T = np.asarray(T, float)
        if T.shape != (4, 4):
            raise ValueError("goal pose must be a 4x4 matrix in SE(3)")
        R, p = se3.split(T)
        if (not np.allclose(R.T @ R, np.eye(3), atol=1e-6) or not np.isclose(np.linalg.det(R), 1.0, atol=1e-6)):
            raise ValueError("rotation block of goal pose is not in SO(3)")

        p_base = self.skeleton(np.zeros(self.dof))[0]
        if np.linalg.norm(p - p_base) > self.max_reach() + pos_tol: # if outside of radius of sum of lengths
            return (False, None) if return_solution else False

        for attempt in range(restarts):
            theta, ok = self.numerical_ik(T, seed=attempt)
            if not ok:
                continue
            # re check fk on theta,
            # _finish_ik can move solution off the target, verifies accuracy and feasability
            if (se3.geodesic_distance(self.fk(theta), T) < pos_tol  
                    and self._within_limits(theta)):
                return (True, theta) if return_solution else True
            
        return (False, None) if return_solution else False

    

    # ---------------------------------------------------------------
    # Jacobians (velocities)
    # ---------------------------------------------------------------
    def space_jacobian(self, thetas : np.ndarray) -> np.ndarray:
        # V = J*theta' is the linear+angular velocity (theta' joint velocity vector)
        # i-th column in J is the i-th screw axis transported through the joints ahead of it 
        # J_i = Ad_{exp([S1]t1) ... exp([S_{i-1}]t_{i-1})} S_i
        thetas = self._theta(thetas)
        J = np.zeros((6, self.dof))
        T = np.eye(4) 
        for i, (S, th) in enumerate(zip(self.screws, thetas)):  # "analytic" construction
            J[:, i] = se3.adjoint(T) @ S
            T = T @ se3.exp(th * S)
        return J
 
    def body_jacobian(self, thetas : np.ndarray) -> np.ndarray:
        # Ad_{T^{-1}} * space jacobian 
        T = self.fk(thetas) 
        return se3.adjoint(se3.inverse(T)) @ self.space_jacobian(thetas)

    # (self evaluation)
    def finite_diff_space_jacobian(self, thetas : np.ndarray, eps=1e-6) -> np.ndarray: 
        # second way of computing Jacobian
        # log( Tn * T0^{-1} ) / eps  +  small wiggle
        thetas = self._theta(thetas)
        T = self.fk(thetas)  
        J = np.zeros((6, self.dof))
        for i in range(self.dof):
            dth = np.zeros(self.dof)
            dth[i] = eps
            J[:, i] = se3.log(self.fk(thetas + dth) @ se3.inverse(T)) / eps  # infinitessimal movement 
        return J

    def validate(self, trials=200, seed=0, tol=1e-6) -> float:
        # analytic vs finite difference Jacobian, returns max_error (validating se3.exp, se3.log and se3.adjoint)
        rng = np.random.default_rng(seed)
        max_err = 0.0
        for _ in range(trials):
            th = rng.uniform(-np.pi, np.pi, self.dof)
            err = float(np.max(np.abs(self.space_jacobian(th) - self.finite_diff_space_jacobian(th)))) # comparison
            max_err = max(max_err, err)
        assert max_err < tol, f"{self.name}: Jacobian mismatch {max_err:.2e} > {tol:.0e}"
        return max_err


    # ---------------------------------------------------------------
    # inverse kinematics (EE pose => thetas)
    # ---------------------------------------------------------------
    # analytic inverse kinematics
    def _wrap(self, a : float) -> float:
        # wrap to [-pi, pi]
        return (a + np.pi) % (2 * np.pi) - np.pi

    def _wrist_zyx(self, Rw :np.ndarray) -> list:
            # computes (alpha, beta, gamma) from Rw = R_z(alpa) R_y(beta) R_x(gamma)
            # see matrix R in B.1 (Modern Robotics)
            if abs(Rw[2,0]) < 1- 1e-9: # check that sin(beta) != ±1
                beta = np.arctan2(-Rw[2,0], np.sqrt(Rw[0,0]**2 + Rw[1,0]**2))
                out = []
                for betas in (beta, -beta):
                    out.append([(np.arctan2(Rw[1,0], Rw[0,0]),
                                betas,
                                np.arctan2(Rw[2,1], Rw[2,2]))])
                    return out
            if Rw[2,0] > 0: # gimbal lock 
                return [(0.0, -np.pi/2, -np.arctan2(Rw[0,1], Rw[1,1]))]
            return [(0.0, np.pi/2, np.arctan2(Rw[0,1], Rw[1,1]))]
    
    def _wrist_zyz(self, Rw : np.ndarray) -> list:
        #
        # computes (alpha, beta, gamma) from Rw = R_z(alpha) R_y(beta) R_z(gamma)
        if abs(Rw[2,2]) < 1 - 1e-9: # if cos(beta) != ±1 then sin(beta) != 0     (recall numpy.array)
            beta = np.arctan2(np.sqrt(Rw[0,2]**2 + Rw[1,2]**2), Rw[2,2])
            out = []
            for betas in (beta, -beta):
                c = np.sign(np.sin(beta))
                out.append((np.arctan2(c * Rw[1,2], c * Rw[0,2]),
                          betas,
                          np.arctan2(c * Rw[2,1], -c * Rw[2,0]))) 
                return out
        if  Rw[2,2] > 0 : # gimbal lock (beta = 0)
            return [(0.0, 0.0, np.arctan2(Rw[1,0], Rw[0,0]))]
        return [(0,0, np.pi, np.arctan2(-Rw[1,0], -Rw[0,0]))]

    def puma_6r(self, T : np.ndarray, L1 : float, a2 : float, a3 : float, tool_offset=0.0) -> list: # with no shoulder offset
        # T is ee pose (in se3) from that we got a set of solutions for angles of joints
        # at most 2^3 solutions 
        # with a zyz_wrist (pinch and rotate)
        R, p = se3.split(T)
        R_M = se3.split(self.M)[0]
        p_wc = p - tool_offset * R[:, 2] # position wrist center
        px, py, pz = p_wc
        r_xy = np.hypot(px,py) # horizontal length
        s = pz - L1 # height above shoulder
        D = (r_xy**2 + s**2 - a2**2 - a3**2) / (2*a2*a3)   # law of cosines (for cos theta_3)
        if abs(D) > 1: # target out of reach
            return []
        base_angle = np.arctan2(px,py)
        solutions = []
        # 2 base solutions (left and right)
        for theta_1, r in [(base_angle, r_xy), (self._wrap(base_angle + np.pi), -r_xy)]:
            # 2 elbow solutions (up and down)
            for elbow in (+1.0, -1.0):
                theta_3 = np.arctan2(elbow * np.sqrt(1 - D**2), D)
                theta_2 = np.arctan2(r,s) - np.arctan2(a3 * np.sin(theta_3), a2 + a3 * np.cos(theta_3))
                # PoE equations
                R03 = se3.split(se3.exp(theta_1 * self.screws[0]) @ 
                                se3.exp(theta_2 * self.screws[1]) @
                                se3.exp(theta_3 * self.screws[2]))[0] 
                R36 = R03.T @ R @ R_M.T
                for theta_4, theta_5, theta_6 in self._wrist_zyz(R36) : # 2 wrist solutions (except gimbal lock)
                    solutions.append(np.array([theta_1, theta_2, theta_3, 
                                               theta_4, theta_5, theta_6]))
        return solutions
    
    # numerical inverse kinematics ----------------------------------------------

    def numerical_ik(self, T : np.ndarray, theta_init=None, max_iter=200,
                    eps_omega=1e-6, eps_v=1e-6, max_step=0.5, seed=0) -> list:
        # uses Newton-Raphson method for nonlinear root-finding (check 6.2.2)
        # follows step by step the algo for end-effector configuration in SE(3)
        T = np.asarray(T, float)
        if theta_init is None:
            theta = np.random.default_rng(seed).uniform(-np.pi, np.pi, self.dof)
        else:
            theta = self._theta(theta_init).astype(float).copy()

        for _ in range(max_iter):
            Vb = se3.log(se3.inverse(self.fk(theta)) @ T) # find body twist
            if np.linalg.norm(Vb[:3]) < eps_omega and np.linalg.norm(Vb[3:]) < eps_v: # while norms > epsilons
                return self._finish_ik(theta), True

            theta_next = np.linalg.pinv(self.body_jacobian(theta)) @ Vb # Moore-Pensrose psudo-inverse
            step = np.linalg.norm(theta_next)
            if step > max_step: # scalling 
                theta_next *= max_step/step
            theta = theta + theta_next
        return self._finish_ik(theta), False # maxed out iterations, probably not accurate

    def _finish_ik(self, theta : np.ndarray) -> np.ndarray:
        # check if joint angles are possible regarding joint_limits
        theta = self._wrap(theta)
        if self.joint_limits is not None:
            theta = np.clip(theta, self.joint_limits[:,0], self.joint_limits[:,1])
        return theta

    

    # ---------------------------------------------------------------
    # interpolation (! take thetas as argument, not poses !) 
    # ---------------------------------------------------------------
    def geodesic_distance(self, theta_a : float, theta_b : float, length_scale=1.0) -> float:
        # computes geodesic distance (see se3.py)
        geod_dist = se3.geodesic_distance(self.fk(theta_a), self.fk(theta_b), w_rot=length_scale ** 2, w_trans=1.0)
        return geod_dist

    def cartesian_trajectory(self, theta_a : float, theta_b : float, n = 20) -> np.ndarray:
        # computes screw interpolate (continuous evolution of pose, see se3)
        Ta , Tb = self.fk(theta_a), self.fk(theta_b)
        traj_pts = [se3.screw_interpolate(Ta, Tb, t) for t in np.linspace(0, 1, n)]
        return traj_pts
    

