"""
animations.py -- animations for : - SO(3)/SE(3) interpolation (screwed and decoupeled)
                                  - inverse/forward kinematics tracking 
"""

import os
import numpy as np
import so3
import se3
import kinematics
import matplotlib
matplotlib.use("Agg") # saves to file
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter



class Poses:
    # Class to animate interpolation (geodesic) from one pose T0 to another pose T1
    def __init__(self, T1_pos : np.ndarray, T1_ori : np.ndarray, T0_pos = [0.0, 0.0, 0.0], T0_ori = np.eye(3), N=60):
        self.T0 = se3.make(T0_ori, T0_pos) # starting pose
        self.T1 = se3.make(T1_ori, T1_pos) # end pose
        self.N = N
        self.geod_dist = round(se3.geodesic_distance(self.T0, self.T1, 1.0), 4) 
        self.half_time_pos =  np.round(se3.split(se3.screw_interpolate(self.T0, self.T1, 0.5))[1], 4)
        self.screws, self.decpl, self.screw_pos, self.decpl_pos = self._make_screws(self.N)
        self.fig, self.ax  = self._geometry()
        self.screw_triad, self.decpl_triad = self.triads()
    def __repr__(self):
        return f"Start : {self.T0} \n Finish : {self.T1} \n Geodesic distance : {self.geod_dist} \n Half time position : {self.half_time_pos}"

    def _make_screws(self):
        ts = np.linspace(0.0, 1.0, self.N) # number of frames in animation
        screws = [se3.screw_interpolate(self.T0, self.T1, t) for t in ts] # list of screws
        decpl =  [se3.decoupled_interpolate(self.T0, self.T1, t) for t in ts] # list of decoupled screws
        screw_pos = np.array([se3.split(T)[1] for T in screws])  
        decpl_pos = np.array([se3.split(T)[1] for T in decpl])   
        return screws, decpl, screw_pos, decpl_pos

    def _geometry(self):
        # set up the 3d space
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d") 
        # set up all plot stuff
        allpts = np.vstack([self.screw_pos, self.decpl_pos])
        lo, hi = allpts.min(0) - 1.2, allpts.max(0) + 1.2
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect((hi - lo))
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.view_init(elev=22, azim=-60) 
        # faint full reference paths (so the arc vs chord shapes read even at rest)
        ax.plot(self.screw_pos[:, 0], self.screw_pos[:, 1], self.screw_pos[:, 2], color="#D85A30", lw=1, alpha=0.25)
        ax.plot(self.decpl_pos[:, 0], self.decpl_pos[:, 1], self.decpl_pos[:, 2], color="#1D9E75", lw=1, alpha=0.25, ls="--")
        # endpoint markers
        ax.scatter(*self.screw_pos[0], color="k", s=40)
        ax.scatter(*self.screw_pos[-1], color="k", s=40)
        ax.text(*self.screw_pos[0], "  T0", fontsize=11)
        ax.text(*self.screw_pos[-1], "  T1", fontsize=11)
        ax.legend(loc="upper left")
        ax.set_title("Moving a rigid body from T0 to T1: screw vs decoupled")
        return (fig, ax) 

    def triads(self):
        screw_triad = self.draw_triad("-", 1.0)     # solid  = screw motion
        decpl_triad = self.draw_triad("--", 0.65)   # dashed = decoupled motion
        return (screw_triad, decpl_triad)

    def draw_triad(self, style, alpha, AXIS_COLORS = ("#D14520", "#3B8BD4", "#1D9E75")):
        # draws 3 Line3D for a body frame's (x,y,z axes)
        return [self.ax.plot([], [], [], style, color=c, lw=2.2, alpha=alpha)[0] for c in AXIS_COLORS]

    def set_triad(self, handles, T, AXIS_LEN = 0.9):
        R, p = se3.split(T)
        for i, h in enumerate(handles):
            tip = p + AXIS_LEN * R[:, i]
            h.set_data([p[0], tip[0]], [p[1], tip[1]])
            h.set_3d_properties([p[2], tip[2]])

    def update(self, frame):
        # the function called in FuncAnimation
        # trails from origin
        screw_trail, = self.ax.plot([], [], [], color="#D85A30", lw=2.5, label = "screw (geodesic)")
        decpl_trail, = self.ax.plot([], [], [], color="#1D9E75", lw=2.5, ls="--", label = "decoupled (slerp)")
        self.set_triad(self.screw_triad, self.screws[frame])
        self.set_triad(self.decpl_triad, self.decpl[frame])
        screw_trail.set_data(self.screw_pos[:frame + 1, 0], self.screw_pos[:frame + 1, 1])
        screw_trail.set_3d_properties(self.screw_pos[:frame + 1, 2])
        decpl_trail.set_data(self.decpl_pos[:frame + 1, 0], self.decpl_pos[:frame + 1, 1])
        decpl_trail.set_3d_properties(self.decpl_pos[:frame + 1, 2])
        self.ax.view_init(elev=22, azim=-60 + 0.35 * frame)   # slow orbit for depth
        return (*self.screw_triad, *self.decpl_triad, screw_trail, decpl_trail)

    def anim_pose_inerpolation(self):
        anim = FuncAnimation(self.fig, self.update, frames=self.N, interval=55, blit=False)  # timed animation
        out = os.path.join(os.path.dirname(__file__), "pose_interpolation.gif") # output gif file
        anim.save(out, writer=PillowWriter(fps=20))
        print("wrote in", out)












class Kins:
    # Class used to animate joints of a Manipulator going from a starting pose to an end pose
    def __init__(self, manipulator, goal, N = 120, starting_pose = None):
        self.manipulator = manipulator # in starting position
        self.N = N
        self.goal = goal # element in SE(3)
        self.start = self._initial_pose() #element in SE(3)
        self.geod_dist = round(se3.geodesic_distance(self.start, self.goal, 1.0), 4)
        if starting_pose is not None:
            self.start = starting_pose
            # and make call (still have to figure thetas)
        # check if goal pose is reachable
        reachable = self.manipulator.is_reachable(self.goal)
        if not reachable:
            raise ValueError(f"{self.manipulator.name}: pose is not reachable")
        self.screw_path = self._make_screws()
        self.skeletons, self.tips, self.hi, self.lo = self.skeleton()
        self.fig, self.ax, self.arm_line = self._geometry()
        
    def __repr__(self):
        return f"{self.manipulator.__repr__()}, gedodesic distance between start and finish is {self.geod_dist}"
    
    def _initial_pose(self):
        # returns initial ee pose
        thetas = np.zeros(self.manipulator.dof)
        return self.manipulator.fk(thetas) # an element in SE(3)

    def _geometry(self):
        fig = plt.figure(figsize=(6.5, 6.5))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_xlim(self.lo[0], self.hi[0]); ax.set_ylim(self.lo[1], self.hi[1]); ax.set_zlim(0, self.hi[2])
        ax.set_box_aspect((self.hi[0] - self.lo[0], self.hi[1] - self.lo[1], self.hi[2]))
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.view_init(elev=18, azim=-70)
        ax.set_title(f"IK tracking of {self.manipulator.dof}R arm from start to goal")
        arm_line, = ax.plot([], [], [], "-o", color="#3B6FD4", lw=4, markersize=6,
                            markerfacecolor="#D85A30", zorder=5)
        return fig, ax, arm_line

    def _make_screws(self):
        ts = np.linspace(0.0, 1.0, self.N) # number of frames in animation
        screws = [se3.screw_interpolate(self.start, self.goal, t) for t in ts] # list of screws
        return screws

    def skeleton(self):
        joint_path = []
        theta_prev = None

        for i,S in enumerate(self.screw_path):
            if theta_prev is None:
                sol, ok = self.manipulator.numerical_ik(S, theta_init=theta_prev)
                if not ok: 
                    for s in range(6): # try 6 more times to find good enough thetas with Newton-Raphson
                        sol, ok = self.manipulator.numerical_ik(S, theta_init=theta_prev, seed=(s+1))
                        if ok:
                            break
            else: # to make numerical IK smooth, better to use previous theta computed as initial theta (finds closer solution)
                sol, ok = self.manipulator.numerical_ik(S, theta_init=theta_prev)
                if not ok: 
                    for s in range(6): # try 6 more times to find good enough thetas with Newton-Raphson
                        sol, ok = self.manipulator.numerical_ik(S, theta_init=theta_prev, seed=(s+1))
                        if ok:
                            break
                else:
                    d = sol - theta_prev
                    d = (d + np.pi) % (2 * np.pi) - np.pi
                    sol = theta_prev + d

            if not ok:
                print(f"warning, frame {i} did not converge")
            print("iter", i)

            theta_prev = sol.copy()
            joint_path.append(sol.copy())

        skeletons = [self.manipulator.skeleton(th) for th in joint_path]
        tips = np.array([sk[-1] for sk in skeletons]) # end effector of skeleton
        allpts = np.vstack(skeletons)
        lo, hi = allpts.min(0) - 0.15, allpts.max(0) + 0.15
        return skeletons, tips, hi, lo

    def update(self, frame):
        sk = self.skeletons[frame]
        self.arm_line.set_data(sk[:, 0], sk[:, 1]); self.arm_line.set_3d_properties(sk[:, 2])
        #tip_trail.set_data(self.tips[:frame + 1, 0], self.tips[:frame + 1, 1])
        #tip_trail.set_3d_properties(self.tips[:frame + 1, 2])
        return self.arm_line

    def anim_inverse_kinematics(self):
        anim = FuncAnimation(self.fig, self.update, frames=self.N, interval=40, blit=False)
        out = os.path.join(os.path.dirname(__file__), "inverse_kinematics_tracking.gif")
        anim.save(out, writer=PillowWriter(fps=25))
        print("wrote", out)














# examples
'''
# 1.
T0_pos = [0, 0, 0]
T0_ori = np.eye(3)
#T1_ori = so3.exp(np.array([0.3, 1.0, 0.5]) / np.linalg.norm([0.3, 1.0, 0.5]) * 2.4)
T1_ori = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
T1_pos = [3.0, 1.5, 2.0]
my_poses = Poses(T1_pos, T1_ori, T0_pos, T0_ori)
print(my_poses)
my_poses.anim_pose_inerpolation()
'''


def practice_6dof():
    # test encoding of 6r arm (no offset)
    L1, L2, L3, L4 = 0.40, 0.40, 0.20, 0.10 # lengts
    joints = [
    ([0, 0, 1], [0, 0, 0]),                       # joint 1 at the base
    ([0, 1, 0], [0, 0, L1]),                       # joint 2, L1 up
    ([0, 1, 0], [0, 0, L1 + L2]),                  # joint 3, L1+L2 up
    ([0, 0, 1], [0, 0, L1 + L2 + L3]),             # wrist, L1+L2+L3 up
    ([0, 1, 0], [0, 0, L1 + L2 + L3]),
    ([0, 0, 1], [0, 0, L1 + L2 + L3]),]
    M = se3.make(np.eye(3), [0, 0, L1 + L2 + L3 + L4])  # tip = total length
    # skeleton for drawing: base, shoulder, elbow, wrist, tool tip.
    draw_points = [
        (0, [0, 0, 0]),                        # base (fixed)
        (1, [0, 0, L1]),                       # after joint 1
        (2, [0, 0, L1 + L2]),                  # after joint 2
        (3, [0, 0, L1 + L2 + L3])]             # wrist (after joint 3)
    return kinematics.Manipulator.from_revolute(joints, M, name="didactic_6dof", draw_points=draw_points)

arm = practice_6dof()
print(arm)
T_goal = se3.make(np.eye(3), [0.3, 0.3, 0.7]) # this pose is reachable

animate_arm = Kins(arm, T_goal, N=60)
animate_arm.anim_inverse_kinematics()
