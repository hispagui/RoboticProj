This repositery consists of two parts; the first one is a topological Motion Planner and the second is a Rigid Body Motion tracker. 
The main documentation used for this project is the book _Modern Robotics, Mechanics, Planning, and Control_ by Kevin M. Lynch and Frank C. Park.

Both parts of this project rely on the group structures (and methods) for the Lie groups $SO(3)$ and $SE(3)$

# Topological Motion Planner

    This part studies the geometry and topology of the configuration space, how the free space arround obstacles splits trajectories into distinct homotopy classes, and wetheer that topological structure actually helps a planner.

# Rigid Body Motion

    This section focuses on the kinematics of manipulators (with $n$ rotational joints), it computes _forward kinematics_, _inverse kinematics_ and jacobians  
    
 it detects homotopy classes of paths between two points 





In this project, we study the geometry and the topology of the configuration space : how orientations and poses live on curved manifolds,

how the free space around obstacles splits trajectories into distinct homotopy classes, and whether that topological structure actually helps a planner.

In this 



## Thread A - The topology of the configuration space itself (no obstacles)
    Here a configuration is an orientation $C = SO(3)$ or a full pose $C = SE(3)$.
    These spaces are curved and topologically non trivial, $SO(3) \cong \mathbb{RP}^3$ and $S^3$ sits above it as its universal cover, already $\pi_1(SO(3)) = \mathbb{Z}/2$.

    In this thread we represent interpolation between poses correctly *(see so3.py and se3.py)*

## Thread B - The topology of the free space (with obstacles)
    This is coming up next