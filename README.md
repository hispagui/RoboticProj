This repositery consists of two parts; the first one is a topological Motion Planner and the second one is a Rigid Body Motion tracker. 
The main documentation used for this project is the book _Modern Robotics, Mechanics, Planning, and Control_ by Kevin M. Lynch and Frank C. Park.

Both parts of this project rely on the group structures (and methods) of the Lie groups $SO(3)$ and $SE(3)$

## Topological Motion Planner

This part studies the geometry and topology of the configuration space, how the free space arround obstacles splits trajectories into distinct homotopy classes, and wetheer that topological structure actually helps a planner.

Here a configuration space is an orientation $\mathcal{C} = SO(3)$ or a full pose $\mathcal{C} = SE(3)$.
These spaces are curved and topologically non trivial, $SO(3) \cong \mathbb{RP}^3$ and $S^3$ sits above it as its universal cover, already $\pi_1(SO(3)) = \mathbb{Z}/2$.

MISSING : implementation of homotopy class of paths for graphs (useful for later comparisons)


## Rigid Body Motion

This section focuses on the kinematics of manipulators (with $n$ rotational joints), it computes _forward kinematics_, _inverse kinematics_ and _jacobians_.
