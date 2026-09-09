These notes are made using the incredibly weel written _Modern Robotics : Mechanics, Plannning and Control_ by Kevin M. Lynch and Frank C. Park.
(freely available on ...)

Use of $SE(3)$ and $SO(3)$ Lie groups to encode some of the key concepts of rotations and rigid-body motions in 3D.

__Notations :__ Suppose we have a robot consisting of $n$ links. 

$\cdot$ The link reference frames are respectively labeled $\{1\}, \{2\}, ..., \{n\}$, with $\{0\} = \{b\}$ being the base frame and $\{e\}$ being the end-effector's frame attached to the last link. The link lengths are denoted by $L_1, L_2, ..., L_n$, and the joint angles are given by $\theta_1, \theta_2, ..., \theta_n$.

$\cdot$ The matrices $T_{bm} \in SE(3)$ denotes the configuration of the (body) frame $\{m\}$ as seen from the base frame $\{b\}$. 
The pose $T$ (often $T = T_{0n}$) of the end-effector frame is a function of all of $L_i$ and $\theta_i$ for $i \in [1,n]$. 
We denote $M \in SE(3)$ the configuration of the end-effector relative to the base frame when the robot is in its zero position (all the joint angles have value zero).

$\cdot$ Links are often described as _screws_ and denoted $\mathcal{S} = \{\omega, v\} \in \mathbb{R}^6$ with $\omega$ a position and $v$ a unit vector in direction of the axis.

$\cdot$ We will denote the _skew-symmetric matrix_ associated with $\omega\in\mathbb{R}^3$ and $\mathcal{S}$ respectively by $[w]$ and $[\mathcal{S}]$.



In the following we will use the exponential coordinate representation : 

```math
\exp :  [\mathcal{S}]\theta \longrightarrow T \in SE(3)
```

```math
\log : T \longrightarrow [\mathcal{S}]\theta \in se(3).
```

And we have $\exp([\mathcal{S}]\theta) = e^{[\mathcal{S}]\theta} = I + [\mathcal{S}]\theta + [\mathcal{S}]^2\frac{\theta^2}{2!} + [\mathcal{S}]^3\frac{\theta^3}{3!} + ...$


# Forward Kinematics
Answers the question, what is the end effector's pose $T$ given the joint angles $\theta_i$ and the link lengths $L_i$ ?

Idea : One can combine consecutive transformations, $T_{01} T_{12} T_{23} = T_{03}$

## The Product of Exponentials Formula
There exist two formulations of this formula; in the base frame and in the end-effector frame.

In the base frame, we call this equation the _space form_ :

$\cdot$ Suppose joint $n$ is displaced by to some joint value $\theta_n$. The end-effector frame $M$ then undergoes a displacement of the form 
```math 
T = e^{[\mathcal{S}_n]\theta_n}M$, where $T \in SE(3)
``` 
is the new configuration of the end-effector frame.

$\cdot$ If we assume that the joint $n-1$ is also allowed to move, then the end-effector undergoes a displacement of the form 
```math
T = e^{[\mathcal{S}_{n-1}]\theta_{n-1}}(e^{[\mathcal{S}_n]\theta_n}M)
```

$\cdot$ Continuing this reasonning and allowing all joints to move, it follows that 
```math
T(\theta) =  e^{[\mathcal{S}_{1}]\theta_{1}} ... e^{[\mathcal{S}_{n-1}]\theta_{n-1}} e^{[\mathcal{S}_n]\theta_n}M.
```


In the end-effector frame, this representation is called the _body form_:
Using $e^{M^{-1}PM} = M^{-1}e^P M$, the above formula yields
```math
\begin{aligned}
T(\theta)
    &= e^{[\mathcal{S}_{1}]\theta_{1}} ... e^{[\mathcal{S}_n]\theta_n}M \\
    &= e^{[\mathcal{S}_{1}]\theta_{1}} ... Me^{M^{-1}[\mathcal{S}_n]M\theta_n} \\
    &= e^{[\mathcal{S}_{1}]\theta_{1}} ... Me^{M^{-1}[\mathcal{S}_{n-1}]M\theta_{n-1}} e^{M^{-1}[\mathcal{S}_n]M\theta_n}\\
    &= M e^{M^{-1}[\mathcal{S}_{1}]M\theta_{1}} ... e^{M^{-1}[\mathcal{S}_n]M\theta_n}\\
    &= M e^{[\mathcal{B}_{1}]\theta_{1}} ... e^{[\mathcal{B}_{n}]\theta_{n}}
\end{aligned}
```

where each 
```math 
[\mathcal{B}_i] = [Ad_{M^{-1}}]\mathcal{S}_i
``` 
represents the joint axes as screws axes in the end-effector frame.



# Inverse Kinematics
Given the end-effector's pose $X \in SE(3)$, what are the valid joints angles that such that $T(\theta) = X$.
## Analytic Inverse Kinematics
For the famous PUMA 6R arm (see _Modern Robotics_, 6.1.1)

[...]

## Numerical Inverse Kinematics
Understand Newton-Raphson method; 
make an initial guess, tweak the guess through iterations and get closer to actual solution.
The tweaking makes use of the (pseudo-)inverse of the Jacobian (so we have to assume that our fk is differentiable)

[...]

see _Modern Robotics_, 6.2.2

## Jacobians
We want to find the twist $V \in se(3)$ of the end-effector.
We obtain the joint velocities from joint positions.
And the twist is obtained through the Jacobian; $V_{ee} = J(q)\dot{q}$ for $q$ the coordinates of the joints.

[...]

see _Modern Robotics_ chapter 5


