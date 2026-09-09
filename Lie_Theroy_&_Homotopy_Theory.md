

## Lie Theory and Riemannian Geometry


### SO(3) and SE(3) are Lie groups

__Definition__ :
A _Lie group_ is a subset $G$ of $\mathbb{R}^n$ such that $G$ is a group and a manifold in $\mathbb{R}^n$ and both the group operation $\cdot : G\rightarrow G$ and the inverse operation $.^{-1}: G\rightarrow G$ are smooth functions.

__Definition__ :
A _Lie algebra_ is an algebra $A$ togther with a _Lie bracket_ operator $\[.,.\] : A \times A \rightarrow A$.
An important fact is that the Lie algebra $A$ associated to a Lie group $M$ happens to be the tangent space
at the identity element $1$, that is: $A = T_1M$.

__Deffinition__ :
Lie bracket.

__Examples__ :
The $3$D rotation group $SO(3)$ is a Lie group (homeomorphic to $\mathbb{RP}^3$).
The group of rigid motions in $3$D or the special Euclidean group $SE(3)$ is also a Lie group (homeomorphic to $\mathbb{RP}^3 \times \mathbb{R}^3$, not as a group tho !).
The


__Representations of coordinates (in SO(3))__ :
- With 3 angles (_Euler angles_) $\alpha$, $\beta$ and $\gamma$ $\in \[ -\pi, \pi )$, rotation matrices can be composed : $R(\alpha, \beta, \gamma) = R_3(\gamma)R_2(\beta)R_1(\alpha) := f(r)$ 
    where $r = [\alpha, \beta, \gamma]$.

- By _exponential representation_ let $\omega \in \mathbb{R}^3$ be a vector, a rotation in 3D can be expressed by a rotation axis $\omega$ and a rotation angle about that axis $||\omega||$.
The matrix exponential of $\omega$ yields the $3\times 3$ skew-symmetric matrix in $SO(3)$.

- _Axis_, let $q = (\omega, \alpha)\in \mathbb{R}^4$ with $||\omega|| = 1$, its basically the same idea as above, a rotation of angle of $\alpha$ arround the $\omega$-axis.
And Rodrigues' rotation formula maps rotation vector to rotation matrix in $SO(3)$.

- _Quaternions_, let $q = (w,x,y,z)\in S^3$ (so with $||q||=1$), we set $w = \cos\frac{\theta}{2}$ and $(x,y,z) = \alpha \sin\frac{\theta}{2}$, then as above, we have a rotation of angle $\theta$ arround the axis $\alpha$.



### Riemannian Geometry
__Definition__ : 
Let $M$ be an $n$-dimensional manifold.
A _Riemannian metric_ on $M$ is a smoothly varying positive-definite inner-product $g_p$ on each tangent space $T_pM$.
We say that the pair $(M,g)$ is a _Riemannian manifold_.

We now define what is the analogue of "straight lines" in a Riemannian manifold $(M,g)$ :

__Definition__ : 
Let $\gamma : I \rightarrow M$ be a regular path.
Consider $D_t: \mathcal{X}(\gamma) \rightarrow \mathcal{X}(\gamma)$ the _covariant derivative along_ $\gamma$.
The velocity of $\gamma$ (its derivation $t \mapsto \gamma'(t)$) defines a vector field along $\gamma$, we call $D_t\gamma'$ the _acceleration_ of $\gamma$. 
When $D_t\gamma' = 0$, we say that $\gamma$ is a _geodesic_.

The idea is now to use geodesics to "explore" $M$, imagine sending probes with zero acceleration with all different velocities, after one second they report back their positions, giving you a "map" of $M$. Formally :

__Definition__ : 
Let $p\in M$ and $v\in T_pM$. We write $\gamma_v:I\rightarrow M$ for the (maximal) geodesic with $\gamma(0) = p$ and $\gamma'(0) = v$. We denote by $U_p \subset T_pM$ the set of those $v$ for which $\gamma_v(1)$ is well-defined. 
The _exponential map_ at $p$, written $exp_p:U_p\rightarrow M$, is defined by $exp_p(v) = \gamma_v(1)$.

__Remark__ : 
In the context of Lie groups, the exponential map is map from the Lie algebra (tangent space) to the corresponding Lie group (manifold).
There is an "inverse" computation that we call the _logarithm map_ which sends an element from the Lie group (manifold) to the Lie algebra (tangent space).
These are the exact concepts we used when writting the interpolation functions found in se3.py and so3.py











## Homotopy and Homology theory

__Definition__ : A _homotopy_ between two continuous functions $f$ and $g$ from $X$ to $Y$ is a continuous function $H : X \times [0,1]  \rightarrow Y$ such that $H(x,0) = f(x)$ and $H(x,1) = g(x)$ for all $x\in X$.

__Definition__ : Two paths $\gamma_1$ and $\gamma_2$ are said to be _homotopic_ if there exists a homotopy taking $\gamma_1$ to $\gamma_2$.

__Definition__ : The _fundamental group_ of a pointed topological space $(X,x)$ denoted $\pi_1(X,x)$ is the group of equivalence classes under homotopy of the loops based at $x$ in $X$.

__Example__ : The fundamental group of $\mathbb{R}^2$ is trivial, $\pi_1(\mathbb{R}^2 - \{p\}) = \mathbb{Z}$ and $\pi_1(\mathbb{R}^2 - \left\{p_1, p_2, ..., p_n \right\}) = F_n$ the free group on $n$ generators.


__Definition__ : Let $X$ be a topological space and $(C_{\star}, d_{\star})$ a _chain complexe_, with $d_n : C_n \rightarrow C_{n-1}, each $C_n$ is an abelian group and $d_{n-1} \circ d_n = 0 \forall n$.
We define $B_n := \im d_{n+1}$ and $Z_n := \ker\: d_n$ to be respectivelly the group of _boundaries_ and the group of _cycles_.
Finally $H_n(X) := Z_n /B_n$ is the n$^{th}$ _homology group_.


__Theorem (Hurewicz)__ : For $X$ connected, there exists a homomorphism $h : \pi_1(X) \rightarrow H_1(X)$ which is surjective and $\ker \: h = [\pi_1(X), \pi_1(X)]$ the commutator subgroup. This implies that the abelianisation of the first fundamental group is isomorphic to the first homology group, $\pi_1(X)^{ab} \cong H_1(X)$.