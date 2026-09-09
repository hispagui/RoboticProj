"""
homotopy.py : topological signatures of paths in R^2 minus a set of points (punctured plane)
    use homotopy theory to compute winding arround holes, and path homotopy classes in punctured plane.
    main functions are reduce_word and winding_signature, both compute going arround obstacles (in different ways)
    reduce_word :       homotopy idea, a word in the free group F_n (for n obstacles),
                        distinguishes paths that are homologous but not homotopic (obtained from ray-crossing construction).
    winding_signature : produces an n-vector of net windings per obstacles, abelianisation, it forgets the order of windings,
                        (cheaper than reduce_word)

                             
                        
    To add : A* planners and comparisons, and homotopy classes on graphs or cell-complexes

"""

import numpy as np


# HOMOTOPY (free groups in R^2)

class Paths:
    def __init__(self, top_space, path : list):
        self.top_space = top_space
        self.path = path

    # could be interesting to have paths objects and compare/add,sub them and make some sort of homotopy arithmetic

        

def _segment_crossings(p0 : list, p1 : list, obstacles : list) -> list:
    # crossings regarding vertical rays (like laser beams, with base in obstacles) when going from p0 to p1
    # obstacles = [[a1,b1], [a2,b2], [a3,b3], ...]
    # return [(t, +/-(k+1)), ...] in the order they occur, + is crossing left to right and - right to left
    x0, y0 = p0
    x1, y1 = p1
    out = []
    if x0 == x1:                       # p0-p1 is a vertical line (no crossing)
        return out
    for k, (xk, yk) in enumerate(obstacles):
        t = (xk - x0) / (x1 - x0)      # "when" the segment meets the line x = x_k 
        if not (0.0 <= t < 1.0):
            continue
        y_cross = y0 + t * (y1 - y0)
        if y_cross > yk:                    # above the base of "laser" => crosses the upward ray
            sign = 1 if x1 > x0 else -1     # 1 if crossed from the left and -1 if crossed from the right
            out.append((t, sign * (k + 1))) # time of crossing with signed index 
    out.sort(key=lambda c: c[0])
    return [label for _, label in out]

def _reduce(word : list) -> tuple:
    # cancels back-and-forth winding arround nothing (should disapear)
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)

def reduced_word(path : list, obstacles : list) -> tuple:
    # path as a list of points
    # returns tuple of signed int : 1 "arround obstacles[0] one way", -2 "arround obstacles[1] another way"
    # "constructs" homotopy classes
    path = np.asarray(path, dtype=float)
    obstacles = np.asarray(obstacles, dtype=float).reshape(-1, 2)
    letters = []
    for a, b in zip(path[:-1], path[1:]):        # iterates consecutive pair of points
        letters.extend(_segment_crossings(a, b, obstacles))
    return _reduce(letters)

def are_homotopic(path_a : list, path_b : list, obstacles: list) -> bool:
    # returns True if the two paths belong to the same homotopy class (need to have same endpoints)
    return reduced_word(path_a, obstacles) == reduced_word(path_b, obstacles)

def word_to_str(word : str, names = None) -> str:
    # returns usual notation for computation in free groups : (1, -2, 1) -> 'a b^-1 a
    if not word:
        return "e"                      # identity / trivial class
    parts = []
    for letter in word:
        k = abs(letter) - 1
        name = names[k] if names else chr(ord('a') + k)
        parts.append(name if letter > 0 else name + "^-1")
    return " ".join(parts)


# "HOMOLOGY" (abelian)

def winding_signature(path : list, obstacles : list) -> np.ndarray:
    # more topological / less geometrical reduce_word
    # H_1 = pi_1^{ab} (first homology group is realised as the abelianisation of the fundamental group)
    path = np.asarray(path, dtype=float)
    obstacles = np.asarray(obstacles, dtype=float).reshape(-1, 2)
    out = np.zeros(len(obstacles))
    for k, (xk, yk) in enumerate(obstacles):
        total = 0.0
        for a, b in zip(path[:-1], path[1:]):      # iterates consecutive pair of points
            a0 = np.arctan2(a[1] - yk, a[0] - xk)
            a1 = np.arctan2(b[1] - yk, b[0] - xk)
            d = a1 - a0
            d = (d + np.pi) % (2 * np.pi) - np.pi   # wrap into (-pi, pi]
            total += d
        out[k] = total / (2 * np.pi)     
    return out                                      # total of 

def are_homologous(path_a, path_b, obstacles, atol=1e-6) -> bool:
    # return True if the two paths have the same winding about every obstacle
    return np.allclose(winding_signature(path_a, obstacles),
                       winding_signature(path_b, obstacles), atol=atol)

def abelianize(word : str, n_obstacles : int) -> np.ndarray:
    # counts net signed ray-crossings
    # maps free group onto Z_n 
    # two words with the same abelianisation are homologous : in the same homology class (may be in different homotopy class)
    v = np.zeros(n_obstacles, dtype=int)
    for letter in word:
        v[abs(letter) - 1] += 1 if letter > 0 else -1
    return v

# for close loops, abelianize and winding_signature agree up to a sign

