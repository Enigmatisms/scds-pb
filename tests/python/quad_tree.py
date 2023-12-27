import sys
import time
sys.path.append("../../")
import numpy as np
import matplotlib.pyplot as plt
from build.tests import smt

def full_test():
    n_nearest  = 8
    point_size = 2
    num_pts = 2000000           # use 7 to visualize a small-scale test
    radius  = 0.01

    query_pt = np.float32([0.23, 0.23])
    
    range = np.float32([
        [0.5, 0.5],
        [0.5, 0.5]
    ])
    quad_tree = smt.QuadTreef(range, 32, 4)
    pts = np.random.rand(num_pts, 2).astype(np.float32)

    start_time = time.time()
    quad_tree.insert(pts)
    print(f"Inserting {num_pts} points took: {(time.time() - start_time) * 1000.:.3f} ms")
    
    start_time = time.time()
    nns = quad_tree.search_nn(query_pt, n_nearest, radius)
    print(f"Tree query point took: {(time.time() - start_time) * 1000.:.3f} ms")
    print(f"{nns.shape[0]} nearest neighbor found.")
    
    start_time = time.time()
    nns_bf = quad_tree.search_nn_bf(query_pt, n_nearest, radius)
    print(f"Brute force took: {(time.time() - start_time) * 1000.:.3f} ms")
    print(f"{nns_bf.shape[0]} nearest neighbor found.")
    
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.scatter(pts[:, 0], pts[:, 1], s = point_size, color = 'k', label = 'Quad Tree Points')
    ax.scatter([query_pt[0]], [query_pt[1]], s = point_size, color = 'b', label = 'Query')

    square = plt.Rectangle(query_pt - radius, radius * 2, radius * 2, color='blue', fill=True, label='Square', alpha = 0.2)
    ax.add_patch(square)

    ax.scatter(nns[:, 0], nns[:, 1], s = point_size << 1, color = 'r', label = 'Nearest neighbors')
    ax.scatter(nns_bf[:, 0], nns_bf[:, 1], s = point_size, color = 'g', label = 'Nearest neighbors (BF)')

    plt.xticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
    plt.yticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])

    ax.grid(axis = 'both')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    np.random.seed(2)
    full_test()