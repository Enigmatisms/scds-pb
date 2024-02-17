import sys
import tqdm
import time
sys.path.append("../../")
import numpy as np
import matplotlib.pyplot as plt
from build.tests import kdt, smt

"""
200k points, depth 32, 4 nn for leaf nodes

unordered_set implementation
Insert time: 364.75483 ms
Tree search time: 0.01971 ms
BF search time: 10.15662 ms

vector implementation
Insert time: 83.24949 ms
Tree search time: 0.01561 ms
BF search time: 0.19610 ms
"""

def repeated_test():
    max_depth     = 10
    n_nearest     = 8
    leaf_node_max = 64
    num_pts       = 600000           # use 7 to visualize a small-scale test
    radius        = 0.01
    query_pt = np.float32([0.23, 0.23])
    repeat_epochs = 30
    
    tree_range = np.float32([
        [0.5, 0.5],
        [0.5, 0.5]
    ])
    
    time_sums = np.float64([0, 0, 0])
    for _ in tqdm.tqdm(range(repeat_epochs)):
        kdtree = kdt.KDTree2f(tree_range, max_depth, leaf_node_max)
        pts = np.random.rand(num_pts, 2).astype(np.float32)
        
        start_time = time.time()
        kdtree.insert(pts)
        time_sums[0] += (time.time() - start_time) * 1000
        
        start_time = time.time()
        _ = kdtree.search_nn(query_pt, n_nearest, radius)
        time_sums[1] += (time.time() - start_time) * 1000

    time_sums /= repeat_epochs
    print(f"Insert time: {time_sums[0]:.5f} ms")
    print(f"Tree search time: {time_sums[1]:.5f} ms")
    
def draw_rectangle(ax, node_range, color='blue', alpha=0.5, line_width=1.5, label=None):
    x, y, hw, hh = node_range
    rectangle = plt.Rectangle((x - hw, y - hh), hw * 2, hh * 2, linewidth=line_width, 
                    edgecolor=color, facecolor=color, alpha=alpha, fill = False, label=label)
    ax.add_patch(rectangle)

def full_test():
    n_nearest  = 8
    point_size = 4
    point_per_node = 1    # large scale: 4, small scale: 1
    num_pts = 8192      # use 7 to visualize a small-scale test
    radius  = 0.05         # large scale: 0.01

    query_pt = np.float32([0.23, 0.23])
    
    range = np.float32([
        [0.5, 0.5],
        [0.5, 0.5]
    ])
    
    pts = np.random.rand(num_pts, 2).astype(np.float32)

    start_time = time.time()
    quad_tree = smt.QuadTreef(range, 32, point_per_node)
    quad_tree.insert(pts)
    print(f"Quad-Tree: Inserting {num_pts} points took: {(time.time() - start_time) * 1000.:.3f} ms. Tree size = {quad_tree.size()}")

    kdtree = kdt.KDTree2f(range, 32, point_per_node)
    start_time = time.time()
    kdtree.build_tree(pts)
    print(f"K-D tree building tree with {num_pts} point(s) took: {(time.time() - start_time) * 1000.:.3f} ms. Tree size = {kdtree.size()}")
    
    start_time = time.time()
    nns = kdtree.search_nn(query_pt, n_nearest, radius)
    print(f"(KDTree) Tree query point took: {(time.time() - start_time) * 1000.:.3f} ms")
    print(f"{nns.shape[0]} nearest neighbor found, depth: {kdtree.depth()}")

    start_time = time.time()
    nns_bf = quad_tree.search_nn(query_pt, n_nearest, radius)
    print(f"(QuadTree) Tree query point took: {(time.time() - start_time) * 1000.:.3f} ms")
    print(f"{nns_bf.shape[0]} nearest neighbor found.")

    start_time = time.time()
    bf_points = quad_tree.search_nn_bf(query_pt, n_nearest, radius)
    print(f"Brute force took: {(time.time() - start_time) * 1000.:.3f} ms")
    print(f"{bf_points.shape[0]} nearest neighbor found.")
    
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.scatter(pts[:, 0], pts[:, 1], s = point_size, color = 'k', label = 'Quad Tree Points')
    ax.scatter([query_pt[0]], [query_pt[1]], s = point_size, color = 'b', label = 'Query')

    square = plt.Rectangle(query_pt - radius, radius * 2, radius * 2, color='blue', fill=True, label='Search range', alpha = 0.2)
    ax.add_patch(square)

    ax.scatter(nns[:, 0], nns[:, 1], s = point_size << 1, color = 'r', label = 'Nearest neighbors (k-d tree)')
    ax.scatter(nns_bf[:, 0], nns_bf[:, 1], s = point_size, color = 'g', label = 'Nearest neighbors (quad tree)')

    plt.xticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
    plt.yticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])

    ax.grid(axis = 'both')
    ax.legend()
    plt.show()
    
def tree_structure_viz():
    point_size = 2
    point_per_node = 5    # large scale: 4, small scale: 1
    num_pts = 1000          # use 7 to visualize a small-scale test
    radius  = 0.1         # large scale: 0.01

    query_pt = np.float32([0.23, 0.23])
    square = plt.Rectangle(query_pt - radius, radius * 2, radius * 2, color='blue', fill=True, label='Search range', alpha = 0.2)
    
    range = np.float32([
        [0.5, 0.5],
        [0.5, 0.5]
    ])
    
    kdtree = kdt.KDTree2f(range, 12, point_per_node)
    pts = np.random.rand(num_pts, 2).astype(np.float32)

    start_time = time.time()
    kdtree.build_tree(pts)
    print(f"Inserting {num_pts} points took: {(time.time() - start_time) * 1000.:.3f} ms. Tree size = {kdtree.size()}")
    
    non_leaves, leaves = kdtree.tree_structure()
    
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.scatter(pts[:, 0], pts[:, 1], s = point_size, color = 'k', label = 'K-D Tree Points')
    ax.scatter([query_pt[0]], [query_pt[1]], s = point_size, color = 'b', label = 'Query')
    ax.add_patch(square)

    num_non_leaves = non_leaves.shape[0]
    num_leaves     = leaves.shape[0]
    for i, node_range in enumerate(non_leaves):
        draw_rectangle(ax, node_range, color = 'black', line_width = 2.5)
        if i == num_non_leaves - 1:
            draw_rectangle(ax, node_range, color = 'black', line_width = 2.5, label = 'Non-leaf Node Range')
    for i, node_range in enumerate(leaves):
        draw_rectangle(ax, node_range, color = 'red')
        if i == num_leaves - 1:
            draw_rectangle(ax, node_range, color = 'red', label = 'Leaf Node Range')

    plt.xticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
    plt.yticks([0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])

    ax.legend()
    plt.show()

if __name__ == "__main__":
    np.random.seed(2)
    # full_test()
    tree_structure_viz()