A KD-tree (K-dimensional tree) is a data structure used for efficient multidimensional search operations, particularly in spatial search applications. It divides space into regions to organize and optimize the search for nearest neighbors or points within a certain distance threshold.

**Working of KD-Tree:**
1. **Building the Tree:**
   - Choose a splitting dimension and median value based on the points' coordinates.
   - Partition the points into two subsets, one with coordinates less than the median and another with coordinates greater.
   - Recursively apply the process to each subset until a leaf node is reached.

2. **Querying:**
   - Start at the root and traverse down the tree based on the splitting dimensions.
   - At each node, decide which subtree to explore based on the query point's coordinates.
   - Update the nearest neighbor and search the other subtree if there's a chance of finding closer points.

**Example:**
Consider a 2D space with points (1, 3), (2, 5), (3, 8), (4, 2), (5, 7). Let's build a KD-tree:

1. **Build the Tree:**
   - Choose the x-coordinate as the splitting dimension. Median value: 3.
   - Partition into (1, 3), (2, 5), (4, 2) and (3, 8), (5, 7).
   - Recursive split on each subset, creating leaf nodes for individual points.

   ```
   Tree:
            (3, 8)
           /      \
     (1, 3)  (2, 5)
                  \
                (4, 2) 
   ```

2. **Querying:**
   - For a query point (3, 4), start at the root (3, 8).
   - Traverse left since 3 < 3.
   - Reach leaf node (2, 5).
   - Update nearest neighbor (2, 5).
   - Backtrack, traverse right to (4, 2).
   - Update nearest neighbor (4, 2).

   The nearest neighbor to (3, 4) is (4, 2).

KD-trees are efficient for applications like spatial databases, machine learning, and computational geometry. They reduce the search space, optimizing nearest neighbor searches in multidimensional datasets.
