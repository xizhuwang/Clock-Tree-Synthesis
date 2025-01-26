1. Programming Language: C++  
2. Compilation Environment: Linux GNU C++  
   - Compiler Installation: `sudo apt-get install g++`  
   - Environment Setup: `sudo apt-get install build-essential`  
   - Editor: `sudo apt-get install vim`  
   - Compilation: `g++ M11307409.cpp -o cts -Wall`  
   - Execution: `./cts input.cts output.cts`  
   - Output: `cat output.cts`  

### Program Workflow:
1. Read the file.  
2. Find the total center of mass (centroid) in the entire coordinate system.  
3. Divide the entire coordinate system into four equal parts (excluding the source point).  
4. Calculate the center of mass for each partition.  
5. Connect the centroids of each partition, starting vertically and then horizontally (perform collision checks; stop the path if a collision occurs).  
6. Further divide each quadrant into four equal sub-regions (16 smaller areas in total).  
7. Assign each receiving point to a sub-region and calculate the center of mass.  
8. Connect the centroid of each sub-region to its smaller areas, starting horizontally and then vertically (perform collision checks; stop the path if a collision occurs).  
9. Connect the centroid of each small region to the receiving points, starting vertically and then horizontally (perform collision checks; stop the path if a collision occurs).  
10. Connect the source point to the total center of mass, starting horizontally and then vertically (perform collision checks; stop the path if a collision occurs).  
11. Generate the output in the specified format.  

### Reflection:
Initially, the program used MST (Minimum Spanning Tree) for connections, followed by a grid-based overlap detection. However, unresolved segments were tackled using the A* algorithm, which frequently encountered routing failures. To address this, the program was restructured to calculate centroids combined with collision detection. Paths were computed based on points and later merged into line segments for output. This approach resulted in an overall connection strategy resembling an H-tree structure.
