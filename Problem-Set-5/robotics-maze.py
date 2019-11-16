############################## DO NOT MODIFY (but please read) ##############################################
# These are helper functions 

class Vertex():
    """A class for vertices"""

    def __init__(self, parent=None, position=None):
        self.parent = parent # Keeps track of parent of a vertex
        self.position = position # Keeps track of position of a vertex

        self.C = float('inf') # Estimate of cost-to-come
        self.H = 0 # Heuristic
        self.F = float('inf') # F = C + H (as in class)

    def __eq__(self, other): # Allows you to check if two vertices are the same by doing "vertex_1 == vertex_2"
        return self.position == other.position   
    
def GetBestVertex(Q):
    """Get the vertex "x_vertex" in Q that has the lowest value of F"""
    """Returns x_vertex and the index of x_vertex in Q"""
    x_vertex = Q[0]
    x_index = 0
    for index, item in enumerate(Q):
        if item.F < x_vertex.F:
            x_vertex = item
            x_index = index
            
    return x_vertex, x_index   

def getNeighbors(x_vertex, maze):
    """Get neighbors of x_vertex"""
    neighbors = []
    for new_position in [(-1, 0), (0, 1), (1,0), (0,-1)]: # Adjacent vertices

        # Get node position
        vertex_position = (x_vertex.position[0] + new_position[0], x_vertex.position[1] + new_position[1])

        # Make sure it is within range
        if vertex_position[0] > (len(maze) - 1) or vertex_position[0] < 0 or vertex_position[1] > (len(maze[len(maze)-1]) -1) or vertex_position[1] < 0:
            continue

        # Make sure it is not occupied by an obstacle
        if maze[len(maze)-vertex_position[1]-1][vertex_position[0]] != 0:
            continue

        # Create new vertex
        new_vertex = Vertex(None, vertex_position)

        # Append
        neighbors.append(new_vertex)
        
    return neighbors
############################################################################################################




################ TO DO: FILL THIS IN #####################################################
def computeH(vertex, B_vertex):
    """Function for computing heuristic H(vertex). To compute H(vertex), use computeH(vertex,B_vertex)"""
    """Recall that the heuristic is an underestimate of the cost-to-go from a vertex to the goal"""
    """You can use the heuristic we discussed in class. You can try other heuristics too if you are curious."""
    x_dist = abs(B_vertex.position[0] - vertex.position[0])
    y_dist = abs(B_vertex.position[1] - vertex.position[1])
    return x_dist + y_dist
##########################################################################################


######### This is the function that implements A star; you will fill in parts of this #########
def astar(maze, A, B):
    """Returns a list of tuples as a path from A to B in the given maze"""
    
    ############################## DO NOT MODIFY ##############################################
    # Create start and end vertices
    A_vertex = Vertex(None, A) 
    B_vertex = Vertex(None, B)
    A_vertex.C = 0
    A_vertex.H = computeH(A_vertex, B_vertex)
    A_vertex.F = A_vertex.H

    # Initialize Q and "dead" state list
    Q = []
    Q.append(A_vertex)
    DeadSet = []
    ############################################################################################

    # Loop until you get to the goal
    while len(Q) > 0:

        ################ TO DO: IMPLEMENT GetBestVertex(Q) #####################################
        # Get the current vertex: "x_vertex", i.e., the one that has the lowest value of F
        x_vertex, x_index = GetBestVertex(Q) # Implement this function
        ########################################################################################

        ################### DO NOT MODIFY THIS ##############################################
        # Check if we are at the goal vertex B
        if x_vertex == B_vertex:
            # If we are, backtrack to get the path
            path = []
            current = x_vertex
            while current is not None:
                path.append(current.position)
                current = current.parent
            return path[::-1] # Return reversed path
        #####################################################################################
        
        ########################### TO DO ############################################
        # Remove x_vertex from Q, add to dead list
        Q.remove(x_vertex)
        DeadSet.append(x_vertex)
        #####################################################################################

        ################### DO NOT MODIFY THIS ##############################################
        # Generate neighbors
        neighbors = getNeighbors(x_vertex, maze)
        #####################################################################################
        
        # TO DO: FILL THIS IN #############################
        # Loop through neighbors, update costs, etc.
        for x_prime in neighbors:
            if x_prime in DeadSet:
                continue
            
            length = 1 #abs(x_vertex.position[0] - x_prime.position[0]) + abs(x_vertex.position[1] - x_prime.position[1])
            tentative_C = x_vertex.C + length
            
            if tentative_C < x_prime.C:
                x_prime.parent = x_vertex
                x_prime.C = tentative_C
                x_prime.H = x_prime.C + x_prime.H
                
                if x_prime not in Q:
                    Q.append(x_prime)
                
                    
        ##########################################################
        
##############################################################################################        
        


############################## DO NOT MODIFY ##############################################
maze = [[0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 1, 1, 0],
        [0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
        [1, 1, 1, 1, 1, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]]
    
A = (0, 0)
B = (7, 6)
    
path = astar(maze, A, B)
print(path)

    
############################################################################################    


# This cell visualizes the path you found above
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np

maze_modified = 1.0 - np.array(maze)
for v in path: 
    maze_modified[len(maze)-v[1]-1][v[0]] = 0.8
maze_modified[len(maze)-A[1]-1][A[0]] = 0.7
maze_modified[len(maze)-B[1]-1][B[0]] = 0.7
    
plt.imshow(maze_modified); # , cmap='hot');
plt.axis('off');
plt.show()
