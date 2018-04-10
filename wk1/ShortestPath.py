#!/usr/bin/python

import numpy as np
import yaml
import heapq
import sys
import math
        
class PriorityQueue(object):
    def __init__(self, data = []):
        self.heapq = data
        heapq.heapify(self.heapq)
        self.locationmap = dict({item[-1] : item for item in data})
        self.counter = len(data)
        self.REMOVED = "<remove marker>"

    def insert(self, node, priority = 0):
        if node in self.locationmap:
            self.remove(node)
        entry = [priority, self.counter, node]
        self.counter += 1
        self.locationmap[node] = entry
        heapq.heappush(self.heapq, entry)

    def getPriority(self, node):
        return self.locationmap[node][0]

    def remove(self, node):
        entry = self.locationmap.pop(node)
        entry[-1] = self.REMOVED
        return entry[0]

    def pop(self):
        while self.heapq:
            priority, counter, node = heapq.heappop(self.heapq)
            if node is not self.REMOVED:
                del self.locationmap[node]
                return priority, node
        raise KeyError('pop from empty prority queue')

    def size(self):
        return len(self.locationmap)


def indexToLoc(index, xspacing, yspacing):
    return ((index[1]+0.5)*xspacing, (index[0]+0.5)*yspacing)

def locToIndex(loc, xspacing, yspacing):
    return (int(np.rint(loc[1]/yspacing - 0.5)), int(np.rint(loc[0]/xspacing - 0.5)))

def findValidNode(node, loc, occupancy_map, xspacing, yspacing):
    if occupancy_map[node] == 0:
        return node

    maxrow, maxcol = occupancy_map.shape
    closest = None
    maxdist = sys.float_info.max
    for delta in [(1,0), (0,1), (-1,0), (0,-1)]:
        row = node[0] + delta[0]
        col = node[1] + delta[1]
        if row < 0 or row >= maxrow or col < 0 or col >= maxcol or occupancy_map[(row,col)] == 1:
            continue
        pos = indexToLoc((row, col), xspacing, yspacing)
        d = math.sqrt((pos[0] - loc[0])**2 + (pos[1] - loc[1])**2)
        if d < maxdist:
            maxdist = d
            closest = (row, col)
    if closest is None:
        raise ValueError("can not find proper grid index for ", pos)
    return closest

def dijkstras(occupancy_map,x_spacing,y_spacing,start,goal):
    """
    Implements Dijkstra's shortest path algorithm
    Input:
    occupancy_map - an N by M numpy array of boolean values (represented
        as integers 0 and 1) that represents the locations of the obstacles
        in the world
    x_spacing - parameter representing spacing between adjacent columns
    y_spacing - parameter representing spacing between adjacent rows
    start - a 3 by 1 numpy array of (x,y,theta) for the starting position 
    goal - a 3 by 1 numpy array of (x,y,theta) for the finishing position 
    Output: 
    path: list of the indices of the nodes on the shortest path found
        starting with "start" and ending with "end" (each node is in
        metric coordinates)
    """
    ROWS, COLS = occupancy_map.shape
    #convert physical location to index in the grid
    startNode = locToIndex(start, x_spacing, y_spacing)
    startingNodeLoc = indexToLoc(startNode, x_spacing, y_spacing)
    initialcost = math.sqrt((startingNodeLoc[0] - start[0])**2 + (startingNodeLoc[1] - start[1])**2)
    goalNode   = locToIndex(goal, x_spacing, y_spacing)
    
    freelist = np.where(occupancy_map == 0)
    if occupancy_map[startNode[0], startNode[1]] != 0:
        #raise ValueError("start : ({}, {}) invalid, is an obstacle".format(startNode[0], startNode[1]))
        startNode = findValidNode(startNode, start, occupancy_map, x_spacing, y_spacing)
    if occupancy_map[goalNode[0], goalNode[1]] != 0:
        #raise ValueError("goal: ({}, {}) invalid, is an obstacle".format(goalNode[0], goalNode[1]))
        goalNode = findValidNode(goalNode, goal, occupancy_map, x_spacing, y_spacing)
    candidate = [ [sys.float_info.max, 
                   i, (freelist[0][i], freelist[1][i])] for i in range(len(freelist[0]))] 
    visited = set([])
    queue = PriorityQueue(candidate)
    paths = {}
    found = False

    #update initial cost
    queue.remove(startNode)
    queue.insert(startNode, initialcost)
    paths[startNode] = None
    updateInitial(occupancy_map, ROWS, COLS, start, startNode, 0, 1, queue, paths, x_spacing, y_spacing, initialcost)
    updateInitial(occupancy_map, ROWS, COLS, start, startNode, 0, -1, queue, paths, x_spacing, y_spacing, initialcost)
    updateInitial(occupancy_map, ROWS, COLS, start, startNode, 1, 0, queue, paths, x_spacing, y_spacing, initialcost)
    updateInitial(occupancy_map, ROWS, COLS, start, startNode, -1, 0, queue, paths, x_spacing, y_spacing, initialcost)
    while queue.size() > 0:
        priority, current = queue.pop()
        if current == goalNode:
            found = True
            break
        #not reaching goal node yet, for each of its neighbor, update the weight
        visited.add(current)
        update(occupancy_map, ROWS, COLS, current, 0, 1, priority, queue, paths, visited, x_spacing, y_spacing)
        update(occupancy_map, ROWS, COLS, current, 0, -1, priority, queue, paths, visited, x_spacing, y_spacing)
        update(occupancy_map, ROWS, COLS, current, 1, 0, priority, queue, paths, visited, x_spacing, y_spacing)
        update(occupancy_map, ROWS, COLS, current, -1, 0, priority, queue, paths, visited, x_spacing, y_spacing)
        
    if not found:
        raise ValueError("fail to find shortest path")
    node = goalNode
    shortestpath = []
    while node is not None:
        shortestpath.append(node)
        node = paths[node]
    #shortestpath.append(startNode)
    #print (startNode)
    #print ('*', list(reversed(shortestpath)))
    #print (goalNode)
    p = list(reversed([ indexToLoc(n, x_spacing, y_spacing) for n in shortestpath]))
    #start and final position may not fall on center of the grid
    if abs(p[0][0] - start[0]) > 0.0005 or abs(p[0][1] - start[1]) > 0.0005:
        p.insert(0, [start[0][0], start[1][0]])
    if abs(p[-1][0] - goal[0]) > 0.0005 or abs(p[-1][1] - goal[1]) > 0.0005:
        p.append([goal[0][0], goal[1][0]])
    res = np.array(p)
    print (res)
    return res

def updateInitial(occupancy_map, maxrow, maxcol, start, startNode, deltax, deltay, queue, paths, 
                  x_spacing, y_spacing, initialcost):
    x = startNode[0] + deltax
    y = startNode[1] + deltay
    #check if updated position is out of bound
    if x < 0 or x > maxrow or y < 0 or y > maxcol:
        return
    #check if it is an obstacle
    if occupancy_map[x, y] == 1:
        return

    startpos = indexToLoc((startNode[0], startNode[1]), x_spacing, y_spacing)
    deltaloc = indexToLoc((x,y), x_spacing, y_spacing)

    if deltay == 0 and abs(start[0][0] - deltaloc[0]) > 1e-3:
        return
    elif deltax == 0 and abs(start[1][0] - deltaloc[1]) > 1e-3:
        return
    delta = math.sqrt((deltaloc[0] - start[0][0])**2 + (deltaloc[1] - start[1][0])**2) 
    delta1 = math.sqrt((deltaloc[0] - startpos[0])**2 + (deltaloc[1] - startpos[1])**2) + initialcost 
    queue.remove((x,y))
    if delta < delta1:
        queue.insert((x,y), delta) 
        paths[(x,y)] = None
    else:
        queue.insert((x,y), delta1)
        paths[(x,y)] = startNode

def update(occupancy_map, maxrow, maxcol, loc, deltax, deltay, priority, queue, paths, visited,
            x_spacing, y_spacing):
    x = loc[0] + deltax
    y = loc[1] + deltay
    #check if updated position is out of bound
    if x < 0 or x >= maxrow or y < 0 or y >= maxcol:
        return
    #check if it is an obstacle
    if occupancy_map[x, y] == 1:
        return
    #have we visit this node already
    if (x,y) in visited:
        return
    itsPriority = queue.getPriority((x,y))
    deltaloc = indexToLoc((x,y), x_spacing, y_spacing)
    delta = math.sqrt((deltaloc[0] - loc[0])**2 + (deltaloc[1] - loc[1])**2)
    if priority + delta < itsPriority:
        #to update priority queue, have to first remove it, then add it back with updated information
        queue.remove((x,y))
        queue.insert((x,y), priority+delta)
        paths[(x,y)] = loc
    
def test():
    """
    Function that provides a few examples of maps and their solution paths
    """
    test_map1 = np.array([
              [1, 1, 1, 1, 1, 1, 1, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 0, 0, 0, 0, 0, 0, 1],
              [1, 1, 1, 1, 1, 1, 1, 1]])
    x_spacing1 = 0.13
    y_spacing1 = 0.2
    start1 = np.array([[0.3], [0.3], [0]])
    goal1 = np.array([[0.6], [1], [0]])
    path1 = dijkstras(test_map1,x_spacing1,y_spacing1,start1,goal1)
    true_path1 = np.array([
        [ 0.3  ,  0.3  ],
        [ 0.325,  0.3  ],
        [ 0.325,  0.5  ],
        [ 0.325,  0.7  ],
        [ 0.455,  0.7  ],
        [ 0.455,  0.9  ],
        [ 0.585,  0.9  ],
        [ 0.600,  1.0  ]
        ])
    if np.array_equal(path1,true_path1):
      print("Path 1 passes")

    test_map2 = np.array([
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0],
             [1, 1, 1, 1, 1, 1, 1, 1],
             [1, 0, 0, 1, 1, 0, 0, 1],
             [1, 0, 0, 1, 1, 0, 0, 1],
             [1, 0, 0, 1, 1, 0, 0, 1],
             [1, 0, 0, 0, 0, 0, 0, 1],
             [1, 0, 0, 0, 0, 0, 0, 1],
             [1, 1, 1, 1, 1, 1, 1, 1]])
    start2 = np.array([[0.5], [1.0], [1.5707963267948966]])
    goal2 = np.array([[1.1], [0.9], [-1.5707963267948966]])
    x_spacing2 = 0.2
    y_spacing2 = 0.2
    path2 = dijkstras(test_map2,x_spacing2,y_spacing2,start2,goal2)
    true_path2 = np.array([[ 0.5,  1.0],
                           [ 0.5,  1.1],
                           [ 0.5,  1.3],
                           [ 0.5,  1.5],
                           [ 0.7,  1.5],
                           [ 0.9,  1.5],
                           [ 1.1,  1.5],
                           [ 1.1,  1.3],
                           [ 1.1,  1.1],
                           [ 1.1,  0.9]])
    if np.array_equal(path2,true_path2):
      print("Path 2 passes")

def test_for_grader():
    """
    Function that provides the test paths for submission
    """
    test_map1 = np.array([
              [1, 1, 1, 1, 1, 1, 1, 1, 1],
              [1, 0, 1, 0, 0, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 1, 0, 1, 0, 1, 0, 1],
              [1, 0, 0, 0, 1, 0, 0, 0, 1],
              [1, 1, 1, 1, 1, 1, 1, 1, 1]])
    x_spacing1 = 1
    y_spacing1 = 1
    start1 = np.array([[1.5], [1.5], [0]])
    goal1 = np.array([[7.5], [1], [0]])
    path1 = dijkstras(test_map1,x_spacing1,y_spacing1,start1,goal1)
    s = 0
    for i in range(len(path1)-1):
      s += np.sqrt((path1[i][0]-path1[i+1][0])**2 + (path1[i][1]-path1[i+1][1])**2)
    print("Path 1 length:")
    print(s)


    test_map2 = np.array([
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 1],
            [1, 1, 1, 1, 1, 1, 1, 1]])
    start2 = np.array([[0.4], [0.4], [1.5707963267948966]])
    goal2 = np.array([[0.4], [1.8], [-1.5707963267948966]])
    x_spacing2 = 0.2
    y_spacing2 = 0.2
    path2 = dijkstras(test_map2,x_spacing2,y_spacing2,start2,goal2)
    s = 0
    for i in range(len(path2)-1):
      s += np.sqrt((path2[i][0]-path2[i+1][0])**2 + (path2[i][1]-path2[i+1][1])**2)
    print("Path 2 length:")
    print(s)



def main():
    # Load parameters from yaml
    param_path = 'params.yaml' # rospy.get_param("~param_path")
    f = open(param_path,'r')
    params_raw = f.read()
    f.close()
    params = yaml.load(params_raw)
    # Get params we need
    occupancy_map = np.array(params['occupancy_map'])
    pos_init = np.array(params['pos_init'])
    pos_goal = np.array(params['pos_goal'])
    x_spacing = params['x_spacing']
    y_spacing = params['y_spacing']
    path = dijkstras(occupancy_map,x_spacing,y_spacing,pos_init,pos_goal)
    print(path)

if __name__ == '__main__':
    #main()
    #test()
    test_for_grader()

