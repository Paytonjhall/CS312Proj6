#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import heapq
import random


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number of solutions found during search, the
        best solution found.  You may use the other three field however you like.
        algorithm</returns>
    '''

    def fancy(self, time_allowance=60):
	# Get/Set Variables
        results = {}
        cities = self._scenario.getCities()
        start_time = time.time()
        sortedCities = sorted(cities, key=lambda city: city._x)
	# Call divide and conquer function
        bssf = self.divideAndConquerRec(sortedCities)
        if bssf.cost == float('inf'):
            sortedCities = sorted(cities, key=lambda city: city._y)
            bssf = self.divideAndConquerRec(sortedCities)
		# If it doens't work, shuffle cities and try again.
            if bssf.cost == float('inf'):
                while time.time() - start_time < time_allowance:
                    random.shuffle(sortedCities)
                    bssf = self.divideAndConquerRec(sortedCities)
                    if bssf.cost < float('inf'):
                        break
	# Return time and route
        end_time = time.time()
        print(bssf.route)
        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = 0
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    def divideAndConquerRec(self, cities):
	# Get Number of Cities, if 1, return, if not, split them up and do it again
        cityNum = len(cities)
        if cityNum == 1:
            return TSPSolution(cities)
        cityOne = self.divideAndConquerRec(cities[:cityNum // 2])
        cityTwo = self.divideAndConquerRec(cities[cityNum // 2:])
        return self.combineCities(cityOne, cityTwo)

    def combineCities(self, cityOne, cityTwo):
	# Make a best route, Test each city distance within that city, if its better than our best, set it.
        bestRoute = []
        cityOneArr = cityOne.route
        cityTwoArr = cityTwo.route
        bestRoute = TSPSolution(cityOneArr + cityTwoArr)
        bestCost = bestRoute.cost
        for i in range(len(cityOneArr)):
            for j in range(len(cityTwoArr)):
                tempRoute = cityOneArr[:i]
                tempRoute += cityTwoArr[j:]
                tempRoute += cityTwoArr[:j]
                tempRoute += cityOneArr[i:]
                tempSol = TSPSolution(tempRoute)
                if tempSol.cost < bestCost:
                    bestCost = tempSol.cost
                    bestRoute = tempSol
        return bestRoute


    ''' <summary>
        This is the entry point for the default solver
        which just finds a valid random tour.  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of solution,
        time spent to find solution, number of permutations tried during search, the
        solution found, and three null values for fields not used for this
        algorithm</returns>
    '''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True
        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
        This is the entry point for the greedy solver, which you must implement for
        the group project (but it is probably a good idea to just do it for the branch-and
        bound project as a way to get your feet wet).  Note this could be used to find your
        initial BSSF.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number of solutions found, the best
        solution found, and three null values for fields not used for this
        algorithm</returns>
    '''


    # Maybe do divide and conquer for this one

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution,
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints:
        max queue size, total number of states created, and number of pruned states.</returns>
    '''

    def branchAndBound(self, time_allowance=60.0):
        # I just copied all these variables from the randomTour method
        # Variables:
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 1
        max_queue = 0
        num_pruned = 0
        total_states = 0
        start_time = time.time()

        # Copy the random tour bssf and try to beat it!
        # Time and space complexity: O(n)
        random_bssf = self.defaultRandomTour().get('soln')
        self.bssf = self.TSPSolution(random_bssf.cost, [], [], random_bssf.route, 0)

        # Create matrix for cost analysis
        # Time and space complexity: O(n^2)
        self.init_cost = np.zeros((ncities, ncities))
        for i in range(ncities):
            for j in range(ncities):
                if i == j:
                    # if duplicates, set to inf
                    self.init_cost[i, j] = np.inf
                    continue
                # else, set to real value
                self.init_cost[i, j] = cities[i].costTo(cities[j])
        cur_path = [cities[0]._name]
        cur_city = [cities[0]]
        # set up first matrix using helper function
        first_matrix = self._reduceMatrix(self.TSPSolution(0, self.init_cost, [], [], 0), 0, True)
        # get costs
        self.init_cost = first_matrix._matrix.copy()
        self.min_heap = [first_matrix]
        heapq.heapify(self.min_heap)

        # Run branch and bound algorithm
        # Time and space complexity: O(n^2)
        while time.time() - start_time < time_allowance:
            if len(self.min_heap) == 0:
                break
            if len(self.min_heap) > max_queue:
                max_queue = len(self.min_heap)
            # Heap pop complexity: O(nlogn)
            cur_obj = heapq.heappop(self.min_heap)
            for i in range(len(cities)):
                # Skip already visited cities
                if i in cur_obj._path:
                    continue

                # Create new matrix for new paths
                if self.init_cost[cur_obj._index, i] != np.inf:
                    total_states += 1
                    red_obj = self._reduceMatrix(cur_obj, i)

                    # If lower cost add it to bssf
                    # else prune the path
                    if red_obj._cost < self.bssf._cost:
                        # Heap push complexity: O(nlogn)
                        heapq.heappush(self.min_heap, red_obj)

                        # Set new bssf if better
                        if len(red_obj._path) == ncities:
                            count += 1
                            self.bssf = self.TSPSolution(red_obj._cost, red_obj._matrix.copy(), red_obj._path.copy(),
                                                         red_obj._city_path.copy(), red_obj._index)
                    else:
                        # prune
                        num_pruned += 1
        end_time = time.time()
        # set result for returning
        results['cost'] = self.bssf._cost if self.bssf is not None else np.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = TSPSolution(self.bssf._city_path)
        results['max'] = max_queue
        results['total'] = total_states
        results['pruned'] = num_pruned

        return results

# Helper function for reducing matrix, sets up matrix and sets inf to un-'goable' areas
    def _reduceMatrix(self, cur_obj, dest, first=False):
        matrix = cur_obj._matrix.copy()
        cost = cur_obj._matrix.copy()
        a = cur_obj._index
        b = dest
        reduction_cost = 0

        # Set values to inf if not first matrix, add parent cost
        if not first:
            matrix[a, b] = np.inf
            matrix[b, a] = np.inf
            matrix[a] = np.inf
            matrix[:, b] = np.inf
            reduction_cost += cur_obj._cost
            reduction_cost += cost[cur_obj._index, dest]

        # Create reduction cost matrix
        for i in range(len(matrix)):
            min_row = matrix[i].min()
            if min_row != np.inf:
                reduction_cost += min_row
                matrix[i] = matrix[i] - min_row
        for i in range(len(matrix)):
            min_col = matrix[:, i].min()
            if min_col != np.inf:
                reduction_cost += min_col
                matrix[:, i] = matrix[:, i] - min_col
        # add path to solution
        path = cur_obj._path.copy()
        path.append(dest)
        city_path = cur_obj._city_path.copy()
        city_path.append(self._scenario.getCities()[dest])
        return self.TSPSolution(reduction_cost, matrix, path, city_path, dest)


    def greedy(self, time_allowance=60.0):
        start_time = time.time()
		results = {}
		ncities = len(self._scenario._cities)
		matrix, lowerBound = self.createMatrix(False)
		total = 1
		
		currentIndex = 0
		path = []
		path.append(self._scenario._cities[currentIndex])
		# While not all cities have been visited, this algorithm finds the next shortest path,
		# and sets unavailable paths to infinity.
		while len(path) < ncities:
			matrix[currentIndex, 0] = np.inf
			nextMin = np.amin(matrix[currentIndex, :])
			lowerBound += nextMin
			nextIndex = np.where(matrix[currentIndex, :] == nextMin)[0][0]
			path.append(self._scenario._cities[nextIndex])
			matrix[:, nextIndex] += np.repeat(np.inf, ncities)
			matrix[currentIndex, :] += np.repeat(np.inf, ncities)
			currentIndex = nextIndex
			total += 1
		
		if matrix[currentIndex, 0] == np.inf: # Check if the path is viable or not.
			self.bssf = TSPSolution([self._scenario._cities[0], self._scenario._cities[0]])
		else:
			self.bssf = TSPSolution(path)
		
		end_time = time.time()
		
		results['cost'] = self.bssf.cost
		results['time'] = end_time - start_time
		results['count'] = None
		results['soln'] = self.bssf
		results['max'] = None
		results['total'] = total
		results['pruned'] = None
		return results

    # Class used to organize data for returning to GUI
    class TSPSolution:

        def __init__(self, cost, matrix, path, city_path, index):
            self._matrix = matrix
            self._cost = cost
            self._path = path
            self._city_path = city_path
            self._index = index

        def __lt__(self, cmp):
            if len(self._path) > len(cmp._path):
                return True
            if len(self._path) == len(cmp._path) and self._cost < cmp._cost:
                return True
            return False
