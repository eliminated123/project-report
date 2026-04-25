import numpy as np
from numpy import random
import matplotlib.pyplot as plt

class Lattice:
    def __init__(self, size = 100, seed = 1234, sigma = 0.1, gamma = 0.005, lockdown_start= (10**100) , movement_prob=1.0, lockdown_strength=0.0): ##initialising variables 
        self.size = size
        self.grid = np.zeros((size,size), dtype=int) ## creates grid using arrays of size x size with only vacant sites represented by 0
        self.rng = np.random.default_rng(seed)
        
        self.sigma = sigma
        self.gamma = gamma

        self.sus_pop = [] ##empty list to store population of each agent
        self.exp_pop = []
        self.inf_pop = []
        self.rec_pop = []

        self.lockdown_start = lockdown_start
        self.movement_prob = movement_prob ##movement_prob set to 1 to ensure that the model initially allows movement
        self.lockdown_strength = lockdown_strength
        self.current_step = 0 ##starts the step count at 0

    def record_counts(self):
            self.sus_pop.append(np.sum(self.grid == 1))
            self.exp_pop.append(np.sum(self.grid == 2))
            self.inf_pop.append(np.sum(self.grid == 3))
            self.rec_pop.append(np.sum(self.grid == 4)) ##records total number of each agents in lattice and adds it to their corresponding list to store population

    def initialise_grid(self, total_particles = 250, ratio_exp = 0.05, ratio_rec = 0.0): 
        init_exp = int(ratio_exp * total_particles) ##calculates number of intial exposed agents
        init_rec = int(ratio_rec * total_particles) ##calculates number of intial recovered/vaccinated agents
        init_sus = total_particles - init_exp - init_rec ##calculates number of initial susceptible agents

        total_cells = self.size * self.size ##total number of cells
        
        chosen_index= self.rng.choice(total_cells, total_particles, replace = False) ##randomly selects total particle numbers from 0 to total cells, with no duplicates

        coordinates = np.unravel_index(chosen_index, (self.size,self.size)) ## converts indices to grid positions (co-ordinates)
        
        self.grid[coordinates[0][:init_exp], coordinates[1][:init_exp]] = 2 ##takes first exposed positions and changes grid co-ordinates to 2 to represent exposed agents
        self.grid[coordinates[0][init_exp:init_exp + init_rec], coordinates[1][init_exp:init_exp + init_rec]] = 4  ##takes the positions after the last exposed until to the number == init_exp + init_red and changed grid value to 4 
        self.grid[coordinates[0][init_exp + init_rec:], coordinates[1][init_exp + init_rec:]] = 1 ##takes remaining positions and changes grid co-ordinates to 1 to represent susceptible agents
        
        self.record_counts()
        
        return(self.grid) ##returns updated grid

    def get_neighbours(self, i , j):  
        neighbours = [] ##empty list to store the neighbouring cells of a co-ordinate
        if i > 0:
            neighbours.append((i - 1, j))
        if i < self.size - 1:
            neighbours.append((i + 1, j))
        if j > 0:
            neighbours.append((i, j - 1))
        if j < self.size - 1:
            neighbours.append((i, j + 1)) ## finds and stores cells above and next to the co-ordinate (parameters chosen so the grid is modelled as a hard wall)

        if i > 0 and j > 0:
            neighbours.append((i - 1, j - 1))
        if i > 0 and j < self.size - 1:
            neighbours.append((i - 1, j + 1))
        if i < self.size - 1 and j > 0:
            neighbours.append((i + 1, j - 1))
        if i < self.size -1  and j < self.size - 1 :
            neighbours.append((i + 1, j + 1)) ##find and stores diagonal cells

        return neighbours ##returns list of neighbours


    def move_agents(self):
        agents = np.argwhere(self.grid !=0) ##identifies the non-vacant sites in grid
        self.rng.shuffle(agents) ##randomly reorganises agents in array

        if self.current_step > self.lockdown_start:
            movement_prob = self.lockdown_strength ##when lock down starts strength will reduce probability of movement 
        else:
            movement_prob = 1

        for i, j in agents:
            if self.grid[i,j] == 0:
                continue ##identifies if the agent has already moved by checking if co-ordinate is empty
            
            if self.rng.random() > movement_prob:
                continue ##when lockdown starts the movement of agents will stop as it skips over the code below
            
            state = self.grid[i, j] ##records value (agent) in cell 
            neighbours = self.get_neighbours(i, j) ##find neighbours using function above

            ni, nj = neighbours[self.rng.integers(len(neighbours))] ##randomises the neightbours in list and selects one

            if self.grid[ni, nj] == 0: ##checks if the chosen neighbouring cell is empty
                self.grid[ni, nj] = state ##updates chosen cell to the value of the agent (moves agent)
                self.grid[i, j] = 0 ##updates original position of agent back to 0 (empty)

    def update_agents(self):
            agents = np.argwhere(self.grid !=0) ##locates all filled sites
            for i, j in agents:
                if self.grid[i,j] == 1:
                    neighbours = self.get_neighbours(i, j) ##records neighbouring cells to the susceptible agent
                    values = [self.grid[ni,nj] for ni, nj in neighbours] ##finds the values of the neigbouring cells in grid
                    prob_inf = 1 *  (values.count(3)/(len(neighbours))) ##calculates probability of infection depending on number of infected agents surrounding 
                    
                    if self.rng.random() < prob_inf: ##compares random value (0-1) to calculated probability value above
                        self.grid[i, j] = 2 ##if random value is smaller than agent moves from susceptible to exposed

                elif self.grid[i,j] == 2:
                    if self.rng.random() < self.sigma:  
                        self.grid[i,j] = 3 ##for exposed agent if the random value is smaller than chosen value of sigma agent becomes infected
                
                elif self.grid[i,j] == 3:
                    if self.rng.random() < self.gamma:
                        self.grid[i,j] = 4 ##for infected agent if the random value is smaller than chosen value of gamma agent recovers
    
    
    def step(self):
        self.move_agents() ##runs the moving agent function
        self.update_agents() ##calulates whether agents moves down SEIR
        self.record_counts() ##records the population of each agent at each step
        self.current_step += 1 ##updates the step counter for lockdown model 

    def run(self, MC_steps):
        for i in range (MC_steps):
            self.step() ##runs monte carlo simulation over a given number of steps
    
    def plot_population(self):
        time = range(len(self.sus_pop)) ##creates a list which is the number of monte carlo steps
        fig, ax = plt.subplots()
        plt.plot(time, self.sus_pop, label = 'Susceptible')
        plt.plot(time, self.exp_pop, label = 'Exposed')
        plt.plot(time, self.inf_pop, label = 'Infected')
        plt.plot(time, self.rec_pop, label = 'Recovered') ##plots of each agents population at each step against step
        plt.xlabel("Monte Carlo Step")
        plt.ylabel("number of agents")
        plt.legend()
        plt.grid(True) ##adds grid to plot for better visualisation

        plt.show()