# -*- coding: utf-8 -*-
"""
Simple Power Plant Economic Dispatch using Linear Programming in Python

Setup: Power Co. operates a 200MW power plant that consists of four 
gas-fired turbines. The cost to operate each generator/turbine (in $/hr) is a
quadratic function of the power generation (MW). To solve with an LP, we 
will linearize the quadratic cost function. Each generator has a min
and max generation capacity (MW) limits.

Question: What is optimal generation level for each generator/turbine? 
The goal of this code to find the optimal generation level for each turbine 
so as to minimize costs while satisfying demand. 

LP problem statement:
 
   min SUM_ij (cost_coeff(i,j) * x(i,j))

   s.t. SUM_ij (x(i,j)) = demand - SUM_i (p_tot_min_i)
        x_lb(i,j) <= x(i,j) <= x_ub(i,j) for every i,j

 notes:
     generator i, production interval j
     x(i,j) is choice variable- generator i's incremental production within 
         production interval j. Example: if production interval j=1 is 
         (20MW,40MW) and generator i produces 35MW, x(i,j)=15MW 
     cost_coeff(i,j): vector of coefficients for lineraized cost function
         these params are the slopes of the linearized cost function
     demand: total demand requirement (MW)
     p_tot_min_i: min power that generator i must produce. note that x is 
         incremental production above the minimum amount of power a generator must
         produce
     x_lb(i,j): lower bound power generation for generator i and in production
         interval j (MW)
      x_ub(i,j): upper bound power generation for generator i and in production
         interval j (MW)  

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog

# number of line segments to make when linearizing quadartic cost functions
n_seg = 4

class generator:
    # each generator has a name, cost function, lower, upper bounds for power gen
    def __init__(self, name, quad_cost, p_min, p_max):
        self.name = name
        # parameters [a,b,c] of the quadartic cost function a + b*x +cx**2
        # cost in $/hr and x= MW power generation
        self.quad_cost = quad_cost
        self.p_max = p_max    # min power level (MW)
        self.p_min = p_min    # max power level (MW)
        # lower and upper bound for each choice variable is the start and end
        # for each line segment
        self.p_lb_ub = [(0, (self.p_max - self.p_min) / n_seg)]*n_seg
        # min cost to run generator at min power level
        self.p_min_cost = np.dot(self.quad_cost, np.array([1, p_min, p_min**2]))
        
    # linearized cost function
    def lin_cost(self):
        c_params=[]
        # iterate through each line segment of the cost function
        for i in range(1, n_seg + 1):
            # starting point of line segment
            p_strt = (1/n_seg)*(self.p_max-self.p_min)*(i-1)
            # ending point of line segment
            p_end = (1/n_seg)*(self.p_max-self.p_min)*i
            # cost at the starting point
            c_strt = np.dot(self.quad_cost, np.array([1, p_strt, p_strt**2])) 
            # cost at the ending point
            c_end = np.dot(self.quad_cost, np.array([1, p_end, p_end**2])) 
            # slope of the linearized cost function over the segment
            c_delta = ((self.p_max-self.p_min)/n_seg)*(c_end - c_strt)
            c_params.append(c_delta)
        return c_params

    # plot the generator's quadratic cost function
    def plot_cost(self):
        power = np.linspace(self.p_min, self.p_max, 50)
        plt.plot(power, self.quad_cost[0] + self.quad_cost[1]*power + self.quad_cost[2]*power**2,
                 label=self.name)
        

# number of generators
n_gen = 4
# instances of the generator class
gen_a = generator('gen_a', [30, 0.5, 0.05], 10, 45)
gen_b = generator('gen_b', [40, 0.4, 0.04], 15, 60)
gen_c = generator('gen_c', [50, 0.3, 0.03], 20, 55)
gen_d = generator('gen_d', [60, 0.2, 0.02], 15, 40)

# list of generators
gen_lst = [gen_a, gen_b, gen_c, gen_d]
# cost parameters from linearized cost functions
cost_params=[]
# upper and lower bounds for each choice variable
p_bounds=[]
# total min power generation summed across all generators
p_tot_min = 0
# total min cost summed across all generators
p_tot_min_cost = 0
for gen in gen_lst:
    cost_params = cost_params + gen.lin_cost()
    p_bounds = p_bounds + gen.p_lb_ub
    p_tot_min += gen.p_min
    p_tot_min_cost += gen.p_min_cost
    # plot the cost functions
    gen.plot_cost()

# label and legend for plot
plt.ylabel('cost ($ per hour)')
plt.xlabel('power (MW)')
plt.legend
plt.show()

# total demand requirement is 150MW
demand = 150
# matrix of coeffs
A = np.ones((1, 16))


# linear prog 
econ_dis = linprog(cost_params, A_eq=A, b_eq=demand-p_tot_min, 
                   bounds=p_bounds, options={"disp":True})
print(econ_dis)

min_cost = econ_dis.fun + p_tot_min_cost
print('The minimum cost to meet %.0f MW demand is $%.0f' % (demand, min_cost))

print('Optimal power generation:')
j=0
for gen in gen_lst:
    gen_x = econ_dis.x[j:j+n_seg].sum(axis=0) + gen.p_min
    print('Generator %s: %.1f MW' % (gen.name, gen_x))
    j += 4
