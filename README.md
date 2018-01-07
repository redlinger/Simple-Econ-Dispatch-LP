# Simple-Econ-Dispatch-LP-Python
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
