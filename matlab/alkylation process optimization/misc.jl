
# Initial values
c = [0.063, 5.04, 0.035, 10, 3.36] # Cost
d4l  = 99/100
d4u  = 100/99
d7l  = 99/100
d7u  = 100/99
d9l  = 9/10
d9u  = 10/9
d10l = 99/100
d10u = 100/99

# Bounds for variables
bounds = [
    (1, 2000),
    (1, 16000),
    (1, 120),
    (1, 5000),
    (1, 2000),
    (85, 93),
    (90, 95),
    (3, 12),
    (1.2, 4),
    (145, 162)
]

# Initial guess
global x_k = [1745, 12000, 110, 3048, 1974, 89.2, 92.8, 8, 3.6, 145]
# Optimum
global x_star = [1698, 15818, 54.1, 3031, 2000, 90.1, 95, 10.5, 1.6, 154]

# Objective function
function profit_function(x)
    return c[1]*x[4]*x[7] - c[2]*x[1] - c[3]*x[2] - c[4]*x[3] - c[5]*x[5]
end

grad_pf(x) = [
    -c[2], -c[3], -c[4], c[1]*x[7], -c[5], 0, c[1]*x[4], 0, 0, 0
]

# Equality constraints
function Ineq_constraints_and_gradients()
    constraints = Dict{Symbol, Tuple{Function, Function}}()

    constraints[:con1] = (
        x -> x[1] * (1.12 + 0.13167 * x[8] - 0.00667 * x[8]^2) - d4l * x[4],
        x -> [(1.12 + 0.13167 * x[8] - 0.00667 * x[8]^2), 0, 0, -d4l, 0, 0, 0, (0.13167 - 2 * 0.00667 * x[8]), 0, 0]
    )

    constraints[:con2] = (
        x -> -x[1] * (1.12 + 0.13167 * x[8] - 0.00667 * x[8]^2) + d4u * x[4],
        x -> [-(1.12 + 0.13167 * x[8] - 0.00667 * x[8]^2), 0, 0, d4u, 0, 0, 0, -(0.13167 - 2 * 0.00667 * x[8]), 0, 0]
    )

    constraints[:con3] = (
        x -> (86.35 + 1.098 * x[8] - 0.038 * x[8]^2 + 0.325 * (x[6] - 89)) - d7l * x[7],
        x -> [0, 0, 0, 0, 0, 0.325, -d7l, (1.098 - 2 * 0.038 * x[8]), 0, 0]
    )

    constraints[:con4] = (
        x -> -(86.35 + 1.098 * x[8] - 0.038 * x[8]^2 + 0.325 * (x[6] - 89)) + d7u * x[7],
        x -> [0, 0, 0, 0, 0, -0.325, d7u, -(1.098 - 2 * 0.038 * x[8]), 0, 0]
    )

    constraints[:con5] = (
        x -> (35.82 - 0.222 * x[10]) - d9l * x[9],
        x -> [0, 0, 0, 0, 0, 0, 0, 0, -d9l, -0.222]
    )

    constraints[:con6] = (
        x -> -(35.82 - 0.222 * x[10]) + d9u * x[9],
        x -> [0, 0, 0, 0, 0, 0, 0, 0, d9u, 0.222]
    )

    constraints[:con7] = (
        x -> (-133 + 3 * x[7]) - d10l * x[10],
        x -> [0, 0, 0, 0, 0, 0, 3, 0, 0, -d10l]
    )

    constraints[:con8] = (
        x -> -(-133 + 3 * x[7]) + d10u * x[10],
        x -> [0, 0, 0, 0, 0, 0, -3, 0, 0, d10u]
    )

    return constraints
end

# Equality constraints
function Eq_constraints_and_gradients()
    constraints = Dict{Symbol, Tuple{Function, Function}}()
    
    constraints[:con9] = (
        x -> 1.22 * x[4] - x[1] - x[5],
        x -> [-1, 0, 0, 1.22, -1, 0, 0, 0, 0, 0]
    )

    constraints[:con10] = (
        x -> 98000 * x[3]/( x[4] * x[9] + 1000 * x[3]) - x[6] ,
        x -> [0, 0, 98000 * x[4]*x[9] / (x[4]*x[9] + 1000*x[3])^2, 
              -98000 * x[3]*x[9] / (x[4]*x[9] + 1000*x[3])^2, 0, -1, 0, 0, 
              -98000 * x[3]*x[4] / (x[4]*x[9] + 1000*x[3])^2, 0]
    )

    constraints[:con11] = (
        x -> (x[2] + x[5])/x[1] - x[8],
        x -> [-(x[2]+x[5])/x[1]^2 , 1/x[1], 0, 0, 1/x[1], 0, 0, -1, 0, 0]
    )
    return constraints
end

Eqconstraints = Eq_constraints_and_gradients()
Ineqconstraints = Ineq_constraints_and_gradients()

# Feasibility check function
function check_feasibility(x, constraints)
    println("Checking feasibility...")
    for (con_name, (con_func, _)) in constraints
        println("$con_name: $(con_func(x))")
    end
end

function check_feasibility_true(x, constraints)
    feas_list = []
    for (con_name, (con_func, grad_func)) in constraints
        if con_func(x) >= -1e-1
            append!(feas_list, true)
        else
            append!(feas_list, false)
        end
    end
    return feas_list
end
