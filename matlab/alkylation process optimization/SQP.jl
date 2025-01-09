## Alkylation Process Optimization in Julia

using JuMP, LinearAlgebra, Gurobi, Plots
include("misc.jl")

# Optimization Setup
max_iter = 100
tol = 1e-3
global mu = 1e0
slack_weight = 100000
epsilon = 1e1
step_size = [(bounds[i][2]-bounds[i][1])/2 for i=1:length(bounds)]
f_star = profit_function(x_star)
f_plot = [profit_function(x_k)]
for iter in 1:max_iter
    println("\nIteration -------------------------------------------------------------------- $iter\n")

    # Formulate LP problem
    model = Model(Gurobi.Optimizer)

    @variable(model, -step_size[i] <= d[i=1:10] <= step_size[i]) # Search Directions
    @variable(model, slack[i=1:11]) # Slack variables

    # Bounds
    for i in 1:10
        @constraint(model, d[i] + x_k[i] <= bounds[i][2] ) # Upper
        @constraint(model, d[i] + x_k[i] >= bounds[i][1] ) # Lower
    end

    iter = 1
    # Constraints (Inequalities)
    for (con_name, (con_func, grad_func)) in Ineqconstraints
        grad_k = grad_func(x_k)
        f_k = con_func(x_k)
        @constraint(model, f_k + dot(grad_k, d)  >= 0)
        #iter = iter + 1
    end

    # Constraints (Equalities)
    
    for (con_name, (con_func, grad_func)) in Eqconstraints
        grad_k = grad_func(x_k)
        f_k = con_func(x_k)
        @constraint(model, f_k + dot(grad_k, d) + slack[iter] == 0)
        iter = iter + 1
    end

    L = mu * I
    # Objective function
    @objective(model, Max, profit_function(x_k) + dot(grad_pf(x_k), d) +
                            sum(sum(d[i] * L[i, j] * d[j] for i in 1:10) for j in 1:10) -
                            slack_weight*sum(slack.^2))

    # Solve the problem
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        println("Optimization did not converge.")
        global step_size = step_size * 3
        global mu = 0.9 * mu
        continue
    end

    x_new = value.(d) + x_k
    f_new = objective_value(model)

    println("\nUpdated Variables:")
    println("x_new:\n$x_new")
    println("Search Direction:\n$(value.(d))")

    append!(f_plot, profit_function(x_k))
    global mu = 0.9 * mu

    #Convergence check
    if abs( (profit_function(x_k) - f_star) / f_star) <= tol
        println("Converged: Relative error is within tolerance.")
        break
    end

    # Check prediction
    diff = profit_function(x_new) - profit_function(x_k)
    better = diff > 0
    feas_list = check_feasibility_true(x_k, Ineqconstraints)
    if better
        global x_k = x_new
        global step_size = step_size * 2
    else
        global step_size = step_size*0.5
    end

end

println("Optimal solution:")
println(x_k)
println("With value:")
println(profit_function(x_k))
println("\nFeasibility check:")
check_feasibility(x_k, Ineqconstraints)
check_feasibility(x_k, Eqconstraints)

L = length(f_plot)
p = plot(f_plot, xlabel="Iteration", ylabel="Objective Value", title="Profit for Linearized model", ylims=(500,2000))
plot!([f_star for i in 1:L])
savefig(p,"earnings.png")