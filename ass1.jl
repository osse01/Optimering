using JuMP, HiGHS

m = 10
L = 100
M = 3
l = [11 13 17 23 37 41 51 61 71 79]
  
A = [0 1 0 0 0
     0 1 0 0 0
     1 0 0 0 0
     0 0 1 0 0
     0 0 0 1 0
     0 0 0 0 1
     0 0 0 1 0
     0 0 1 0 0
     0 1 0 0 0
     1 0 0 0 0]

b = ones(1,m) # Demand

epsilon = 0.1
maxIter = 20

function columnGeneration(dualVar)
    print("\n", dualVar)
    columnModel = Model(HiGHS.Optimizer)
    @variable(columnModel, a[1:m]>=0, Int)
    @constraint(columnModel, column, sum(l[j] * a[j] for j in 1:m) <= L)
    @objective(columnModel, Max, sum(dualVar[i] * a[i] for i in 1:m))
    print("\nColumn Generation!\n")
    optimize!(columnModel)
    print("\n", value.(a))
    return vec(value.(a))
end

for iter in 1:maxIter
    n = size(A,2) 
    LPmodel = Model(HiGHS.Optimizer)
    @variable(LPmodel, x[1:n]>=0)
    @constraint(LPmodel, demand[j in 1:m], sum( A[j,i]*x[i] for i in 1:n ) >= b[j])
    @constraint(LPmodel, maximum[j in 1:m], sum( A[j,i]*x[i] for i in 1:n ) <= M )
    @objective(LPmodel, Min, sum(x) )
    print("\nLP Problem!\n")
    optimize!(LPmodel)
    dualVar = [dual(demand[i]) for i in 1:m]
    column = columnGeneration(dualVar)
    reducedCost = 1 - dualVar' * column
    # Stopping condition
    if reducedCost >= epsilon
        break
    end
    global A = [A column]
    print(size(A))
end

# Run last time
n = size(A,2) 
LPmodel = Model(HiGHS.Optimizer)
@variable(LPmodel, x[1:n]>=0)
@constraint(LPmodel, demand[j in 1:m], sum( A[j,i]*x[i] for i in 1:n ) >= b[j])
@constraint(LPmodel, maximum[j in 1:m], sum( A[j,i]*x[i] for i in 1:n ) <= M )
@objective(LPmodel, Min, sum(x) )
optimize!(LPmodel)

print("\nx values:\n", value.(x), "\n")
print("Objective value: ", objective_value(LPmodel), "\n")
