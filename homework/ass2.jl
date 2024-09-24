using JuMP, HiGHS

epsilon = 0.1 # Stopping condition

c = [10 12 11 10 9 5 5 11]
A = [13 8 10 5 11 13 15 7
     11 14 11 12 8 13 9 10
     7 12 13 11 10 7 8 12]
D1 = [10 6 12 10
      9 11 11 10
      9 10 7 8]
D2 = [13 10 10 9
      11 11 11 7
      9 12 13 13]
b = [126 123 114]
e1 = [74 71 69]
e2 = [70 73 73]

N = size(A,2)
rowA = size(A,1)

function SUB()
    SUB = Model(HiHGS.Optimizer)
    # Code...
end

function RMP(x)
    RMP = Model(HiGHS.Optimizer)

    p = size(x,2)
    @variable(RMP, lambda[1:iter]>=0)
    @objective(RMP, Max, sum( (c[i]*x[i])*lambda[i] for i in 1:p))
    @constraint(RMP, coupling[1:rowA], sum(A*x[i]*lambda[i] for i in 1:p) <= b[1:rowA])
    @constraint(RMP, convex, sum(lambda) == 1)
    optimize!(RMP)
    
end

x1 = zeros(1,4) # Start point
x2 = zeros(1,4) # Start point
for iter in 1:100
    RMP()
    print(value.(lambda))
    SUB()
    print(value.(x))

end
