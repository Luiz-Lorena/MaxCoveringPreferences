using JuMP #, Gurobi
import HiGHS, Gurobi, CPLEX
import MultiObjectiveAlgorithms as MOA
import MathOptInterface as MOI

mutable struct Solucao
    objetivos
    centers
    selected
    tempo
end

function execute_model_LibGEOS(circles, population, C, I, p)

    ids = sortperm(population,rev=true)
    best_ids = ids[1:p]
    best_populations = population[best_ids]

    n = size(C,1)
    model = Model()
    #set_silent(model)
    set_optimizer(model, () -> MOA.Optimizer(Gurobi.Optimizer))
    set_attribute(model, MOA.Algorithm(), MOA.DominguezRios())
    @variable(model, y[i in 1:n], Bin)
    @variable(model, x[i in 1:n], Bin)
    @variable(model, z[i in 1:n,j in 1:n], Bin)

    @expression(model, coverage_expr, sum(y[i] for i in 1:n))
    @expression(model, preference_expr, sum(y[i]*population[i] for i in 1:n))
    @expression(model, intersection_expr, sum(z[i,j]*I[i,j] for i in 1:n, j in 1:n))

    @objective(model, Max, [coverage_expr, preference_expr, intersection_expr])
    
    @constraint(model, sum(x[i] for i in 1:n) == p)
    for i in 1:n
        @constraint(model, sum(C[i,j]*x[j] for j in 1:n) >= y[i])
    end
    for i in 1:n
        for j in 1:n
            @constraint(model, z[i,j] <= x[i])
            @constraint(model, z[i,j] <= x[j])
            @constraint(model, z[i,j] >= x[i] + x[j] - 1)
        end
    end

    cplex_time = time()
    optimize!(model)
    cplex_time = time() - cplex_time
    solution_summary(model)
    println("Time: $cplex_time")
    println("Total soluções: $(result_count(model))")

    objetivos = zeros(result_count(model),3)
    
    for id in 1:result_count(model)
        println("##############")
        println("Solução: $id")
        println("Objetivos: ", objective_value(model; result = id))

        objetivos[id,1] = value(coverage_expr; result = id)
        objetivos[id,2] = value(preference_expr; result = id)
        objetivos[id,3] = value(intersection_expr; result = id)

        centers = [i for i in 1:n if value(x[i]; result = id) > 0.5 ]
        #plot_solution_3(id, points, population, centers, radius, p)

        println("Centers: $centers")
        println("Population centers: ",population[centers])

        println("Best ids: $best_ids")
        println("Best population: $best_populations")

        #plot_solution_LibGEOS(circles, centers, best_ids)
    end

    # #plot_pareto_front(objetivos)
    # #plot_pareto_front_percentual(objetivos_percentual)
    # return objetivos
end

function execute_model_new(population, C, p)
    n = size(C,1)
    model = Model()
    #set_silent(model)
    set_optimizer(model, () -> MOA.Optimizer(Gurobi.Optimizer))
    set_attribute(model, MOA.Algorithm(), MOA.DominguezRios())

    #MOI.set(model, MOI.TimeLimitSec(), 120)
    #MOI.set(model, MOI.RelativeGapTolerance(), 0.005)

    @variable(model, y[i in 1:n], Bin)
    @variable(model, x[i in 1:n], Bin)
    @expression(model, coverage_expr, sum(y[i] for i in 1:n))
    @expression(model, preference_expr, sum(y[i]*population[i] for i in 1:n))
    @objective(model, Max, [coverage_expr, preference_expr])
    @constraint(model, sum(x[i] for i in 1:n) == p)
    for i in 1:n
        @constraint(model, sum(C[i,j]*x[j] for j in 1:n) >= y[i])
    end
    cplex_time = time()
    optimize!(model)
    cplex_time = time() - cplex_time
    solution_summary(model)
    println("Time: $cplex_time")
    println("Total soluções: $(result_count(model))")

    objetivos = zeros(result_count(model),2)

    solucoes = []

    for id in 1:result_count(model)
        println("##############")
        println("Solução: $id")
        println("Objetivos: ", objective_value(model; result = id))

        objetivos[id,1] = value(coverage_expr; result = id)
        objetivos[id,2] = value(preference_expr; result = id) 
        
        centers = [i for i in 1:n if value(x[i]; result = id) > 0.6]
        selected = [i for i in 1:n if value(y[i]; result = id) > 0.6]

        solucao = Solucao(objetivos,centers,selected,cplex_time)
        push!(solucoes,solucao)
    end

    plot_pareto_front(objetivos)

    return solucoes
end

function execute_model(B,C,p,π)
    # =================================== #
    # Max Covering with preferences model #
    # =================================== #
    n = size(C,1)
    model = Model()
    #set_time_limit_sec(model, 120)
    #set_silent(model)
    #set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
    set_optimizer(model, () -> MOA.Optimizer(Gurobi.Optimizer))
    #set_optimizer(model, () -> MOA.Optimizer(CPLEX.Optimizer))
    #set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    set_attribute(model, MOA.Algorithm(), MOA.DominguezRios())
    #set_attribute(model, MOA.Algorithm(), MOA.TambyVanderpooten())
    #set_attribute(model, MOA.SolutionLimit(), 10)
    #set_attribute(model, MOA.TimeLimitSec(), 120)
    #set_optimizer_attribute(model, "MIPGap", 0.001)
    # MOI.set(model, MOI.RelativeGapTolerance(), 0.01)

    @variable(model, x[i in 1:n], Bin)
    @variable(model, y[i in 1:n], Bin)

    @expression(model, coverage_expr, sum(y[i] for i in 1:n))
    @expression(model, preference_expr, sum((x[i]+x[j]-1) * abs(B[i,j]) for i in 1:n, j in 1:n))
    @objective(model, Max, [coverage_expr, preference_expr])

    @constraint(model, sum(x[i] for i in 1:n) == p)
    for i in 1:n
        @constraint(model, sum(C[i,j]*x[j] for j in 1:n) >= y[i])
    end
    cplex_time = time()
    optimize!(model)
    cplex_time = time() - cplex_time
    solution_summary(model)

    for id in 1:result_count(model)
        println("##############")
        println("Solução: $id")
        println("Objetivos: ", objective_value(model; result = id))

        valores_y = []
        for i in 1:n
            push!(valores_y,(value(y[i]; result = id)))
        end
        #valores_y = [(value(y[i]; result = id) > 0.6)]
        println("Y: ", valores_y)

        # points_chosen = [i for i in 1:n if value(x[i]; result = id) > 0.6]
        # println("Centros: $points_chosen")

        # Get centers from solution
        total = 0  # total centers chosen from preference list
        centers = Int64[]
        for i in 1:n
            if value(x[i]; result = id) > 0.5
                push!(centers,i)
                # count if it is in the preference list
                if i in π[1:p]
                    total += 1          
                end
            end
        end
        total_percentual = round((total/p)*100,digits=2)  # percentual of centers chosen

        println("Centros: $centers")
        println("Ranking: $(π[1:p])")

        total_coverage = value(coverage_expr; result = id)
        total_coverage_percentual = round( (total_coverage/n)*100,digits=2)

        println("Total Cobertura: $total_coverage ($total_coverage_percentual %)")
        
        println("Total Restrições: $total ($total_percentual %)")     
    end

end