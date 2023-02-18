# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# TESTUJ NOWA F-KCJE, RANDOM ZMIANA 

#@load "C:\\Users\\User\\Documents\\PhDWorkspace.jld2"

"""using JLD2
using FileIO
jldsave("C:\\Users\\User\\Documents\\PhDWorkspace.jld2"; ex3_v200_total_surplus,
ex3_v200_producer_surplus,
ex3_v200_consumer_surplus,
ex3_v200_price,
ex3_v200_quantity_produced,
ex3_v200_quantity_sold,
ex3_v200_reselling,
ex3_v200_producer_surplus_singleton,
ex3_v200_quality,
ex3_v200_durability,
ex3_v200_margin,
ex3_v200_quality_exp,
ex3_v200_durability_exp,
ex3_v200_sv, 

ex3_v201_total_surplus,
ex3_v201_producer_surplus,
ex3_v201_consumer_surplus,
ex3_v201_price,
ex3_v201_quantity_produced,
ex3_v201_quantity_sold,
ex3_v201_reselling,
ex3_v201_producer_surplus_singleton,
ex3_v201_quality,
ex3_v201_durability,
ex3_v201_margin,
ex3_v201_quality_exp,
ex3_v201_durability_exp,
ex3_v201_sv,

ex3_v200201_H,
ex3_v200201_Li,
ex3_v200201_Lw,
ex3_v200201_sd,
ex3_v200201_cu,
ex3_v200201_sm,
ex3_v200201_nl,
ex3_v200201_ρm,
ex3_v200201_cc,
ex3_v200201_m)"""

# Impact on social communication of simulation dynamics

ex3_v200_total_surplus = Vector{Vector{Float64}}()
ex3_v200_producer_surplus = Vector{Vector{Float64}}()
ex3_v200_consumer_surplus = Vector{Vector{Float64}}()
ex3_v200_price = Vector{Vector{Vector{Float64}}}()
ex3_v200_quantity_produced = Vector{Vector{Vector{Float64}}}()
ex3_v200_quantity_sold = Vector{Vector{Vector{Float64}}}()
ex3_v200_reselling = Vector{Vector{Vector{Float64}}}()
ex3_v200_producer_surplus_singleton = Vector{Vector{Vector{Float64}}}()
ex3_v200_buying_history = []
ex3_v200_quality = Vector{Vector{Vector{Float64}}}()
ex3_v200_durability = Vector{Vector{Vector{Float64}}}()
ex3_v200_margin = Vector{Vector{Vector{Float64}}}()
ex3_v200_quality_exp = Vector{Vector{Vector{Float64}}}()
ex3_v200_durability_exp = Vector{Vector{Vector{Float64}}}()
ex3_v200_sv = Vector{Vector{Float64}}()

ex3_v201_total_surplus = Vector{Vector{Float64}}()
ex3_v201_producer_surplus = Vector{Vector{Float64}}()
ex3_v201_consumer_surplus = Vector{Vector{Float64}}()
ex3_v201_price = Vector{Vector{Vector{Float64}}}()
ex3_v201_quantity_produced = Vector{Vector{Vector{Float64}}}()
ex3_v201_quantity_sold = Vector{Vector{Vector{Float64}}}()
ex3_v201_reselling = Vector{Vector{Vector{Float64}}}()
ex3_v201_producer_surplus_singleton = Vector{Vector{Vector{Float64}}}()
ex3_v201_buying_history = []
ex3_v201_quality = Vector{Vector{Vector{Float64}}}()
ex3_v201_durability = Vector{Vector{Vector{Float64}}}()
ex3_v201_margin = Vector{Vector{Vector{Float64}}}()
ex3_v201_quality_exp = Vector{Vector{Vector{Float64}}}()
ex3_v201_durability_exp = Vector{Vector{Vector{Float64}}}()
ex3_v201_sv = Vector{Vector{Float64}}()

ex3_v202_total_surplus = Vector{Vector{Float64}}()
ex3_v202_producer_surplus = Vector{Vector{Float64}}()
ex3_v202_consumer_surplus = Vector{Vector{Float64}}()
ex3_v202_price = Vector{Vector{Vector{Float64}}}()
ex3_v202_quantity_produced = Vector{Vector{Vector{Float64}}}()
ex3_v202_quantity_sold = Vector{Vector{Vector{Float64}}}()
ex3_v202_reselling = Vector{Vector{Vector{Float64}}}()
ex3_v202_producer_surplus_singleton = Vector{Vector{Vector{Float64}}}()
ex3_v202_buying_history = []
ex3_v202_quality = Vector{Vector{Vector{Float64}}}()
ex3_v202_durability = Vector{Vector{Vector{Float64}}}()
ex3_v202_margin = Vector{Vector{Vector{Float64}}}()
ex3_v202_quality_exp = Vector{Vector{Vector{Float64}}}()
ex3_v202_durability_exp = Vector{Vector{Vector{Float64}}}()
ex3_v202_sv = Vector{Vector{Float64}}()

ex3_v200201_H = Vector{Int64}()
ex3_v200201_Li = Vector{Float64}()
ex3_v200201_Lw = Vector{Float64}()
ex3_v200201_sd = Vector{Int64}()
ex3_v200201_cu = Vector{Float64}()
ex3_v200201_sm = Vector{Float64}()
ex3_v200201_nl = Vector{Int64}()
ex3_v200201_ρm = Vector{Float64}()
ex3_v200201_cc = Vector{Float64}()
ex3_v200201_m = Vector{Float64}()

for i in 1:50

    println(i)

    N = 200

    li = sample(0.5:0.1:1.0)
    lw = sample(0.5:0.1:1.0)
    h = sample(3:8)
    sd = sample(1:100000)
    cu = sample(0.0:0.1:0.5)
    sm = sample(0.10:0.05:0.30)
    nl = sample(N:N:(5*N))
    ρm = rand(Uniform(0.5, 0.95))
    cc = sample(0.4:0.1:0.6)
    m = sample(1.1:0.1:1.3)

    Random.seed!(sd)

    ex3_v200_sim = TO_GO(500, 2, N, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], cu, true, 0, [ρm, 1.], "softmax", ["internal knowledge", "internal knowledge"], [0., 0.], h, "dist", TriangularDist(0,1,0.5))

    Random.seed!(sd)

    ex3_v201_sim = TO_GO(500, 2, N, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], cu, true, 0, [ρm, 1.], "softmax", ["market research", "internal knowledge"], [sm, 0.], h, "dist", TriangularDist(0,1,0.5))

    #Random.seed!(sd)

    #ex3_v202_sim = TO_GO(500, 2, N, nl, [cc, cc], [m, m], "random", li, lw, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], cu, true, 0, [ρm, 1.], "softmax", ["random", "internal knowledge"], [0., 0.], h, "dist", TriangularDist(0,1,0.5))

    push!(ex3_v200201_H, h)
    push!(ex3_v200201_Li, li)
    push!(ex3_v200201_Lw, lw)
    push!(ex3_v200201_sd, sd)
    push!(ex3_v200201_cu, cu)
    push!(ex3_v200201_sm, sm)
    push!(ex3_v200201_nl, nl)
    push!(ex3_v200201_ρm, ρm)
    push!(ex3_v200201_cc, cc)
    push!(ex3_v200201_m, m)

    push!(ex3_v200_total_surplus, calculate_surplus(ex3_v200_sim, "total", false))
    push!(ex3_v200_producer_surplus, calculate_surplus(ex3_v200_sim, "producer", false))
    push!(ex3_v200_consumer_surplus, calculate_surplus(ex3_v200_sim, "consumer,total",false))
    push!(ex3_v200_quality, getfield.(ex3_v200_sim.sellers, :quality_history))
    push!(ex3_v200_durability, getfield.(ex3_v200_sim.sellers, :durability_history))
    push!(ex3_v200_margin, getfield.(ex3_v200_sim.sellers, :margin_history))
    push!(ex3_v200_price, calculate_price_history.(ex3_v200_sim.sellers; product_life = h))
    push!(ex3_v200_quantity_produced, getfield.(ex3_v200_sim.sellers, :quantity_produced_history))
    push!(ex3_v200_quantity_sold, getfield.(ex3_v200_sim.sellers, :quantity_sold_history))
    push!(ex3_v200_producer_surplus_singleton, calculate_profit_history.(ex3_v200_sim.sellers))
    push!(ex3_v200_reselling, getfield.(ex3_v200_sim.sellers, :reselling_history))
    push!(ex3_v200_buying_history, getfield.(ex3_v200_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v200_quality_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v200_sim.buyers, :quality_expectation_history)]) for y in 1:2])
    push!(ex3_v200_durability_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v200_sim.buyers, :durability_expectation_history)]) for y in 1:2])
    push!(ex3_v200_sv, getfield.(ex3_v200_sim.buyers, :signal_volume))

    push!(ex3_v201_total_surplus, calculate_surplus(ex3_v201_sim, "total", false))
    push!(ex3_v201_producer_surplus, calculate_surplus(ex3_v201_sim, "producer", false))
    push!(ex3_v201_consumer_surplus, calculate_surplus(ex3_v201_sim, "consumer,total",false))
    push!(ex3_v201_quality, getfield.(ex3_v201_sim.sellers, :quality_history))
    push!(ex3_v201_durability, getfield.(ex3_v201_sim.sellers, :durability_history))
    push!(ex3_v201_margin, getfield.(ex3_v201_sim.sellers, :margin_history))
    push!(ex3_v201_price, calculate_price_history.(ex3_v201_sim.sellers; product_life = h))
    push!(ex3_v201_quantity_produced, getfield.(ex3_v201_sim.sellers, :quantity_produced_history))
    push!(ex3_v201_quantity_sold, getfield.(ex3_v201_sim.sellers, :quantity_sold_history))
    push!(ex3_v201_producer_surplus_singleton, calculate_profit_history.(ex3_v201_sim.sellers))
    push!(ex3_v201_reselling, getfield.(ex3_v201_sim.sellers, :reselling_history))
    push!(ex3_v201_buying_history, getfield.(ex3_v201_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v201_quality_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v201_sim.buyers, :quality_expectation_history)]) for y in 1:2])
    push!(ex3_v201_durability_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v201_sim.buyers, :durability_expectation_history)]) for y in 1:2])
    push!(ex3_v201_sv, getfield.(ex3_v201_sim.buyers, :signal_volume))

    """push!(ex3_v202_total_surplus, calculate_surplus(ex3_v202_sim, "total", false))
    push!(ex3_v202_producer_surplus, calculate_surplus(ex3_v202_sim, "producer", false))
    push!(ex3_v202_consumer_surplus, calculate_surplus(ex3_v202_sim, "consumer,total",false))
    push!(ex3_v202_quality, getfield.(ex3_v202_sim.sellers, :quality_history))
    push!(ex3_v202_durability, getfield.(ex3_v202_sim.sellers, :durability_history))
    push!(ex3_v202_margin, getfield.(ex3_v202_sim.sellers, :margin_history))
    push!(ex3_v202_price, calculate_price_history.(ex3_v202_sim.sellers; product_life = h))
    push!(ex3_v202_quantity_produced, getfield.(ex3_v202_sim.sellers, :quantity_produced_history))
    push!(ex3_v202_quantity_sold, getfield.(ex3_v202_sim.sellers, :quantity_sold_history))
    push!(ex3_v202_producer_surplus_singleton, calculate_profit_history.(ex3_v202_sim.sellers))
    push!(ex3_v202_reselling, getfield.(ex3_v202_sim.sellers, :reselling_history))
    push!(ex3_v202_buying_history, getfield.(ex3_v202_sim.buyers, :unit_buying_selling_history))
    push!(ex3_v202_quality_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v202_sim.buyers, :quality_expectation_history)]) for y in 1:2])
    push!(ex3_v202_durability_exp, [mean([getindex.(x,y) for x in getfield.(ex3_v202_sim.buyers, :durability_expectation_history)]) for y in 1:2])
    push!(ex3_v202_sv, getfield.(ex3_v202_sim.buyers, :signal_volume))"""

end

##### HISTORIA ############

# Czy zyski producentów prowadzących badania i nie się różnią?

ex3_4_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów")
plot_ecdf(false, sum.(getindex.(ex3_v202_producer_surplus_singleton, 1)), "Producent działa losowo")

ex3_4_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta", xlim = (-500, 0), legend = :topleft)
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów", xlim = (-500, 0))

ApproximateTwoSampleKSTest(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))) # producent robiący badania osiąga statystycznie istotnie wyższe zyski

# Czy jeśli nie robi się badań, to jest większe zagrożenie stratą?

transition_percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
transition_percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

# Tak, jest - czy oczekiwana strata jest wyższa?

mean(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)) .< 0])
mean(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .< 0])

median(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)) .< 0])
median(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .< 0])

# Tak, w przypadku średniej i mediany oczekiwanych strat, jeśli nie prowadzi się badań to są one wyższe

### Czy zmienność zysków w czasie również jest niższa, jeśli prowadzi się badania?

ex3_4_pl = plot_ecdf(true, std.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, std.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów")

UnequalVarianceTTest(std.(getindex.(ex3_v200_producer_surplus_singleton, 1)), std.(getindex.(ex3_v201_producer_surplus_singleton, 1)))
ApproximateTwoSampleKSTest(std.(getindex.(ex3_v200_producer_surplus_singleton, 1)), std.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

# Z czego wynikają wyższe zyski producentów?

selling_income_200 = sum.(multi.(getindex.(ex3_v200_price, 1), getindex.(ex3_v200_quantity_sold, 1)))
utilization_cost_200 = (1 .- ex3_v200201_cu) .* sum.(multi.(divide.(getindex.(ex3_v200_price, 1), getindex.(ex3_v200_margin, 1)), getindex.(ex3_v200_quantity_produced, 1) .- getindex.(ex3_v200_quantity_sold, 1)))

TR200 = selling_income_200 .+ utilization_cost_200

selling_income_201 = sum.(multi.(getindex.(ex3_v201_price, 1), getindex.(ex3_v201_quantity_sold, 1)))
utilization_cost_201 = (1 .- ex3_v200201_cu) .* sum.(multi.(divide.(getindex.(ex3_v201_price, 1), getindex.(ex3_v201_margin, 1)), getindex.(ex3_v201_quantity_produced, 1) .- getindex.(ex3_v201_quantity_sold, 1)))

TR201 = selling_income_201 .+ utilization_cost_201

ex3_4_pl = plot_ecdf(true, TR200, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, TR201, "Producent bada oczekiwania konsumentów")

ex3_4_pl = plot_ecdf(true, selling_income_200, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, selling_income_201, "Producent bada oczekiwania konsumentów")

ex3_4_pl = plot_ecdf(true, utilization_cost_200, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, utilization_cost_201, "Producent bada oczekiwania konsumentów")

# Nie z przychodu, przychody firm kierujących się swoim przeczuciem są wyższe

fTC(Q, K, D, H, c) = Q .* c .* K .* (1 .- D .^ H) ./ (1 .- D)
fATC(K, D, H, c) = c .* K .* (1 .- D .^ H) ./ (1 .- D)


CS200 = sum.(fTC.(getindex.(ex3_v200_quantity_produced,1), getindex.(ex3_v200_quality,1), getindex.(ex3_v200_durability,1), ex3_v200201_H, ex3_v200201_cc))
CS201 = sum.(fTC.(getindex.(ex3_v201_quantity_produced,1), getindex.(ex3_v201_quality,1), getindex.(ex3_v201_durability,1), ex3_v200201_H, ex3_v200201_cc))

ATC200 = mean.(fATC.(getindex.(ex3_v200_quality,1), getindex.(ex3_v200_durability,1), ex3_v200201_H, ex3_v200201_cc))
ATC201 = mean.(fATC.(getindex.(ex3_v201_quality,1), getindex.(ex3_v201_durability,1), ex3_v200201_H, ex3_v200201_cc))

ex3_4_pl = plot_ecdf(true, CS200, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, CS201, "Producent bada oczekiwania konsumentów")

ex3_4_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_quantity_produced,1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_quantity_produced,1)), "Producent bada oczekiwania konsumentów")

StatsPlots.density((mean.(getindex.(ex3_v200_quantity_produced,1)) .- mean.(getindex.(ex3_v201_quantity_produced,1))) ./ mean.(getindex.(ex3_v200_quantity_produced,1)))

ex3_4_pl = plot_ecdf(true, ATC200, "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, ATC201, "Producent bada oczekiwania konsumentów")



####

ex3_4_pl = plot_ecdf(true, std.(getindex.(ex3_v200_quantity_sold, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, std.(getindex.(ex3_v201_quantity_sold, 1)), "Producent bada oczekiwania konsumentów")

UnequalVarianceTTest(std.(getindex.(ex3_v200_quantity_sold, 1)), std.(getindex.(ex3_v201_quantity_sold, 1)))
ApproximateTwoSampleKSTest(std.(getindex.(ex3_v200_quantity_sold, 1)), std.(getindex.(ex3_v201_quantity_sold, 1)))



####

ex3_4_pl = plot_ecdf(true, mean.(getindex.(ex3_v200_quality, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_quality, 1)), "Producent bada oczekiwania konsumentów")

ApproximateTwoSampleKSTest(std.(getindex.(ex3_v200_quality, 1)), std.(getindex.(ex3_v201_quality, 1)))

ex3_4_pl = plot_ecdf(true, mean.(getindex.(ex3_v200_durability, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_durability, 1)), "Producent bada oczekiwania konsumentów")

mean(mean.(getindex.(ex3_v200_durability, 1)))
mean(mean.(getindex.(ex3_v201_durability, 1)))

UnequalVarianceTTest(mean.(getindex.(ex3_v200_durability, 1)), mean.(getindex.(ex3_v201_durability, 1)))
ApproximateTwoSampleKSTest(mean.(getindex.(ex3_v200_durability, 1)), mean.(getindex.(ex3_v201_durability, 1)))

ex3_4_pl = plot_ecdf(true, mean.(getindex.(ex3_v200_price, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_price, 1)), "Producent bada oczekiwania konsumentów")

ApproximateTwoSampleKSTest(mean.(getindex.(ex3_v200_price, 1)), std.(getindex.(ex3_v201_price, 1)))


Plots.plot(mean(getindex.(ex3_v200_quality,1)))
Plots.plot!(mean(getindex.(ex3_v201_quality,1)))
Plots.plot!(mean(getindex.(ex3_v200_quality,2)))
Plots.plot!(mean(getindex.(ex3_v201_quality,2)))

Plots.plot(mean(getindex.(ex3_v200_durability,1)))
Plots.plot!(mean(getindex.(ex3_v201_durability,1)))
Plots.plot!(mean(getindex.(ex3_v200_durability,2)))
Plots.plot!(mean(getindex.(ex3_v201_durability,2)))

Plots.plot(mean(getindex.(ex3_v200_margin,1)))
Plots.plot!(mean(getindex.(ex3_v201_margin,1)))
Plots.plot!(mean(getindex.(ex3_v200_margin,2)))
Plots.plot!(mean(getindex.(ex3_v201_margin,2)))

Plots.plot(mean(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.plot!(mean(getindex.(ex3_v201_producer_surplus_singleton, 1)))
Plots.plot!(mean(getindex.(ex3_v200_producer_surplus_singleton, 2)))
Plots.plot!(mean(getindex.(ex3_v201_producer_surplus_singleton, 2)))

Plots.plot(mean(getindex.(ex3_v200_quality,1)))
Plots.plot!(mean(getindex.(ex3_v200_quality_exp, 1)))

Plots.plot(mean(getindex.(ex3_v201_quality,1)))
Plots.plot!(mean(getindex.(ex3_v201_quality_exp, 1)))

k=4
Plots.plot(getindex.(ex3_v201_quality,1)[k])
Plots.plot!(getindex.(ex3_v201_quality_exp, 1)[k])

Plots.plot(getindex.(ex3_v200_quality,2)[k])
Plots.plot!(getindex.(ex3_v200_quality_exp, 2)[k])

Plots.plot(mean(getindex.(ex3_v200_quality,2)))
Plots.plot!(mean([getindex.(x,2) for x in ex3_v200_quality_exp]))

Plots.plot(mean(getindex.(ex3_v201_quality,2)))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v201_quality_exp]))

###############################

Plots.scatter(xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)), mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

mod201 = GLM.lm(@formula(y~1/x), DataFrame(x = xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)), y = mean.(getindex.(ex3_v201_producer_surplus_singleton, 1))))

Plots.plot!(xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)), coef(mod201)[1] .+ coef(mod201)[2] ./ xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)))

Plots.scatter(xydiff.(getindex.(ex3_v200_quality, 1), getindex.(ex3_v200_quality_exp,1)), mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

#87


getindex.(ex3_v201_producer_surplus_singleton, 1)[BitVector([all(x .<= 100) for x in getindex.(ex3_v201_quantity_produced, 1)])]

ex3_4_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)[BitVector([all(x .<= 40) for x in getindex.(ex3_v200_quantity_produced, 1)])]), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)[BitVector([all(x .<= 40) for x in getindex.(ex3_v201_quantity_produced, 1)])]), "Producent bada oczekiwania konsumentów")

collect(1:475)[BitVector([any(x .> 100) for x in getindex.(ex3_v200_quantity_produced, 1)])]
collect(1:475)[BitVector([any(x .> 100) for x in getindex.(ex3_v201_quantity_produced, 1)])]

k = 45

Plots.plot(ex3_v200_quantity_produced[k][1])
Plots.plot!(ex3_v200_quantity_sold[k][1])

Plots.plot(ex3_v201_quantity_produced[k][1])
Plots.plot!(ex3_v201_quantity_sold[k][1])

Plots.plot(ex3_v201_quality[k][1])
Plots.plot!(ex3_v201_quality_exp[k][1])
Plots.plot!(twinx(), ex3_v201_producer_surplus_singleton[k][1], color = "black", linealpha = 0.5)

Plots.plot!(ex3_v200_quality[k][1])
Plots.plot!(ex3_v200_quality_exp[k][1])

Plots.plot(ex3_v201_quantity_produced[k][1])
Plots.plot!(ex3_v201_quantity_sold[k][1])

U(k,d,h,β,p) = β .* k .* (1 .-(0.85 .* d) .^ h) ./ (1 .- (0.85 .* d))  - p

p = Plots.plot()

for i in 0.1:0.1:0.9

    p = Plots.plot!(mean(U.(getindex.(ex3_v201_quality, 1),
        getindex.(ex3_v201_durability, 1),
        ex3_v200201_H[1:178],
        i,
        getindex.(ex3_v201_price, 1)
        )))

end

p = Plots.plot()

for i in 0.1:0.1:0.9

    p = Plots.plot!(mean(U.(getindex.(ex3_v200_quality, 1),
        getindex.(ex3_v200_durability, 1),
        ex3_v200201_H,
        i,
        getindex.(ex3_v200_price, 1)
        )))

end

p

k = 0.7

Plots.plot(mean(U.(getindex.(ex3_v201_quality, 1),
getindex.(ex3_v201_durability, 1),
ex3_v200201_H,
k,
getindex.(ex3_v201_price, 1)
)))

Plots.plot!(mean(U.(getindex.(ex3_v200_quality, 1),
getindex.(ex3_v200_durability, 1),
ex3_v200201_H,
k,
getindex.(ex3_v200_price, 1)
)))

Plots.plot( U(ex3_v201_quality_exp[k][1], ex3_v201_durability_exp[k][1], ex3_v200201_H[k], 0.9, ex3_v201_price[k][1]), color = "red")
Plots.plot!()


Plots.plot!(twinx(), ex3_v201_margin[k][1])
Plots.plot!(twinx(), ex3_v201_quality[k][1])


Plots.plot!(twinx(), U(ex3_v201_quality[k][1], ex3_v201_durability[k][1], ex3_v200201_H[k], 0.9, ex3_v201_price[k][1]), color = "blue")



Plots.plot(U(ex3_v201_quality_exp[k][1], ex3_v201_durability_exp[k][1], ex3_v200201_H[k], 0.9, ex3_v201_price[k][1]), color = "red")
Plots.plot!(twinx(), ex3_v201_price[k][1])

Plots.scatter(U(ex3_v201_quality_exp[k][1], ex3_v201_durability_exp[k][1], ex3_v200201_H[k], 0.9, ex3_v201_price[k][1]), ex3_v201_quantity_sold[k][1])

##############################

# Gęstość, różnice między oczekiwaniami a realnymi charakterystykami

ex3_p1_pl = StatsPlots.density(xydiff.(getindex.(ex3_v200_quality, 1), getindex.(ex3_v200_quality_exp,1)), label = "Producent nie prowadzi badań", xlabel = "Odchylenie parametrów produktu [jakość] od oczekiwań konsumentów", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Funkcja gęstości odchyleń parametrów produktu od oczekiwań konsumentów", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)), label = "Producent prowadzi badania")

StatsPlots.density!(xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)), label = "Producent prowadzi badania")
Plots.plot!([0,0], [0,10], color = "black")

ex3_v201_quality[1][1]
ex3_v201_quality_exp[1][1]

savefigs(ex3_p1_pl, "\\plots\\ex2\\PL density diff expectations vs average")

mean(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp]))

StatsPlots.density(mean.(getindex.(ex3_v200_quality, 1)))


ex3_p1_eng = StatsPlots.density(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp]), label = "No research", xlabel = "Average Quality deviation from consumer expectation", ylabel = "f(x)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Density function", legendfontsize = 6, legend = :topleft, normalize = true)
StatsPlots.density!(xydiff.(getindex.(ex3_v201_quality, 1), [getindex.(x,1) for x in ex3_v201_quality_exp]), label = "Research")
Plots.plot!([0,0], [0,10], color = "black")

savefigs(ex3_p1_eng, "\\plots\\ex2\\ENG density diff expectations vs average")

####################################################################

p_ex3_2_pl = Plots.boxplot(xydiff.(getindex.(ex3_v200_quality, 1), getindex.(ex3_v200_quality_exp,1)), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Jakość")
Plots.boxplot!(xydiff.(getindex.(ex3_v201_quality, 1), getindex.(ex3_v201_quality_exp,1)))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

savefigs(p_ex3_2_pl, "\\plots\\ex2\\PL boxplot diff expectations vs average quality")

p_ex3_2_eng = Plots.boxplot(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp]), ylabel = "Deviation of average parameters from expectations", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Quality")
Plots.boxplot!(xydiff.(getindex.(ex3_v201_quality, 1), [getindex.(x,1) for x in ex3_v201_quality_exp]))
Plots.plot!(xticks = ([1,2], ["No research", "Research"]))

savefigs(p_ex3_2_eng, "\\plots\\ex2\\ENG boxplot diff expectations vs average quality")

####################################################################

p_ex3_3_pl = Plots.boxplot(xydiff.(getindex.(ex3_v200_durability, 1), [getindex.(x,1) for x in ex3_v200_durability_exp]), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Trwałość")
Plots.boxplot!(xydiff.(getindex.(ex3_v201_durability, 1), [getindex.(x,1) for x in ex3_v201_durability_exp]))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

savefigs(p_ex3_3_pl, "\\plots\\ex2\\PL boxplot diff expectations vs average durability")

p_ex3_3_eng = Plots.boxplot(xydiff.(getindex.(ex3_v200_durability, 1), [getindex.(x,1) for x in ex3_v200_durability_exp]), ylabel = "Deviation of average parameters from expectations", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Durability")
Plots.boxplot!(xydiff.(getindex.(ex3_v201_durability, 1), [getindex.(x,1) for x in ex3_v201_durability_exp]))
Plots.plot!(xticks = ([1,2], ["No research", "Research"]))

savefigs(p_ex3_3_eng, "\\plots\\ex2\\ENG boxplot diff expectations vs average durability")

#############################################################################################################


percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), 36.5)
percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), 26)

median(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
median(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

savefigs(ex3_4_pl, "\\plots\\ex1\\PL ECDF producer profit")

ex3_4_eng = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "No research", xlabel = "Profit", ylabel = "F(x)", title = "ECDF - Profit")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Research")

savefigs(ex3_4_eng, "\\plots\\ex1\\ENG ECDF producer profit")

1 - percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), 5) / percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), 5)

percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), 36.12)
percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), 25.65)



mean(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .< 0]) / mean(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)) .< 0])


###########################################################################################

gbp_df = DataFrame(x = repeat(ex3_v200201_H, 2), y = vcat(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))), g = repeat( ["Firma nie prowadzi badań konsumenckich", "Firma prowadzi badania konsumenckie"],inner = lastindex(ex3_v200201_H)))

gbp_df = sort!(gbp_df, :x)

ex3_5_pl = @df gbp_df groupedboxplot(:x, :y, group = :g, xlabel = "Liczba okresów przydatności dobra", ylabel = "Zysk firmy")


percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), 5)
percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), 5)

###########################################################################################

ex3_6_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], "Producent bada oczekiwania konsumentów")

ex3_7_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], "Producent bada oczekiwania konsumentów")

percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], 5)
percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], 5)

percentile(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], 5)
percentile(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], 5)

Plots.plot(ex3_6_pl, ex3_7_pl, layout = (1,2), plot_title = "Dystybuanta empiryczna - zysk", plot_titlefontsize = 12)

############################################################################################

k1 = 2000
k2 = 5000
k3 = 7000

ex3_8_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "<300")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], "")

ex3_9_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "400-k1")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], "")

ex3_10_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "700-900")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], "")

ex3_11_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "1000-1200")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], "")

Plots.plot(ex3_8_pl, ex3_9_pl, ex3_10_pl, ex3_11_pl)

############################################################################################

Plots.scatter(ex3_v200201_Li, mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.scatter!(ex3_v200201_Li, mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

λ = sort(unique(ex3_v200201_Li))
θ = sort(unique(ex3_v200201_Lw))

mean_surplus_no_research = [mean(mean.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_Li .== li) .& (ex3_v200201_Lw .== lw)]) for li in λ, lw in θ]

ex4_p3 = StatsPlots.heatmap(λ, θ, mean_surplus_no_research', xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta badającego konsumentów", titlefontsize = 8)


mean_surplus_research = [mean(mean.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_Li .== li) .& (ex3_v200201_Lw .== lw)]) for li in λ, lw in θ]

ex4_p3 = StatsPlots.heatmap(λ, θ, mean_surplus_research', xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta badającego konsumentów", titlefontsize = 8)

mean_surplus_difference = [mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[(ex3_v200201_Li .== li) .& (ex3_v200201_Lw .== lw)]) for li in λ, lw in θ]

ex4_p3 = StatsPlots.heatmap(λ, θ, mean_surplus_difference', xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta badającego konsumentów", titlefontsize = 8)
nl=200

ex3_4_eng = Plots.plot(sort(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))), (1:length(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))))./length(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))), xlabel = "Profit of firms", ylabel = "Empirical Cumulative Distribution Function", title = "Emprical Cumulative Distribution Function of firms' profits", label = "Firm doesn't invest in market research", legend=:bottomright, legendfontsize = 12)
Plots.plot!(sort(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))), (1:length(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))))./length(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))), label = "Firm invests in market research", titlefontsize = 14, xlabelfontsize = 14, ylabelfontsize = 14, size = (750, 550))

ex4_p5 = Plots.scatter(sort(unique(ex3_v200201_nl[1:1020])) / 200, [mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[ex3_v200201_nl[1:1020] .== nl]) for nl in sort(unique(ex3_v200201_nl))], xlabel = "Number of connections per single consumer", ylabel = "Additional profit from market research", title = "Market research incremental effect", color = "red", legend = false, titlefontsize = 14, xlabelfontsize = 14, ylabelfontsize = 14, size = (750, 550))

mod = GLM.lm(@formula(y~x), DataFrame(x = 1 ./ (sort(unique(ex3_v200201_nl[1:1020]))/200), y = [mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[ex3_v200201_nl[1:1020] .== nl]) for nl in sort(unique(ex3_v200201_nl))]))

using Measures

Plots.plot!(coef(mod)[1] .+ coef(mod)[2] ./ (sort(unique(ex3_v200201_nl[1:1020]))/200))

p12 = Plots.plot(ex3_4_eng, ex4_p5, layout=(1,2), size = (1500, 550), margin = 15mm)

savefigs(p12, "\\plots\\ex1\\SATU")

# największa różnica pomiędzy zyskami producentów robiących badania i nie jest wtedy, 

########################################################################

Plots.scatter()

#########################################################################



TR200 = mean.(multi.(getindex.(ex3_v200_price, 1), getindex.(ex3_v200_quantity_sold, 1)))
TR201 = mean.(multi.(getindex.(ex3_v201_price, 1), getindex.(ex3_v201_quantity_sold, 1)))

Plots.plot(sort(mean.(getindex.(ex3_v200_price, 1))), (1:length(mean.(getindex.(ex3_v200_price, 1))))./length(mean.(getindex.(ex3_v200_price, 1))), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)
Plots.plot!(sort(mean.(getindex.(ex3_v201_price, 1))), (1:length(mean.(getindex.(ex3_v201_price, 1))))./length(mean.(getindex.(ex3_v201_price, 1))), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)

Plots.plot(sort(mean.(getindex.(ex3_v200_quantity_sold, 1))), (1:length(mean.(getindex.(ex3_v200_quantity_sold, 1))))./length(mean.(getindex.(ex3_v200_quantity_sold, 1))), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)
Plots.plot!(sort(mean.(getindex.(ex3_v201_quantity_sold, 1))), (1:length(mean.(getindex.(ex3_v201_quantity_sold, 1))))./length(mean.(getindex.(ex3_v201_quantity_sold, 1))), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)

Plots.plot(sort(TR200), (1:length(TR200))./length(TR200), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)
Plots.plot!(sort(TR201), (1:length(TR201))./length(TR201), xlabel = "Profit", ylabel = "ECDF", title = "ECDF - Profit", label = "No research", legend=:bottomright, legendfontsize = 8)



#########################################################################

bounds = collect(-50:5:100)

function cut_integer_bounds(x::Vector{Float64},bounds::Vector)

    value_counts = Int64[]

    for item in eachindex(bounds)[Not(length(bounds))]
        push!(value_counts, sum((x .>= bounds[item]) .& (x .< bounds[item+1])))
    end

    return value_counts

end

function earthmoverdistance(a::Vector, b::Vector) 
    return sum( abs, cumsum( a ) .- cumsum( b ) )
end

emd1 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], bounds)
) / sum(ex3_v200201_nl .<= k1)

emd2 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], bounds)
) / sum((ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2))

emd3 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], bounds)
) / sum((ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3))

emd4 = earthmoverdistance(
    cut_integer_bounds(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], bounds),
    cut_integer_bounds(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], bounds)
) / sum(ex3_v200201_nl .>= k3)

median_distance(x,y) = median(x) - median(y)
mean_distnace(x,y) = mean(x) - mean(y)

emd_df = DataFrame(
    Stan = ["k1", "k1+200,1200", "k2 + 200, 1k1+200", "k3"],
    EMD = round.([emd1, emd2, emd3, emd4], digits = 2),
    Mediana = round.([median_distance(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1]),
    median_distance(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)]),
    median_distance(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)]),
    median_distance(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3])], digits = 2),
    Średnia = round.([mean_distnace(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .<= k1]),
    mean_distnace(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k1+200) .& (ex3_v200201_nl .<= k2)]),
    mean_distnace(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= k2 + 200) .& (ex3_v200201_nl .<= k3)]),
    mean_distnace(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .>= k3])], digits = 2)
)

latexify(emd_df, env = :table) |> print

############################################################################################


plot_ecdf(true, getindex.([var.(x) for x in ex3_v200_producer_surplus_singleton], 1), "", xlabel = "Var (" * L"\pi" * ")", ylabel = "F(x)", title = "Wariancja zysku " * L"\pi")
plot_ecdf(false, getindex.([var.(x) for x in ex3_v201_producer_surplus_singleton], 1), "")

###########################################################################################

ex3_11_pl = plot_ecdf(true, sum.(getindex.(ex3_v200_quantity_produced, 1)) .- sum.(getindex.(ex3_v200_quantity_sold, 1)), "", xlabel = "Zysk producenta", ylabel = "F(x)", title = "1000-1200")
plot_ecdf(false, sum.(getindex.(ex3_v201_quantity_produced, 1)) .- sum.(getindex.(ex3_v201_quantity_sold, 1)), "", xlim = (-1000, 2000))

##################################################################



#################################################################3

Plots.boxplot([sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_sm .== sm] for sm in sort(unique(ex3_v200201_sm))])
Plots.boxplot([sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_sm .== sm] for sm in sort(unique(ex3_v200201_sm))])

Plots.scatter(ex3_v200201_ρm, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.scatter(ex3_v200201_ρm, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v200201_ρm, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v200201_cu, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v200201_cu, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(ex3_v200201_cu, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v200201_Li, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v200201_Li, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(ex3_v200201_Li, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

Plots.scatter(ex3_v200201_Lw, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter(ex3_v200201_Lw, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), smooth = true)

##########################################################################################


round(iqr(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp])), digits = 3)
round(iqr(xydiff.(getindex.(ex3_v201_quality, 2), [getindex.(x,2) for x in ex3_v201_quality_exp])), digits = 3)
round(cv(xydiff.(getindex.(ex3_v201_quality, 2), [getindex.(x,2) for x in ex3_v201_quality_exp])), digits = 3)

p_ex3_10 = Plots.boxplot(xydiff.(getindex.(ex3_v200_quality, 1), [getindex.(x,1) for x in ex3_v200_quality_exp]), ylabel = "Odchylenie parametrów produktu od oczekiwań", legend = nothing, ylabelfontsize = 8, dpi = 300, title = "Jakość")
Plots.boxplot!(xydiff.(getindex.(ex3_v201_quality, 1), [getindex.(x,1) for x in ex3_v201_quality_exp]))
Plots.plot!(xticks = ([1,2], ["Producent nie prowadzi badań", "Producent prowadzi badania"]))

#%%%

Plots.plot(mean(getindex.(ex3_v200_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v201_quality, 1)))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v200_quality_exp]))
Plots.plot!(mean([getindex.(x,1) for x in ex3_v201_quality_exp]))

Plots.plot(ex3_v200_quality[1][1])

getindex.(ex3_v200_quality_exp,1)


# Trend zmian jakości i trwałości

Plots.plot(mean(getindex.(ex3_v200_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v201_quality, 1)))
Plots.plot!(mean(getindex.(ex3_v200_quality, 2)))
Plots.plot!(mean(getindex.(ex3_v201_quality, 2)))

Plots.plot(mean(getindex.(ex3_v200_durability, 1)))
Plots.plot!(mean(getindex.(ex3_v201_durability, 1)))
Plots.plot!(mean(getindex.(ex3_v200_durability, 2)))
Plots.plot!(mean(getindex.(ex3_v201_durability, 2)))

#%%%%

Plots.plot((ex3_v200_quality)[1][1])
Plots.plot!((ex3_v201_quality)[1][1])

StatsPlots.density(mean.(ex3_v200_sv))
StatsPlots.density!(mean.(ex3_v201_sv))

ApproximateTwoSampleKSTest(mean.(ex3_v200_sv), mean.(ex3_v201_sv))

Plots.scatter(mean.(ex3_v200_sv), sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), smooth = true)
Plots.scatter!(mean.(ex3_v201_sv), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), smooth = true)

Plots.scatter(mean.(ex3_v200_sv)[ex3_v200201_H .< 3], sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .< 3], smooth = true)
Plots.scatter!(mean.(ex3_v201_sv), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), smooth = true)

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .== x] for x in sort(unique(ex3_v200201_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 3000))

ex1_p10 = StatsPlots.boxplot([sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .== x] for x in sort(unique(ex3_v200201_H))], legend = false, xlabel = "Maksymalna przydatność dobra [liczba okresów]", ylabel = "Zysk firmy", title = "Zysk firmy, która nie bada konsumentów", ylim = (-1500, 3000))

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v200_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 2)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)), "Producent bada oczekiwania konsumentów")

UnequalVarianceTTest(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .< 5], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_H .>= 5], "Producent bada oczekiwania konsumentów")

#########################################################################################

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .<= 300], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .<= 300], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= 400) .& (ex3_v200201_nl .<= k1)], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= 400) .& (ex3_v200201_nl .<= k1)], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= 700) .& (ex3_v200201_nl .<= 900)], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[(ex3_v200201_nl .>= 700) .& (ex3_v200201_nl .<= 900)], "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))[ex3_v200201_nl .>= 1000], "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))[ex3_v200201_nl .>= 1000], "Producent bada oczekiwania konsumentów")

####

mean(sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))
mean(sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

mean(sum.(getindex.(ex3_v201_quantity_sold, 1)))
mean(sum.(getindex.(ex3_v200_quantity_sold, 1)))

mean(sum.(getindex.(ex3_v201_quantity_produced, 1)))
mean(sum.(getindex.(ex3_v200_quantity_produced, 1)))

mean(mean.(getindex.(ex3_v201_price, 1)))
mean(mean.(getindex.(ex3_v200_price, 1)))

mean(mean.(getindex.(ex3_v201_price, 1))) / mean(mean.(getindex.(ex3_v201_margin, 1)))
mean(mean.(getindex.(ex3_v200_price, 1))) / mean(mean.(getindex.(ex3_v200_margin, 1)))

[sum(ex3_v201_quantity_sold[k][1] .* ex3_v201_price[k][1]) - sum(ex3_v201_quantity_produced[k][1] .* ex3_v201_price[k][1] ./ ex3_v201_margin[k][1]) .+ sum((ex3_v201_quantity_produced[k][1] .- ex3_v201_quantity_sold[k][1]) * 0.5 .* ex3_v201_price[k][1] ./ ex3_v201_margin[k][1]) for k in 1:400] .- sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))

mean()

k=3
sum(ex3_v201_quantity_sold[k][1] .* ex3_v201_price[k][1]) - sum(ex3_v201_quantity_produced[k][1] .* ex3_v201_price[k][1] ./ ex3_v201_margin[k][1]) .+ sum((ex3_v201_quantity_produced[k][1] .- ex3_v201_quantity_sold[k][1]) * 0.5 .* ex3_v201_price[k][1] ./ ex3_v201_margin[k][1])

sum(ex3_v200_quantity_sold[k][1] .* ex3_v200_price[k][1]) - sum(ex3_v200_quantity_produced[k][1] .* ex3_v200_price[k][1] ./ ex3_v200_margin[k][1]) .+ sum((ex3_v200_quantity_produced[k][1] .- ex3_v200_quantity_sold[k][1]) * 0.5 .* ex3_v200_price[k][1] ./ ex3_v200_margin[k][1])

calculate_utility(k,d,h) = k .* (1 .- d .^ h) ./ (1 .- d)

Plots.plot(calculate_utility(ex3_v201_quality[k][1], ex3_v201_durability[k][1], ex3_v200201_H[k]))
Plots.plot!(calculate_utility(ex3_v200_quality[k][1], ex3_v200_durability[k][1], ex3_v200201_H[k]))

Plots.plot!(calculate_utility(getindex.(ex3_v201_quality_exp[k], 1), getindex.(ex3_v201_durability_exp[k], 1), ex3_v200201_H[k]))
Plots.plot!(calculate_utility(getindex.(ex3_v200_quality_exp[k], 1), getindex.(ex3_v200_durability_exp[k], 1), ex3_v200201_H[k]))

Plots.plot(ex3_v201_price[k][1])
Plots.plot!(ex3_v200_price[k][1])

Plots.plot(calculate_utility(ex3_v201_quality[k][1], ex3_v201_durability[k][1], ex3_v200201_H[k]) ./ ex3_v201_price[k][1])
Plots.plot!(calculate_utility(ex3_v200_quality[k][1], ex3_v200_durability[k][1], ex3_v200201_H[k]) ./ ex3_v200_price[k][1])

Plots.plot!(calculate_utility(getindex.(ex3_v201_quality_exp[k], 1), getindex.(ex3_v201_durability_exp[k], 1), ex3_v200201_H[k]) ./ ex3_v201_price[k][1])
Plots.plot!(calculate_utility(getindex.(ex3_v200_quality_exp[k], 1), getindex.(ex3_v200_durability_exp[k], 1), ex3_v200201_H[k]) ./ ex3_v200_price[k][1])

CSV.write(pwd() * "\\data_test.csv", DataFrame(b_qs = ex3_v201_quantity_sold[k][1], b_p = ex3_v201_price[k][1], b_qp = ex3_v201_quantity_produced[k][1], b_m = ex3_v201_margin[k][1], nb_qs = ex3_v200_quantity_sold[k][1], nb_p = ex3_v200_price[k][1], nb_qp = ex3_v200_quantity_produced[k][1], nb_m = ex3_v200_margin[k][1]))

#####

function consecutive_negative(x)
    current_counter = 0
    best_counter = 0
    for xi in x
        if xi < 0
            current_counter += 1
        else
            if current_counter > best_counter
                best_counter = current_counter
            end

            current_counter = 0

        end
    end
    return best_counter
end

StatsPlots.density(consecutive_negative.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
StatsPlots.density!(consecutive_negative.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

####

ex3_p1 = plot_ecdf(true, mean.(getindex.(ex3_v200_quality, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Jakość producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_quality, 1)), "Producent bada oczekiwania konsumentów")

ex3_p1 = plot_ecdf(true, mean.(getindex.(ex3_v200_durability, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Jakość producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, mean.(getindex.(ex3_v201_durability, 1)), "Producent bada oczekiwania konsumentów")



plot_ecdf(true, sum.(ex3_v200_total_surplus), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(ex3_v201_total_surplus), "Producent bada oczekiwania konsumentów")

plot_ecdf(true, sum.(ex3_v200_consumer_surplus), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(ex3_v201_consumer_surplus), "Producent bada oczekiwania konsumentów")

#%%%

function speed_of_learing(x, xe)
    dx = diff(x)
    dxe = diff(xe)
    return sum(sign.(dx) .== sign.(dxe))
end

RMSE.([getindex.(x,1) for x in ex3_v200_quality_exp], getindex.(ex3_v200_quality, 1))

Plots.scatter(RMSE.([getindex.(x,1) for x in ex3_v200_quality_exp], getindex.(ex3_v200_quality, 1)), sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.scatter!(RMSE.([getindex.(x,1) for x in ex3_v201_quality_exp], getindex.(ex3_v201_quality, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

Plots.scatter(speed_of_learing.([getindex.(x,1) for x in ex3_v200_quality_exp], getindex.(ex3_v200_quality, 1)), sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.scatter!(speed_of_learing.([getindex.(x,1) for x in ex3_v201_quality_exp], getindex.(ex3_v201_quality, 1)), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

Plots.scatter(std.([getindex.(x,1) for x in ex3_v200_quality_exp]), sum.(getindex.(ex3_v200_producer_surplus_singleton, 1)))

Plots.scatter!(std.([getindex.(x,1) for x in ex3_v201_quality_exp]), sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

sv200 = cut_integer(mean.(ex3_v200_sv), 10)[1]
sv201 = cut_integer(mean.(ex3_v201_sv), 10)[1]

[std(mean.(getindex.(ex3_v200_producer_surplus_singleton, 1))[sv200 .== x]) for x in sort(unique(sv200))]

[std(mean.(getindex.(ex3_v201_producer_surplus_singleton, 1))[sv201 .== x]) for x in sort(unique(sv201))]

Plots.scatter(mean.(ex3_v200_sv), mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))
Plots.scatter!(mean.(ex3_v201_sv), mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

mean.(ex3_v200_sv) .- mean.(ex3_v201_sv)

Plots.scatter(mean.(ex3_v200_sv), mean.(ex3_v201_sv), smooth = true, xlim = (0,150), ylim = (0,150))


Plots.scatter(mean.(ex3_v200201_Lw), mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)))

λ = cut_integer(float.(ex3_v200201_Li), 5)[1]
θ = cut_integer(float.(ex3_v200201_Lw), 5)[1]

hm_rs = [mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[(λ .== l) .& (θ .== t)]) for l in sort(unique(λ)), t in sort(unique(θ))]

mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[(λ .== 2) .& (θ .== 3)])

mean((mean.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- mean.(getindex.(ex3_v200_producer_surplus_singleton, 1)))[(λ .== 3) .& (θ .== 2)])

Plots.plot([mean(x) for x in eachrow(hm_rs)]) # by λ
Plots.plot([mean(x) for x in eachcol(hm_rs)]) # by θ

ex4_p1 = StatsPlots.heatmap(sort(unique(λ)), sort(unique(θ)), hm_rs', xlabel = "λ, produkty oceniane osobiście", ylabel = "θ, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita", titlefontsize = 8)

####

Y = sum.(getindex.(ex3_v201_producer_surplus_singleton, 1)) .- sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))

Y1 = sum.(getindex.(ex3_v200_producer_surplus_singleton, 1))
Y2 = sum.(getindex.(ex3_v201_producer_surplus_singleton, 1))

H = float.(ex3_v200201_H)
λ = float.(ex3_v200201_Li)
θ = float.(ex3_v200201_Lw)
sd = float.(ex3_v200201_sd)
cu = float.(ex3_v200201_cu)
sm = float.(ex3_v200201_sm)
nl = float.(ex3_v200201_nl)
rm = float.(ex3_v200201_ρm)


k1 = float.(mean.(getindex.(ex3_v200_quality, 1)))
k2 = float.(mean.(getindex.(ex3_v201_quality, 1)))

Y = k2.-k1

d1 = float.(mean.(getindex.(ex3_v200_durability, 1)))
d2 = float.(mean.(getindex.(ex3_v201_durability, 1)))

m1 = float.(mean.(getindex.(ex3_v200_margin, 1)))
m2 = float.(mean.(getindex.(ex3_v201_margin, 1)))

r1 = mean.(getindex.(ex3_v200_reselling, 1))
r2 = mean.(getindex.(ex3_v201_reselling, 1))

ek1 = float.(mean.([getindex.(x,1) for x in ex3_v200_quality_exp]))
ek2 = float.(mean.([getindex.(x,1) for x in ex3_v201_quality_exp]))

ed1 = float.(mean.([getindex.(x,1) for x in ex3_v200_durability_exp]))
ed2 = float.(mean.([getindex.(x,1) for x in ex3_v201_durability_exp]))

optimal_research_100 = GLM.lm(@formula(Y1 ~ λ + θ + k1 + d1 + r1 + ek1 + ed1), DataFrame(Y1 = Y1, H = H, λ = λ, θ = θ, sd = sd, k1 = k1, k2 = k2, d1 = d1, d2 = d2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, ek1 = ek1, ek2 = ek2, ed1 = ed1, ed2 = ed2))

round(adjr2(optimal_research_100), digits = 2)

optimal_research_101 = GLM.lm(@formula(Y2 ~ λ + θ + k2 + d2 + r2 + ek2 + ed2), DataFrame(Y1 = Y1, Y2 = Y2, H = H, λ = λ, θ = θ, sd = sd, k1 = k1, k2 = k2, d1 = d1, d2 = d2, m1 = m1, m2 = m2, r1 = r1, r2 = r2, ek1 = ek1, ek2 = ek2, ed1 = ed1, ed2 = ed2))

round(adjr2(optimal_research_101), digits = 2)



####