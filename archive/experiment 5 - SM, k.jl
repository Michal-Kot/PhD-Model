# Experiment 5

include(pwd() * "\\methods\\methods.jl")

# Impact of secondary market on consumer beliefs

# Dlaczego wyższe durability powoduje spadek consumer surplus?

ex5_v1_total_surplus = []
ex5_v1_producer_surplus = []
ex5_v1_consumer_surplus = []
ex5_v1_price = []
ex5_v1_quality = []
ex5_v1_durability = []
ex5_v1_margin = []
ex5_v1_quantity_produced = []
ex5_v1_quantity_sold = []
ex5_v1_quantity_leased = []
ex5_v1_producer_surplus_singleton = []
ex5_v1_reselling = []
ex5_v1_consumer_surplus_pm = []
ex5_v1_consumer_surplus_sm_s = []
ex5_v1_consumer_surplus_sm_b = []
ex5_v1_quality_exp = []
ex5_v1_durability_exp = []

ex5_v2_total_surplus = []
ex5_v2_producer_surplus = []
ex5_v2_consumer_surplus = []
ex5_v2_price = []
ex5_v2_quality = []
ex5_v2_durability = []
ex5_v2_margin = []
ex5_v2_quantity_produced = []
ex5_v2_quantity_sold = []
ex5_v2_quantity_leased = []
ex5_v2_producer_surplus_singleton = []
ex5_v2_reselling = []
ex5_v2_consumer_surplus_pm = []
ex5_v2_consumer_surplus_sm_s = []
ex5_v2_consumer_surplus_sm_b = []
ex5_v2_quality_exp = []
ex5_v2_durability_exp = []
ex5_v1_H = []
ex5_v1_Li = []
ex5_v1_Lw = []


for i in 1:500

    if (mod(i,10) == 0) | (i == 1) | (i == 2)
        println(i)
    end

    li = rand()
    lw = rand()
    h = sample(2:8)

    push!(ex5_v1_H, h)
    push!(ex5_v1_Li, li)
    push!(ex5_v1_Lw, lw)

    ex5_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.0, 1.0], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, true, 1, [0.7, 1.], "softmax", [false, true], [0.0, 0.1], h, false)

    push!(ex5_v1_total_surplus, calculate_surplus(ex5_v1_sim, "total", false))
    push!(ex5_v1_producer_surplus, calculate_surplus(ex5_v1_sim, "producer", false))
    push!(ex5_v1_consumer_surplus, calculate_surplus(ex5_v1_sim, "consumer,total", false))
    push!(ex5_v1_price, calculate_price_history.(ex5_v1_sim.sellers; product_life = h))
    push!(ex5_v1_quality, getfield.(ex5_v1_sim.sellers, :quality_history))
    push!(ex5_v1_durability, getfield.(ex5_v1_sim.sellers, :durability_history))
    push!(ex5_v1_margin, getfield.(ex5_v1_sim.sellers, :margin_history))
    push!(ex5_v1_quantity_produced, getfield.(ex5_v1_sim.sellers, :quantity_produced_history))
    push!(ex5_v1_quantity_sold, getfield.(ex5_v1_sim.sellers, :quantity_sold_history))
    push!(ex5_v1_producer_surplus_singleton, calculate_profit_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_reselling, getfield.(ex5_v1_sim.sellers, :reselling_history))
    push!(ex5_v1_consumer_surplus_pm, calculate_surplus(ex5_v1_sim, "consumer,pm", false))
    push!(ex5_v1_consumer_surplus_sm_b, calculate_surplus(ex5_v1_sim, "consumer,sm,b", false))
    push!(ex5_v1_consumer_surplus_sm_s, calculate_surplus(ex5_v1_sim, "consumer,sm,s", false))
    push!(ex5_v1_quality_exp, mean([b.quality_expectation_history for b in ex5_v1_sim.buyers]))
    push!(ex5_v1_durability_exp, mean([b.durability_expectation_history for b in ex5_v1_sim.buyers]))

    ex5_v2_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.0, 1.0], "random", li, lw, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[.8, 2.], [.8, 2.]], 0.10, false, 1, [0.7, 1.], "softmax", [false, true], [0.0, 0.1], h, false)
    push!(ex5_v2_total_surplus, calculate_surplus(ex5_v2_sim, "total", false))
    push!(ex5_v2_producer_surplus, calculate_surplus(ex5_v2_sim, "producer", false))
    push!(ex5_v2_consumer_surplus, calculate_surplus(ex5_v2_sim, "consumer,total",false))
    push!(ex5_v2_price, calculate_price_history.(ex5_v2_sim.sellers; product_life = h))
    push!(ex5_v2_quality, getfield.(ex5_v2_sim.sellers, :quality_history))
    push!(ex5_v2_durability, getfield.(ex5_v2_sim.sellers, :durability_history))
    push!(ex5_v2_margin, getfield.(ex5_v2_sim.sellers, :margin_history))
    push!(ex5_v2_quantity_produced, getfield.(ex5_v2_sim.sellers, :quantity_produced_history))
    push!(ex5_v2_quantity_sold, getfield.(ex5_v2_sim.sellers, :quantity_sold_history))
    push!(ex5_v2_producer_surplus_singleton, calculate_profit_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_reselling, getfield.(ex5_v2_sim.sellers, :reselling_history))
    push!(ex5_v2_consumer_surplus_pm, calculate_surplus(ex5_v2_sim, "consumer,pm", false))
    push!(ex5_v2_consumer_surplus_sm_b, calculate_surplus(ex5_v2_sim, "consumer,sm,b", false))
    push!(ex5_v2_consumer_surplus_sm_s, calculate_surplus(ex5_v2_sim, "consumer,sm,s", false))
    push!(ex5_v2_quality_exp, mean([b.quality_expectation_history for b in ex5_v2_sim.buyers]))
    push!(ex5_v2_durability_exp, mean([b.durability_expectation_history for b in ex5_v2_sim.buyers]))

end

ex3_p1 = plot_ecdf(true, sum.(getindex.(ex5_v1_producer_surplus_singleton, 1)), "Producent nie bada oczekiwań konsumentów", xlabel = "Zysk producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - zysk producenta")
plot_ecdf(false, sum.(getindex.(ex5_v1_producer_surplus_singleton, 2)) .- 400 * 0.05 * 200 * 0.001, "Producent bada oczekiwania konsumentów")

calc_u(k,d,h) = k * (1-d^h)/(1-d)

u_v1_1 = calc_u.(mean.(getindex.(ex5_v1_quality, 1)), mean.(getindex.(ex5_v1_durability, 1)), ex5_v1_H)
u_v1_2 = calc_u.(mean.(getindex.(ex5_v1_quality, 2)), mean.(getindex.(ex5_v1_durability, 2)), ex5_v1_H)
u_v2_1 = calc_u.(mean.(getindex.(ex5_v2_quality, 1)), mean.(getindex.(ex5_v2_durability, 1)), ex5_v1_H)
u_v2_2 = calc_u.(mean.(getindex.(ex5_v2_quality, 2)), mean.(getindex.(ex5_v2_durability, 2)), ex5_v1_H)

# Wysoka użyteczność

v1_best_rn = u_v1_1 .== max.(u_v1_1, u_v1_2)
v1_worst_rn = .!(v1_best_rn)
plot_ecdf(true, mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)[v1_best_rn]), "RW tak, WJ tak, CR nie"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
v1_best_ry = u_v1_2 .== max.(u_v1_1, u_v1_2)
v1_worst_ry = .!(v1_best_ry)
plot_ecdf(false, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2)[v1_best_ry]), "RW tak, WJ tak, CR tak")
v2_best_rn = u_v2_1 .== max.(u_v2_1, u_v2_2)
v2_worst_rn = .!(v2_best_rn)
plot_ecdf(false, mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)[v2_best_rn]), "RW nie, WJ tak, CR nie")
v2_best_ry = u_v2_2 .== max.(u_v2_1, u_v2_2)
v2_worst_ry = .!(v2_best_ry)
plot_ecdf(false, mean.(getindex.(ex5_v2_producer_surplus_singleton, 2)[v2_best_ry]), "RW nie, WJ tak, CR tak")

# Niska użyteczność

plot_ecdf(true, mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)[v1_worst_rn]), "RW tak, WJ nie, CR nie"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2)[v1_worst_ry]), "RW tak, WJ nie, CR tak")
plot_ecdf(false, mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)[v2_worst_rn]), "RW nie, WJ nie, CR nie")
plot_ecdf(false, mean.(getindex.(ex5_v2_producer_surplus_singleton, 2)[v2_worst_ry]), "RW nie, WJ nie, CR tak")


#####

high_lw = ex5_v1_Lw .> 0.5

plot_ecdf(true, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_best_rn)[high_lw]), "wysoka λ, wysoka jakość"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_best_rn)[.!high_lw]), "niska λ, wysoka jakość")

UnequalVarianceTTest(mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_best_rn)[high_lw]), mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_best_rn)[.!high_lw]))

plot_ecdf(true, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_worst_rn)[high_lw]), "wysoka λ, niska jakość"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_worst_rn)[.!high_lw]), "niska λ, niska jakość")

UnequalVarianceTTest(mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_worst_rn)[high_lw]), mean.(getindex.(ex5_v1_producer_surplus_singleton, 2 .- v1_worst_rn)[.!high_lw]))

plot_ecdf(true, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2)[high_lw]), "RW tak, WJ tak, CR nie"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex5_v1_producer_surplus_singleton, 2)[.!high_lw]), "RW tak, WJ tak, CR tak")

plot_ecdf(true, mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)[high_lw]), "RW tak, WJ tak, CR nie"; xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)[.!high_lw]), "RW tak, WJ tak, CR tak")

