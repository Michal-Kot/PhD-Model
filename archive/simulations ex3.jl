include(pwd() * "\\methods\\methods.jl")

ex3_v1_total_surplus = []
ex3_v1_producer_surplus = []
ex3_v1_consumer_surplus = []
ex3_v1_price = []
ex3_v1_quantity = []
ex3_v1_advertising = []
ex3_v1_producer_surplus_singleton = []
ex3_v1_quality_expectation = []
ex3_v1_durability_expectation = []
ex3_v1_margin = []

ex3_v6_total_surplus = []
ex3_v6_producer_surplus = []
ex3_v6_consumer_surplus = []
ex3_v6_price = []
ex3_v6_quantity = []
ex3_v6_advertising = []
ex3_v6_advertising_highest = []
ex3_v6_producer_surplus_singleton = []
ex3_v6_quality_expectation = []
ex3_v6_durability_expectation = []
ex3_v6_margin = []

ex3_v16_d = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = fill(rand(m_d), 2)

    d_init = sample(1:5,2)

    push!(ex3_v16_d, d_init)

    # 1.0

    ex3_v1_sim = TO_GO(2, 500, 250, 0.25, 0.25; q = [1.0, 1.0], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = nl, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex3_v1_total_surplus, calculate_total_surplus(ex3_v1_sim, "total"))
    push!(ex3_v1_producer_surplus, calculate_total_surplus(ex3_v1_sim, "producer"))
    push!(ex3_v1_consumer_surplus, calculate_total_surplus(ex3_v1_sim, "consumer"))
    push!(ex3_v1_price, mean(calculate_price_history.(ex3_v1_sim.sellers)))
    push!(ex3_v1_quantity, mean(getfield.(ex3_v1_sim.sellers, :quantity_history)))
    push!(ex3_v1_advertising, mean(getfield.(ex3_v1_sim.sellers, :advertising_history)))
    push!(ex3_v1_producer_surplus_singleton, sum.(calculate_profit_history.(ex3_v1_sim.sellers)))
    push!(ex3_v1_quality_expectation, calculate_expectation(ex3_v1_sim, "quality", true))
    push!(ex3_v1_durability_expectation, calculate_expectation(ex3_v1_sim, "durability", true))
    push!(ex3_v1_margin, getfield.(ex3_v1_sim.sellers, :margin_history))     

    # 1.5

    ex3_v6_sim = TO_GO(2, 500, 250, 0.25, 0.25; q = [1.5, 0.5], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = nl, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex3_v6_total_surplus, calculate_total_surplus(ex3_v6_sim, "total"))
    push!(ex3_v6_producer_surplus, calculate_total_surplus(ex3_v6_sim, "producer"))
    push!(ex3_v6_consumer_surplus, calculate_total_surplus(ex3_v6_sim, "consumer"))
    push!(ex3_v6_price, mean(calculate_price_history.(ex3_v6_sim.sellers)))
    push!(ex3_v6_quantity, mean(getfield.(ex3_v6_sim.sellers, :quantity_history)))
    push!(ex3_v6_advertising, mean(getfield.(ex3_v6_sim.sellers, :advertising_history)))
    push!(ex3_v6_producer_surplus_singleton, sum.(calculate_profit_history.(ex3_v6_sim.sellers)))
    push!(ex3_v6_quality_expectation, calculate_expectation(ex3_v6_sim, "quality", true))
    push!(ex3_v6_durability_expectation, calculate_expectation(ex3_v6_sim, "durability", true))
    push!(ex3_v6_margin, getfield.(ex3_v6_sim.sellers, :margin_history))   
  

end

scatter(mean.(ex3_v16_d), ex3_v1_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex3_v16_lambda, 1), Float64.(ex3_v1_total_surplus), "blue",0.50)
scatter!(mean.(ex3_v16_d), ex3_v6_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex3_v16_lambda, 1), Float64.(ex3_v6_total_surplus), "green")

scatter(ex3_v16_nl, getindex.(ex3_v6_producer_surplus_singleton,1), markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex3_v16_lambda, 1), Float64.(ex3_v1_total_surplus), "blue",0.50)
scatter!(ex3_v16_nl, getindex.(ex3_v6_producer_surplus_singleton,2), markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex3_v16_lambda, 1), Float64.(ex3_v6_total_surplus), "green")


L1 = getindex.(ex3_v16_lambda, 1) 
L2 = getindex.(ex3_v16_lambda, 2)

L1, L1u = cut_integer(L1, 20)
L2, L2u = cut_integer(L2, 20)

mean_c(x) = length(x) > 0 ? mean(x) : NaN

mean_c(ex3_v1_total_surplus[(L1 .== 1) .& (L2 .== 20)])

ex4_p1 = heatmap(L1u, L2u, [mean_c(ex3_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]', xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość", titlefontsize = 10, clim=(25000, 65000))

savefig(ex4_p1, pwd() * "\\plots\\ex4_equal quality heatmap lammbdas.svg")

ex4_p2 = heatmap(L1u, L2u, [mean(ex3_v6_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]', xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, różna jakość", titlefontsize = 10, clim=(25000,65000))

savefig(ex4_p2, pwd() * "\\plots\\ex4_different quality heatmap lammbdas.svg")

scatter(getindex.(ex3_v16_lambda, 2), ex3_v1_total_surplus)
scatter!(getindex.(ex3_v16_lambda, 2), ex3_v6_total_surplus)

scatter(getindex.(ex3_v16_lambda, 1), getindex.(ex3_v6_durability_expectation, 1), markeralpha = 0.25)

scatter(getindex.(ex3_v16_lambda, 2), getindex.(ex3_v1_producer_surplus_singleton, 2))

