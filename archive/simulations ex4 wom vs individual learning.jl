include(pwd() * "\\methods\\methods.jl")


ex4_v1_total_surplus = []
ex4_v1_producer_surplus = []
ex4_v1_consumer_surplus = []
ex4_v1_price = []
ex4_v1_quantity = []
ex4_v1_advertising = []
ex4_v1_producer_surplus_singleton = []
ex4_v1_quality_expectation = []
ex4_v1_durability_expectation = []
ex4_v1_margin = []

ex4_v2_total_surplus = []
ex4_v2_producer_surplus = []
ex4_v2_consumer_surplus = []
ex4_v2_price = []
ex4_v2_quantity = []
ex4_v2_advertising = []
ex4_v2_advertising_highest = []
ex4_v2_producer_surplus_singleton = []
ex4_v2_quality_expectation = []
ex4_v2_durability_expectation = []
ex4_v2_margin = []

ex4_v16_lambda = []

for i in 1:10000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = fill(rand(m_d), 2)

    d_init = [2,2]

    λ_ind = rand()
    λ_wom = rand()

    push!(ex4_v16_lambda, [λ_ind, λ_wom])

    # 1.0

    ex4_v1_sim = TO_GO(2, 500, 250, λ_ind, λ_wom; q = [1.0, 1.0], m = m_init, c = [0.5,0.5], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex4_v1_total_surplus, calculate_total_surplus(ex4_v1_sim, "total"))
    push!(ex4_v1_producer_surplus, calculate_total_surplus(ex4_v1_sim, "producer"))
    push!(ex4_v1_consumer_surplus, calculate_total_surplus(ex4_v1_sim, "consumer"))
    push!(ex4_v1_price, mean(calculate_price_history.(ex4_v1_sim.sellers)))
    push!(ex4_v1_quantity, mean(getfield.(ex4_v1_sim.sellers, :quantity_history)))
    push!(ex4_v1_advertising, mean(getfield.(ex4_v1_sim.sellers, :advertising_history)))
    push!(ex4_v1_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v1_sim.sellers)))
    push!(ex4_v1_quality_expectation, calculate_expectation(ex4_v1_sim, "quality", true))
    push!(ex4_v1_durability_expectation, calculate_expectation(ex4_v1_sim, "durability", true))
    push!(ex4_v1_margin, getfield.(ex4_v1_sim.sellers, :margin_history))     

    # 1.5

    ex4_v2_sim = TO_GO(2, 500, 250, λ_ind, λ_wom; q = [1.5, 0.5], m = m_init, c = [0.5,0.5], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex4_v2_total_surplus, calculate_total_surplus(ex4_v2_sim, "total"))
    push!(ex4_v2_producer_surplus, calculate_total_surplus(ex4_v2_sim, "producer"))
    push!(ex4_v2_consumer_surplus, calculate_total_surplus(ex4_v2_sim, "consumer"))
    push!(ex4_v2_price, mean(calculate_price_history.(ex4_v2_sim.sellers)))
    push!(ex4_v2_quantity, mean(getfield.(ex4_v2_sim.sellers, :quantity_history)))
    push!(ex4_v2_advertising, mean(getfield.(ex4_v2_sim.sellers, :advertising_history)))
    push!(ex4_v2_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v2_sim.sellers)))
    push!(ex4_v2_quality_expectation, calculate_expectation(ex4_v2_sim, "quality", true))
    push!(ex4_v2_durability_expectation, calculate_expectation(ex4_v2_sim, "durability", true))
    push!(ex4_v2_margin, getfield.(ex4_v2_sim.sellers, :margin_history))   
  

end

scatter(getindex.(ex4_v16_lambda, 1), ex4_v1_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v1_total_surplus), "blue",0.50)
scatter!(getindex.(ex4_v16_lambda, 1), ex4_v2_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v2_total_surplus), "green")

scatter(getindex.(ex4_v16_lambda, 2), ex4_v1_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex4_v16_lambda, 2), Float64.(ex4_v1_total_surplus), "blue")
scatter!(getindex.(ex4_v16_lambda, 2), ex4_v2_total_surplus, markeralpha = 0.25, smooth = true)
#add_smoothing_spline(getindex.(ex4_v16_lambda, 2), Float64.(ex4_v2_total_surplus), "green")

L1 = getindex.(ex4_v16_lambda, 1) 
L2 = getindex.(ex4_v16_lambda, 2)

L1, L1u = cut_integer(L1, 20)
L2, L2u = cut_integer(L2, 20)

hm_eq = [mean_c(ex4_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'
hm_dq = [mean_c(ex4_v2_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

max_plot = maximum([maximum(hm_eq), maximum(hm_dq)])
min_plot = minimum([minimum(hm_eq), minimum(hm_dq)])

ex4_p1 = heatmap(L1u, L2u, hm_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość", titlefontsize = 10, clim = (min_plot, max_plot))

savefig(ex4_p1, pwd() * "\\plots\\ex4_equal quality heatmap lambdas.svg")

ex4_p2 = heatmap(L1u, L2u, hm_dq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, różna jakość", titlefontsize = 10, clim = (min_plot, max_plot))

savefig(ex4_p2, pwd() * "\\plots\\ex4_different quality heatmap lambdas.svg")

######## Effect of number of links on surplus


ex4_v3_total_surplus = []
ex4_v3_producer_surplus = []
ex4_v3_consumer_surplus = []
ex4_v3_price = []
ex4_v3_quantity = []
ex4_v3_advertising = []
ex4_v3_producer_surplus_singleton = []
ex4_v3_quality_expectation = []
ex4_v3_durability_expectation = []
ex4_v3_margin = []

ex4_v4_total_surplus = []
ex4_v4_producer_surplus = []
ex4_v4_consumer_surplus = []
ex4_v4_price = []
ex4_v4_quantity = []
ex4_v4_advertising = []
ex4_v4_advertising_highest = []
ex4_v4_producer_surplus_singleton = []
ex4_v4_quality_expectation = []
ex4_v4_durability_expectation = []
ex4_v4_margin = []

ex4_v16_nl = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = fill(rand(m_d), 2)

    d_init = [2,2]

    nl = sample(1:2000)

    push!(ex4_v16_nl, nl)

    # 1.0

    ex4_v3_sim = TO_GO(2, 500, 250, 0.25, 0.25; q = [1.0, 1.0], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = nl, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex4_v3_total_surplus, calculate_total_surplus(ex4_v3_sim, "total"))
    push!(ex4_v3_producer_surplus, calculate_total_surplus(ex4_v3_sim, "producer"))
    push!(ex4_v3_consumer_surplus, calculate_total_surplus(ex4_v3_sim, "consumer"))
    push!(ex4_v3_price, mean(calculate_price_history.(ex4_v3_sim.sellers)))
    push!(ex4_v3_quantity, mean(getfield.(ex4_v3_sim.sellers, :quantity_history)))
    push!(ex4_v3_advertising, mean(getfield.(ex4_v3_sim.sellers, :advertising_history)))
    push!(ex4_v3_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v3_sim.sellers)))
    push!(ex4_v3_quality_expectation, calculate_expectation(ex4_v3_sim, "quality", true))
    push!(ex4_v3_durability_expectation, calculate_expectation(ex4_v3_sim, "durability", true))
    push!(ex4_v3_margin, getfield.(ex4_v3_sim.sellers, :margin_history))     

    # 1.5

    ex4_v4_sim = TO_GO(2, 500, 250, 0.25, 0.25; q = [1.5, 0.5], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, num_links = nl, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex4_v4_total_surplus, calculate_total_surplus(ex4_v4_sim, "total"))
    push!(ex4_v4_producer_surplus, calculate_total_surplus(ex4_v4_sim, "producer"))
    push!(ex4_v4_consumer_surplus, calculate_total_surplus(ex4_v4_sim, "consumer"))
    push!(ex4_v4_price, mean(calculate_price_history.(ex4_v4_sim.sellers)))
    push!(ex4_v4_quantity, mean(getfield.(ex4_v4_sim.sellers, :quantity_history)))
    push!(ex4_v4_advertising, mean(getfield.(ex4_v4_sim.sellers, :advertising_history)))
    push!(ex4_v4_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v4_sim.sellers)))
    push!(ex4_v4_quality_expectation, calculate_expectation(ex4_v4_sim, "quality", true))
    push!(ex4_v4_durability_expectation, calculate_expectation(ex4_v4_sim, "durability", true))
    push!(ex4_v4_margin, getfield.(ex4_v4_sim.sellers, :margin_history))   
  

end

ex4_p3 = scatter(ex4_v16_nl, ex4_v3_total_surplus, markeralpha = 0.25, smooth = true, xlabel = "Liczba połączeń w grafie", ylabel = "Nadwyżka całkowita", label = "Identyczna jakość")
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v3_total_surplus), "blue",0.50)
scatter!(ex4_v16_nl, ex4_v4_total_surplus, markeralpha = 0.25, smooth = true, label = "Różna jakość")
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v4_total_surplus), "green")

savefig(ex4_p3, pwd() * "\\plots\\ex4_p3_graph connections on total surplus.svg")

GLM.lm(@formula(eq_s ~ nl), DataFrame(nl = Int64.(ex4_v16_nl), eq_s = Float64.(ex4_v3_total_surplus)))
GLM.lm(@formula(dq_s ~ nl), DataFrame(nl = Int64.(ex4_v16_nl), dq_s = Float64.(ex4_v4_total_surplus)))


ex4_p4 = scatter(ex4_v16_nl, getindex.(ex4_v4_producer_surplus_singleton,1), markeralpha = 0.25, smooth = true, xlabel = "Liczba połączeń w grafie", ylabel = "Nadwyżka producenta", label = "Producent dobra o wysokiej jakości")
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v3_total_surplus), "blue",0.50)
scatter!(ex4_v16_nl, getindex.(ex4_v4_producer_surplus_singleton,2), markeralpha = 0.25, smooth = true, label = "Producent dobra o niskiej jakości")
#add_smoothing_spline(getindex.(ex4_v16_lambda, 1), Float64.(ex4_v4_total_surplus), "green")

savefig(ex4_p4, pwd() * "\\plots\\ex4_p3_graph connections on producer surplus.svg")

GLM.lm(@formula(dq_s1 ~ nl), DataFrame(nl = Int64.(ex4_v16_nl), dq_s1 = Float64.(getindex.(ex4_v4_producer_surplus_singleton,1))))
GLM.lm(@formula(dq_s2 ~ nl), DataFrame(nl = Int64.(ex4_v16_nl), dq_s2 = Float64.(getindex.(ex4_v4_producer_surplus_singleton,2))))