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

ex4_v6_total_surplus = []
ex4_v6_producer_surplus = []
ex4_v6_consumer_surplus = []
ex4_v6_price = []
ex4_v6_quantity = []
ex4_v6_advertising = []
ex4_v6_advertising_highest = []
ex4_v6_producer_surplus_singleton = []
ex4_v6_quality_expectation = []
ex4_v6_durability_expectation = []
ex4_v6_margin = []

ex4_v16_nn = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = fill(rand(m_d), 2)

    d_init = [2,2]

    n_neigh = sample(0:10:2000)

    push!(ex4_v16_nn, n_neigh)

    # 1.0

    ex4_v1_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1., 1.], num_links = n_neigh)
    push!(ex4_v1_total_surplus, calculate_total_surplus(ex4_v1_sim, "total"))
    push!(ex4_v1_producer_surplus, calculate_total_surplus(ex4_v1_sim, "producer"))
    push!(ex4_v1_consumer_surplus, calculate_total_surplus(ex4_v1_sim, "consumer"))
    push!(ex4_v1_price, mean(calculate_price_history.(ex4_v1_sim.sellers)))
    push!(ex4_v1_quantity, mean(getfield.(ex4_v1_sim.sellers, :quantity_history)))
    push!(ex4_v1_advertising, mean(getfield.(ex4_v1_sim.sellers, :advertising_history)))
    push!(ex4_v1_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v1_sim.sellers, ex4_v1_sim.function_args.γ, ex4_v1_sim.function_args.num_buyers)))
    push!(ex4_v1_quality_expectation, calculate_expectation(ex4_v1_sim, "quality", true))
    push!(ex4_v1_durability_expectation, calculate_expectation(ex4_v1_sim, "durability", true))
    push!(ex4_v1_margin, getfield.(ex4_v1_sim.sellers, :margin_history))     

    # 1.5

    ex4_v6_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 0.5], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.5, 0.5], num_links = n_neigh)
    push!(ex4_v6_total_surplus, calculate_total_surplus(ex4_v6_sim, "total"))
    push!(ex4_v6_producer_surplus, calculate_total_surplus(ex4_v6_sim, "producer"))
    push!(ex4_v6_consumer_surplus, calculate_total_surplus(ex4_v6_sim, "consumer"))
    push!(ex4_v6_price, mean(calculate_price_history.(ex4_v6_sim.sellers)))
    push!(ex4_v6_quantity, mean(getfield.(ex4_v6_sim.sellers, :quantity_history)))
    push!(ex4_v6_advertising, mean(getfield.(ex4_v6_sim.sellers, :advertising_history)))
    push!(ex4_v6_producer_surplus_singleton, sum.(calculate_profit_history.(ex4_v6_sim.sellers, ex4_v6_sim.function_args.γ, ex4_v6_sim.function_args.num_buyers)))
    push!(ex4_v6_quality_expectation, calculate_expectation(ex4_v6_sim, "quality", true))
    push!(ex4_v6_durability_expectation, calculate_expectation(ex4_v6_sim, "durability", true))
    push!(ex4_v6_margin, getfield.(ex4_v6_sim.sellers, :margin_history))   
  

end


scatter(ex4_v16_nn, ex4_v1_total_surplus)
scatter!(ex4_v16_nn, ex4_v6_total_surplus)

scatter(ex4_v16_nn, ex4_v1_consumer_surplus)
scatter!(ex4_v16_nn, ex4_v6_consumer_surplus)

scatter(ex4_v16_nn, ex4_v1_producer_surplus)
scatter!(ex4_v16_nn, ex4_v6_producer_surplus)

scatter(ex4_v16_nn, getindex.(ex4_v1_producer_surplus_singleton,1))
scatter!(ex4_v16_nn, getindex.(ex4_v1_producer_surplus_singleton,2))

scatter(ex4_v16_nn, getindex.(ex4_v6_producer_surplus_singleton,1))
scatter!(ex4_v16_nn, getindex.(ex4_v6_producer_surplus_singleton,2))