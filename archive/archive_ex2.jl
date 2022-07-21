
#################################################################

ex1_v5_total_surplus = []
ex1_v5_producer_surplus = []
ex1_v5_consumer_surplus = []
ex1_v5_price = []
ex1_v5_quantity = []
ex1_v5_advertising = []
ex1_v5_producer_surplus_singleton = []

ex1_v6_total_surplus = []
ex1_v6_producer_surplus = []
ex1_v6_consumer_surplus = []
ex1_v6_price = []
ex1_v6_quantity = []
ex1_v6_advertising = []
ex1_v6_advertising_highest = []
ex1_v6_producer_surplus_singleton = []

ex1_v56_r = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = rand(m_d,4)

    r_d = Uniform(0.2,0.5)
    r = rand(r_d, 4)

    push!(ex1_v56_r, r)

    ex1_v5_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = r, d = [1,1,1,1], q_init = [1., 1., 1., 1.], num_links = 1000)
    push!(ex1_v5_total_surplus, calculate_total_surplus(ex1_v5_sim, "total"))
    push!(ex1_v5_producer_surplus, calculate_total_surplus(ex1_v5_sim, "producer"))
    push!(ex1_v5_consumer_surplus, calculate_total_surplus(ex1_v5_sim, "consumer"))
    push!(ex1_v5_price, mean(calculate_price_history.(ex1_v5_sim.sellers)))
    push!(ex1_v5_quantity, mean(getfield.(ex1_v5_sim.sellers, :quantity_history)))
    push!(ex1_v5_advertising, mean(getfield.(ex1_v5_sim.sellers, :advertising_history)))
    push!(ex1_v5_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v5_sim.sellers, ex1_v5_sim.function_args.γ, ex1_v5_sim.function_args.num_buyers)))

    ex1_v6_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = r, d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex1_v6_total_surplus, calculate_total_surplus(ex1_v6_sim, "total"))
    push!(ex1_v6_producer_surplus, calculate_total_surplus(ex1_v6_sim, "producer"))
    push!(ex1_v6_consumer_surplus, calculate_total_surplus(ex1_v6_sim, "consumer"))
    push!(ex1_v6_price, mean(calculate_price_history.(ex1_v6_sim.sellers)))
    push!(ex1_v6_quantity, mean(getfield.(ex1_v6_sim.sellers, :quantity_history)))
    push!(ex1_v6_advertising, mean(getfield.(ex1_v6_sim.sellers, :advertising_history)))
    push!(ex1_v6_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v6_sim.sellers, ex1_v6_sim.function_args.γ, ex1_v6_sim.function_args.num_buyers)))

end

plot_ecdf(ex1_v5_total_surplus, "equal Q", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v6_total_surplus, "not equal Q", "Total Surplus", "Probability", "ECDF", false)

scatter(getindex.(ex1_v56_r,1), getindex.(ex1_v5_producer_surplus_singleton,1),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,2), getindex.(ex1_v5_producer_surplus_singleton,2),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,3), getindex.(ex1_v5_producer_surplus_singleton,3),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,4), getindex.(ex1_v5_producer_surplus_singleton,4),smooth=true,markeralpha=0.25)

scatter(getindex.(ex1_v56_r,1), getindex.(ex1_v6_producer_surplus_singleton,1),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,2), getindex.(ex1_v6_producer_surplus_singleton,2),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,3), getindex.(ex1_v6_producer_surplus_singleton,3),smooth=true,markeralpha=0.25)
scatter!(getindex.(ex1_v56_r,4), getindex.(ex1_v6_producer_surplus_singleton,4),smooth=true,markeralpha=0.25)

plot_ecdf(ex1_v5_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v6_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex1_v5_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v6_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v5_price), "equal Q", "Price", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v6_price), "not equal Q", "Price", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v5_quantity), "equal Q", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v6_quantity), "not equal Q", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v5_advertising), "equal Q", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v6_advertising), "not equal Q", "Advertising", "Probability", "ECDF", false)