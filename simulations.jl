include(pwd() * "\\methods.jl")

####################### EXPERIMENT 1, IMPACT OF QUALITY DIFFERENCES ON MARKET EFFECTIVENESS ############################################

# Cel 1: Weryfikacja hipotezy, że rynki gdzie występuje zróżnicowanie jakości produktów są bardziej efektywne niż pozostałe. 

# Efekt 1: różnicowanie jakości produktu powoduje, że rynek jest bardziej efektywny - suma nadwyżek producentów i konsumentów jest wyższa
# Efekt 2: w przypadku rynku z różnicowaniem produktu, komunikacja reklamowa pozwala istotnie zwiększyć efektywność, podczas gdy na rynku z jednakową jakością produktu efekt nie jest obserwowany

ex1_v1_total_surplus = []
ex1_v1_producer_surplus = []
ex1_v1_consumer_surplus = []
ex1_v1_price = []
ex1_v1_quantity = []
ex1_v1_advertising = []

ex1_v2_total_surplus = []
ex1_v2_producer_surplus = []
ex1_v2_consumer_surplus = []
ex1_v2_price = []
ex1_v2_quantity = []
ex1_v2_advertising = []
ex1_v2_advertising_highest = []

ex1_v12_init_margin = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_init = rand(Uniform(0,0.2))

    push!(ex1_v12_init_margin, m_init)

    ex1_v1_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.], q_init = [1., 1., 1., 1.], num_links = 1000)
    push!(ex1_v1_total_surplus, calculate_total_surplus(ex1_v1_sim, "total"))
    push!(ex1_v1_producer_surplus, calculate_total_surplus(ex1_v1_sim, "producer"))
    push!(ex1_v1_consumer_surplus, calculate_total_surplus(ex1_v1_sim, "consumer"))
    push!(ex1_v1_price, mean(calculate_price_history.(ex1_v1_sim.sellers)))
    push!(ex1_v1_quantity, mean(getfield.(ex1_v1_sim.sellers, :quantity_history)))
    push!(ex1_v1_advertising, mean(getfield.(ex1_v1_sim.sellers, :advertising_history)))

    ex1_v2_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex1_v2_total_surplus, calculate_total_surplus(ex1_v2_sim, "total"))
    push!(ex1_v2_producer_surplus, calculate_total_surplus(ex1_v2_sim, "producer"))
    push!(ex1_v2_consumer_surplus, calculate_total_surplus(ex1_v2_sim, "consumer"))
    push!(ex1_v2_price, mean(calculate_price_history.(ex1_v2_sim.sellers)))
    push!(ex1_v2_quantity, mean(getfield.(ex1_v2_sim.sellers, :quantity_history)))
    push!(ex1_v2_advertising, mean(getfield.(ex1_v2_sim.sellers, :advertising_history)))

end

plot_ecdf(ex1_v1_total_surplus, "equal Q", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v2_total_surplus, "not equal Q", "Total Surplus", "Probability", "ECDF", false)

plot_ecdf(ex1_v1_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v2_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex1_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v1_price), "equal Q", "Price", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v2_price), "not equal Q", "Price", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v1_quantity), "equal Q", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v2_quantity), "not equal Q", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v1_advertising), "equal Q", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v2_advertising), "not equal Q", "Advertising", "Probability", "ECDF", false)

p = plot(mean(ex1_v1_price), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "equal quality", legend=:topright)

plot!(mean(ex1_v2_price), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "not equal quality", legend=:topright)

EqualVarianceTTest(mean(ex1_price_eq), mean(ex2_price_dq))

p = plot(mean(ex1_v1_quantity), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex1_v2_quantity), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_quantity_eq), mean(ex2_quantity_dq))

p = plot(mean(ex1_v1_advertising), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex1_v2_advertising), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_advertising_eq), mean(ex2_advertising_dq))

####################### SPECIAL CASE: NOT EQUAL INITIAL MARGINS

ex1_v3_total_surplus = []
ex1_v3_producer_surplus = []
ex1_v3_consumer_surplus = []
ex1_v3_price = []
ex1_v3_quantity = []
ex1_v3_advertising = []
ex1_v3_producer_surplus_singleton = []

ex1_v4_total_surplus = []
ex1_v4_producer_surplus = []
ex1_v4_consumer_surplus = []
ex1_v4_price = []
ex1_v4_quantity = []
ex1_v4_advertising = []
ex1_v4_advertising_highest = []
ex1_v4_producer_surplus_singleton = []

ex1_v34_init_margin = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)

    m_init = rand(m_d,4)

    push!(ex1_v34_init_margin, m_init)

    ex1_v3_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.], q_init = [1., 1., 1., 1.], num_links = 1000)
    push!(ex1_v3_total_surplus, calculate_total_surplus(ex1_v3_sim, "total"))
    push!(ex1_v3_producer_surplus, calculate_total_surplus(ex1_v3_sim, "producer"))
    push!(ex1_v3_consumer_surplus, calculate_total_surplus(ex1_v3_sim, "consumer"))
    push!(ex1_v3_price, mean(calculate_price_history.(ex1_v3_sim.sellers)))
    push!(ex1_v3_quantity, mean(getfield.(ex1_v3_sim.sellers, :quantity_history)))
    push!(ex1_v3_advertising, mean(getfield.(ex1_v3_sim.sellers, :advertising_history)))
    push!(ex1_v3_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v3_sim.sellers, ex1_v3_sim.function_args.γ, ex1_v3_sim.function_args.num_buyers)))

    ex1_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex1_v4_total_surplus, calculate_total_surplus(ex1_v4_sim, "total"))
    push!(ex1_v4_producer_surplus, calculate_total_surplus(ex1_v4_sim, "producer"))
    push!(ex1_v4_consumer_surplus, calculate_total_surplus(ex1_v4_sim, "consumer"))
    push!(ex1_v4_price, mean(calculate_price_history.(ex1_v4_sim.sellers)))
    push!(ex1_v4_quantity, mean(getfield.(ex1_v4_sim.sellers, :quantity_history)))
    push!(ex1_v4_advertising, mean(getfield.(ex1_v4_sim.sellers, :advertising_history)))
    push!(ex1_v4_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v4_sim.sellers, ex1_v4_sim.function_args.γ, ex1_v4_sim.function_args.num_buyers)))

end

plot_ecdf(ex1_v3_total_surplus, "equal Q", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v4_total_surplus, "not equal Q", "Total Surplus", "Probability", "ECDF", false)

margin_diff_12_34 = (getindex.(ex1_v34_init_margin, 1) .+ getindex.(ex1_v34_init_margin, 2)) / 2 - (getindex.(ex1_v34_init_margin, 3) .+ getindex.(ex1_v34_init_margin, 4)) / 2
scatter(margin_diff_12_34, ex1_v4_total_surplus, smooth = true)

scatter(getindex.(ex1_v34_init_margin, 1), getindex.(ex1_v3_producer_surplus_singleton, 1), markeralpha = 0.25, smooth = true, ylim = (0, 3000))
scatter!(getindex.(ex1_v34_init_margin, 1), getindex.(ex1_v4_producer_surplus_singleton, 1), markeralpha = 0.25, smooth = true, ylim = (0,3000))

plot_ecdf(ex1_v3_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v4_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex1_v3_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v4_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v3_price), "equal Q", "Price", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v4_price), "not equal Q", "Price", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v3_quantity), "equal Q", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v4_quantity), "not equal Q", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean(ex1_v3_advertising), "equal Q", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean(ex1_v4_advertising), "not equal Q", "Advertising", "Probability", "ECDF", false)

####################### EXPERIMENT 2, IMPACT OF QUALITY DIFFERENCES ON ADVERTISING & MARKET EFFECTIVENESS ############################################

# Cel 1: Weryfikacja hipotezy, że rynki gdzie występuje zróżnicowanie jakości produktów są bardziej efektywne niż pozostałe. 

# Efekt 1: różnicowanie jakości produktu powoduje, że rynek jest bardziej efektywny - suma nadwyżek producentów i konsumentów jest wyższa
# Efekt 2: w przypadku rynku z różnicowaniem produktu, komunikacja reklamowa pozwala istotnie zwiększyć efektywność, podczas gdy na rynku z jednakową jakością produktu efekt nie jest obserwowany

# v1 - equal quality, advertising allowed

ex2_v1_total_surplus = []
ex2_v1_producer_surplus = []
ex2_v1_consumer_surplus = []
ex2_v1_price = []
ex2_v1_quantity = []
ex2_v1_advertising = []

# v2 - equal quality, advertising NOT allowed

ex2_v2_total_surplus = []
ex2_v2_producer_surplus = []
ex2_v2_consumer_surplus = []
ex2_v2_price = []
ex2_v2_quantity = []
ex2_v2_advertising = []
ex2_v2_advertising_highest = []

# v3 - NOT equal quality, advertising allowed

ex2_v3_total_surplus = []
ex2_v3_producer_surplus = []
ex2_v3_consumer_surplus = []
ex2_v3_price = []
ex2_v3_quantity = []
ex2_v3_advertising = []
ex2_v3_advertising_highest = []

# v4 - NOT equal quality, advertising NOT allowed

ex2_v4_total_surplus = []
ex2_v4_producer_surplus = []
ex2_v4_consumer_surplus = []
ex2_v4_price = []
ex2_v4_quantity = []
ex2_v4_advertising = []
ex2_v4_advertising_highest = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex2_v1_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1., 1., 1., 1.])
    push!(ex2_v1_total_surplus, calculate_total_surplus(ex2_v1_sim, "total"))
    push!(ex2_v1_producer_surplus, calculate_total_surplus(ex2_v1_sim, "producer"))
    push!(ex2_v1_consumer_surplus, calculate_total_surplus(ex2_v1_sim, "consumer"))
    push!(ex2_v1_price, mean(calculate_price_history.(ex2_v1_sim.sellers)))
    push!(ex2_v1_quantity, mean(getfield.(ex2_v1_sim.sellers, :quantity_history)))
    push!(ex2_v1_advertising, mean(getfield.(ex2_v1_sim.sellers, :advertising_history)))

    ex2_v2_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1., 1., 1., 1.], variant_advertising = false)
    push!(ex2_v2_total_surplus, calculate_total_surplus(ex2_v2_sim, "total"))
    push!(ex2_v2_producer_surplus, calculate_total_surplus(ex2_v2_sim, "producer"))
    push!(ex2_v2_consumer_surplus, calculate_total_surplus(ex2_v2_sim, "consumer"))
    push!(ex2_v2_price, mean(calculate_price_history.(ex2_v2_sim.sellers)))
    push!(ex2_v2_quantity, mean(getfield.(ex2_v2_sim.sellers, :quantity_history)))
    push!(ex2_v2_advertising, mean(getfield.(ex2_v2_sim.sellers, :advertising_history)))

    ex2_v3_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex2_v3_total_surplus, calculate_total_surplus(ex2_v3_sim, "total"))
    push!(ex2_v3_producer_surplus, calculate_total_surplus(ex2_v3_sim, "producer"))
    push!(ex2_v3_consumer_surplus, calculate_total_surplus(ex2_v3_sim, "consumer"))
    push!(ex2_v3_price, mean(calculate_price_history.(ex2_v3_sim.sellers)))
    push!(ex2_v3_quantity, mean(getfield.(ex2_v3_sim.sellers, :quantity_history)))
    push!(ex2_v3_advertising, mean(getfield.(ex2_v3_sim.sellers, :advertising_history)))

    ex2_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], d_init = [1.,1.,1.,1.,], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50], variant_advertising = false)
    push!(ex2_v4_total_surplus, calculate_total_surplus(ex2_v4_sim, "total"))
    push!(ex2_v4_producer_surplus, calculate_total_surplus(ex2_v4_sim, "producer"))
    push!(ex2_v4_consumer_surplus, calculate_total_surplus(ex2_v4_sim, "consumer"))
    push!(ex2_v4_price, mean(calculate_price_history.(ex2_v4_sim.sellers)))
    push!(ex2_v4_quantity, mean(getfield.(ex2_v4_sim.sellers, :quantity_history)))
    push!(ex2_v4_advertising, mean(getfield.(ex2_v4_sim.sellers, :advertising_history)))

end

plot_ecdf(ex2_v1_total_surplus, "equal Q, advertising", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v2_total_surplus, "equal Q, no advertising", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v3_total_surplus, "not equal Q, advertising", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v4_total_surplus, "not equal Q, no advertising", "Total Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v1_consumer_surplus, "equal Q, advertising", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v2_consumer_surplus, "equal Q, no advertising", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v3_consumer_surplus, "not equal Q, advertising", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v4_consumer_surplus, "not equal Q, no advertising", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v1_producer_surplus, "equal Q, advertising", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v2_producer_surplus, "equal Q, no advertising", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v3_producer_surplus, "not equal Q, advertising", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v4_producer_surplus, "not equal Q, no advertising", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean(ex2_v1_price), "equal Q, advertising", "Price", "Probability", "ECDF", true)
plot_ecdf(mean(ex2_v2_price), "equal Q, no advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v3_price), "not equal Q, advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v4_price), "not equal Q, no advertising", "Price", "Probability", "ECDF", false)

plot_ecdf(mean(ex2_v1_quantity), "equal Q, advertising", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean(ex2_v2_quantity), "equal Q, no advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v3_quantity), "not equal Q, advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v4_quantity), "not equal Q, no advertising", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean(ex2_v1_advertising), "equal Q, advertising", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean(ex2_v2_advertising), "equal Q, no advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v3_advertising), "not equal Q, advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean(ex2_v4_advertising), "not equal Q, no advertising", "Advertising", "Probability", "ECDF", false)


####################### EXPERIMENT 3, IMPACT OF INITIAL MARGIN ON END RESULTS ###########################################

# Cel 1: Weryfikacja hipotezy, że rynki gdzie występuje zróżnicowanie jakości produktów są bardziej efektywne niż pozostałe. 

# Efekt 1: różnicowanie jakości produktu powoduje, że rynek jest bardziej efektywny - suma nadwyżek producentów i konsumentów jest wyższa
# Efekt 2: w przypadku rynku z różnicowaniem produktu, komunikacja reklamowa pozwala istotnie zwiększyć efektywność, podczas gdy na rynku z jednakową jakością produktu efekt nie jest obserwowany

margin_test = collect(0:0.05:0.50)
max_iter = 1000
initial_margin_sensitivity_eq = []

for iter in 1:max_iter
    println(iter)
    m1 = sample(margin_test)
    m2 = sample(margin_test)
    sim_with_obs_13 = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0], m = [m1, m2], c = [0.6, 0.6], ϵ = [0.33, 0.33], a = [0.0, 0.0], r = [0.4, 0.4], q_init = [1.0, 1.0], d = [1,1], d_init = [1.,1.], num_links = 1000)
    push!(initial_margin_sensitivity_eq, (m = getfield(getfield(sim_with_obs_13, :function_args), :m), profit = sum.(calculate_profit_history.(sim_with_obs_13.sellers, sim_with_obs_13.function_args.γ, sim_with_obs_13.function_args.num_buyers)), price = mean(calculate_price_history.(sim_with_obs_13.sellers)), quantity = mean(getfield.(sim_with_obs_13.sellers, :quantity_history))))
end

m = getfield.(initial_margin_sensitivity_eq, :m)
profit = getfield.(initial_margin_sensitivity_eq, :profit)
price = getfield.(initial_margin_sensitivity_eq, :price)
quantity = getfield.(initial_margin_sensitivity_eq, :quantity)

profit_margin = [mean(sum.(profit)[(getindex.(m,1) .== m1) .& (getindex.(m,2) .== m2)]) for m1 in margin_test, m2 in margin_test]
heatmap(profit_margin)

quantity_margin = [mean(sum.(quantity)[(getindex.(m,1) .== m1) .& (getindex.(m,2) .== m2)]) for m1 in margin_test, m2 in margin_test]
heatmap(quantity_margin)

price_margin = [mean(mean.(price)[(getindex.(m,1) .== m1) .& (getindex.(m,2) .== m2)]) for m1 in margin_test, m2 in margin_test]
heatmap(price_margin)
