include(pwd() * "\\methods\\methods.jl")

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

    ex1_v1_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], q_init = [1., 1., 1., 1.], num_links = 1000)
    push!(ex1_v1_total_surplus, calculate_total_surplus(ex1_v1_sim, "total"))
    push!(ex1_v1_producer_surplus, calculate_total_surplus(ex1_v1_sim, "producer"))
    push!(ex1_v1_consumer_surplus, calculate_total_surplus(ex1_v1_sim, "consumer"))
    push!(ex1_v1_price, mean(calculate_price_history.(ex1_v1_sim.sellers)))
    push!(ex1_v1_quantity, mean(getfield.(ex1_v1_sim.sellers, :quantity_history)))
    push!(ex1_v1_advertising, mean(getfield.(ex1_v1_sim.sellers, :advertising_history)))

    ex1_v2_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
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

plot_ecdf(mean.(ex1_v1_price), "equal Q", "Price", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v2_price), "not equal Q", "Price", "Probability", "ECDF", false)

plot_ecdf(sum.(ex1_v1_quantity), "equal Q", "Quantity", "Probability", "ECDF", true)
plot_ecdf(sum.(ex1_v2_quantity), "not equal Q", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean.(ex1_v1_advertising), "equal Q", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v2_advertising), "not equal Q", "Advertising", "Probability", "ECDF", false)

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

####################### EXPERIMENT 2, IMPACT OF QUALITY DIFFERENCES ON ADVERTISING & MARKET EFFECTIVENESS ############################################

# Cel 1: Sprawdzenie, w jaki sposób zmienia się efektywność komunikacji reklamowej, gdy na rynku występuje (oraz nie) zróżnicowanie jakości produktów

# Wynik:

# v3 - equal quality, advertising allowed

ex2_v3_total_surplus = []
ex2_v3_producer_surplus = []
ex2_v3_consumer_surplus = []
ex2_v3_price = []
ex2_v3_quantity = []
ex2_v3_advertising = []

# v4 - equal quality, advertising NOT allowed

ex2_v4_total_surplus = []
ex2_v4_producer_surplus = []
ex2_v4_consumer_surplus = []
ex2_v4_price = []
ex2_v4_quantity = []
ex2_v4_advertising = []
ex2_v4_advertising_highest = []

# v5 - NOT equal quality, advertising allowed

ex2_v5_total_surplus = []
ex2_v5_producer_surplus = []
ex2_v5_consumer_surplus = []
ex2_v5_price = []
ex2_v5_quantity = []
ex2_v5_advertising = []
ex2_v5_advertising_highest = []

# v6 - NOT equal quality, advertising NOT allowed

ex2_v6_total_surplus = []
ex2_v6_producer_surplus = []
ex2_v6_consumer_surplus = []
ex2_v6_price = []
ex2_v6_quantity = []
ex2_v6_advertising = []
ex2_v6_advertising_highest = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_init = rand(Uniform(0,0.2), 4)

    ex2_v3_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1., 1., 1., 1.])
    push!(ex2_v3_total_surplus, calculate_total_surplus(ex2_v3_sim, "total"))
    push!(ex2_v3_producer_surplus, calculate_total_surplus(ex2_v3_sim, "producer"))
    push!(ex2_v3_consumer_surplus, calculate_total_surplus(ex2_v3_sim, "consumer"))
    push!(ex2_v3_price, mean(calculate_price_history.(ex2_v3_sim.sellers)))
    push!(ex2_v3_quantity, mean(getfield.(ex2_v3_sim.sellers, :quantity_history)))
    push!(ex2_v3_advertising, mean(getfield.(ex2_v3_sim.sellers, :advertising_history)))

    ex2_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1., 1., 1., 1.], variant_advertising = false)
    push!(ex2_v4_total_surplus, calculate_total_surplus(ex2_v4_sim, "total"))
    push!(ex2_v4_producer_surplus, calculate_total_surplus(ex2_v4_sim, "producer"))
    push!(ex2_v4_consumer_surplus, calculate_total_surplus(ex2_v4_sim, "consumer"))
    push!(ex2_v4_price, mean(calculate_price_history.(ex2_v4_sim.sellers)))
    push!(ex2_v4_quantity, mean(getfield.(ex2_v4_sim.sellers, :quantity_history)))
    push!(ex2_v4_advertising, mean(getfield.(ex2_v4_sim.sellers, :advertising_history)))

    ex2_v5_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex2_v5_total_surplus, calculate_total_surplus(ex2_v5_sim, "total"))
    push!(ex2_v5_producer_surplus, calculate_total_surplus(ex2_v5_sim, "producer"))
    push!(ex2_v5_consumer_surplus, calculate_total_surplus(ex2_v5_sim, "consumer"))
    push!(ex2_v5_price, mean(calculate_price_history.(ex2_v5_sim.sellers)))
    push!(ex2_v5_quantity, mean(getfield.(ex2_v5_sim.sellers, :quantity_history)))
    push!(ex2_v5_advertising, mean(getfield.(ex2_v5_sim.sellers, :advertising_history)))

    ex2_v6_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50], variant_advertising = false)
    push!(ex2_v6_total_surplus, calculate_total_surplus(ex2_v6_sim, "total"))
    push!(ex2_v6_producer_surplus, calculate_total_surplus(ex2_v6_sim, "producer"))
    push!(ex2_v6_consumer_surplus, calculate_total_surplus(ex2_v6_sim, "consumer"))
    push!(ex2_v6_price, mean(calculate_price_history.(ex2_v6_sim.sellers)))
    push!(ex2_v6_quantity, mean(getfield.(ex2_v6_sim.sellers, :quantity_history)))
    push!(ex2_v6_advertising, mean(getfield.(ex2_v6_sim.sellers, :advertising_history)))

end

plot_ecdf(ex2_v3_total_surplus, "equal Q, advertising", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v4_total_surplus, "equal Q, no advertising", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v5_total_surplus, "not equal Q, advertising", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v6_total_surplus, "not equal Q, no advertising", "Total Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v3_consumer_surplus, "equal Q, advertising", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v4_consumer_surplus, "equal Q, no advertising", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v5_consumer_surplus, "not equal Q, advertising", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v6_consumer_surplus, "not equal Q, no advertising", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v3_producer_surplus, "equal Q, advertising", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v4_producer_surplus, "equal Q, no advertising", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v5_producer_surplus, "not equal Q, advertising", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v6_producer_surplus, "not equal Q, no advertising", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v3_price), "equal Q, advertising", "Price", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v4_price), "equal Q, no advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v5_price), "not equal Q, advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v6_price), "not equal Q, no advertising", "Price", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v3_quantity), "equal Q, advertising", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v4_quantity), "equal Q, no advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v5_quantity), "not equal Q, advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v6_quantity), "not equal Q, no advertising", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v3_advertising), "equal Q, advertising", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v4_advertising), "equal Q, no advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v5_advertising), "not equal Q, advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v6_advertising), "not equal Q, no advertising", "Advertising", "Probability", "ECDF", false)


####################### EXPERIMENT 3, IMPACT OF QUALITY AND DURABILITY DIFFERENCES & MARKET EFFECTIVENESS ############################################

# Cel 1: Sprawdzenie, w jaki sposób zmienia się efektywność komunikacji reklamowej, gdy na rynku występuje (oraz nie) zróżnicowanie jakości produktów

# Wynik:

# v7 - equal quality, equal durability

ex2_v7_total_surplus = []
ex2_v7_producer_surplus = []
ex2_v7_consumer_surplus = []
ex2_v7_price = []
ex2_v7_quantity = []
ex2_v7_advertising = []

# v8 - equal quality, NOT equal durability

ex2_v8_total_surplus = []
ex2_v8_producer_surplus = []
ex2_v8_consumer_surplus = []
ex2_v8_price = []
ex2_v8_quantity = []
ex2_v8_advertising = []
ex2_v8_advertising_highest = []

# v9 - NOT equal quality, equal durability

ex2_v9_total_surplus = []
ex2_v9_producer_surplus = []
ex2_v9_consumer_surplus = []
ex2_v9_price = []
ex2_v9_quantity = []
ex2_v9_advertising = []
ex2_v9_advertising_highest = []

# v10 - NOT equal quality, NOT equal durability

ex2_v10_total_surplus = []
ex2_v10_producer_surplus = []
ex2_v10_consumer_surplus = []
ex2_v10_price = []
ex2_v10_quantity = []
ex2_v10_advertising = []
ex2_v10_advertising_highest = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_init = rand(Uniform(0,0.2))
    m_init = fill(m_init, 4)

    d_init = sample(1:4,4) # 🍎 tu poprawić

    ex2_v7_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, q_init = [1., 1., 1., 1.])
    push!(ex2_v7_total_surplus, calculate_total_surplus(ex2_v7_sim, "total"))
    push!(ex2_v7_producer_surplus, calculate_total_surplus(ex2_v7_sim, "producer"))
    push!(ex2_v7_consumer_surplus, calculate_total_surplus(ex2_v7_sim, "consumer"))
    push!(ex2_v7_price, mean(calculate_price_history.(ex2_v7_sim.sellers)))
    push!(ex2_v7_quantity, mean(getfield.(ex2_v7_sim.sellers, :quantity_history)))
    push!(ex2_v7_advertising, mean(getfield.(ex2_v7_sim.sellers, :advertising_history)))

    ex2_v8_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, q_init = [1., 1., 1., 1.])
    push!(ex2_v8_total_surplus, calculate_total_surplus(ex2_v8_sim, "total"))
    push!(ex2_v8_producer_surplus, calculate_total_surplus(ex2_v8_sim, "producer"))
    push!(ex2_v8_consumer_surplus, calculate_total_surplus(ex2_v8_sim, "consumer"))
    push!(ex2_v8_price, mean(calculate_price_history.(ex2_v8_sim.sellers)))
    push!(ex2_v8_quantity, mean(getfield.(ex2_v8_sim.sellers, :quantity_history)))
    push!(ex2_v8_advertising, mean(getfield.(ex2_v8_sim.sellers, :advertising_history)))

    ex2_v9_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex2_v9_total_surplus, calculate_total_surplus(ex2_v9_sim, "total"))
    push!(ex2_v9_producer_surplus, calculate_total_surplus(ex2_v9_sim, "producer"))
    push!(ex2_v9_consumer_surplus, calculate_total_surplus(ex2_v9_sim, "consumer"))
    push!(ex2_v9_price, mean(calculate_price_history.(ex2_v9_sim.sellers)))
    push!(ex2_v9_quantity, mean(getfield.(ex2_v9_sim.sellers, :quantity_history)))
    push!(ex2_v9_advertising, mean(getfield.(ex2_v9_sim.sellers, :advertising_history)))

    ex2_v10_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex2_v10_total_surplus, calculate_total_surplus(ex2_v10_sim, "total"))
    push!(ex2_v10_producer_surplus, calculate_total_surplus(ex2_v10_sim, "producer"))
    push!(ex2_v10_consumer_surplus, calculate_total_surplus(ex2_v10_sim, "consumer"))
    push!(ex2_v10_price, mean(calculate_price_history.(ex2_v10_sim.sellers)))
    push!(ex2_v10_quantity, mean(getfield.(ex2_v10_sim.sellers, :quantity_history)))
    push!(ex2_v10_advertising, mean(getfield.(ex2_v10_sim.sellers, :advertising_history)))

end

plot_ecdf(ex2_v7_total_surplus, "equal Q, equal D", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v8_total_surplus, "equal Q, ~equal D", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v9_total_surplus, "~equal Q, equal D", "Total Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v10_total_surplus, "~equal Q, ~equal D", "Total Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v7_consumer_surplus, "equal Q, equal D", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v8_consumer_surplus, "equal Q, ~equal D", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v9_consumer_surplus, "~equal Q, equal D", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v10_consumer_surplus, "~equal Q, ~equal D", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex2_v7_producer_surplus, "equal Q, equal D", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v8_producer_surplus, "equal Q, ~equal D", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v9_producer_surplus, "~equal Q, equal D", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex2_v10_producer_surplus, "~equal Q, ~equal D", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v7_price), "equal Q, equal D", "Price", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v8_price), "equal Q, ~equal D", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v9_price), "~equal Q, equal D", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v10_price), "~equal Q, ~equal D", "Price", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v7_quantity), "equal Q, equal D", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v8_quantity), "equal Q, ~equal D", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v9_quantity), "~equal Q, equal D", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v10_quantity), "~equal Q, ~equal D", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v7_advertising), "equal Q, equal D", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v8_advertising), "equal Q, ~equal D", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v9_advertising), "~equal Q, equal D", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v10_advertising), "~equal Q, ~equal D", "Advertising", "Probability", "ECDF", false)

