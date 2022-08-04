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

    ex1_v1_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1.0, 1.0, 1.0, 1.0], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [2,2,2,2], num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex1_v1_total_surplus, calculate_total_surplus(ex1_v1_sim, "total"))
    push!(ex1_v1_producer_surplus, calculate_total_surplus(ex1_v1_sim, "producer"))
    push!(ex1_v1_consumer_surplus, calculate_total_surplus(ex1_v1_sim, "consumer"))
    push!(ex1_v1_price, calculate_price_history.(ex1_v1_sim.sellers))
    push!(ex1_v1_quantity, getfield.(ex1_v1_sim.sellers, :quantity_history))
    push!(ex1_v1_advertising, mean(getfield.(ex1_v1_sim.sellers, :advertising_history)))

    ex1_v2_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1.3, 1.15, 0.85, 0.7], m = fill(m_init, 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [2,2,2,2], num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex1_v2_total_surplus, calculate_total_surplus(ex1_v2_sim, "total"))
    push!(ex1_v2_producer_surplus, calculate_total_surplus(ex1_v2_sim, "producer"))
    push!(ex1_v2_consumer_surplus, calculate_total_surplus(ex1_v2_sim, "consumer"))
    push!(ex1_v2_price, calculate_price_history.(ex1_v2_sim.sellers))
    push!(ex1_v2_quantity, getfield.(ex1_v2_sim.sellers, :quantity_history))
    push!(ex1_v2_advertising, mean(getfield.(ex1_v2_sim.sellers, :advertising_history)))

end

ex1_p1 = plot_ecdf(ex1_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex1_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

savefig(ex1_p1, pwd() * "\\plots\\ex1_surplus_total.svg")

ex1_p2 = plot_ecdf(ex1_v1_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v2_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

savefig(ex1_p2, pwd() * "\\plots\\ex1_surplus_consumer.pdf")

ex1_p3 = plot_ecdf(ex1_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

savefig(ex1_p3, pwd() * "\\plots\\ex1_surplus_producer.pdf")

ex1_p4 = plot_ecdf(mean.(mean.(ex1_v1_price)), "Identyczna jakość", "Średnia cena", "F(x)", "Dystrybuanta empiryczna - średnia cena", true)
plot_ecdf(mean.(mean.(ex1_v2_price)), "Różna jakość", "Średnia cena", "F(x)", "Dystrybuanta empiryczna - średnia cena", false)

savefig(ex1_p4, pwd() * "\\plots\\ex1_price.svg")

ex1_p5 = plot_ecdf(sum.(sum.(ex1_v1_quantity)), "Identyczna jakość", "Wolumen sprzedaży", "F(x)", "Dystrybuanta empiryczna - wolumen sprzedaży", true)
plot_ecdf(sum.(sum.(ex1_v2_quantity)), "Różna jakość", "Wolumen sprzedaży", "F(x)", "Dystrybuanta empiryczna - wolumen sprzedaży", false)

savefig(ex1_p5, pwd() * "\\plots\\ex1_quantity.svg")

plot_ecdf(mean.(ex1_v1_advertising), "equal Q", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v2_advertising), "not equal Q", "Advertising", "Probability", "ECDF", false)

p = plot(mean(ex1_v1_price), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "equal quality", legend=:topright)

plot!(mean(ex1_v2_price), 
    xlabel = "T", ylabel = "Price", 
    title = "Average price", label = "not equal quality", legend=:topright)

EqualVarianceTTest(mean(ex1_price_eq), mean(ex1_price_dq))

p = plot(mean(ex1_v1_quantity), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex1_v2_quantity), 
    xlabel = "T", ylabel = "Quantity", 
    title = "Average quantity", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_quantity_eq), mean(ex1_quantity_dq))

p = plot(mean(ex1_v1_advertising), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "equal quality", legend=:bottomright)
    
plot!(mean(ex1_v2_advertising), 
    xlabel = "T", ylabel = "Advertising", 
    title = "Average advertising", label = "non-equal quality")

EqualVarianceTTest(mean(ex1_advertising_eq), mean(ex1_advertising_dq))

boxplot([trim_outliers(calculate_average_elasticity.(getindex.(ex1_v1_quantity, 1), getindex.(ex1_v1_price, 1))), trim_outliers(calculate_average_elasticity.(getindex.(ex1_v1_quantity, 2), getindex.(ex1_v1_price, 2)))], ylim = (-150,150))

boxplot([trim_outliers(calculate_average_elasticity.(getindex.(ex1_v2_quantity, 1), getindex.(ex1_v2_price, 1))), trim_outliers(calculate_average_elasticity.(getindex.(ex1_v2_quantity, 2), getindex.(ex1_v2_price, 2)))], ylim = (-150,150))

mean(trim_outliers(calculate_average_elasticity.(getindex.(ex1_v1_quantity, 1), getindex.(ex1_v1_price, 1)))) - mean(trim_outliers(calculate_average_elasticity.(getindex.(ex1_v2_quantity, 1), getindex.(ex1_v2_price, 1))))

UnequalVarianceTTest(trim_outliers(calculate_average_elasticity.(getindex.(ex1_v1_quantity, 2), getindex.(ex1_v1_price, 2))), trim_outliers(calculate_average_elasticity.(getindex.(ex1_v2_quantity, 2), getindex.(ex1_v2_price, 2))))

UnequalVarianceTTest(trim_outliers(calculate_average_elasticity.(getindex.(ex1_v1_quantity, 1), getindex.(ex1_v1_price, 1))), trim_outliers(calculate_average_elasticity.(getindex.(ex1_v2_quantity, 1), getindex.(ex1_v2_price, 1))))

####################### EXPERIMENT 3, IMPACT OF QUALITY AND DURABILITY DIFFERENCES & MARKET EFFECTIVENESS ############################################

# Cel 1: Sprawdzenie, w jaki sposób zmienia się efektywność komunikacji reklamowej, gdy na rynku występuje (oraz nie) zróżnicowanie jakości produktów

# Wynik:

# v7 - equal quality, equal durability

ex1_v7_total_surplus = []
ex1_v7_producer_surplus = []
ex1_v7_consumer_surplus = []
ex1_v7_price = []
ex1_v7_quantity = []
ex1_v7_advertising = []
ex1_v7_producer_surplus_singleton = []

# v8 - equal quality, NOT equal durability

ex1_v8_total_surplus = []
ex1_v8_producer_surplus = []
ex1_v8_consumer_surplus = []
ex1_v8_price = []
ex1_v8_quantity = []
ex1_v8_advertising = []
ex1_v8_advertising_highest = []
ex1_v8_producer_surplus_singleton = []

# v9 - NOT equal quality, equal durability

ex1_v9_total_surplus = []
ex1_v9_producer_surplus = []
ex1_v9_consumer_surplus = []
ex1_v9_price = []
ex1_v9_quantity = []
ex1_v9_advertising = []
ex1_v9_advertising_highest = []
ex1_v9_producer_surplus_singleton = []

# v10 - NOT equal quality, NOT equal durability

ex1_v10_total_surplus = []
ex1_v10_producer_surplus = []
ex1_v10_consumer_surplus = []
ex1_v10_price = []
ex1_v10_quantity = []
ex1_v10_advertising = []
ex1_v10_advertising_highest = []
ex1_v10_producer_surplus_singleton = []

ex1_v710_m = []
ex1_v710_d = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_init = rand(Uniform(0,0.2))
    m_init = fill(m_init, 4)

    push!(ex1_v710_m, m_init)

    d_init = shuffle([3,3,1,1])

    push!(ex1_v710_d, d_init)

    ex1_v7_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [2,2,2,2], num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")

    push!(ex1_v7_total_surplus, calculate_total_surplus(ex1_v7_sim, "total"))
    push!(ex1_v7_producer_surplus, calculate_total_surplus(ex1_v7_sim, "producer"))
    push!(ex1_v7_consumer_surplus, calculate_total_surplus(ex1_v7_sim, "consumer"))
    push!(ex1_v7_price, mean(calculate_price_history.(ex1_v7_sim.sellers)))
    push!(ex1_v7_quantity, mean(getfield.(ex1_v7_sim.sellers, :quantity_history)))
    push!(ex1_v7_advertising, mean(getfield.(ex1_v7_sim.sellers, :advertising_history)))
    push!(ex1_v7_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v7_sim.sellers)))

    ex1_v8_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex1_v8_total_surplus, calculate_total_surplus(ex1_v8_sim, "total"))
    push!(ex1_v8_producer_surplus, calculate_total_surplus(ex1_v8_sim, "producer"))
    push!(ex1_v8_consumer_surplus, calculate_total_surplus(ex1_v8_sim, "consumer"))
    push!(ex1_v8_price, mean(calculate_price_history.(ex1_v8_sim.sellers)))
    push!(ex1_v8_quantity, mean(getfield.(ex1_v8_sim.sellers, :quantity_history)))
    push!(ex1_v8_advertising, mean(getfield.(ex1_v8_sim.sellers, :advertising_history)))
    push!(ex1_v8_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v8_sim.sellers)))

    ex1_v9_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1.3, 1.15, 0.85, 0.7], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [2,2,2,2], num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex1_v9_total_surplus, calculate_total_surplus(ex1_v9_sim, "total"))
    push!(ex1_v9_producer_surplus, calculate_total_surplus(ex1_v9_sim, "producer"))
    push!(ex1_v9_consumer_surplus, calculate_total_surplus(ex1_v9_sim, "consumer"))
    push!(ex1_v9_price, mean(calculate_price_history.(ex1_v9_sim.sellers)))
    push!(ex1_v9_quantity, mean(getfield.(ex1_v9_sim.sellers, :quantity_history)))
    push!(ex1_v9_advertising, mean(getfield.(ex1_v9_sim.sellers, :advertising_history)))
    push!(ex1_v9_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v9_sim.sellers)))

    ex1_v10_sim = TO_GO(4, 500, 250, 0.25, 0.25; q = [1.3, 1.15, 0.85, 0.7], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = d_init, num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")
    push!(ex1_v10_total_surplus, calculate_total_surplus(ex1_v10_sim, "total"))
    push!(ex1_v10_producer_surplus, calculate_total_surplus(ex1_v10_sim, "producer"))
    push!(ex1_v10_consumer_surplus, calculate_total_surplus(ex1_v10_sim, "consumer"))
    push!(ex1_v10_price, mean(calculate_price_history.(ex1_v10_sim.sellers)))
    push!(ex1_v10_quantity, mean(getfield.(ex1_v10_sim.sellers, :quantity_history)))
    push!(ex1_v10_advertising, mean(getfield.(ex1_v10_sim.sellers, :advertising_history)))
    push!(ex1_v10_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v10_sim.sellers)))

end

# Rynek, na którym występuje zróżnicowanie jakości i trwałości produktów dominuje, pod względem nadwyżek, pozostałe rynki.

ex1_p5 = plot_ecdf(ex1_v7_total_surplus, "Identyczna jakość, identyczna trwałość", "Nadwyżka całkowita", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex1_v8_total_surplus, "Identyczna jakość, różna trwałość", "Nadwyżka całkowita", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(ex1_v9_total_surplus, "Różna jakość, identyczna trwałość", "Nadwyżka całkowita", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(ex1_v10_total_surplus, "Różna jakość, różna trwałość", "Nadwyżka całkowita", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

savefig(ex1_p5, pwd() * "\\plots\\ex1_durability quality surplus.svg")

UnequalVarianceTTest(Float64.(ex1_v7_total_surplus), Float64.(ex1_v9_total_surplus))
UnequalVarianceTTest(Float64.(ex1_v8_total_surplus), Float64.(ex1_v10_total_surplus))

UnequalVarianceTTest(Float64.(sort(ex1_v8_total_surplus) .- sort(ex1_v10_total_surplus)), Float64.(sort(ex1_v7_total_surplus) .- sort(ex1_v9_total_surplus)))



plot_ecdf(ex1_v7_consumer_surplus, "equal Q, equal D", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v8_consumer_surplus, "equal Q, ~equal D", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex1_v9_consumer_surplus, "~equal Q, equal D", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex1_v10_consumer_surplus, "~equal Q, ~equal D", "Consumer Surplus", "Probability", "ECDF", false)

plot_ecdf(ex1_v7_producer_surplus, "equal Q, equal D", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v8_producer_surplus, "equal Q, ~equal D", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex1_v9_producer_surplus, "~equal Q, equal D", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex1_v10_producer_surplus, "~equal Q, ~equal D", "Producer Surplus", "Probability", "ECDF", false)

####

ex1_v710_d_u = sort(unique(ex1_v710_d))

# Całkowita nadwyżka rynkowa jest maksymalna, gdy produkty o najwyższej jakości mają najniższą trwałość.

total_surplus_durability_eq = [ex1_v8_total_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]

ex3_p1 = boxplot(total_surplus_durability_eq, legend = nothing, xlabel = "Trwałość dobra [liczba iteracji]", ylabel = "Nadwyżka całkowita", title = "Nadwyżka całkowita, dobra o identycznej jakości")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

savefig(ex3_p1, pwd() * "\\plots\\ex3_equal quality durability optim.svg")

total_surplus_durability_dq = [ex1_v10_total_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
ex3_p2 = boxplot(total_surplus_durability_dq, legend = nothing, xlabel = "Trwałość dobra [liczba iteracji]", ylabel = "Nadwyżka całkowita", title = "Nadwyżka całkowita, dobra o różnej jakości")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

savefig(ex3_p2, pwd() * "\\plots\\ex3_diff quality durability optim.svg")

# Consumer surplus

consumer_surplus_durability_eq = [ex1_v8_consumer_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]

boxplot(consumer_surplus_durability_eq, legend = nothing, xlabel = "Durability combination", ylabel = "Consumer surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

consumer_surplus_durability_dq = [ex1_v10_consumer_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(consumer_surplus_durability_dq, legend = nothing, xlabel = "Durability combination", ylabel = "Consumer surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

# Producer surplus

producer_surplus_durability_eq = [ex1_v8_producer_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_surplus_durability_eq, legend = nothing, xlabel = "Durability combination", ylabel = "Producer surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

producer_surplus_durability_dq = [ex1_v10_producer_surplus[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_surplus_durability_dq, legend = nothing, xlabel = "Durability combination", ylabel = "Producer surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

#### Producenci, niezależnie od oferowanej jakości produktów, chcą maksymalizować trwałość swoich produktów vs pozostali producenci

#### Producent 1 osiąga najwyższy zysk, gdy ma wysoką trwałość, a drugi najlepszy produkt niską

producer_1_durability_surplus = [getindex.(ex1_v10_producer_surplus_singleton,1)[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_1_durability_surplus, legend = nothing, xlabel = "Durability combination", ylabel = "Producer 1 surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

producer_2_durability_surplus = [getindex.(ex1_v10_producer_surplus_singleton,2)[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_2_durability_surplus, legend = nothing, xlabel = "Durability combination", ylabel = "Producer 2 surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

producer_3_durability_surplus = [getindex.(ex1_v10_producer_surplus_singleton,3)[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_3_durability_surplus, legend = nothing, xlabel = "Durability combination", ylabel = "Producer 2 surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

producer_4_durability_surplus = [getindex.(ex1_v10_producer_surplus_singleton,4)[[all(x .== y) for x in ex1_v710_d]] for y in ex1_v710_d_u]
boxplot(producer_4_durability_surplus, legend = nothing, xlabel = "Durability combination", ylabel = "Producer 2 surplus")
xticks!(1:length(ex1_v710_d_u), [join(string.(x)) for x in ex1_v710_d_u], xrotation = 90)

###



###

plot_ecdf(mean.(ex1_v7_price), "equal Q, equal D", "Price", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v8_price), "equal Q, ~equal D", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v9_price), "~equal Q, equal D", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v10_price), "~equal Q, ~equal D", "Price", "Probability", "ECDF", false)

plot_ecdf(mean.(ex1_v7_quantity), "equal Q, equal D", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v8_quantity), "equal Q, ~equal D", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v9_quantity), "~equal Q, equal D", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v10_quantity), "~equal Q, ~equal D", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean.(ex1_v7_advertising), "equal Q, equal D", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex1_v8_advertising), "equal Q, ~equal D", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v9_advertising), "~equal Q, equal D", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex1_v10_advertising), "~equal Q, ~equal D", "Advertising", "Probability", "ECDF", false)

