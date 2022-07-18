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

####################### SPECIAL CASE: NOT EQUAL INITIAL MARGINS

# Cel: sprawdzić, w jaki sposób wstępny poziom marż wpływa na końcowe rezultaty symulacji (jeśli w ogóle)
# Wynik: jeśli wstępny poziom marż jest losowy, ale równy dla wszystkich producentów, to nie ma istotnego wpływu na końcowy rezultat niezależnie od różnic w poziomie jakości.
# Jeśli występują różnice w poziomie jakości i wstępne marże nie są równe dla wszystkich producentów, to obserwowana jest relacja: im wyższa marża początkowa i jakość produktu, tym niższy końcowy wynik. Im różnica między producentami jest wyższa (tj. producent najlepszego produktu ma najwyższą marżę, reszta producentów ma je odnowiednio niższe), tym wyższa różnica w wyniku końcowym.

sim_with_obs_12.buyers[1].durability_expectation

function calculate_expectation(sim_res, metric, cumulated = false)

    if metric == "quality"
        if cumulated
            return mean.([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers])
        else
            return [getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers]
        end
    elseif metric == "durability"
        if cumulated
            return mean.([getindex.(mean(getfield.(sim_res.buyers, :durability_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers])
        else
            return [getindex.(mean(getfield.(sim_res.buyers, :durability_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers]
        end
    end

end

ex1_v3_total_surplus = []
ex1_v3_producer_surplus = []
ex1_v3_consumer_surplus = []
ex1_v3_price = []
ex1_v3_quantity = []
ex1_v3_advertising = []
ex1_v3_producer_surplus_singleton = []
ex1_v3_quality_expectation = []
ex1_v3_durability_expectation = []

ex1_v4_total_surplus = []
ex1_v4_producer_surplus = []
ex1_v4_consumer_surplus = []
ex1_v4_price = []
ex1_v4_quantity = []
ex1_v4_advertising = []
ex1_v4_advertising_highest = []
ex1_v4_producer_surplus_singleton = []
ex1_v4_quality_expectation = []
ex1_v4_durability_expectation = []

ex1_v34_init_margin = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.2)
    m_init = rand(m_d,4)

    push!(ex1_v34_init_margin, m_init)

    ex1_v3_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0, 1.0, 1.0], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], q_init = [1., 1., 1., 1.], num_links = 1000)
    push!(ex1_v3_total_surplus, calculate_total_surplus(ex1_v3_sim, "total"))
    push!(ex1_v3_producer_surplus, calculate_total_surplus(ex1_v3_sim, "producer"))
    push!(ex1_v3_consumer_surplus, calculate_total_surplus(ex1_v3_sim, "consumer"))
    push!(ex1_v3_price, mean(calculate_price_history.(ex1_v3_sim.sellers)))
    push!(ex1_v3_quantity, mean(getfield.(ex1_v3_sim.sellers, :quantity_history)))
    push!(ex1_v3_advertising, mean(getfield.(ex1_v3_sim.sellers, :advertising_history)))
    push!(ex1_v3_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v3_sim.sellers, ex1_v3_sim.function_args.γ, ex1_v3_sim.function_args.num_buyers)))
    push!(ex1_v3_quality_expectation, calculate_expectation(ex1_v3_sim, "quality", true))
    push!(ex1_v3_durability_expectation, calculate_expectation(ex1_v3_sim, "durability", true))   

    ex1_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex1_v4_total_surplus, calculate_total_surplus(ex1_v4_sim, "total"))
    push!(ex1_v4_producer_surplus, calculate_total_surplus(ex1_v4_sim, "producer"))
    push!(ex1_v4_consumer_surplus, calculate_total_surplus(ex1_v4_sim, "consumer"))
    push!(ex1_v4_price, mean(calculate_price_history.(ex1_v4_sim.sellers)))
    push!(ex1_v4_quantity, mean(getfield.(ex1_v4_sim.sellers, :quantity_history)))
    push!(ex1_v4_advertising, mean(getfield.(ex1_v4_sim.sellers, :advertising_history)))
    push!(ex1_v4_producer_surplus_singleton, sum.(calculate_profit_history.(ex1_v4_sim.sellers, ex1_v4_sim.function_args.γ, ex1_v4_sim.function_args.num_buyers)))
    push!(ex1_v4_quality_expectation, calculate_expectation(ex1_v4_sim, "quality", true))
    push!(ex1_v4_durability_expectation, calculate_expectation(ex1_v4_sim, "durability", true))   

end

# W jaki sposób różnica marż pomiędzy producentami wysokiej jakości i niskiej jakości wpływa na efektwyność rynku?

# Różnica marż nie ma istotnego wpływu na nadwyżkę na rynku, gdzie nie występuje zróżnicowanie jakości.

margin_diff_12_34 = (getindex.(ex1_v34_init_margin, 1) .+ getindex.(ex1_v34_init_margin, 2)) / 2 - (getindex.(ex1_v34_init_margin, 3) .+ getindex.(ex1_v34_init_margin, 4)) / 2

scatter(margin_diff_12_34, ex1_v3_total_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v3_producer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v3_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)

# Na rynku gdzie występuje różnica jakości, różnica marż wpływa na nadwyżkę konsumenta, nie wpływając istotnie na nadwyżkę producentów.

scatter(margin_diff_12_34, ex1_v4_total_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v4_producer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v4_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)

scatter(margin_diff_12_34, mean.(ex1_v4_price), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Price", legend=nothing)

GLM.lm(@formula(price~margin), DataFrame(margin = margin_diff_12_34, price = mean.(ex1_v4_price)))
GLM.lm(@formula(price~margin), DataFrame(margin = margin_diff_12_34, price = mean.(ex1_v3_price)))

scatter(margin_diff_12_34, getindex.(ex1_v4_durability_expectation, 1), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Expected durability, high-end producer", legend=nothing)
scatter(margin_diff_12_34, getindex.(ex1_v4_quality_expectation, 1), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Expected quality, high-end producer", legend=nothing)

GLM.lm(@formula(qualitye~margin), DataFrame(margin = margin_diff_12_34, qualitye = mean.(ex1_v4_quality_expectation)))
GLM.lm(@formula(qualitye~margin), DataFrame(margin = margin_diff_12_34, qualitye = mean.(ex1_v3_quality_expectation)))

GLM.lm(@formula(durabilitye~margin), DataFrame(margin = margin_diff_12_34, durabilitye = mean.(ex1_v4_durability_expectation)))
GLM.lm(@formula(durabilitye~margin), DataFrame(margin = margin_diff_12_34, durabilitye = mean.(ex1_v3_durability_expectation)))

scatter(margin_diff_12_34, mean.([x .* y for (x,y) in zip(ex1_v4_quality_expectation, ex1_v4_durability_expectation)]),smooth=true)
scatter(margin_diff_12_34, mean.(ex1_v4_quality_expectation), smooth = true, ylim=(0.9,1.0))
scatter(margin_diff_12_34, mean.(ex1_v4_durability_expectation), smooth = true)

scatter(margin_diff_12_34, mean.(ex1_v3_price), smooth = true)
scatter!(margin_diff_12_34, mean.(ex1_v4_price), smooth = true)

scatter(margin_diff_12_34, 0.5 * getindex.(ex1_v4_durability_expectation, 1) .* getindex.(ex1_v4_quality_expectation, 1), smooth = true)

scatter(margin_diff_12_34, getindex.(ex1_v4_durability_expectation, 4))

###############

fill(rand()*0.2,4)

ex1_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = rand(Uniform(0,0.2), 4), c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [2,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])



ex1_v4_sim.function_args.m
plot_margin(ex1_v4_sim)
plot_quantity(ex1_v4_sim)
plot(calculate_expectation(ex1_v4_sim, "durability", false))



###############

scatter(getindex.(ex1_v34_init_margin, 1), getindex.(ex1_v3_producer_surplus_singleton, 1), markeralpha = 0.25, smooth = true, ylim = (0, 3000), label = "equal quality")
scatter!(getindex.(ex1_v34_init_margin, 1), getindex.(ex1_v4_producer_surplus_singleton, 1), markeralpha = 0.25, smooth = true, ylim = (0,3000), label = "not equal quality", xlabel = "Initial margin", ylabel = "Profit")

lm_df = DataFrame(x = margin_diff_12_34, y_eq = getindex.(ex1_v3_producer_surplus_singleton, 1), y_df = getindex.(ex1_v4_producer_surplus_singleton, 1))

model_eq = GLM.lm(@formula(y_eq~x), lm_df)
model_df = GLM.lm(@formula(y_df~x), lm_df)

coef_diff_Paternoster(β1,β2,SEβ1,SEβ2) = (β1 - β2)/sqrt(SEβ1^2+SEβ2^2)

coef_diff_Paternoster(coef(model_eq)[2], coef(model_df)[2], stderror(model_eq)[2], stderror(model_df)[2]) # różnica nie jest istotna statystyczne nawet dla α = 0.10

plot_ecdf(ex1_v3_total_surplus, "equal Q", "Total Surplus", "Probability", "ECDF", true)
plot_ecdf(ex1_v4_total_surplus, "not equal Q", "Total Surplus", "Probability", "ECDF", false)

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

# Cel 1: Sprawdzenie, w jaki sposób zmienia się efektywność komunikacji reklamowej, gdy na rynku występuje (oraz nie) zróżnicowanie jakości produktów

# Wynik:

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

    m_init = rand(Uniform(0,0.2), 4)

    ex2_v1_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1., 1., 1., 1.])
    push!(ex2_v1_total_surplus, calculate_total_surplus(ex2_v1_sim, "total"))
    push!(ex2_v1_producer_surplus, calculate_total_surplus(ex2_v1_sim, "producer"))
    push!(ex2_v1_consumer_surplus, calculate_total_surplus(ex2_v1_sim, "consumer"))
    push!(ex2_v1_price, mean(calculate_price_history.(ex2_v1_sim.sellers)))
    push!(ex2_v1_quantity, mean(getfield.(ex2_v1_sim.sellers, :quantity_history)))
    push!(ex2_v1_advertising, mean(getfield.(ex2_v1_sim.sellers, :advertising_history)))

    ex2_v2_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1., 1., 1. ,1. ], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1., 1., 1., 1.], variant_advertising = false)
    push!(ex2_v2_total_surplus, calculate_total_surplus(ex2_v2_sim, "total"))
    push!(ex2_v2_producer_surplus, calculate_total_surplus(ex2_v2_sim, "producer"))
    push!(ex2_v2_consumer_surplus, calculate_total_surplus(ex2_v2_sim, "consumer"))
    push!(ex2_v2_price, mean(calculate_price_history.(ex2_v2_sim.sellers)))
    push!(ex2_v2_quantity, mean(getfield.(ex2_v2_sim.sellers, :quantity_history)))
    push!(ex2_v2_advertising, mean(getfield.(ex2_v2_sim.sellers, :advertising_history)))

    ex2_v3_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50])
    push!(ex2_v3_total_surplus, calculate_total_surplus(ex2_v3_sim, "total"))
    push!(ex2_v3_producer_surplus, calculate_total_surplus(ex2_v3_sim, "producer"))
    push!(ex2_v3_consumer_surplus, calculate_total_surplus(ex2_v3_sim, "consumer"))
    push!(ex2_v3_price, mean(calculate_price_history.(ex2_v3_sim.sellers)))
    push!(ex2_v3_quantity, mean(getfield.(ex2_v3_sim.sellers, :quantity_history)))
    push!(ex2_v3_advertising, mean(getfield.(ex2_v3_sim.sellers, :advertising_history)))

    ex2_v4_sim = TO_GO(4, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 1.25, 0.75, 0.50], m = m_init, c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [1,1,1,1], num_links = 1000, q_init = [1.5, 1.25, 0.75, 0.50], variant_advertising = false)
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

plot_ecdf(mean.(ex2_v1_price), "equal Q, advertising", "Price", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v2_price), "equal Q, no advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v3_price), "not equal Q, advertising", "Price", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v4_price), "not equal Q, no advertising", "Price", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v1_quantity), "equal Q, advertising", "Quantity", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v2_quantity), "equal Q, no advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v3_quantity), "not equal Q, advertising", "Quantity", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v4_quantity), "not equal Q, no advertising", "Quantity", "Probability", "ECDF", false)

plot_ecdf(mean.(ex2_v1_advertising), "equal Q, advertising", "Advertising", "Probability", "ECDF", true)
plot_ecdf(mean.(ex2_v2_advertising), "equal Q, no advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v3_advertising), "not equal Q, advertising", "Advertising", "Probability", "ECDF", false)
plot_ecdf(mean.(ex2_v4_advertising), "not equal Q, no advertising", "Advertising", "Probability", "ECDF", false)

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