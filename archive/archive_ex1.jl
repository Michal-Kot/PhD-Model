####################### SPECIAL CASE: NOT EQUAL INITIAL MARGINS

# Cel: sprawdzić, w jaki sposób wstępny poziom marż wpływa na końcowe rezultaty symulacji (jeśli w ogóle)
# Wynik: jeśli wstępny poziom marż jest losowy, ale równy dla wszystkich producentów, to nie ma istotnego wpływu na końcowy rezultat niezależnie od różnic w poziomie jakości.
# Jeśli występują różnice w poziomie jakości i wstępne marże nie są równe dla wszystkich producentów, to obserwowana jest relacja: im wyższa marża początkowa i jakość produktu, tym niższy końcowy wynik. Im różnica między producentami jest wyższa (tj. producent najlepszego produktu ma najwyższą marżę, reszta producentów ma je odnowiednio niższe), tym wyższa różnica w wyniku końcowym.

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



margin_diff_12_34 = (getindex.(ex1_v34_init_margin, 1) .+ getindex.(ex1_v34_init_margin, 2)) / 2 - (getindex.(ex1_v34_init_margin, 3) .+ getindex.(ex1_v34_init_margin, 4)) / 2

scatter(margin_diff_12_34, ex1_v3_total_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v3_producer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v3_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)

# Na rynku gdzie występuje różnica jakości, różnica marż wpływa na nadwyżkę konsumenta, nie wpływając istotnie na nadwyżkę producentów.

scatter(margin_diff_12_34, ex1_v4_total_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v4_producer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)
scatter(margin_diff_12_34, ex1_v4_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing)

# Nadwyżka konsumenta, obliczana per capita jako std_reservation_price * quality_expectation * duration_expectation * wtp_advertising - price
# Różnica marż podnosi średnią cenę produktów kupowanych (nie oferowanych) na rynku

scatter(margin_diff_12_34, mean.(ex1_v4_price), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Price", legend=nothing)

GLM.lm(@formula(price~margin), DataFrame(margin = margin_diff_12_34, price = mean.(ex1_v4_price))) # efekt w przypadku rynku ze zróżnicowaniem produktu jest istotny stat.
GLM.lm(@formula(price~margin), DataFrame(margin = margin_diff_12_34, price = mean.(ex1_v3_price))) # efekty w przypadku równych jakości nie jest istotny statystycznie

scatter(margin_diff_12_34, getindex.(ex1_v4_durability_expectation, 1), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Expected durability, high-end producer", legend=nothing)

GLM.lm(@formula(durabilitye~margin), DataFrame(margin = margin_diff_12_34, durabilitye = mean.(ex1_v4_durability_expectation)))
GLM.lm(@formula(durabilitye~margin), DataFrame(margin = margin_diff_12_34, durabilitye = mean.(ex1_v3_durability_expectation)))

scatter(margin_diff_12_34, getindex.(ex1_v4_quality_expectation, 1), smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Expected quality, high-end producer", legend=nothing)

GLM.lm(@formula(qualitye~margin), DataFrame(margin = margin_diff_12_34, qualitye = mean.(ex1_v4_quality_expectation)))
GLM.lm(@formula(qualitye~margin), DataFrame(margin = margin_diff_12_34, qualitye = mean.(ex1_v3_quality_expectation)))



scatter(margin_diff_12_34, mean.([x .* y for (x,y) in zip(ex1_v4_quality_expectation, ex1_v4_durability_expectation)]),smooth=true)
scatter(margin_diff_12_34, mean.(ex1_v4_quality_expectation), smooth = true, ylim=(0.9,1.0))
scatter(margin_diff_12_34, mean.(ex1_v4_durability_expectation), smooth = true)

scatter(margin_diff_12_34, mean.(ex1_v3_price), smooth = true)
scatter!(margin_diff_12_34, mean.(ex1_v4_price), smooth = true)

scatter(margin_diff_12_34, 0.5 * getindex.(ex1_v4_durability_expectation, 1) .* getindex.(ex1_v4_quality_expectation, 1), smooth = true)

scatter(margin_diff_12_34, getindex.(ex1_v4_durability_expectation, 4))

###############

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