# Experiment 5

include(pwd() * "\\methods\\methods.jl")

# Impact of leasing on performance of a company

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

ex5_v3_total_surplus = []
ex5_v3_producer_surplus = []
ex5_v3_consumer_surplus = []
ex5_v3_price = []
ex5_v3_quality = []
ex5_v3_durability = []
ex5_v3_margin = []
ex5_v3_quantity_produced = []
ex5_v3_quantity_sold = []
ex5_v3_quantity_leased = []
ex5_v3_producer_surplus_singleton = []
ex5_v3_reselling = []
ex5_v3_consumer_surplus_pm = []
ex5_v3_consumer_surplus_sm_s = []
ex5_v3_consumer_surplus_sm_b = []

ex5_v4_total_surplus = []
ex5_v4_producer_surplus = []
ex5_v4_consumer_surplus = []
ex5_v4_price = []
ex5_v4_quality = []
ex5_v4_durability = []
ex5_v4_margin = []
ex5_v4_quantity_produced = []
ex5_v4_quantity_sold = []
ex5_v4_quantity_leased = []
ex5_v4_producer_surplus_singleton = []
ex5_v4_reselling = []
ex5_v4_consumer_surplus_pm = []
ex5_v4_consumer_surplus_sm_s = []
ex5_v4_consumer_surplus_sm_b = []

for i in 1:800

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex5_v1_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.75, 1.25], [0.75, 1.25]], [[0.2, 0.6], [0.2, 0.6]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 0)
    push!(ex5_v1_total_surplus, calculate_surplus(ex5_v1_sim, "total", false))
    push!(ex5_v1_producer_surplus, calculate_surplus(ex5_v1_sim, "producer", false))
    push!(ex5_v1_consumer_surplus, calculate_surplus(ex5_v1_sim, "consumer,total", false))
    push!(ex5_v1_price, calculate_price_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_quality, getfield.(ex5_v1_sim.sellers, :quality_history))
    push!(ex5_v1_durability, getfield.(ex5_v1_sim.sellers, :durability_history))
    push!(ex5_v1_margin, getfield.(ex5_v1_sim.sellers, :margin_history))
    push!(ex5_v1_quantity_produced, getfield.(ex5_v1_sim.sellers, :quantity_produced_history))
    push!(ex5_v1_quantity_sold, getfield.(ex5_v1_sim.sellers, :quantity_sold_history))
    push!(ex5_v1_quantity_leased, getfield.(ex5_v1_sim.sellers, :quantity_leased_history))
    push!(ex5_v1_producer_surplus_singleton, calculate_profit_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_reselling, getfield.(ex5_v1_sim.sellers, :reselling_history))
    push!(ex5_v1_consumer_surplus_pm, calculate_surplus(ex5_v1_sim, "consumer,pm", false))
    push!(ex5_v1_consumer_surplus_sm_b, calculate_surplus(ex5_v1_sim, "consumer,sm,b", false))
    push!(ex5_v1_consumer_surplus_sm_s, calculate_surplus(ex5_v1_sim, "consumer,sm,s", false))

    ex5_v2_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[1.0, 1.5], [0.5, 1.0]], [[0.3, 0.6], [0.2, 0.5]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 0)
    push!(ex5_v2_total_surplus, calculate_surplus(ex5_v2_sim, "total", false))
    push!(ex5_v2_producer_surplus, calculate_surplus(ex5_v2_sim, "producer", false))
    push!(ex5_v2_consumer_surplus, calculate_surplus(ex5_v2_sim, "consumer,total",false))
    push!(ex5_v2_price, calculate_price_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_quality, getfield.(ex5_v2_sim.sellers, :quality_history))
    push!(ex5_v2_durability, getfield.(ex5_v2_sim.sellers, :durability_history))
    push!(ex5_v2_margin, getfield.(ex5_v2_sim.sellers, :margin_history))
    push!(ex5_v2_quantity_produced, getfield.(ex5_v2_sim.sellers, :quantity_produced_history))
    push!(ex5_v2_quantity_sold, getfield.(ex5_v2_sim.sellers, :quantity_sold_history))
    push!(ex5_v2_quantity_leased, getfield.(ex5_v2_sim.sellers, :quantity_leased_history))
    push!(ex5_v2_producer_surplus_singleton, calculate_profit_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_reselling, getfield.(ex5_v2_sim.sellers, :reselling_history))
    push!(ex5_v2_consumer_surplus_pm, calculate_surplus(ex5_v2_sim, "consumer,pm", false))
    push!(ex5_v2_consumer_surplus_sm_b, calculate_surplus(ex5_v2_sim, "consumer,sm,b", false))
    push!(ex5_v2_consumer_surplus_sm_s, calculate_surplus(ex5_v2_sim, "consumer,sm,s", false))

    ex5_v3_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1],  [[0.75, 1.25], [0.75, 1.25]], [[0.2, 0.6], [0.2, 0.6]], [[.8, 2.], [.8, 2.]], 0.10, true, false, 0)
    push!(ex5_v3_total_surplus, calculate_surplus(ex5_v3_sim, "total", false))
    push!(ex5_v3_producer_surplus, calculate_surplus(ex5_v3_sim, "producer", false))
    push!(ex5_v3_consumer_surplus, calculate_surplus(ex5_v3_sim, "consumer,total",false))
    push!(ex5_v3_price, calculate_price_history.(ex5_v3_sim.sellers))
    push!(ex5_v3_quality, getfield.(ex5_v3_sim.sellers, :quality_history))
    push!(ex5_v3_durability, getfield.(ex5_v3_sim.sellers, :durability_history))
    push!(ex5_v3_margin, getfield.(ex5_v3_sim.sellers, :margin_history))
    push!(ex5_v3_quantity_produced, getfield.(ex5_v3_sim.sellers, :quantity_produced_history))
    push!(ex5_v3_quantity_sold, getfield.(ex5_v3_sim.sellers, :quantity_sold_history))
    push!(ex5_v3_quantity_leased, getfield.(ex5_v3_sim.sellers, :quantity_leased_history))
    push!(ex5_v3_producer_surplus_singleton, calculate_profit_history.(ex5_v3_sim.sellers))
    push!(ex5_v3_reselling, getfield.(ex5_v3_sim.sellers, :reselling_history))
    push!(ex5_v3_consumer_surplus_pm, calculate_surplus(ex5_v3_sim, "consumer,pm", false))
    push!(ex5_v3_consumer_surplus_sm_b, calculate_surplus(ex5_v3_sim, "consumer,sm,b", false))
    push!(ex5_v3_consumer_surplus_sm_s, calculate_surplus(ex5_v3_sim, "consumer,sm,s", false))


    ex5_v4_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[1.0, 1.5], [0.5, 1.0]], [[0.3, 0.6], [0.2, 0.5]], [[.8, 2.], [.8, 2.]], 0.10, true, false, 0)
    push!(ex5_v4_total_surplus, calculate_surplus(ex5_v4_sim, "total", false))
    push!(ex5_v4_producer_surplus, calculate_surplus(ex5_v4_sim, "producer", false))
    push!(ex5_v4_consumer_surplus, calculate_surplus(ex5_v4_sim, "consumer,total",false))
    push!(ex5_v4_price, calculate_price_history.(ex5_v4_sim.sellers))
    push!(ex5_v4_quality, getfield.(ex5_v4_sim.sellers, :quality_history))
    push!(ex5_v4_durability, getfield.(ex5_v4_sim.sellers, :durability_history))
    push!(ex5_v4_margin, getfield.(ex5_v4_sim.sellers, :margin_history))
    push!(ex5_v4_quantity_produced, getfield.(ex5_v4_sim.sellers, :quantity_produced_history))
    push!(ex5_v4_quantity_sold, getfield.(ex5_v4_sim.sellers, :quantity_sold_history))
    push!(ex5_v4_quantity_leased, getfield.(ex5_v4_sim.sellers, :quantity_leased_history))
    push!(ex5_v4_producer_surplus_singleton, calculate_profit_history.(ex5_v4_sim.sellers))
    push!(ex5_v4_reselling, getfield.(ex5_v4_sim.sellers, :reselling_history))
    push!(ex5_v4_consumer_surplus_pm, calculate_surplus(ex5_v4_sim, "consumer,pm", false))
    push!(ex5_v4_consumer_surplus_sm_b, calculate_surplus(ex5_v4_sim, "consumer,sm,b", false))
    push!(ex5_v4_consumer_surplus_sm_s, calculate_surplus(ex5_v4_sim, "consumer,sm,s", false))


end

# Różnice w jakości nie mają istotnego wpływu na kształtowanie się trwałości produktów. Leasing ma istotny wpływ

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_durability)), "Identyczna jakość, leasing możliwy"; xlabel = "Średnia trwałość dóbr", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Trwałość")
plot_ecdf(false, mean.(mean.(ex5_v2_durability)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_durability)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_durability)), "Różna jakość, leasing niemożliwy")

# Ma to zastosowanie do średniej rynkowej trwałości, jak również trwałości wszystkich graczy rynkowych

plot_ecdf(true, mean.(getindex.(ex5_v1_durability, 1)), "Identyczna jakość, leasing możliwy"; xlabel = "Średnia trwałość dóbr", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(getindex.(ex5_v2_durability, 1)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(getindex.(ex5_v3_durability, 1)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(getindex.(ex5_v4_durability, 1)), "Różna jakość, leasing niemożliwy")

# Różnica w poziomie średniej rynkowej trwałości pomiędzy rynkiem z leasingiem i bez leasingu, jest wyższa na rynku gdzie dobra są homogeniczne

UnequalVarianceTTest(sort(mean.(getindex.(ex5_v1_durability, 1))) .- sort(mean.(getindex.(ex5_v3_durability, 1))), sort(mean.(getindex.(ex5_v2_durability, 1))) .- sort(mean.(getindex.(ex5_v4_durability, 1))))

mean(sort(mean.(getindex.(ex5_v1_durability, 1))) .- sort(mean.(getindex.(ex5_v3_durability, 1))))
mean(sort(mean.(getindex.(ex5_v2_durability, 1))) .- sort(mean.(getindex.(ex5_v4_durability, 1))))

# Leasing nie ma istotnego wpływu na różnice w decyzjach odnośnie jakości produktów

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quality)), "Identyczna jakość, leasing możliwy"; xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(mean.(ex5_v2_quality)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quality)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quality)), "Różna jakość, leasing niemożliwy")

UnequalVarianceTTest(mean.(mean.(ex5_v1_quality)), mean.(mean.(ex5_v3_quality)))
UnequalVarianceTTest(mean.(mean.(ex5_v2_quality)), mean.(mean.(ex5_v4_quality)))

# Na rynku gdzie możliwy jest leasing produkowana jest wyższa liczba produktów

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quantity_produced)), "Identyczna jakość, leasing możliwy", xlabel = "Średni wolumen produkcji", ylabel = "F(x)", title = "Dystrybuanta empiryczna - Wielkość produkcji")
plot_ecdf(false, mean.(mean.(ex5_v2_quantity_produced)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quantity_produced)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quantity_produced)), "Różna jakość, leasing niemożliwy")

UnequalVarianceTTest(mean.(mean.(ex5_v1_quantity_produced)), mean.(mean.(ex5_v1_quantity_produced)))
UnequalVarianceTTest(mean.(mean.(ex5_v1_quantity_produced)), mean.(mean.(ex5_v1_quantity_produced)))

ex5_p1 = plot_ecdf(true, mean.(mean.(ex5_v1_quantity_produced)) .- mean.(mean.(ex5_v1_quantity_sold)) .- mean.(mean.(ex5_v1_quantity_leased)), "Identyczna jakość, leasing możliwy", xlabel = "Liczba niesprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba niesprzedanych sztuk")
plot_ecdf(false, mean.(mean.(ex5_v2_quantity_produced)).- mean.(mean.(ex5_v2_quantity_sold)) .- mean.(mean.(ex5_v2_quantity_leased)), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(ex5_v3_quantity_produced)).- mean.(mean.(ex5_v3_quantity_sold)) .- mean.(mean.(ex5_v3_quantity_leased)), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(ex5_v4_quantity_produced)).- mean.(mean.(ex5_v4_quantity_sold)) .- mean.(mean.(ex5_v4_quantity_leased)), "Różna jakość, leasing niemożliwy")

# W rezultacie, leasing prowadzi do wyższej nadwyżki producenta

ex5_p1 = plot_ecdf(true, mean.(mean.(getindex.(ex5_v1_producer_surplus_singleton, 1))), "Identyczna jakość, leasing możliwy"; xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v2_producer_surplus_singleton, 1))), "Różna jakość, leasing możliwy")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v3_producer_surplus_singleton, 1))), "Identyczna jakość, leasing niemożliwy")
plot_ecdf(false, mean.(mean.(getindex.(ex5_v4_producer_surplus_singleton, 1))), "Różna jakość, leasing niemożliwy")

# Wyższy % leasingowanych produktów, prowadzi do wyższych trwałości
# Czy odwrotnie? Wyższe trwałości prowadzą konsumentów do chętniejszego wybierania leasingu
# Czy to nie jest zaprzeczenie teorii Coase-Swana? Albo uszczegółowienie o kierunek przyczynowości?

ex5_v1_perc_leased = mean.(getindex.(ex5_v1_quantity_leased, 1)) ./ mean.(getindex.(ex5_v1_quantity_produced, 1))
ex5_v2_perc_leased = mean.(getindex.(ex5_v2_quantity_leased, 1)) ./ mean.(getindex.(ex5_v2_quantity_produced, 1))
ex5_v3_perc_leased = mean.(getindex.(ex5_v3_quantity_leased, 1)) ./ mean.(getindex.(ex5_v3_quantity_produced, 1))
ex5_v4_perc_leased = mean.(getindex.(ex5_v4_quantity_leased, 1)) ./ mean.(getindex.(ex5_v4_quantity_produced, 1))

scatter(ex5_v1_perc_leased[mean.(getindex.(ex5_v1_durability,1)) .> 0.32], mean.(getindex.(ex5_v1_durability,1))[mean.(getindex.(ex5_v1_durability,1)) .> 0.32], smooth = true)
scatter!(ex5_v2_perc_leased[mean.(getindex.(ex5_v2_durability,1)) .> 0.42], mean.(getindex.(ex5_v2_durability,1))[mean.(getindex.(ex5_v2_durability,1)) .> 0.42], smooth = true)

# WAŻNE: Im wyższa jakość produktu, tym silnejszy wzrost nadwyżki producenta wraz ze wzrostem trwałości produktu. Producentom dóbr wysokiej jakości opłacalne jest dbanie o ich wysoką trwałość.

scatter(mean.(getindex.(ex5_v1_durability,1)), mean.(getindex.(ex5_v1_producer_surplus_singleton,1)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing możliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v2_durability,1)), mean.(getindex.(ex5_v2_producer_surplus_singleton,1)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing możliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v3_durability,1)), mean.(getindex.(ex5_v3_producer_surplus_singleton,1)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing niemożliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v4_durability,1)), mean.(getindex.(ex5_v4_producer_surplus_singleton,1)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing niemożliwy", markerstrokewidth = 0)


scatter(mean.(getindex.(ex5_v1_durability,2)), mean.(getindex.(ex5_v1_producer_surplus_singleton,2)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing możliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v2_durability,2)), mean.(getindex.(ex5_v2_producer_surplus_singleton,2)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing możliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v3_durability,2)), mean.(getindex.(ex5_v3_producer_surplus_singleton,2)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing niemożliwy", markerstrokewidth = 0)
scatter!(mean.(getindex.(ex5_v4_durability,2)), mean.(getindex.(ex5_v4_producer_surplus_singleton,2)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing niemożliwy", markerstrokewidth = 0)


# Wyższa jakość dóbr podnosi poziom dobrobytu. Zgodne z wynikiem Coase.

scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing możliwy")
scatter!(mean.(mean.(ex5_v2_durability)), sum.(sum.(ex5_v2_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing możliwy")
scatter!(mean.(mean.(ex5_v3_durability)), sum.(sum.(ex5_v3_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość, leasing niemożliwy")
scatter!(mean.(mean.(ex5_v4_durability)), sum.(sum.(ex5_v4_total_surplus)), smooth = true, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, leasing niemożliwy")

# Nie istnieje statystycznie istotna zależność pomiędzy trwałością a nadwyżką konsumenta. 
# Do sprawdzenia - czy ma na to wpływ c? 🍎

scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus)))
scatter!(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_producer_surplus)))

# Ma za to wpływ na nadwyżkę na rynku wtórnym, co wydaje się poprawne.

scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_pm)), smooth=true)
scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_sm_b)), smooth=true)
scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_consumer_surplus_sm_s)), smooth=true)

scatter(mean.(mean.(ex5_v1_durability)), mean.(mean.(ex5_v1_price)))

plot(sort(0.5 * mean.(mean.(ex5_v1_margin))))

scatter!(mean.(mean.(ex5_v1_durability)), 0.5 * sum_of_geom_series.(mean.(mean.(ex5_v1_quality)), mean.(mean.(ex5_v1_durability))))

scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_quantity_sold)), smooth = true)
scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_quantity_leased)))



scatter(mean.(mean.(ex5_v1_durability)), sum.(sum.(ex5_v1_reselling)))
scatter(mean.(mean.(ex5_v2_durability)), sum.(sum.(ex5_v2_reselling)))
scatter(mean.(mean.(ex5_v3_durability)), sum.(sum.(ex5_v3_reselling)))
scatter(mean.(mean.(ex5_v4_durability)), sum.(sum.(ex5_v4_reselling)), smooth = true)