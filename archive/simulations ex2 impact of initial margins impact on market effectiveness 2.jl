include(pwd() * "\\methods\\methods.jl")

sim_single = TO_GO(4, 1000, 100, 0.25, 0.25; q = [1.25, 1.1, 0.9, 0.7], m = [0.2,0.1,0.1,0.1], c = [0.3, 0.3, 0.3, 0.3], ϵ = [0.33, 0.33, 0.33, 0.33], a = [0.05, 0.05, 0.05, 0.05], r = [0.01, 0.01, 0.01, 0.01], d = [7,5,3,1], num_links = 1000, future_discount = (0.3, 0.5), allow_negative_margin = false, consumer_behaviour = "stochastic", seller_behaviour = "deterministic")

sim_single = TO_GO(4, 500, 500, 0.25, 0.25; q = [1., 1., 1. ,1. ], m = [0.2,0.2,0.2,0.2], c = [0.5,0.5,0.5,0.5], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], d = [3,3,3,3], num_links = 1000, consumer_behaviour = "stochastic", seller_behaviour = "deterministic", allow_negative_margin = true)

plot_margin(sim_single, true)
plot_margin(sim_single, false)
plot(calculate_price_history.(sim_single.sellers))
plot_quantity(sim_single)
plot_quality_expectation(sim_single)
plot_profit_history(sim_single)
calculate_total_surplus(sim_single)
plot(calculate_expectation(sim_single, "durability"))

sim_single.buyers[1].durability_expectation_history

plot_advertising(sim_single)

mean(mean.(getfield.(sim_single.buyers, :ad_received_history)))

sim_single.sellers[1].advertising

scatter(sim_single.sellers[4].margin_history, calculate_profit_history(sim_single.sellers[4]))

plot(calculate_profit_history(sim_single.sellers[4]))
plot!(sim_single.sellers[4].quantity_history)
plot!(calculate_price_history(sim_single.sellers[4]))

scatter(calculate_price_history(sim_single.sellers[4]), calculate_profit_history(sim_single.sellers[4]))

sim_single.buyers[1].unit_possessed_history

plot([getindex.(sim_single.buyers[1].unit_possessed_history,x) for x in 1:4])

sim_single.buyers[1].std_reservation_price
sim_single.buyers[1].future_discount


# Cel zweryfikować w jaki sposób losowość marż początkowych wpływa na finalne wyniki symulacji.

ex2_v1_total_surplus = []
ex2_v1_producer_surplus = []
ex2_v1_consumer_surplus = []
ex2_v1_price = []
ex2_v1_quantity = []
ex2_v1_advertising = []
ex2_v1_producer_surplus_singleton = []
ex2_v1_quality_expectation = []
ex2_v1_durability_expectation = []
ex2_v1_margin = []

ex2_v2_total_surplus = []
ex2_v2_producer_surplus = []
ex2_v2_consumer_surplus = []
ex2_v2_price = []
ex2_v2_quantity = []
ex2_v2_advertising = []
ex2_v2_advertising_highest = []
ex2_v2_producer_surplus_singleton = []
ex2_v2_quality_expectation = []
ex2_v2_durability_expectation = []
ex2_v2_margin = []

ex2_v3_total_surplus = []
ex2_v3_producer_surplus = []
ex2_v3_consumer_surplus = []
ex2_v3_price = []
ex2_v3_quantity = []
ex2_v3_advertising = []
ex2_v3_advertising_highest = []
ex2_v3_producer_surplus_singleton = []
ex2_v3_quality_expectation = []
ex2_v3_durability_expectation = []
ex2_v3_margin = []

ex2_v4_total_surplus = []
ex2_v4_producer_surplus = []
ex2_v4_consumer_surplus = []
ex2_v4_price = []
ex2_v4_quantity = []
ex2_v4_advertising = []
ex2_v4_advertising_highest = []
ex2_v4_producer_surplus_singleton = []
ex2_v4_quality_expectation = []
ex2_v4_durability_expectation = []
ex2_v4_margin = []

ex2_v5_total_surplus = []
ex2_v5_producer_surplus = []
ex2_v5_consumer_surplus = []
ex2_v5_price = []
ex2_v5_quantity = []
ex2_v5_advertising = []
ex2_v5_advertising_highest = []
ex2_v5_producer_surplus_singleton = []
ex2_v5_quality_expectation = []
ex2_v5_durability_expectation = []
ex2_v5_margin = []

ex2_v6_total_surplus = []
ex2_v6_producer_surplus = []
ex2_v6_consumer_surplus = []
ex2_v6_price = []
ex2_v6_quantity = []
ex2_v6_advertising = []
ex2_v6_advertising_highest = []
ex2_v6_producer_surplus_singleton = []
ex2_v6_quality_expectation = []
ex2_v6_durability_expectation = []
ex2_v6_margin = []

ex2_v16_init_margin = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    m_d = Uniform(0,0.5)
    m_init = rand(m_d,2)

    push!(ex2_v16_init_margin, m_init)

    d_init = [2,2]

    # 1.0

    ex2_v1_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1., 1.], num_links = 1000)
    push!(ex2_v1_total_surplus, calculate_total_surplus(ex2_v1_sim, "total"))
    push!(ex2_v1_producer_surplus, calculate_total_surplus(ex2_v1_sim, "producer"))
    push!(ex2_v1_consumer_surplus, calculate_total_surplus(ex2_v1_sim, "consumer"))
    push!(ex2_v1_price, mean(calculate_price_history.(ex2_v1_sim.sellers)))
    push!(ex2_v1_quantity, mean(getfield.(ex2_v1_sim.sellers, :quantity_history)))
    push!(ex2_v1_advertising, mean(getfield.(ex2_v1_sim.sellers, :advertising_history)))
    push!(ex2_v1_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v1_sim.sellers, ex2_v1_sim.function_args.γ, ex2_v1_sim.function_args.num_buyers)))
    push!(ex2_v1_quality_expectation, calculate_expectation(ex2_v1_sim, "quality", true))
    push!(ex2_v1_durability_expectation, calculate_expectation(ex2_v1_sim, "durability", true))
    push!(ex2_v1_margin, getfield.(ex2_v1_sim.sellers, :margin_history))   

    # 1.1

    ex2_v2_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.1, 0.9], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.1, 0.9], num_links = 1000)
    push!(ex2_v2_total_surplus, calculate_total_surplus(ex2_v2_sim, "total"))
    push!(ex2_v2_producer_surplus, calculate_total_surplus(ex2_v2_sim, "producer"))
    push!(ex2_v2_consumer_surplus, calculate_total_surplus(ex2_v2_sim, "consumer"))
    push!(ex2_v2_price, mean(calculate_price_history.(ex2_v2_sim.sellers)))
    push!(ex2_v2_quantity, mean(getfield.(ex2_v2_sim.sellers, :quantity_history)))
    push!(ex2_v2_advertising, mean(getfield.(ex2_v2_sim.sellers, :advertising_history)))
    push!(ex2_v2_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v2_sim.sellers, ex2_v2_sim.function_args.γ, ex2_v2_sim.function_args.num_buyers)))
    push!(ex2_v2_quality_expectation, calculate_expectation(ex2_v2_sim, "quality", true))
    push!(ex2_v2_durability_expectation, calculate_expectation(ex2_v2_sim, "durability", true))
    push!(ex2_v2_margin, getfield.(ex2_v2_sim.sellers, :margin_history))   
  

    # 1.2

    ex2_v3_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.2, 0.8], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.2, 0.8], num_links = 1000)
    push!(ex2_v3_total_surplus, calculate_total_surplus(ex2_v3_sim, "total"))
    push!(ex2_v3_producer_surplus, calculate_total_surplus(ex2_v3_sim, "producer"))
    push!(ex2_v3_consumer_surplus, calculate_total_surplus(ex2_v3_sim, "consumer"))
    push!(ex2_v3_price, mean(calculate_price_history.(ex2_v3_sim.sellers)))
    push!(ex2_v3_quantity, mean(getfield.(ex2_v3_sim.sellers, :quantity_history)))
    push!(ex2_v3_advertising, mean(getfield.(ex2_v3_sim.sellers, :advertising_history)))
    push!(ex2_v3_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v3_sim.sellers, ex2_v3_sim.function_args.γ, ex2_v3_sim.function_args.num_buyers)))
    push!(ex2_v3_quality_expectation, calculate_expectation(ex2_v3_sim, "quality", true))
    push!(ex2_v3_durability_expectation, calculate_expectation(ex2_v3_sim, "durability", true))
    push!(ex2_v3_margin, getfield.(ex2_v3_sim.sellers, :margin_history))   
  

    # 1.3

    ex2_v4_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.3, 0.7], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.3, 0.7], num_links = 1000)
    push!(ex2_v4_total_surplus, calculate_total_surplus(ex2_v4_sim, "total"))
    push!(ex2_v4_producer_surplus, calculate_total_surplus(ex2_v4_sim, "producer"))
    push!(ex2_v4_consumer_surplus, calculate_total_surplus(ex2_v4_sim, "consumer"))
    push!(ex2_v4_price, mean(calculate_price_history.(ex2_v4_sim.sellers)))
    push!(ex2_v4_quantity, mean(getfield.(ex2_v4_sim.sellers, :quantity_history)))
    push!(ex2_v4_advertising, mean(getfield.(ex2_v4_sim.sellers, :advertising_history)))
    push!(ex2_v4_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v4_sim.sellers, ex2_v4_sim.function_args.γ, ex2_v4_sim.function_args.num_buyers)))
    push!(ex2_v4_quality_expectation, calculate_expectation(ex2_v4_sim, "quality", true))
    push!(ex2_v4_durability_expectation, calculate_expectation(ex2_v4_sim, "durability", true))
    push!(ex2_v4_margin, getfield.(ex2_v4_sim.sellers, :margin_history))   
  

    # 1.4

    ex2_v5_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.4, 0.6], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.4, 0.6], num_links = 1000)
    push!(ex2_v5_total_surplus, calculate_total_surplus(ex2_v5_sim, "total"))
    push!(ex2_v5_producer_surplus, calculate_total_surplus(ex2_v5_sim, "producer"))
    push!(ex2_v5_consumer_surplus, calculate_total_surplus(ex2_v5_sim, "consumer"))
    push!(ex2_v5_price, mean(calculate_price_history.(ex2_v5_sim.sellers)))
    push!(ex2_v5_quantity, mean(getfield.(ex2_v5_sim.sellers, :quantity_history)))
    push!(ex2_v5_advertising, mean(getfield.(ex2_v5_sim.sellers, :advertising_history)))
    push!(ex2_v5_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v5_sim.sellers, ex2_v5_sim.function_args.γ, ex2_v5_sim.function_args.num_buyers)))
    push!(ex2_v5_quality_expectation, calculate_expectation(ex2_v5_sim, "quality", true))
    push!(ex2_v5_durability_expectation, calculate_expectation(ex2_v5_sim, "durability", true))
    push!(ex2_v5_margin, getfield.(ex2_v5_sim.sellers, :margin_history))   
  

    # 1.5

    ex2_v6_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 0.5], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = d_init, q_init = [1.5, 0.5], num_links = 1000)
    push!(ex2_v6_total_surplus, calculate_total_surplus(ex2_v6_sim, "total"))
    push!(ex2_v6_producer_surplus, calculate_total_surplus(ex2_v6_sim, "producer"))
    push!(ex2_v6_consumer_surplus, calculate_total_surplus(ex2_v6_sim, "consumer"))
    push!(ex2_v6_price, mean(calculate_price_history.(ex2_v6_sim.sellers)))
    push!(ex2_v6_quantity, mean(getfield.(ex2_v6_sim.sellers, :quantity_history)))
    push!(ex2_v6_advertising, mean(getfield.(ex2_v6_sim.sellers, :advertising_history)))
    push!(ex2_v6_producer_surplus_singleton, sum.(calculate_profit_history.(ex2_v6_sim.sellers, ex2_v6_sim.function_args.γ, ex2_v6_sim.function_args.num_buyers)))
    push!(ex2_v6_quality_expectation, calculate_expectation(ex2_v6_sim, "quality", true))
    push!(ex2_v6_durability_expectation, calculate_expectation(ex2_v6_sim, "durability", true))
    push!(ex2_v6_margin, getfield.(ex2_v6_sim.sellers, :margin_history))   
  

end

## Wpływ marż wstępnych na dynamikę marż w czasie symulacji
# Istnieje silna zależność, im wyższa marża początkowa, tym wyższa średnia marża w czasie symulacji

scatter(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v1_margin], 1), markeralpha=0.25, smooth = true, xlabel = "Initial margin", ylabel = "Average margin d. simulation", label = "Q = (1.0, 1.0)")
scatter!(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v2_margin], 1), markeralpha=0.25, smooth = true, label = "Q = (1.1, 0.9)")
scatter!(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v3_margin], 1), markeralpha=0.25, smooth = true)
scatter!(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v4_margin], 1), markeralpha=0.25, smooth = true)
scatter!(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v5_margin], 1), markeralpha=0.25, smooth = true)
scatter!(getindex.(ex2_v16_init_margin, 1), getindex.([mean.(x) for x in ex2_v6_margin], 1), markeralpha=0.25, smooth = true)

### Wyższe poziomy marż pozwalają producentom na osiąganie wyższych zysków, jednak po osiągnięciu poziomu m_π* zysk zaczyna spadać 
### Producentom opłaca się zwiększać marżę, tak by maksymalizować zyski

ex2_p4 = scatter(getindex.([mean.(x) for x in ex2_v1_margin], 1), getindex.(ex2_v1_producer_surplus_singleton, 1), xlabel = "Średnia marża", ylabel = "Nadwyżka pierwszego producenta", label = "K = (1.0, 1.0)", markeralpha = 0.25, title = "Nadwyżka pierwszego producenta")
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v1_margin], 1), getindex.(ex2_v1_producer_surplus_singleton, 1), "blue")
scatter!(getindex.([mean.(x) for x in ex2_v6_margin], 1), getindex.(ex2_v6_producer_surplus_singleton, 1), xlabel = "Średnia marża", ylabel = "Nadwyżka pierwszego producenta", label = "K = (1.5, 0.5)", markeralpha = 0.25)
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v6_margin], 1), getindex.(ex2_v6_producer_surplus_singleton, 1), "green")

savefig(ex2_p4, pwd() * "\\plots\\ex2_p4 surplus 1.svg")

ex2_p5 = scatter(getindex.([mean.(x) for x in ex2_v1_margin], 2), getindex.(ex2_v1_producer_surplus_singleton, 2), xlabel = "Średnia marża", ylabel = "Nadwyżka drugiego producenta", label = "K = (1.0, 1.0)", markeralpha = 0.25, title = "Nadwyżka drugiego producenta")
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v1_margin], 2), getindex.(ex2_v1_producer_surplus_singleton, 2), "blue")
scatter!(getindex.([mean.(x) for x in ex2_v6_margin], 2), getindex.(ex2_v6_producer_surplus_singleton, 2), xlabel = "Średnia marża", ylabel = "Nadwyżka drugiego producenta", label = "K = (1.5, 0.5)", markeralpha = 0.25)
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v6_margin], 2), getindex.(ex2_v6_producer_surplus_singleton, 2), "green")

savefig(ex2_p5, pwd() * "\\plots\\ex2_p4 surplus 2.svg")

### Całkowita nadwyżka producenta na rynku bez zróżnicowania produktu nie zależy od różnic pomiędzy marżami producentów - z uwagi na substytucyjny charakter produktów

d(x) = x[1] - x[2]

scatter(d.([mean.(x) for x in ex2_v1_margin]), sum.(ex2_v1_producer_surplus_singleton), xlabel = "Margin difference", ylabel = "Producer surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v1_margin]), sum.(ex2_v1_producer_surplus_singleton), "blue")
scatter!(d.([mean.(x) for x in ex2_v1_margin]), getindex.(ex2_v1_producer_surplus_singleton, 1), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v1_margin]), getindex.(ex2_v1_producer_surplus_singleton, 1), "green")
scatter!(d.([mean.(x) for x in ex2_v1_margin]), getindex.(ex2_v1_producer_surplus_singleton, 2), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v1_margin]), getindex.(ex2_v1_producer_surplus_singleton, 2), "yellow")

scatter(d.([mean.(x) for x in ex2_v6_margin]), sum.(ex2_v6_producer_surplus_singleton), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v6_margin]), sum.(ex2_v6_producer_surplus_singleton), "blue")

scatter!(d.([mean.(x) for x in ex2_v6_margin]), getindex.(ex2_v6_producer_surplus_singleton, 1), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v6_margin]), getindex.(ex2_v6_producer_surplus_singleton, 1), "green")

scatter!(d.([mean.(x) for x in ex2_v6_margin]), getindex.(ex2_v6_producer_surplus_singleton, 2), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25)
add_smoothing_spline(d.([mean.(x) for x in ex2_v6_margin]), getindex.(ex2_v6_producer_surplus_singleton, 2), "yellow")

###



#margin_diff_12_34 = getindex.(ex2_v16_init_margin, 1) .- getindex.(ex2_v16_init_margin, 2)
margin_diff_12_34 = mean.(getindex.(ex2_v1_margin, 1)) .- mean.(getindex.(ex2_v1_margin, 2)) 

# Całkowita nadwyżka w zależności od marży początkowej: im wyższa różnica pomiędzy marżą firmy najlepszej vs. inne firmy, to efektywność rynku spada

ex2_p1 = scatter(margin_diff_12_34, ex2_v1_total_surplus, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Całkowita nadwyżka", label="Identyczna jakość", markeralpha = 0.25, title = "Nadwyżka całkowita a jakość dóbr")
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v1_total_surplus), "blue")
scatter!(margin_diff_12_34, ex2_v6_total_surplus, label="Rożna jakość", markeralpha = 0.25)
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v6_total_surplus), "green")

savefig(ex2_p1, pwd() * "\\plots\\ex2_surplus vs margin diff.svg")

ex2_p2 = scatter(margin_diff_12_34, ex2_v1_producer_surplus, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka producenta", label="Identyczna jakość", markeralpha = 0.25, title = "Nadwyżka całkowita")
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v1_producer_surplus), "blue")
scatter!(margin_diff_12_34, ex2_v6_producer_surplus, label="Rożna jakość", markeralpha = 0.25)
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v6_producer_surplus), "green")

# Nadwyżka konsumenta spada wraz ze wzrostem różnicy marż. Przyczyna wynika z faktu, że jeśli marża firmy o najwyższej jakości jest zbyt wysoka, to nawet konsumenci poszukujący jakości zdecydują się na zakup dobra o niższej jakości, tym samym ograniczając swoją nadwyżkę. Cel social plannera to zmusić firmę najlepszą jakościowo do obniżenia marży, a dla pozostałych firm to podnieść marże. Cele firmy najlepszej są sprzeczne z celami social plannera.

ex2_p3 = scatter(margin_diff_12_34, ex2_v1_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.0,1.0)", markeralpha = 0.25, title = "Nadwyżka konsumenta a jakość dóbr")
scatter!(margin_diff_12_34, ex2_v2_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.1,0.9)", markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v3_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.2,0.8)", markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v4_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.3,0.7)", markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v5_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.4,0.6)", markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v6_consumer_surplus, smooth = true, xlabel = "Różnica marży, pomiędzy producentem pierwszym a drugim", ylabel = "Nadwyżka konsumenta", label = "K = (1.5,0.5)", markeralpha = 0.25)

savefig(ex2_p3, pwd() * "\\plots\\ex2_consumer surplus vs margin diff.svg")