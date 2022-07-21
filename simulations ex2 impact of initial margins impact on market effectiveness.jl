include(pwd() * "\\methods\\methods.jl")

sim_single = TO_GO(4, 500, 5000, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6,0.6,0.6,0.6], ϵ = [0.33,0.33,0.33,0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 1000)

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

    # 1.0

    ex2_v1_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.0, 1.0], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1., 1.], num_links = 1000)
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

    ex2_v2_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.1, 0.9], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1.1, 0.9], num_links = 1000)
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

    ex2_v3_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.2, 0.8], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1.2, 0.8], num_links = 1000)
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

    ex2_v4_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.3, 0.7], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1.3, 0.7], num_links = 1000)
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

    ex2_v5_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.4, 0.6], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1.4, 0.6], num_links = 1000)
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

    ex2_v6_sim = TO_GO(2, 500, 250, 0.25, 0.25, "deterministic"; q = [1.5, 0.5], m = m_init, c = [0.6,0.6], ϵ = [0.33,0.33], a = [0.0, 0.0], r = [0.4, 0.4], d = [1,1], q_init = [1.5, 0.5], num_links = 1000)
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

scatter(getindex.([mean.(x) for x in ex2_v1_margin], 1), getindex.(ex2_v1_producer_surplus_singleton, 1), xlabel = "Average margin d. simulation", ylabel = "Producer 1 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25, ylim = (0,15000))
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v1_margin], 1), getindex.(ex2_v1_producer_surplus_singleton, 1), "blue")
scatter!(getindex.([mean.(x) for x in ex2_v6_margin], 1), getindex.(ex2_v6_producer_surplus_singleton, 1), xlabel = "Average margin d. simulation", ylabel = "Producer 1 surplus", label = "Q = (1.5, 0.5)", markeralpha = 0.25)
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v6_margin], 1), getindex.(ex2_v6_producer_surplus_singleton, 1), "green")

scatter(getindex.([mean.(x) for x in ex2_v1_margin], 2), getindex.(ex2_v1_producer_surplus_singleton, 2), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.0, 1.0)", markeralpha = 0.25, ylim = (0, 15000))
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v1_margin], 2), getindex.(ex2_v1_producer_surplus_singleton, 2), "blue")
scatter!(getindex.([mean.(x) for x in ex2_v6_margin], 2), getindex.(ex2_v6_producer_surplus_singleton, 2), xlabel = "Average margin d. simulation", ylabel = "Producer 2 surplus", label = "Q = (1.5, 0.5)", markeralpha = 0.25)
add_smoothing_spline(getindex.([mean.(x) for x in ex2_v6_margin], 2), getindex.(ex2_v6_producer_surplus_singleton, 2), "green")

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

margin_diff_12_34 = getindex.(ex2_v16_init_margin, 1) .- getindex.(ex2_v16_init_margin, 2)

# Całkowita nadwyżka w zależności od marży początkowej: im wyższa różnica pomiędzy marżą firmy najlepszej vs. inne firmy, to efektywność rynku spada

scatter(margin_diff_12_34, ex2_v1_total_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Total market surplus", legend=nothing, markeralpha = 0.25)
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v1_total_surplus), "blue")
scatter!(margin_diff_12_34, ex2_v6_total_surplus, smooth = true, legend=nothing, markeralpha = 0.25)
add_smoothing_spline(margin_diff_12_34, Float64.(ex2_v6_total_surplus), "green")

# Nadwyżka konsumenta spada wraz ze wzrostem różnicy marż. Przyczyna wynika z faktu, że jeśli marża firmy o najwyższej jakości jest zbyt wysoka, to nawet konsumenci poszukujący jakości zdecydują się na zakup dobra o niższej jakości, tym samym ograniczając swoją nadwyżkę. Cel social plannera to zmusić firmę najlepszą jakościowo do obniżenia marży, a dla pozostałych firm to podnieść marże. Cele firmy najlepszej są sprzeczne z celami social plannera.

scatter(margin_diff_12_34, ex2_v1_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v2_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v3_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v4_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v5_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)
scatter!(margin_diff_12_34, ex2_v6_consumer_surplus, smooth = true, xlabel = "Margin difference between high quality and low-quality producers", ylabel = "Profit", legend=nothing, markeralpha = 0.25)



########## HypothesisTests
sim_single_1111 = TO_GO(4, 250, 500, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6, 0.6, 0.6, 0.6], ϵ = [0.33, 0.33, 0.33, 0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 1000, d = [2, 2, 2, 2])

sim_single_4321 = TO_GO(4, 250, 500, 0.25, 0.25, "deterministic"; q = [1.6, 1.3, 1.0, 0.7], m = [0.2, 0.2, 0.2, 0.2], c = [0.6, 0.6, 0.6, 0.6], ϵ = [0.33, 0.33, 0.33, 0.33], a = [0.0, 0.0, 0.0, 0.0], r = [0.4, 0.4, 0.4, 0.4], q_init = [1.6, 1.3, 1.0, 0.7], num_links = 1000, d = [2, 3, 2, 2])

plot(calculate_expectation(sim_single_1111, "durability"))
plot(calculate_expectation(sim_single_4321, "durability"))

prod1111 = getindex.(sim_single_1111.dur_his,1)
dur1111 = getindex.(sim_single_1111.dur_his,2)

map(x -> mean(dur1111[prod1111 .== x]), 1:4)

prod4321 = getindex.(sim_single_4321.dur_his,1)
dur4321 = getindex.(sim_single_4321.dur_his,2)

map(x -> mean(dur4321[prod4321 .== x]), 1:4)

(calculate_total_surplus(sim_single_1111,"consumer"),calculate_total_surplus(sim_single_4321,"consumer"))

plot(sim_single_4321.buyers[1].surplus_history[1:10])

StatsPlots.density(sum.(getfield.(sim_single_1111.buyers, :surplus_history)))
StatsPlots.density!(sum.(getfield.(sim_single_4321.buyers, :surplus_history)))

surplus_1111 = [x[y] for (x,y) in zip((getfield.(sim_single_1111.ut_his, :wtp) .- getfield.(sim_single_1111.ut_his, :p)),argmax.(getfield.(sim_single_1111.ut_his, :u)))]

surplus_4321 = [x[y] for (x,y) in zip((getfield.(sim_single_4321.ut_his, :wtp) .- getfield.(sim_single_4321.ut_his, :p)),argmax.(getfield.(sim_single_4321.ut_his, :u)))]

chp1111 = argmax.(getfield.(sim_single_1111.ut_his, :u))
wtp1111 = [x[y] for (x,y) in zip(getfield.(sim_single_1111.ut_his, :wtp) .- getfield.(sim_single_1111.ut_his, :p), chp1111)]

chp4321 = argmax.(getfield.(sim_single_4321.ut_his, :u))
wtp4321 = [x[y] for (x,y) in zip(getfield.(sim_single_4321.ut_his, :wtp) .- getfield.(sim_single_4321.ut_his, :p), chp4321)]

mean(wtp1111[chp1111 .== 1])
mean(wtp4321[chp4321 .== 1])

mean(wtp1111[chp1111 .== 2])
mean(wtp4321[chp4321 .== 2])

sum(wtp1111[chp1111 .== 3])
sum(wtp4321[chp4321 .== 3])

sum(wtp1111[chp1111 .== 4])
sum(wtp4321[chp4321 .== 4])

countmap(chp1111)
countmap(chp4321)

plot(calculate_expectation(sim_single_1111, "durability"))
plot(calculate_expectation(sim_single_4321, "durability"))

18528 *  / 23259
3.03 / 2.41

sort(countmap([argmax(x.u[x.ap]) for x in sim_single_1111.ut_his]))
sort(countmap([argmax(x.u[x.ap]) for x in sim_single_4321.ut_his]))

sim_single_4321.ut_his[1].wtp .- sim_single_4321.ut_his[1].p

(0.5 + sim_single_4321.ut_his[1].qs) .* sim_single_4321.ut_his[1].qe .* sim_single_4321.ut_his[1].de ./ sim_single_4321.ut_his[1].p

mean(surplus_1111)
mean(surplus_4321)

sim_single_1111.ut_his
sim_single_4321.ut_his

plot(argmax.(sim_single_1111.buyers[4].unit_bought_history))

plot(getindex.(sim_single_1111.buyers[4].durability_expectation_history, 3))


sim_single.buyers[4].quality_expectation_history
