# Experiment 4

    include(pwd() * "\\methods\\methods.jl")

    # Impact of secondary market on model dynamics.

    ex4_v1_total_surplus = []
    ex4_v1_producer_surplus = []
    ex4_v1_consumer_surplus = []
    ex4_v1_price = []
    ex4_v1_quantity_produced = []
    ex4_v1_quantity_sold = []
    ex4_v1_quantity_leased = []
    ex4_v1_reselling = []
    ex4_v1_producer_surplus_singleton = []
    ex4_v1_buying_history = []

    ex4_v2_total_surplus = []
    ex4_v2_producer_surplus = []
    ex4_v2_consumer_surplus = []
    ex4_v2_price = []
    ex4_v2_quantity_produced = []
    ex4_v2_quantity_sold = []
    ex4_v2_quantity_leased = []
    ex4_v2_reselling = []
    ex4_v2_producer_surplus_singleton = []
    ex4_v2_buying_history = []

    ex4_v3_total_surplus = []
    ex4_v3_producer_surplus = []
    ex4_v3_consumer_surplus = []
    ex4_v3_price = []
    ex4_v3_quantity_produced = []
    ex4_v3_quantity_sold = []
    ex4_v3_quantity_leased = []
    ex4_v3_reselling = []
    ex4_v3_producer_surplus_singleton = []
    ex4_v3_buying_history = []

    for i in 1:500

        if (mod(i,10) == 0) | (i == 1)
            println(i)
        end

        ex4_v1_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [true, true], [0.1, 0.1], 5, false)
        push!(ex4_v1_total_surplus, calculate_surplus(ex4_v1_sim, "total", true))
        push!(ex4_v1_producer_surplus, calculate_surplus(ex4_v1_sim, "producer", true))
        push!(ex4_v1_consumer_surplus, calculate_surplus(ex4_v1_sim, "consumer,total",true))
        push!(ex4_v1_price, calculate_price_history.(ex4_v1_sim.sellers; product_life = 5))
        push!(ex4_v1_quantity_produced, getfield.(ex4_v1_sim.sellers, :quantity_produced_history))
        push!(ex4_v1_quantity_sold, getfield.(ex4_v1_sim.sellers, :quantity_sold_history))
        push!(ex4_v1_producer_surplus_singleton, calculate_profit_history.(ex4_v1_sim.sellers))
        push!(ex4_v1_reselling, getfield.(ex4_v1_sim.sellers, :reselling_history))
        push!(ex4_v1_buying_history, getfield.(ex4_v1_sim.buyers, :unit_buying_selling_history))

        ex4_v2_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [false, true], [0, 0.1], 5, false)

        push!(ex4_v2_total_surplus, calculate_surplus(ex4_v2_sim, "total", true))
        push!(ex4_v2_producer_surplus, calculate_surplus(ex4_v2_sim, "producer", true))
        push!(ex4_v2_consumer_surplus, calculate_surplus(ex4_v2_sim, "consumer,total",true))
        push!(ex4_v2_price, calculate_price_history.(ex4_v2_sim.sellers; product_life = 5))
        push!(ex4_v2_quantity_produced, getfield.(ex4_v2_sim.sellers, :quantity_produced_history))
        push!(ex4_v2_quantity_sold, getfield.(ex4_v2_sim.sellers, :quantity_sold_history))
        push!(ex4_v2_producer_surplus_singleton, calculate_profit_history.(ex4_v2_sim.sellers))
        push!(ex4_v2_reselling, getfield.(ex4_v2_sim.sellers, :reselling_history))
        push!(ex4_v2_buying_history, getfield.(ex4_v2_sim.buyers, :unit_buying_selling_history))

        ex4_v3_sim = TO_GO(200, 2, 400, 500, [0.4, 0.4], [1.1, 1.1], "random", 0.25, 0.25, "stochastic", 1.1, [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[0.8, 2.], [0.8, 2.]], 0.1, true, 1, [0.7, 1.0], "softmax", [false, false], [0., 0.], 5, false)
        push!(ex4_v3_total_surplus, calculate_surplus(ex4_v3_sim, "total", true))
        push!(ex4_v3_producer_surplus, calculate_surplus(ex4_v3_sim, "producer", true))
        push!(ex4_v3_consumer_surplus, calculate_surplus(ex4_v3_sim, "consumer,total",true))
        push!(ex4_v3_price, calculate_price_history.(ex4_v3_sim.sellers; product_life = 5))
        push!(ex4_v3_quantity_produced, getfield.(ex4_v3_sim.sellers, :quantity_produced_history))
        push!(ex4_v3_quantity_sold, getfield.(ex4_v3_sim.sellers, :quantity_sold_history))
        push!(ex4_v3_producer_surplus_singleton, calculate_profit_history.(ex4_v3_sim.sellers))
        push!(ex4_v3_reselling, getfield.(ex4_v3_sim.sellers, :reselling_history))
        push!(ex4_v3_buying_history, getfield.(ex4_v3_sim.buyers, :unit_buying_selling_history))

    end

# Istnienie rynku wtórnego podnosi dobrobyt

function trim_outliers(x, p = 1)
    lower_bound = percentile(x, p)
    upper_bound = percentile(x,100-p)
    y = x[(x .>= lower_bound) .& (x .<= upper_bound)]
    return y
end

ex4_p1 = plot_ecdf(true, trim_outliers(ex4_v1_total_surplus,1), "Obaj gracze wykonują badania", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, trim_outliers(ex4_v2_total_surplus), "Jeden gracz wykonuje badania")
plot_ecdf(false, trim_outliers(ex4_v3_total_surplus), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p1, pwd() * "\\plots\\ex2\\total surplus research.svg")

# ... oraz nadwyżkę konsumenta ...

ex4_p2 = plot_ecdf(true, ex4_v1_consumer_surplus, "Obaj gracze wykonują badania", xlabel = "Nadwyżka konsumenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka konsumenta")
plot_ecdf(false, ex4_v2_consumer_surplus, "Jeden gracz wykonuje badania")
plot_ecdf(false, ex4_v3_consumer_surplus, "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p2, pwd() * "\\plots\\ex2\\consumer surplus research.svg")

# ... a także nadwyżkę producenta ...

ex4_p3 = plot_ecdf(true, ex4_v1_producer_surplus, "Obaj gracze wykonują badania", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, ex4_v2_producer_surplus, "Jeden gracz wykonuje badania")
plot_ecdf(false, ex4_v3_producer_surplus, "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p3, pwd() * "\\plots\\ex2\\producer surplus research.svg")

# Efekt istnienia rynku wtórnego nie jest istotny w przypadku rynku dóbr homogenicznych, ale jest istotny dla rynku, na którym występuje pionowe zróżnicowanie

ex4_p4 = plot_ecdf(true, mean.(getindex.(ex4_v1_producer_surplus_singleton,1)), "Obaj producenci badają rynek", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v1_producer_surplus_singleton,2)), "")

Plots.savefig(ex4_p4, pwd() * "\\plots\\ex2\\producer surplus research both.svg")

ex4_p5 = plot_ecdf(true, mean.(getindex.(ex4_v2_producer_surplus_singleton,1)), "Producent nie bada rynku", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v2_producer_surplus_singleton,2)), "Producent bada rynek")

Plots.savefig(ex4_p5, pwd() * "\\plots\\ex2\\producer surplus research one.svg")

ex4_p6 = plot_ecdf(true, mean.(getindex.(ex4_v3_producer_surplus_singleton,1)), "Żaden producent nie bada rynku", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(getindex.(ex4_v3_producer_surplus_singleton,2)), "",)

Plots.savefig(ex4_p6, pwd() * "\\plots\\ex2\\producer surplus research none.svg")

ex4_p7 = plot_ecdf(true, sum.(sum.(ex4_v1_quantity_sold)), "Obaj gracze wykonują badania", xlabel = "Liczba sprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba sprzedanych sztuk")
plot_ecdf(false,  sum.(sum.(ex4_v2_quantity_sold)), "Jeden gracz wykonuje badania")
plot_ecdf(false,  sum.(sum.(ex4_v3_quantity_sold)), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p7, pwd() * "\\plots\\ex2\\quantity sold leased.svg")

ex4_p8 = plot_ecdf(true, sum.(sum.(ex4_v1_quantity_produced)) .- sum.(sum.(ex4_v1_quantity_sold)), "Obaj gracze wykonują badania", xlabel = "Liczba niesprzedanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba niesprzedanych sztuk")
plot_ecdf(false,  sum.(sum.(ex4_v2_quantity_produced)) .- sum.(sum.(ex4_v2_quantity_sold)), "Jeden gracz wykonuje badania")
plot_ecdf(false,  sum.(sum.(ex4_v3_quantity_produced)) .- sum.(sum.(ex4_v3_quantity_sold)), "Żaden gracz nie wykonuje badań")

Plots.savefig(ex4_p8, pwd() * "\\plots\\ex2\\quantity unsold.svg")

sum(sum.(sum.(ex4_v1_quantity_leased)) .+ sum.(sum.(ex4_v1_quantity_sold))) / sum( sum.(sum.(ex4_v1_quantity_produced)))

getindex.(ex4_v1_quality)

sum(sum.(sum.(ex4_v2_quantity_leased)) .+ sum.(sum.(ex4_v2_quantity_sold))) / sum( sum.(sum.(ex4_v2_quantity_produced)))

sum(sum.(sum.(ex4_v3_quantity_leased)) .+ sum.(sum.(ex4_v3_quantity_sold))) / sum( sum.(sum.(ex4_v3_quantity_produced)))

Plots.savefig(ex4_p8, pwd() * "\\plots\\quantity unsold unleased.svg")

UnequalVarianceTTest(float.(mean.(getindex.(ex4_v1_producer_surplus_singleton,1))), float.(mean.(getindex.(ex4_v3_producer_surplus_singleton,1))))
UnequalVarianceTTest(float.(mean.(getindex.(ex4_v2_producer_surplus_singleton,1))), float.(mean.(getindex.(ex4_v4_producer_surplus_singleton,1))))

count_x(y;x) = sum(y .== x)

# Dlaczego istnienie rynku wtórnego podnosi nadwyżkę producenta?

ex4_p5 = StatsPlots.histogram(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v1_buying_history]), alpha = 0.25, normalize = true, xlabel = "Częstotliwość zakupów", ylabel = "f(x)", label = "Rynek wtórny istnieje", title = "Częstotliwość zakupów, identyczna jakość", bins = 0:5:70)
StatsPlots.histogram!(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v3_buying_history]), alpha = 0.25, normalize = true, label = "Rynek wtórny nie istnieje", bins = 0:5:70, ylim = (0, 0.08))

Plots.savefig(ex4_p5, pwd() * "\\plots\\ex3_primar freq equal.svg")


UnequalVarianceTTest(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v1_buying_history]), mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v3_buying_history]))

ex4_p6 = StatsPlots.histogram(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v2_buying_history]), alpha = 0.25, normalize = true, xlabel = "Częstotliwość zakupów", ylabel = "f(x)", label = "Rynek wtórny istnieje", title = "Częstotliwość zakupów, różna jakość", bins = 0:5:70)
StatsPlots.histogram!(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v4_buying_history]), alpha = 0.25, normalize = true, label = "Rynek wtórny nie istnieje", bins=0:5:70, ylim = (0, 0.08))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v2_buying_history]))
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v4_buying_history]))
UnequalVarianceTTest(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v2_buying_history]), mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v4_buying_history]))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v1_buying_history])) / 
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v3_buying_history]))

mean(mean.([[count_x(string.(getfield.(y,:d)); x = "b") for y in x] for x in ex4_v2_buying_history])) / 
mean(mean.([[count_x(string.(getfield.(y,:d)); x="b") for y in x] for x in ex4_v4_buying_history]))

Plots.savefig(ex4_p6, pwd() * "\\plots\\ex3_primar freq diff.svg")

# Kto najbardziej zyskuje na rynku wtórnym?

mean(mean(getindex.(ex4_v1_producer_surplus_singleton, 1))) - mean(mean(getindex.(ex4_v3_producer_surplus_singleton, 1)))
UnequalVarianceTTest(mean.(getindex.(ex4_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex4_v3_producer_surplus_singleton, 1)))

mean(mean(getindex.(ex4_v2_producer_surplus_singleton, 1))) - mean(mean(getindex.(ex4_v4_producer_surplus_singleton, 1)))
UnequalVarianceTTest(mean.(getindex.(ex4_v2_producer_surplus_singleton, 1)), mean.(getindex.(ex4_v4_producer_surplus_singleton, 1)))

mean(mean(getindex.(ex4_v1_producer_surplus_singleton, 2))) - mean(mean(getindex.(ex4_v3_producer_surplus_singleton, 2)))
UnequalVarianceTTest(mean.(getindex.(ex4_v1_producer_surplus_singleton, 2)), mean.(getindex.(ex4_v3_producer_surplus_singleton, 2)))

mean(mean(getindex.(ex4_v2_producer_surplus_singleton, 2))) - mean(mean(getindex.(ex4_v4_producer_surplus_singleton, 2)))
UnequalVarianceTTest(mean.(getindex.(ex4_v2_producer_surplus_singleton, 2)), mean.(getindex.(ex4_v4_producer_surplus_singleton, 2)))

UnequalVarianceTTest(sort(mean.(getindex.(ex4_v2_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex4_v1_producer_surplus_singleton, 1))), sort(mean.(getindex.(ex4_v4_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex4_v3_producer_surplus_singleton, 1))))

UnequalVarianceTTest(sort(mean.(getindex.(ex4_v2_producer_surplus_singleton, 1))) .- sort(mean.(getindex.(ex4_v1_producer_surplus_singleton, 1))), sort(mean.(getindex.(ex4_v2_producer_surplus_singleton, 2))) .- sort(mean.(getindex.(ex4_v4_producer_surplus_singleton, 2))))

plot_ecdf(true, mean.(getindex.(ex4_v3_producer_surplus_singleton, 2)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")
plot_ecdf(false, mean.(getindex.(ex4_v4_producer_surplus_singleton, 2)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Całkowita nadwyżka", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka całkowita")


Plots.plot(mean(getindex.(ex4_v3_producer_surplus_singleton, 2)))
# Bo podnosi liczbę produkowanych sztuk 

ex4_p5 = plot_ecdf(true, mean.(mean.(ex4_v1_quantity_produced)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Liczba produkowanych sztuk", ylabel = "F(x)", title = "Dystrybuanta empiryczna - liczba produkowanych sztuk")
plot_ecdf(false, mean.(mean.(ex4_v2_quantity_produced)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex4_v3_quantity_produced)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex4_v4_quantity_produced)), "Różna jakość, rynek wtórny nie istnieje")

Plots.savefig(ex4_p4, pwd() * "\\plots\\ex3_secmar prod quantity.svg")

# Produkcja jest większa, gdyż i sprzedaż jest większa -> rynek wtórny i odsprzedaż zwiększają frequency zakupów

ex4_p5 = plot_ecdf(true, mean.(mean.(ex4_v1_quantity_sold)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(ex4_v2_quantity_sold)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex4_v3_quantity_sold)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex4_v4_quantity_sold)), "Różna jakość, rynek wtórny nie istnieje")

ex4_p4 = plot_ecdf(true, mean.(mean.(ex4_v1_quantity_leased)), "Identyczna jakość, rynek wtórny istnieje", xlabel = "Nadwyżka producenta", ylabel = "F(x)", title = "Dystrybuanta empiryczna - nadwyżka producenta")
plot_ecdf(false, mean.(mean.(ex4_v2_quantity_leased)), "Różna jakość, rynek wtórny istnieje")
plot_ecdf(false, mean.(mean.(ex4_v3_quantity_leased)), "Identyczna jakość, rynek wtórny nie istnieje")
plot_ecdf(false, mean.(mean.(ex4_v4_quantity_leased)), "Różna jakość, rynek wtórny nie istnieje")

Plots.scatter(sum.(sum.(ex4_v1_quantity_produced)), sum.(sum.(ex4_v1_reselling)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, xlabel = "Liczba wyprodukowanych dóbr", ylabel = "Liczba odsprzedanych dóbr", label = "Identyczna jakość", legend = :bottomright)
Plots.scatter!(sum.(sum.(ex4_v2_quantity_produced)), sum.(sum.(ex4_v2_reselling)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość")

Plots.scatter(sum.(getindex.(ex4_v1_quantity_produced, 1)), sum.(getindex.(ex4_v1_reselling, 1)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Identyczna jakość", ylabel = "Liczba odsprzedanych dóbr", xlabel = "Liczba wyprodukowanych dóbr")
Plots.scatter!(sum.(getindex.(ex4_v2_quantity_produced, 1)), sum.(getindex.(ex4_v2_reselling, 1)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, lepszy produkt")
Plots.scatter!(sum.(getindex.(ex4_v2_quantity_produced, 2)), sum.(getindex.(ex4_v2_reselling, 2)), smooth = true, markerstrokewidth = 0, markeralpha = 0.25, linewidth = 3, label = "Różna jakość, gorszy produkt")

GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex4_v1_quantity_produced, 1)), y=sum.(getindex.(ex4_v1_reselling, 1)) ))
GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex4_v2_quantity_produced, 1)), y=sum.(getindex.(ex4_v2_reselling, 1)) ))
GLM.lm(@formula(y~x), DataFrame(x= sum.(getindex.(ex4_v2_quantity_produced, 2)), y=sum.(getindex.(ex4_v2_reselling, 2)) ))

