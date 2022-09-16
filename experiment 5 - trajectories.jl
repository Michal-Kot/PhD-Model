# Experiment 4

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex5_v1_total_surplus = []
ex5_v1_producer_surplus = []
ex5_v1_consumer_surplus = []
ex5_v1_price = []
ex5_v1_quality = []
ex5_v1_durability = []
ex5_v1_margin = []
ex5_v1_quantity_produced = []
ex5_v1_quantity_sold = []
ex5_v1_producer_surplus_singleton = []

ex5_v2_total_surplus = []
ex5_v2_producer_surplus = []
ex5_v2_consumer_surplus = []
ex5_v2_price = []
ex5_v2_quality = []
ex5_v2_durability = []
ex5_v2_margin = []
ex5_v2_quantity_produced = []
ex5_v2_quantity_sold = []
ex5_v2_producer_surplus_singleton = []

ex5_v3_total_surplus = []
ex5_v3_producer_surplus = []
ex5_v3_consumer_surplus = []
ex5_v3_price = []
ex5_v3_quality = []
ex5_v3_durability = []
ex5_v3_margin = []
ex5_v3_quantity_produced = []
ex5_v3_quantity_sold = []
ex5_v3_producer_surplus_singleton = []

ex5_v4_total_surplus = []
ex5_v4_producer_surplus = []
ex5_v4_consumer_surplus = []
ex5_v4_price = []
ex5_v4_quality = []
ex5_v4_durability = []
ex5_v4_margin = []
ex5_v4_quantity_produced = []
ex5_v4_quantity_sold = []
ex5_v4_producer_surplus_singleton = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex5_v1_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.7, 0.7], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [0.1, 0.1], [[0.75, 1.25], [0.75, 1.25]], [[0.3, 0.5], [0.3, 0.5]], [[1., 2.], [1., 2.]], 0.10, true, 0)
    push!(ex5_v1_total_surplus, calculate_surplus(ex5_v1_sim, "total", false))
    push!(ex5_v1_producer_surplus, calculate_surplus(ex5_v1_sim, "producer", false))
    push!(ex5_v1_consumer_surplus, calculate_surplus(ex5_v1_sim, "consumer,total", false))
    push!(ex5_v1_price, calculate_price_history.(ex5_v1_sim.sellers))
    push!(ex5_v1_quality, getfield.(ex5_v1_sim.sellers, :quality_history))
    push!(ex5_v1_durability, getfield.(ex5_v1_sim.sellers, :durability_history))
    push!(ex5_v1_margin, getfield.(ex5_v1_sim.sellers, :margin_history))
    push!(ex5_v1_quantity_produced, getfield.(ex5_v1_sim.sellers, :quantity_produced_history))
    push!(ex5_v1_quantity_sold, getfield.(ex5_v1_sim.sellers, :quantity_sold_history))
    push!(ex5_v1_producer_surplus_singleton, calculate_profit_history.(ex5_v1_sim.sellers))

    ex5_v2_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.7, 0.7], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [0.1, 0.1], [[1.0, 1.5], [0.5, 1.0]], [[0.4, 0.6], [0.2, 0.4]], [[1., 2.], [1., 2.]], 0.10, true, 0)
    push!(ex5_v2_total_surplus, calculate_surplus(ex5_v2_sim, "total", false))
    push!(ex5_v2_producer_surplus, calculate_surplus(ex5_v2_sim, "producer", false))
    push!(ex5_v2_consumer_surplus, calculate_surplus(ex5_v2_sim, "consumer,total",false))
    push!(ex5_v2_price, calculate_price_history.(ex5_v2_sim.sellers))
    push!(ex5_v2_quality, getfield.(ex5_v2_sim.sellers, :quality_history))
    push!(ex5_v2_durability, getfield.(ex5_v2_sim.sellers, :durability_history))
    push!(ex5_v2_margin, getfield.(ex5_v2_sim.sellers, :margin_history))
    push!(ex5_v2_quantity_produced, getfield.(ex5_v2_sim.sellers, :quantity_produced_history))
    push!(ex5_v2_quantity_sold, getfield.(ex5_v2_sim.sellers, :quantity_sold_history))
    push!(ex5_v2_producer_surplus_singleton, calculate_profit_history.(ex5_v2_sim.sellers))

    ex5_v3_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.7, 0.7], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [0.1, 0.1],  [[0.75, 1.25], [0.75, 1.25]], [[0.3, 0.5], [0.3, 0.5]], [[1., 2.], [1., 2.]], 0.10, false, 0)
    push!(ex5_v3_total_surplus, calculate_surplus(ex5_v3_sim, "total", false))
    push!(ex5_v3_producer_surplus, calculate_surplus(ex5_v3_sim, "producer", false))
    push!(ex5_v3_consumer_surplus, calculate_surplus(ex5_v3_sim, "consumer,total",false))
    push!(ex5_v3_price, calculate_price_history.(ex5_v3_sim.sellers))
    push!(ex5_v3_quality, getfield.(ex5_v3_sim.sellers, :quality_history))
    push!(ex5_v3_durability, getfield.(ex5_v3_sim.sellers, :durability_history))
    push!(ex5_v3_margin, getfield.(ex5_v3_sim.sellers, :margin_history))
    push!(ex5_v3_quantity_produced, getfield.(ex5_v3_sim.sellers, :quantity_produced_history))
    push!(ex5_v3_quantity_sold, getfield.(ex5_v3_sim.sellers, :quantity_sold_history))
    push!(ex5_v3_producer_surplus_singleton, calculate_profit_history.(ex5_v3_sim.sellers))


    ex5_v4_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.7, 0.7], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [0.1, 0.1], [[1.0, 1.5], [0.5, 1.0]], [[0.4, 0.6], [0.2, 0.4]], [[1., 2.], [1., 2.]], 0.10, false, 0)
    push!(ex5_v4_total_surplus, calculate_surplus(ex5_v4_sim, "total", false))
    push!(ex5_v4_producer_surplus, calculate_surplus(ex5_v4_sim, "producer", false))
    push!(ex5_v4_consumer_surplus, calculate_surplus(ex5_v4_sim, "consumer,total",false))
    push!(ex5_v4_price, calculate_price_history.(ex5_v4_sim.sellers))
    push!(ex5_v4_quality, getfield.(ex5_v4_sim.sellers, :quality_history))
    push!(ex5_v4_durability, getfield.(ex5_v4_sim.sellers, :durability_history))
    push!(ex5_v4_margin, getfield.(ex5_v4_sim.sellers, :margin_history))
    push!(ex5_v4_quantity_produced, getfield.(ex5_v4_sim.sellers, :quantity_produced_history))
    push!(ex5_v4_quantity_sold, getfield.(ex5_v4_sim.sellers, :quantity_sold_history))
    push!(ex5_v4_producer_surplus_singleton, calculate_profit_history.(ex5_v4_sim.sellers))


end

ex5_data = [ex5_v1_total_surplus,
ex5_v1_producer_surplus,
ex5_v1_consumer_surplus,
ex5_v1_price,
ex5_v1_quality,
ex5_v1_durability,
ex5_v1_margin,
ex5_v1_quantity_produced,
ex5_v1_quantity_sold,
ex5_v1_producer_surplus_singleton,

ex5_v2_total_surplus,
ex5_v2_producer_surplus,
ex5_v2_consumer_surplus,
ex5_v2_price,
ex5_v2_quality,
ex5_v2_durability,
ex5_v2_margin,
ex5_v2_quantity_produced,
ex5_v2_quantity_sold,
ex5_v2_producer_surplus_singleton,

ex5_v3_total_surplus,
ex5_v3_producer_surplus,
ex5_v3_consumer_surplus,
ex5_v3_price,
ex5_v3_quality,
ex5_v3_durability,
ex5_v3_margin,
ex5_v3_quantity_produced,
ex5_v3_quantity_sold,
ex5_v3_producer_surplus_singleton,

ex5_v4_total_surplus,
ex5_v4_producer_surplus,
ex5_v4_consumer_surplus,
ex5_v4_price,
ex5_v4_quality,
ex5_v4_durability,
ex5_v4_margin,
ex5_v4_quantity_produced,
ex5_v4_quantity_sold,
ex5_v4_producer_surplus_singleton]

serialize("C:\\Users\\User\\Documents\\Model-PhD-LocalData\\ex5_data", ex5_data)

ex5_data = deserialize("C:\\Users\\User\\Documents\\Model-PhD-LocalData\\ex5_data")

ex5_v1_quality = ex5_data[5]
ex5_v1_durability = ex5_data[6]
ex5_v1_margin = ex5_data[7]

ex5_v2_quality = ex5_data[15]
ex5_v2_durability = ex5_data[16]
ex5_v2_margin = ex5_data[17]

ex5_v3_quality = ex5_data[25]
ex5_v3_durability = ex5_data[26]
ex5_v3_margin = ex5_data[27]

ex5_v4_quality = ex5_data[35]
ex5_v4_durability = ex5_data[36]
ex5_v4_margin = ex5_data[37]

euclidean_distance(x1,y1,x2,y2) = sqrt.((x1 .- x2).^2 .+ (y1 .- y2).^2)

ex5_v1_dist = mean.([euclidean_distance(d1d2[1],k1k2[1],d1d2[2],k1k2[2]) for (d1d2, k1k2) in zip(ex5_v1_durability,ex5_v1_quality)])
ex5_v2_dist = mean.([euclidean_distance(d1d2[1],k1k2[1],d1d2[2],k1k2[2]) for (d1d2, k1k2) in zip(ex5_v2_durability,ex5_v2_quality)])
ex5_v3_dist = mean.([euclidean_distance(d1d2[1],k1k2[1],d1d2[2],k1k2[2]) for (d1d2, k1k2) in zip(ex5_v3_durability,ex5_v3_quality)])
ex5_v4_dist = mean.([euclidean_distance(d1d2[1],k1k2[1],d1d2[2],k1k2[2]) for (d1d2, k1k2) in zip(ex5_v4_durability,ex5_v4_quality)])

using CategoricalArrays

ex5_v1_dist_cat = cut(ex5_v1_dist, 5)
ex5_v2_dist_cat = cut(ex5_v2_dist, 4)
ex5_v3_dist_cat = cut(ex5_v3_dist, 5)
ex5_v4_dist_cat = cut(ex5_v4_dist, 5)

countmap(ex5_v1_dist_cat)

ex5_v1_margin_diff = mean.(getindex.(ex5_v1_margin,1)) .- mean.(getindex.(ex5_v1_margin,2))
ex5_v2_margin_diff = mean.(getindex.(ex5_v2_margin,1)) .- mean.(getindex.(ex5_v2_margin,2))
ex5_v3_margin_diff = mean.(getindex.(ex5_v3_margin,1)) .- mean.(getindex.(ex5_v3_margin,2))
ex5_v4_margin_diff = mean.(getindex.(ex5_v4_margin,1)) .- mean.(getindex.(ex5_v4_margin,2))

ex5_v1_kd_diff = mean.(getindex.(ex5_v1_quality, 1)) ./ (1 .- mean.(getindex.(ex5_v1_durability, 1))) .- mean.(getindex.(ex5_v1_quality, 2)) ./ (1 .- mean.(getindex.(ex5_v1_durability, 2)))
ex5_v2_kd_diff = mean.(getindex.(ex5_v2_quality, 1)) ./ (1 .- mean.(getindex.(ex5_v2_durability, 1))) .- mean.(getindex.(ex5_v2_quality, 2)) ./ (1 .- mean.(getindex.(ex5_v2_durability, 2)))

scatter(ex5_v1_kd_diff ./ ex5_v1_margin_diff, mean.(getindex.(ex5_v1_producer_surplus_singleton,1)))
scatter!(ex5_v2_kd_diff ./ ex5_v2_margin_diff, mean.(getindex.(ex5_v2_producer_surplus_singleton,1)) .- mean.(getindex.(ex5_v2_producer_surplus_singleton,2)))

using Contour
using PlotlyJS

x = cut(ex5_v2_kd_diff,10)
y = cut(ex5_v2_margin_diff,10)
zz = [mean(mean.(getindex.(ex5_v2_producer_surplus_singleton,1))[(x .== xx) .& (y .== yy)]) for xx in sort(unique(x)), yy in sort(unique(y))]

PlotlyJS.plot(PlotlyJS.contour(z=zz))

levels()

for cl in levels(c1)
    lvl = level(cl) # the z-value of this contour level
    for line in lines(cl)
        xs, ys = coordinates(line) # coordinates of this line segment
        plot(xs, ys, color=lvl) # pseuod-code; use whatever plotting package you prefer
    end
end

ex5_v1_dist_cat = cut(ex5_v1_kd_diff, 5)

p = plot()
[plot_smoothingspline(ex5_v1_margin_diff[ex5_v1_dist_cat .== x], mean.(getindex.(ex5_v1_producer_surplus_singleton,1))[ex5_v1_dist_cat .== x]) for x in sort(unique(ex5_v1_dist_cat))]
p

p = plot()
[plot_smoothingspline(ex5_v2_margin_diff[ex5_v2_dist_cat .== x], mean.(getindex.(ex5_v2_producer_surplus_singleton,1))[ex5_v2_dist_cat .== x]) for x in sort(unique(ex5_v2_dist_cat))]
p

p = plot()
[plot_smoothingspline(ex5_v3_margin_diff[ex5_v3_dist_cat .== x], mean.(getindex.(ex5_v3_producer_surplus_singleton,1))[ex5_v3_dist_cat .== x]) for x in sort(unique(ex5_v3_dist_cat))]
p

p = plot()
[plot_smoothingspline(ex5_v4_margin_diff[ex5_v4_dist_cat .== x], mean.(getindex.(ex5_v4_producer_surplus_singleton,1))[ex5_v4_dist_cat .== x]) for x in sort(unique(ex5_v4_dist_cat))]
p

using SmoothingSplines

function plot_smoothingspline(X,Y)
    spl = fit(SmoothingSpline, X, Y, 0.25) # λ=250.0
    Ypred = predict(spl) # fitted vector
    plot!(sort(X), Ypred[sortperm(X)])
end

X = mean.(getindex.(ex5_v1_margin,1))[ex5_v1_dist_cat .== ex5_v1_dist_cat[5]]
Y = mean.(getindex.(ex5_v1_producer_surplus_singleton,1))[ex5_v1_dist_cat .== ex5_v1_dist_cat[5]]



scatter(X,Y)


ex5_p1 = plot_ecdf(mean.(ex5_v1_total_surplus), "Identyczna jakość, rynek wtórny istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(mean.(ex5_v2_total_surplus), "Różna jakość, rynek wtórny istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(mean.(ex5_v3_total_surplus), "Identyczna jakość, rynek wtórny nie istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(mean.(ex5_v4_total_surplus), "Różna jakość, rynek wtórny nie istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)


ex5_p2 = plot_ecdf(ex5_v1_consumer_surplus, "Identyczna jakość, rynek wtórny istnieje", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex5_v2_consumer_surplus, "Różna jakość, rynek wtórny istnieje", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex5_v3_consumer_surplus, "Identyczna jakość, rynek wtórny nie istnieje", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex5_v4_consumer_surplus, "Różna jakość, rynek wtórny nie istnieje", "Consumer Surplus", "Probability", "ECDF", false)

ex5_p3 = plot_ecdf(ex5_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex5_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex5_v3_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex5_v4_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex5_p3 = plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v1_producer_surplus_singleton, 1)]), 2.5), "Identyczna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", true)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v2_producer_surplus_singleton, 1)]),2.5), "Różna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v3_producer_surplus_singleton, 1)]),2.5), "Identyczna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v4_producer_surplus_singleton, 1)]),2.5), "Różna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)

ex5_p3 = plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v1_producer_surplus_singleton, 1)]), 2.5), "Identyczna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", true)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v2_producer_surplus_singleton, 2)]),2.5), "Różna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v3_producer_surplus_singleton, 2)]),2.5), "Identyczna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex5_v4_producer_surplus_singleton, 2)]),2.5), "Różna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)

sec_mar_better = mean.([x[51:150] for x in getindex.(ex5_v2_producer_surplus_singleton, 1)]) .- mean.([x[51:150] for x in getindex.(ex5_v4_producer_surplus_singleton, 1)])
sec_mar_worse = mean.([x[51:150] for x in getindex.(ex5_v2_producer_surplus_singleton, 2)]) .- mean.([x[51:150] for x in getindex.(ex5_v4_producer_surplus_singleton, 2)])

UnequalVarianceTTest(sec_mar_better, sec_mar_worse)

ex5_p3 = plot_ecdf(sum.(sum.(ex5_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex5_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex5_p3 = plot_ecdf(sum.(sum.(ex5_v1_quantity_sold)) / 200, "Identyczna jakość, rynek wtórny istnieje", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex5_v2_quantity_sold)) / 200, "Różna jakość, rynek wtórny istnieje", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(sum.(sum.(ex5_v3_quantity_sold)) / 200, "Identyczna jakość, rynek wtórny nie istnieje", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(sum.(sum.(ex5_v4_quantity_sold)) / 200, "Różna jakość, rynek wtórny nie istnieje", "Producer Surplus", "Probability", "ECDF", false)





UnequalVarianceTTest(mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)))

UnequalVarianceTTest(mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex5_v3_producer_surplus_singleton, 1)))

UnequalVarianceTTest(mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)), mean.(getindex.(ex5_v4_producer_surplus_singleton, 1)))

mean.(getindex.(ex5_v2_producer_surplus_singleton, 1))
mean.(getindex.(ex5_v3_producer_surplus_singleton, 1))
mean.(getindex.(ex5_v4_producer_surplus_singleton, 1))
mean.(getindex.(ex5_v4_producer_surplus_singleton, 1))

mean.(mean.(ex5_v1_price))

ex5_p3 = plot_ecdf(mean.(mean.(ex5_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex5_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

#########
#%%



scatter(mean.(getindex.(ex5_v1_margin,1)), mean.(getindex.(ex5_v1_producer_surplus_singleton, 1)))
scatter!(mean.(getindex.(ex5_v2_margin,1)), mean.(getindex.(ex5_v2_producer_surplus_singleton, 1)))
scatter!(mean.(getindex.(ex5_v3_margin,1)), mean.(getindex.(ex5_v3_producer_surplus_singleton, 1)))
scatter!(mean.(getindex.(ex5_v4_margin,1)), mean.(getindex.(ex5_v4_producer_surplus_singleton, 1)))

plot(mean(getindex.(ex5_v1_durability,1)), label = "Identyczna jakość, rynek wtórny istnieje", legend = :right, xlabel = "Iteracje", ylabel = "Średnia trwałość produktów")
plot!(mean(getindex.(ex5_v2_durability,1)), label = "Różna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v3_durability,1)), label = "Identyczna jakość, rynek wtórny nie istnieje")
plot!(mean(getindex.(ex5_v4_durability,1)), label = "Różna jakość, rynek wtórny nie istnieje")
plot!(mean(getindex.(ex5_v2_durability,2)), label = "Różna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v4_durability,2)), label = "Różna jakość, rynek wtórny nie istnieje")

plot(mean(getindex.(ex5_v1_quality,1)), label = "Identyczna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v2_quality,1)), label = "Różna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v3_quality,1)), label = "Identyczna jakość, rynek wtórny nie istnieje")
plot!(mean(getindex.(ex5_v4_quality,1)), label = "Różna jakość, rynek wtórny nie istnieje")
plot!(mean(getindex.(ex5_v2_quality,2)), label = "Różna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v4_quality,2)), label = "Różna jakość, rynek wtórny nie istnieje")

plot(mean(getindex.(ex5_v1_margin,1)), label = "Identyczna jakość, rynek wtórny istnieje", legend = :bottomright)
plot!(mean(getindex.(ex5_v3_margin,1)), label = "Identyczna jakość, rynek wtórny nie istnieje")

plot!(mean(getindex.(ex5_v2_margin,1)), label = "Różna jakość, rynek wtórny istnieje", legend = :bottomright)
plot!(mean(getindex.(ex5_v4_margin,1)), label = "Różna jakość, rynek wtórny nie istnieje")

plot!(mean(getindex.(ex5_v2_margin,2)), label = "Różna jakość, rynek wtórny istnieje", legend = :bottomright)
plot!(mean(getindex.(ex5_v4_margin,2)), label = "Różna jakość, rynek wtórny nie istnieje")

plot(mean(getindex.(ex5_v1_quantity_produced,1)), label = "Identyczna jakość, rynek wtórny istnieje", legend = :topleft)
plot!(mean(getindex.(ex5_v1_quantity_sold,1)), label = "Identyczna jakość, rynek wtórny istnieje", legend = :topleft)
plot!(mean(getindex.(ex5_v2_margin,1)), label = "Różna jakość, rynek wtórny istnieje")
plot!(mean(getindex.(ex5_v3_margin,1)), label = "Identyczna jakość, rynek wtórny nie istnieje")
plot!(mean(getindex.(ex5_v4_margin,1)), label = "Różna jakość, rynek wtórny nie istnieje")