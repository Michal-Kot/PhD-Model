##################################### Unit testing ####################################################

include(pwd() * "\\methods\\methods.jl")

########################################### Functions 

function create_buying(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function create_destroy(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "destroyed"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function create_resell(_buyer, product_id)
    b = []
    for item in _buyer.unit_buying_selling_history
        if item.decision == "sell, secondary market"
            if item.product == product_id
                push!(b, item.t)
            end
        end
    end
    return b
end

function check_impact_of_signals(sim_single, buyer_id, product_id)

    _buyer = sim_single.buyers[buyer_id]

    own_signals = []
    outer_signals = []

    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == product_id
                push!(own_signals, _buyer.quality_expectation_history[item.t][product_id] - _buyer.quality_expectation_history[item.t - 1][product_id])
            end
        end
    end 

    for _buyer in sim_single.buyers[Not(buyer_id)]
        for item in _buyer.unit_buying_selling_history
            if item.decision == "buy, primary market"
                if item.product == product_id
                    push!(outer_signals, _buyer.quality_expectation_history[item.t][product_id] - _buyer.quality_expectation_history[item.t - 1][product_id])
                end
            end
        end 
    end

    return (own_signals, outer_signals)

end

function quality_sampled_products(buyers, product_id)

    quality_of_unit_possessed_history = [[] for _ in 1:lastindex(buyers[1].quality_expectation_history)]


    for _buyer in buyers
        for item in _buyer.unit_buying_selling_history
            if item.decision == "buy, primary market"
                if item.product == product_id
                    push!(quality_of_unit_possessed_history[item.t], _buyer.quality_of_unit_possessed_history[item.t])
                end
            end
        end
    end

    return quality_of_unit_possessed_history

end

############################## Individual level ###############################################

########################################### If signals change expectation

sim_single = TO_GO(25, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.75, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

## Number of signals

sim_single.buyers[2].unit_buying_selling_history
sim_single.buyers[1].signal_volume
countmap(getfield.(sim_single.buyers[2].unit_buying_selling_history, :decision))["buy, primary market"]

sim_single.buyers[1].unit_buying_selling_history
sim_single.buyers[2].signal_volume
countmap(getfield.(sim_single.buyers[1].unit_buying_selling_history, :decision))["buy, primary market"]

## Charts

Plots.plot(sim_single.sellers[1].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1))
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 1), seriestype = :vline)

Plots.plot(sim_single.sellers[2].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 2))
Plots.plot!(create_buying(sim_single.buyers[1], 2) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 2) .- 0.1, seriestype = :vline, 
linestyle = :dash)

########################################### Check if 0 weight of comms results in no difference in expectations

## EXP = 0, WOM = 1

sim_single = TO_GO(50, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

prod(x) = any(x .> 0) ? argmax(x) : 0

Plots.plot(prod.(sim_single.buyers[1].unit_possessed_history), color = :black)
Plots.plot!(create_buying(sim_single.buyers[1], 1), seriestype = :vline, color = :red)
Plots.plot!(create_buying(sim_single.buyers[1], 2), seriestype = :vline, color = :green)
Plots.plot!(create_destroy(sim_single.buyers[1], 1), seriestype = :vline, color = :red, linestyle = :dash)
Plots.plot!(create_destroy(sim_single.buyers[1], 2), seriestype = :vline, color = :green, linestyle = :dash)

[any(x .== 0) for x in check_impact_of_signals(sim_single, 1, 1)]

Plots.plot(sim_single.sellers[1].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1))
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 0.1, seriestype = :vline, linestyle = :dash)

Plots.plot(sim_single.sellers[2].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 2))
Plots.plot!(create_buying(sim_single.buyers[1], 2) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 2) .- 0.1, seriestype = :vline, 
linestyle = :dash)

[any(x .== 0) for x in check_impact_of_signals(sim_single, 1, 2)]

## EXP = 1, WOM = 0

sim_single = TO_GO(50, 2, 2, 1, [0.4, 0.4], [1.0, 1.0], "random", 1.0, 0.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

Plots.plot(prod.(sim_single.buyers[1].unit_possessed_history), color = :black)
Plots.plot!(create_buying(sim_single.buyers[1], 1), seriestype = :vline, color = :red)
Plots.plot!(create_buying(sim_single.buyers[1], 2), seriestype = :vline, color = :green)
Plots.plot!(create_destroy(sim_single.buyers[1], 1), seriestype = :vline, color = :red, linestyle = :dash)
Plots.plot!(create_destroy(sim_single.buyers[1], 2), seriestype = :vline, color = :green, linestyle = :dash)

[any(x .== 0) for x in check_impact_of_signals(sim_single, 1, 1)]

Plots.plot(sim_single.sellers[1].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 1))
Plots.plot!(create_buying(sim_single.buyers[1], 1) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 1) .- 0.1, seriestype = :vline, linestyle = :dash)

Plots.plot(sim_single.sellers[2].quality_history)
Plots.plot!(getindex.(sim_single.buyers[1].quality_expectation_history, 2))
Plots.plot!(create_buying(sim_single.buyers[1], 2) .+ 0.1, seriestype = :vline, 
linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[2], 2) .- 0.1, seriestype = :vline, 
linestyle = :dash)

########################################### Reselling probability changes due to secondary market

sim_single = TO_GO(50, 2, 200, 400, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

Plots.plot(getfield.(sim_single.buyers, :reselling_probability_history)[1], color = :black, linewidth = 2)
Plots.plot!(create_buying(sim_single.buyers[1], 1), linetype = :vline, color = "red", linestyle = :dash)
Plots.plot!(create_buying(sim_single.buyers[1], 2), linetype = :vline, color = "red", linestyle = :dash)
Plots.plot!(create_resell(sim_single.buyers[1], 1), linetype = :vline, color = "green", linestyle = :dash)
Plots.plot!(create_resell(sim_single.buyers[1], 2), linetype = :vline, color = "green", linestyle = :dash)

## No secondary market -> flat probability

sim_single = TO_GO(50, 2, 200, 400, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, false, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

Plots.plot(getfield.(sim_single.buyers, :reselling_probability_history)[1])
Plots.plot!(create_buying(sim_single.buyers[1], 1), linetype = :vline, color = "grey")
Plots.plot!(create_buying(sim_single.buyers[1], 2), linetype = :vline, color = "grey")
Plots.plot!(create_resell(sim_single.buyers[1], 1), linetype = :vline, color = "red")
Plots.plot!(create_resell(sim_single.buyers[1], 2), linetype = :vline, color = "red")

########################################### Max breakage time

sim_single = TO_GO(1000, 2, 1000, 2000, [0.4, 0.4], [1.0, 1.0], "random", 0.0, 1.0, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, false, 0, [0.7, 1.], "softmax", [false, true], [0, 0.5], 4, false)

history_of_trading = getfield.(sim_single.buyers, :unit_buying_selling_history)

max_h = 0

for b in history_of_trading
    for item in eachindex(b)
        if item > 1
            if b[item - 1].decision == "buy, primary market"
                if b[item].decision == "destroyed"
                    if b[item].t - b[item-1].t > max_h
                        max_h = b[item].t - b[item-1].t
                    end
                end
            end
        end
    end
end

max_h

###########################################

########################################### Meso level #############################################

########################################### No links
## -> no signals shall be received

sim_single = TO_GO(200, 2, 400, 0, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

maximum(getfield.(sim_single.buyers, :signal_volume))

########################################### No learning
## -> expectation shall be flat

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0., 0., "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

allequal([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)])
allequal([getindex.(x,2) for x in getfield.(sim_single.buyers, :quality_expectation_history)])

########################################### Sampled product quality
## -> in the mean shall be 'close' to producers' decisions

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0., 0., "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[.8, 2.], [.8, 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

plot_quantity(sim_single.sellers, 1)

p = Plots.plot(sim_single.sellers[1].quality_history, color="black", legend = nothing)
qs = quality_sampled_products(sim_single.buyers, 1)
Plots.plot!(mean_nothing.(qs),color = "blue")
for t in 1:200
    if length(qs[t]) > 0
        Plots.plot!([t,t],[mean(qs[t])-std(qs[t]), mean(qs[t])+std(qs[t])],color="blue")
    end
end
p

########################################### System level ###########################################

########################################### No secondary market
## -> empty reselling

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, false, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

Plots.plot(getfield.(sim_single.sellers, :reselling_history))

########################################### Secondary market
## -> market effectiveness?

sim_single = TO_GO(200, 2, 400, 600, [0.4, 0.4], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

Plots.plot(getfield.(sim_single.sellers, :reselling_history))

########################################### Initial margin 2.0
## -> no initial trading, sellers to decrease

sim_single = TO_GO(200, 2, 400, 600, [1.0, 1.0], [2.0, 2.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, true, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

Plots.plot(getfield.(sim_single.sellers, :quantity_sold_history))
Plots.plot(getfield.(sim_single.sellers, :margin_history))

########################################### Cost coefficient == 1
## -> no trading?

sim_single = TO_GO(200, 2, 400, 600, [1.0, 1.0], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [[0.2, 0.8], [0.2, 0.8]], [[0.2, 0.8], [0.2, 0.8]], [[1., 2.], [1., 2.]], 0.50, false, 0, [0.7, 1.], "softmax", [false, true], [0, 0.1], 8, false)

Plots.plot(getfield.(sim_single.sellers, :quantity_sold_history))
Plots.plot(getfield.(sim_single.sellers, :selling_income_history))

ke = [getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]
de = [getindex.(x,1) for x in getfield.(sim_single.buyers, :durability_expectation_history)]
srp = getfield.(sim_single.buyers, :std_reservation_price)
fd = getfield.(sim_single.buyers, :future_discount)

u(β,ρ,k,d;H) = β .* k .* (1 .- (ρ.*d).^H) ./ (1 .- (ρ.*d))

ue = u.(srp, fd, ke, de; H = 8)

pr = calculate_price_history(sim_single.sellers[1]; product_life = 8)
pr = fill(pr, 400)

_surplus = [x .> 0 for x in ue .- pr]

any(_surplus[1])

srp[5]
fd[5]

getindex.(calculate_price_history.(sim_single.sellers; product_life= 8), 68)

sim_single.sellers[1].cost_coefficient
sim_single.sellers[1].margin_history

dec = sim_single.buyers[5].unit_buying_selling_history[1]
dec.surplus .== dec.utility .- dec.prices
sim_single.buyers[5].quality_expectation_history[4]
Plots.plot(ke[5])
Plots.plot!(sim_single.sellers[1].quality_history)
Plots.plot(sim_single.sellers[1].margin_history)

utility_expectation_buyers = [[] for _ = 1:200] # to T

for _buyer in sim_single.buyers
    for item in _buyer.unit_buying_selling_history
        if item.decision == "buy, primary market"
            if item.product == 1
                push!(utility_expectation_buyers[item.t], item.utility[1])
            end
        end
    end
end

mean_nothing(x) = length(x) == 0 ? missing : mean(x)



kd(k,d;H) = k .* (1 .- d .^ H) ./ (1 .- d)

Plots.plot(kd.(sim_single.sellers[1].quality_history, sim_single.sellers[1].durability_history; H = 8))
Plots.plot!(mean_nothing.(utility_expectation_buyers))