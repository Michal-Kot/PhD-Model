using Distributions, StatsBase
using OrderedCollections
using Plots
using StatsPlots
using LinearAlgebra
using InvertedIndices
using Graphs
using GraphPlot, Compose, Cairo, Fontconfig
using HypothesisTests
using SmoothingSplines
using GLM
using DataFrames

### TO DO: add duration impact on buyers decisions


################################# STRUCTS ####################################################################

"""

Sellers form a supply side of the model

Sellers' main goal is to maximize profit π
π = pq - c(q)p
p - price of product, q - quantity of sales, c - function of cost, linear on q, a - advertising intensity

"""

mutable struct seller
    id::Int64
    quality::Float64
    durability::Int64
    average_cost::Float64
    margin::Float64
    quantity_history::Vector{Float64}
    margin_history::Vector{Float64}
    advertising::Float64
    advertising_history::Vector{Float64}
    risk_decisions::Float64
    persuasiveness::Float64
end

function create_sellers(num_sellers::Int64,q::Vector{Float64},c::Vector{Float64},m::Vector{Float64},a::Vector{Float64},r::Vector{Float64},p::Vector{Float64},d::Vector{Int64})::Vector{seller}
    @assert length(q) == num_sellers
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(a) == num_sellers
    @assert length(r) == num_sellers
    @assert length(p) == num_sellers
    sellers_vector = []
    for s in 1:num_sellers
        new_seller = seller(s, q[s], d[s], c[s], m[s], [], [], a[s], [], r[s], p[s])
        push!(sellers_vector, new_seller)
    end
    return sellers_vector
end

"""

Buyers form a demand side of the model

Buyers main goal is to maximize utility of consumption υ(q) = s*eq - p, subject to budget constraint ρ*eq >= p
s - seeking for quality parameter, eq = expected quality of the product, p - price of product, ρ - willingness to pay for unit of quality

"""

mutable struct buyer
    id::Int64
    neighbours::Vector{Int64}
    broken_product::Bool 
    std_reservation_price::Float64
    quality_seeking::Float64
    quality_expectation::Vector{Float64}
    quality_expectation_history::Vector{Vector{Float64}}
    durability_expectation::Vector{Float64}
    durability_expectation_history::Vector{Vector{Float64}}
    unit_bought::Vector{Bool}
    unit_bought_history::Vector{Vector{Bool}}
    unit_bought_time::Int64
    quality_of_unit_bought::Float64
    quality_of_unit_bought_history::Vector{Vector{Float64}}
    ad_received::Vector{Bool}
    ad_received_history::Vector{Vector{Bool}}
    surplus_history::Vector{Float64}
end

function create_buyers(num_buyers::Int64, num_sellers::Int64, network::SimpleGraph{Int64}, initial_quality_expectation::Vector{Float64}=fill(1.0, num_sellers), initial_duration_expectation::Vector{Float64}=fill(5.0, num_sellers))::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)
        new_buyer = buyer(b, #id
            my_neighbours, #neighbours
            true, #broken_product
            rand(Uniform(0,1)), #std_reservation_price
            rand(Uniform(0,1)), #quality_seeking
            initial_quality_expectation, #quality_expectation
            [], #quality_expectation_history
            initial_duration_expectation, #duration_expectation
            [], #duration_expectation_history
            fill(false, num_sellers), #unit_bought
            [], #unit_bought_history,
            0, #unit_bought_time
            0.0, #quality_of_unit_bought
            [], #quality_of_unit_bought_history
            [], #ad_received
            [], #ad_received_history
            []) #surplus_history
        push!(buyers_vector, new_buyer)
    end
    return buyers_vector
end

"""

Create network of a given type - vertices are buyers

if random type, then use Erdos Renyi algorithm, specify number of vertices and edges
if preferential attachment type , then use Albert Barabasi algorithm, specify number of vertices and number of links for a single vertex (constant)

"""

function create_network(type::String; num_buyers::Int64, num_links::Int64=0, pref_attachment_links::Int64=0)

    if type == "random"
        return Graphs.erdos_renyi(num_buyers, num_links)
    elseif type == "preferential_attachment"
        return Graphs.barabasi_albert(num_buyers, pref_attachment_links)
    end
end

#################################### AUX FUNCTIONS ##############################################################

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64
    return max(min(x,ub),lb)
end

function calculate_cost(_seller::seller)::Float64
    cost = _seller.durability * _seller.quality * _seller.average_cost
    return cost
end

function calculate_price(_seller::seller)::Float64
    price = calculate_cost(_seller) * (1 + _seller.margin)
    return price
end

function calculate_price_history(_seller::seller)::Vector{Float64}
    price = calculate_cost(_seller) .* (1 .+ _seller.margin_history)
    return price
end

function u2w(u::Vector{Float64}, p_min::Float64 = 0.1)::Vector{Float64}
    @assert p_min < 1.0
    ut = u .- minimum(u)
    wgt = ut ./ sum(ut) .* (1 - p_min)
    wgt[argmin(u)] = p_min
    return wgt
end

function calculate_profit_history(_seller::seller, γ::Float64, num_buyers::Int64)::Vector{Float64}
    (calculate_price_history(_seller) .- calculate_cost(_seller)) .* _seller.quantity_history .- γ * num_buyers * _seller.advertising_history
end

function plot_phasediagram(_seller::seller, function_args::Vector)
    x = calculate_price_history(_seller)
    y = calculate_profit_history(_seller, function_args.γ, function_args.num_buyers)
    d_x = diff(x)
    push!(d_x, 0)
    d_y = diff(y)
    push!(d_y,0)
    labs = string.(1:length(x))
    points_labels = ([mod(x,100) == 0 ? labs[x] : "" for x in 1:length(labs)])
    points_labels[1] = "1"
    quiver(x,y,quiver=(d_x, d_y), color = :green)
    scatter!(x,y, xlabel = "Price", ylabel = "Profit", markershape = :none, markercolor = :white, markerstrokecolor = :white, series_annotations = points_labels, markersize = 0)
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}
    x = fill(false, n)
    x[k] = true
    return x
end

function calculate_total_surplus(sim_res, return_type::String = "total", cumulated::Bool = true)

    producer_surplus = sum(calculate_profit_history.(sim_res.sellers, sim_res.function_args.γ, sim_res.function_args.num_buyers))
    consumer_surplus = sum(getfield.(sim_res.buyers, :surplus_history))
    total_surplus = producer_surplus .+ consumer_surplus

    if return_type == "total"

        if cumulated 

            return sum(total_surplus)

        else

            return total_surplus

        end

    elseif return_type == "producer"

        if cumulated 

            return sum(producer_surplus)

        else

            return producer_surplus

        end

    elseif return_type == "consumer"

        if cumulated 

            return sum(consumer_surplus)

        else

            return consumer_surplus

        end

    end

end

function seller_price_advertising_adjustment(sellers, iter, δ, α, γ, num_buyers, variant_advertising, allow_negative_margin)
    if iter == 1

        # no change to price in iter 1

        for _seller in sellers

            push!(_seller.margin_history, _seller.margin)
            push!(_seller.advertising_history, _seller.advertising)

        end

    elseif iter == 2

        # random change to price in iter 2
        # price changes by δ
        # advertising changes by α

        for _seller in sellers

            change_margin = sample([-1,0,1], Weights(fill(1/3,3)))
            new_margin = _seller.margin + change_margin * rand() * δ
            _seller.margin = new_margin
            push!(_seller.margin_history, _seller.margin)

            if variant_advertising
                change_advertising = sample([0,1], Weights(fill(1/2,2)))
                new_advertising = _seller.advertising + change_advertising * rand() * α
                new_advertising = in_boundaries(new_advertising, 0.0, 0.10)
                _seller.advertising = new_advertising
            end

            push!(_seller.advertising_history, _seller.advertising)

        end

    else

        # from iter 3 to the end, firm decides about the price, by comparing price change and resulting profit change

        for _seller in sellers

            # margin has to be in tha range between cost of production & 2
            #margin_ub = 2 - _seller.quality * _seller.average_cost

            # how did margin change?
            margin_change = sign(_seller.margin_history[end] - _seller.margin_history[end-1])

            # how did advertising change?

            advertising_change = sign(_seller.advertising_history[end] - _seller.advertising_history[end-1])

            # how did profit change?

            profit_history = calculate_profit_history(_seller, γ, num_buyers)
            profit_change = profit_history[end] - profit_history[end-1]
            
            if profit_change > 0

                if rand() <= _seller.risk_decisions
                    new_margin = _seller.margin + margin_change * rand() * δ
                else
                    new_margin = _seller.margin
                end

                if rand() <= _seller.risk_decisions
                    new_advertising = _seller.advertising + advertising_change * rand() * α
                else
                    new_advertising = _seller.advertising
                end

            elseif profit_change < 0

                if rand() <= _seller.risk_decisions
                    new_margin = _seller.margin_history[end-1] + (-1) * margin_change * rand() * δ
                else
                    new_margin = _seller.margin_history[end-1]
                end

                if rand() <= _seller.risk_decisions
                    new_advertising = _seller.advertising_history[end-1] + (-1) * advertising_change * rand() * α
                else
                    new_advertising = _seller.advertising_history[end-1]
                end

            elseif profit_change == 0

                if rand() <= _seller.risk_decisions
                    margin_change = sample([-1,0,1], Weights(fill(1/3,3)))
                    new_margin = _seller.margin + margin_change * rand() * δ
                else
                    new_margin = _seller.margin
                end

                if rand() <= _seller.risk_decisions
                    advertising_change = sample([-1,0,1], Weights(fill(1/3,3)))
                    new_advertising = _seller.advertising + advertising_change * rand() * α
                else
                    new_advertising = _seller.advertising
                end

            end

            if iter >= 6

                if all((_seller.margin_history[(end-4):end] .* _seller.quantity_history[(end-4):end]) .<= 0) #& (_seller.margin > 0)
                    new_margin = _seller.margin - rand() * δ
                end

                if (all(calculate_profit_history(_seller, γ, num_buyers)[(end-4):end] .< 0)) & (_seller.margin < 0)
                    new_margin = _seller.margin + δ
                end

            end

            if !allow_negative_margin
                new_margin = in_boundaries(new_margin, 0, 2)
            end

            _seller.margin = new_margin
            push!(_seller.margin_history, _seller.margin)

            if variant_advertising
                new_advertising = in_boundaries(new_advertising, 0.0, 0.10)
                _seller.advertising = new_advertising
            end

            push!(_seller.advertising_history, _seller.advertising)
        
        end

    end

    return sellers

end

function plot_quantity(sim_res)
    groupedbar(reduce(hcat,[getfield.(sim_res.sellers, :quantity_history)[x] for x in 1:length(sim_res.sellers)]),linecolor=nothing, bar_position = :stack, xlabel = "Time", ylabel = "Quantity of sales", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_quality_expectation(sim_res)
    plot([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers], label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_profit_history(sim_res)
    plot(calculate_profit_history.(sim_res.sellers, sim_res.function_args.γ, sim_res.function_args.num_buyers))
end

function plot_advertising(sim_res)
    plot(getfield.(sim_res.sellers, :advertising_history), xlabel = "Time", ylabel = "Advertising intensity [share of population targeted]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_margin(sim_res)
    Plots.plot([getfield.(sim_res.sellers, :margin_history)[x] for x in 1:length(sim_res.sellers)], xlabel = "Time", ylabel = "Margin", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_ecdf(metric, label, xlabel, ylabel, title, is_new)
    if is_new
        Plots.plot(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright)
    else
        Plots.plot!(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright)
    end
end


############################################ MAIN SIM FUNCTION #####################################################

function TO_GO(num_sellers::Int64, num_buyers::Int64, max_iter::Int64, λ_ind::Float64, λ_wom::Float64, consumer_behaviour::String; q::Vector{Float64} = fill(1.0, num_sellers), c::Vector{Float64} = fill(0.8, num_sellers), m::Vector{Float64} = fill(0.2, num_sellers), a::Vector{Float64} = fill(0.05, num_sellers), r::Vector{Float64} = fill(0.05, num_sellers), ϵ::Vector{Float64} = fill(1/3, num_sellers), q_init::Vector{Float64} = fill(1.0, num_sellers), p::Vector{Float64} = fill(0.05, num_sellers), d::Vector{Int64} = fill(5, num_sellers), d_init::Vector{Float64} = fill(5.0, num_sellers), num_links::Int64 = 200, δ::Float64 = 0.05, γ::Float64 = 0.050, α::Float64 = 0.005, variant_advertising::Bool = true, allow_negative_margin::Bool = true)

    function_args = (num_sellers = num_sellers, num_buyers = num_buyers, max_iter = max_iter, λ_ind = λ_ind, λ_wom = λ_wom, q = q, c = c, m = m, a = a, r = r, ϵ = ϵ, q_init = q_init, p=p, d=d, d_init=d_init, num_links = num_links, δ = δ, γ = γ, α=α)

    # create sellers and buyers

    sellers = create_sellers(num_sellers, q, c, m, a, r, p, d)

    buyers_network = create_network("random", num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, q_init, d_init)

    # utlity function definition

    u(k,d,p;β) = ((0.5 + β) .* k .* d) ./ p

    # TEMP 🍎

    utility_history = []

    # sellers are choosing new price

    #margin_lb = -2
    #margin_ub = 2

    for iter in 1:max_iter

        sellers = seller_price_advertising_adjustment(sellers, iter, δ, α, γ, num_buyers, variant_advertising, allow_negative_margin)

        #### RYNEK ####

        demand = Int64[]

        # REKLAMY

        advertising_intensity = getfield.(sellers, :advertising)

        for _buyer in buyers
            # Advertising

            if rand() < sum(advertising_intensity)
                firm_ad_seen = sample(1:num_sellers, Weights(advertising_intensity))
                advertising_bool = fill(false, num_sellers)
                advertising_bool[firm_ad_seen] = true
                _buyer.ad_received = advertising_bool
                push!(_buyer.ad_received_history, _buyer.ad_received)
            else
                _buyer.ad_received = fill(false, num_sellers)
                push!(_buyer.ad_received_history, _buyer.ad_received)
            end
            
        end 

        for _buyer in buyers

            # Product usage -> if broken

            if iter >= 2

                if any(_buyer.unit_bought)

                    product_durability = sellers[argmax(_buyer.unit_bought)].durability
                    product_amortization = iter - _buyer.unit_bought_time
                    probability_to_break = 0.5 / product_durability * product_amortization

                    if rand() < probability_to_break
                        _buyer.broken_product = true
                        _buyer.durability_expectation = _buyer.durability_expectation .+ λ_ind .* _buyer.unit_bought .* (product_amortization .- _buyer.durability_expectation)
                        for idn in _buyer.neighbours
                            buyers[idn].durability_expectation = buyers[idn].durability_expectation .+ λ_wom .* _buyer.unit_bought .* (product_amortization .- buyers[idn].durability_expectation)
                        end
                    end

                end

            end

        end

        prices = calculate_price.(sellers)

        for _buyer in buyers

            # buying new product in iter == 1, or replacing product in iter >= 2

            if _buyer.broken_product | all(_buyer.unit_bought .== false)

                wtp_advertising = 1 .+ (_buyer.ad_received .* p)

                buyer_wtp = _buyer.std_reservation_price .* _buyer.quality_expectation .* _buyer.durability_expectation .* wtp_advertising 

                available_products = prices .<= buyer_wtp # 🍎 rozspisać

                if any(available_products)

                    utility = u.(_buyer.quality_expectation, _buyer.durability_expectation, prices; β = _buyer.quality_seeking)

                    push!(utility_history, (id=_buyer.id, qs=_buyer.quality_seeking, qe=_buyer.quality_expectation, p=prices, ap=available_products, u=utility, wtp = buyer_wtp, ad = _buyer.ad_received, de=_buyer.durability_expectation,srp = _buyer.std_reservation_price))

                    if consumer_behaviour == "deterministic"

                        # consumer choses product that maximizes utility

                        chosen_product = argmax(utility .- 1e10 .* .!available_products)

                    elseif consumer_behaviour == "stochastic"

                        # consumer randomly choses products, accordingly to utility (may choose not the best one)

                    weight = u2w(utility) .* available_products

                        chosen_product = sample(1:num_sellers, Weights(weight))

                    end

                    _buyer.unit_bought = create_bool_purchase(num_sellers, chosen_product)
                    _buyer.broken_product = false
                    _buyer.unit_bought_time = iter

                    push!(_buyer.unit_bought_history, _buyer.unit_bought)
                    push!(demand, chosen_product)
    
                    buyer_surplus = _buyer.quality_expectation[chosen_product] * _buyer.std_reservation_price * _buyer.durability_expectation[chosen_product]  * wtp_advertising[chosen_product] - prices[chosen_product]# 🍎
                    push!(_buyer.surplus_history, buyer_surplus)

                else

                    _buyer.unit_bought = fill(false, num_sellers)
                    _buyer.broken_product = false
                    push!(_buyer.unit_bought_history, _buyer.unit_bought)
                    push!(_buyer.surplus_history, 0.0)

                end

            else

                push!(_buyer.surplus_history, 0.0)

            end

        end

        for _buyer in buyers

            # calculate quality driven by the product

            if any(_buyer.unit_bought) # if product is possessed

                chosen_product = argmax(_buyer.unit_bought)

                _buyer.quality_of_unit_bought = rand(Uniform(sellers[chosen_product].quality - ϵ[chosen_product], sellers[chosen_product].quality + ϵ[chosen_product]))

                push!(_buyer.quality_of_unit_bought_history, _buyer.unit_bought .*  _buyer.quality_of_unit_bought)

                _buyer.quality_expectation = _buyer.quality_expectation .+ λ_ind .* _buyer.unit_bought .* (_buyer.quality_of_unit_bought .- _buyer.quality_expectation)

                for idn in _buyer.neighbours
                    buyers[idn].quality_expectation = buyers[idn].quality_expectation .+ λ_wom .* _buyer.unit_bought .* (_buyer.quality_of_unit_bought .- buyers[idn].quality_expectation)
                end
            
            else

                _buyer.quality_of_unit_bought = 0.0
                push!(_buyer.quality_of_unit_bought_history, fill(0.0, num_sellers))

            end

        end

        # kupujący weryfikują oczekiwaną jakość produktów

        for _buyer in buyers

            push!(_buyer.quality_expectation_history, _buyer.quality_expectation)
            push!(_buyer.durability_expectation_history, _buyer.durability_expectation)

        end

        # firma weryfikuje Q sprzedaży

        for _seller in sellers

            push!(_seller.quantity_history, count(demand .== _seller.id))
        end

    end

    return (sellers = sellers, buyers = buyers, function_args = function_args, ut_his = utility_history)

end

__precompile__()