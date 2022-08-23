function simulate_outcome_of_strategy(my_k, my_d, my_m, my_c, comp_k, comp_d, comp_p)

    s = LinRange(0:0.01:1.0)

    my_p = my_c * my_k / (1 - my_d) * my_m
    
    my_w = s .* my_k / (1 - my_d) .- my_p

    comp_w = [s .* ck / (1 - cd) .- cp for (ck, cd, cp) in zip(comp_k, comp_d, comp_p)] # to napisać wektorowo, gdy length(comp) > 0

    is_my_w_positive = my_w .>= 0

    is_my_w_best = reduce(.*, [my_w .> cw for cw in comp_w])

    my_demand = count(is_my_w_best .& is_my_w_positive)

    my_profit = my_demand * (my_p - my_c * my_k / (1 - my_d))

    return my_profit, my_demand

end

strategy_outcomes = []

for comb in vec(collect(Iterators.product([0.9,1.0,1.1], [0.8,0.85,0.90], [1.1, 1.15, 1.20])))

    my_k = comb[1]
    my_d = comb[2]
    my_m = comb[3]

    y = simulate_outcome_of_strategy(my_k, my_d, my_m, 0.7, [1.1, 0.9, 0.7], [0.7, 0.8, 0.9], [0.6 * 1.1 / (1 - 0.7) * 1.2, 0.5 * 0.9 / (1 - 0.8) * 1.2, 0.6 * 0.7 / (1 - 0.9) * 1.2])

    push!(strategy_outcomes, (comb, y[1], y[2]))
    
end

k, d, m = 1.0, 0.6, 1.15

y_hist = []

for iter in 1:50

    strategy_outcomes = []

    k_range = [k - 0.05, k, k + 0.05]
    d_range = [d - 0.05, d, d + 0.05]
    m_range = [m - 0.05, m, m + 0.05]



    for comb in vec(collect(Iterators.product(k_range, d_range, m_range)))
    
        my_k = comb[1]
        my_d = comb[2]
        my_m = comb[3]
    
        y = simulate_outcome_of_strategy(my_k, my_d, my_m, 0.5, [1.1, 0.9, 0.7], [0.6, 0.6, 0.6], [0.6 * 1.1 / (1 - 0.6) * 1.2, 0.5 * 0.9 / (1 - 0.6) * 1.2, 0.6 * 0.7 / (1 - 0.6) * 1.2])

        push!(strategy_outcomes, (comb, y[1], y[2]))
        
    end

    push!(y_hist, strategy_outcomes)

    k, d, m = getindex.(strategy_outcomes,1)[argmax(getindex.(strategy_outcomes, 2))]

    k = in_boundaries(k, 0.5, 1.5)
    d = in_boundaries(d, 0.2, 0.8)
    m = in_boundaries(m, 0.0, 2.0)

    println((k,d,m))
    println(maximum(getindex.(strategy_outcomes, 2)))

end

s1_b = argmax(getindex.(y_hist[1], 2))
y_hist[1][s1_b]

getindex.(strategy_outcomes,1)[argmax(getindex.(strategy_outcomes,2))]




function calculate_optimal_margin_and_advertising(my_m, my_dq, comp_p, comp_dq, my_c, δ_m, my_a, δ_a, my_i, comp_a, comp_i, seller_behaviour = "stochastic", verbose = false)

    my_m_current = my_m
    my_m_up = my_m + δ_m
    my_m_down = my_m - δ_m

    my_k_current = my_k
    my_k_sup = my_k + δ_k
    my_k_sdown = my_k - δ_k

    my_d_current = my_d
    my_d_up = my_d + δ_d
    my_d_down = my_d - δ_d


    s = LinRange(0:0.01:1.0)

    comp_w = [s .* cdq .* (1 + ca * ci) for (cdq, ca, ci) in zip(comp_dq, comp_a, comp_i)]

    m_a_profit = []

    for a in [my_a_adsdown, my_a_current, my_a_adsup]

        my_w = s .* my_dq .* (1 + a * my_i)

        my_spl_current = my_w .- my_p_current
        my_spl_priceup = my_w .- my_p_priceup
        my_spl_pricedown = my_w .- my_p_pricedown

        comp_spl = [cw .- cp for (cw,cp) in zip(comp_w, comp_p)]

        my_spl_diff_current = reduce(.*,[my_spl_current .> cspl for cspl in comp_spl])
        my_spl_diff_priceup = reduce(.*,[my_spl_priceup .> cspl for cspl in comp_spl])
        my_spl_diff_pricedown = reduce(.*,[my_spl_pricedown .> cspl for cspl in comp_spl])

        my_demand_current = count((my_spl_current .>= 0) .& (my_spl_diff_current))
        my_demand_priceup = count((my_spl_priceup .>= 0) .& (my_spl_diff_priceup))
        my_demand_pricedown = count((my_spl_pricedown .>= 0) .& (my_spl_diff_pricedown))

        my_margin_current = my_p_current - my_dq * my_c
        my_margin_priceup = my_p_priceup - my_dq * my_c
        my_margin_pricedown = my_p_pricedown - my_dq * my_c

        my_π_current = my_demand_current * my_margin_current - 2 * a
        my_π_priceup = my_demand_priceup * my_margin_priceup - 2 * a
        my_π_pricedown = my_demand_pricedown * my_margin_pricedown - 2 * a

        push!(m_a_profit, [my_π_pricedown, my_π_current, my_π_priceup])

    end

    profit_matrix = hcat(m_a_profit...)

    new_m = -1
    new_a = -1

    if seller_behaviour == "deterministic"

        my_decision = Tuple(argmax(profit_matrix))

    elseif seller_behaviour == "stochastic"

        my_decision = sample(vec(reshape(CartesianIndices(profit_matrix), 9, 1)), Weights(vec(reshape(profit_matrix, 9, 1))))

    end

    m_decision = my_decision[1]
    a_decision = my_decision[2]
    
    if m_decision == 1

        new_m = my_m - δ_m

    elseif m_decision == 2

        new_m = my_m

    elseif m_decision == 3

        new_m = my_m + δ_m

    end

    if a_decision == 1

        new_a = my_a - δ_a

    elseif a_decision == 2

        new_a = my_a

    elseif a_decision == 3

        new_a = my_a + δ_a

    end

    if verbose

        if my_decision == (1,1)

            println("M down, A down")

        elseif my_decision == (1,2)

            println("M down, A keep")

        elseif my_decision == (1,3)

            println("M down, A up")

        elseif my_decision == (2,1)

            println("M keep, A down")

        elseif my_decision == (2,2)

            println("M keep, A keep")

        elseif my_decision == (2,3)

            println("M keep, A up")

        elseif my_decision == (3,1)

            println("M up, A down")

        elseif my_decision == (3,2)

            println("M up, A keep")

        elseif my_decision == (3,3)

            println("M up, A up")

        end

    end

    return new_m, new_a

end