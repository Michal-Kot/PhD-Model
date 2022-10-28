############## użyteczność vs. k i d

using Plots
using LaTeXStrings
pyplot()

function utility_stream(k,d;h)
    return k * (1-d^h)/(1-d)
end

[average_cost(k,d;h=1) for k in K, d in D]

K = LinRange(0.010:0.010:0.990)
D = LinRange(0.010:0.010:0.990)

ac_p1 = contour(K,D, (K,D)->utility_stream(K,D;h=1), levels=50, xlabel = "Wstępna jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 1", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
ac_p2 = contour(K,D, (K,D)->utility_stream(K,D;h=2), levels=50, xlabel = "Wstępna jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 2", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
ac_p5 = contour(K,D, (K,D)->utility_stream(K,D;h=5), levels=50, xlabel = "Wstępna jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 5", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)
ac_p10 = contour(K,D, (K,D)->utility_stream(K,D;h=10), levels=50, xlabel = "Wstępna jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 10", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

ac_p = plot(ac_p1, ac_p2, ac_p5, ac_p10, layout = (2,2))

savefig(ac_p, pwd() * "\\Documents\\GitHub\\PhD-Model\\thesis_plots\\contour_utility_stream.pdf")

############## sieć powiązań

using Graphs
using GraphPlot
using Compose, Cairo, Fontconfig

function neighbourhood(g, i, k)
    n_v = []
    for kk in 1:k
        n_v_l = []
        if kk == 1
            n_v_l = Graphs.neighbors(g,i)
        else
            for ii in n_v
                n_v_l = unique(vcat(n_v_l, Graphs.neighbors(g,ii)...))
            end
        end
        n_v = unique(vcat(n_v, n_v_l...))
    end
    n_v = n_v[n_v .!= i]
    return unique(n_v)
end

G = Graphs.erdos_renyi(30, 45)
nlist = Vector{Vector{Int}}(undef, 2) # two shells
nlist[1] = 1:10 # first shell
nlist[2] = 11:Graphs.nv(G) # second shell
locs_x, locs_y = shell_layout(g, nlist)

# n = 1

println(cgrad(:matter, 5, categorical = true)[1])

colors = fill(cgrad(:matter, 5, categorical = true)[5], Graphs.nv(G))
colors[neighbourhood(G, 1, 4)] .= cgrad(:matter, 5, categorical = true)[4]
colors[neighbourhood(G, 1, 3)] .= cgrad(:matter, 5, categorical = true)[3]
colors[neighbourhood(G, 1, 2)] .= cgrad(:matter, 5, categorical = true)[2]
colors[neighbourhood(G, 1, 1)] .= cgrad(:matter, 5, categorical = true)[1]
colors[1] = RGBA(1.,1.,1.)
gr1 = gplot(G, locs_x, locs_y, nodelabel = 1:Graphs.nv(G), nodefillc = colors, nodestrokec = "black", nodestrokelw  = 1)
draw(PDF(pwd() * "\\Documents\\GitHub\\PhD-Model\\thesis_plots\\network_graph.pdf", 16cm, 16cm), gr1)

############## oczekiwana dalsza użyteczność

using Plots
using LaTeXStrings

H = 10
ts = 1
k = 0.5
d = 0.5
β = 0.75
ρ = 0.95

function utility_next(t, ts, H, β, k, d, ρ)
    @assert t >= ts
    @assert t <= ts + H
    β * k * d^(t-ts) * (1 - (ρ * d)^(H+ts-t)) / (1 - ρ * d)
end

[utility_next(t, ts, H, β, k, d, ρ) for t in 1:10]

eu_p = plot([utility_next(t, ts, H, β, k, 0.25, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.25", title = "Oczekiwana użyteczność z konsumpcji a " * L"\hat D_{ijt}")

plot!([utility_next(t, ts, H, β, k, 0.50, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.50")

plot!([utility_next(t, ts, H, β, k, 0.75, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.75")

savefig(eu_p, pwd() * "\\thesis_plots\\expected_further_utility.pdf")


############## cechy dobra a oczekiwania

using Plots
using LaTeXStrings
pyplot()

ρ = 0.95

function q(t, ts, H, k, d)
    @assert t >= ts
    @assert t <= ts + H
    k * d^(t-ts) * (1 - d^(H+ts-t)) / (1 - d)
end

function beta_req(ke,de;t,ts,H,k,d,c,m,ρ)
    return c * m * q(t,ts,H,k,d) / q(t,ts,H, ke*ρ, de*ρ)
end

c = 0.4
m = 1.1

q(1,1,10,0.5,0.5)
q(1,1,10,0.5,0.5)
beta_req(0.5,0.5;t=1,ts=1,H=10,k=0.5,d=0.5,c=c,m=m,ρ=ρ)

K = LinRange(0.30:0.010:0.70)
D = LinRange(0.30:0.010:0.70)

beta_min_p1 = contour(K,D, (K,D)->beta_req(K,D;t=1,ts=1,H=10,k=0.5,d=0.5,c=c,m=m,ρ=ρ), levels=100, xlabel = "Oczekiwana jakość dobra, " * L"K_{ijt}^E", ylabel = "Oczekiwana trwałość dobra, " * L"D_{ijt}^E", title = "Minimalna wartość parametru " * L"\beta_i" * " aby doszło do zakupu", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

savefig(beta_min_p1, pwd() * "\\thesis_plots\\beta_min_to_purchase.pdf")