#### Data Analysis ####

using CSV, DataFrames, DataFramesMeta, Query, Chain, Plots, StatsBase

df_eurostat = CSV.read("C:\\Users\\User\\Downloads\\namq_10_fcs_linear.csv", DataFrame)

#### Porównanie z krajami o podobnym poziomie PKB ####

df_countries_durables = @chain df_eurostat begin
    @subset(:na_item .== "P311_S14")
    @subset(:s_adj .== "SCA")
    @subset(:unit .== "CP_MEUR")
    @subset(parse.(Int, SubString.(:TIME_PERIOD,1,4)) .>= 2010)
    @subset(map(x -> x ∈ ["PL", "NL", "SE", "BE", "IE"], :geo))
    @by([:TIME_PERIOD, :geo], :mean_turnover = mean(:OBS_VALUE))
    unstack(:TIME_PERIOD, :geo, :mean_turnover)
end

countries_names_pl = Dict(["PL", "NL", "SE", "BE", "IE"] .=> ["Polska", "Holandia", "Szwecja", "Belgia", "Irlandia"])

countries_colors = Dict(["PL", "NL", "SE", "BE", "IE"] .=> ["red", "orange", "blue", "black", "green"])

p1 = Plots.plot()

for cc in axes(df_countries_durables,2)
    if cc > 1
        Plots.plot!(df_countries_durables.TIME_PERIOD, df_countries_durables[:,cc], label = countries_names_pl[names(df_countries_durables)[cc]], color = countries_colors[names(df_countries_durables)[cc]])
    end
end

Plots.plot!(xlabel = "Rok, kwartał", ylabel = "Wydatki gospodarstw domowych, w MLN €", ylim = (0,10000), legend = :outertopright, xtickfontsize = 6)

p1

Plots.savefig(p1, pwd() * "\\thesis_plots\\durables_expenditure.pdf")

#### Porównanie z krajami o podobnym poziomie PKB per capita ####

df_countries_durables2 = @chain df_eurostat begin
    @subset(:na_item .== "P311_S14")
    @subset(:s_adj .== "SCA")
    @subset(:unit .== "CP_MEUR")
    @subset(parse.(Int, SubString.(:TIME_PERIOD,1,4)) .>= 2010)
    @subset(map(x -> x ∈ ["PL", "HR", "HU", "EL", "RO"], :geo))
    @by([:TIME_PERIOD, :geo], :mean_turnover = mean(:OBS_VALUE))
    unstack(:TIME_PERIOD, :geo, :mean_turnover)
end

df_countries_total2 = @chain df_eurostat begin
    @subset(:na_item .== "P31_S14")
    @subset(:s_adj .== "SCA")
    @subset(:unit .== "CP_MEUR")
    @subset(parse.(Int, SubString.(:TIME_PERIOD,1,4)) .>= 2010)
    @subset(map(x -> x ∈ ["PL", "HR", "HU", "EL", "RO"], :geo))
    @by([:TIME_PERIOD, :geo], :mean_turnover = mean(:OBS_VALUE))
    unstack(:TIME_PERIOD, :geo, :mean_turnover)
end

countries_names_pl2 = Dict(["PL", "HR", "HU", "EL", "RO"] .=> ["Polska", "Chorwacja", "Węgry", "Grecja", "Rumunia"])

countries_colors2 = Dict(["PL", "HR", "HU", "EL", "RO"] .=> ["red", "grey", "green", "blue", "orange"])

p2 = Plots.plot()

for cc in axes(df_countries_durables,2)
    if cc > 1
        Plots.plot!(df_countries_durables2.TIME_PERIOD, 100*df_countries_durables2[:,cc] ./ df_countries_total2[:,cc], label = countries_names_pl2[names(df_countries_durables2)[cc]], color = countries_colors2[names(df_countries_durables2)[cc]])
    end
end

Plots.plot!(xlabel = "Rok, kwartał", ylabel = "Udział wydatków na dobra trwałe, w %", ylim = (2,10.1), legend = :outertopright, xtickfontsize = 6)

p2

Plots.savefig(p2, pwd() * "\\thesis_plots\\durables_percentage.pdf")