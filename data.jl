df_node = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/node_info.csv", DataFrame)
df_gen = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/gen_info.csv", DataFrame)
df_load = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/load_info.csv", DataFrame)
df_line = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/line_info.csv", DataFrame)
df_DC = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/DC_info.csv", DataFrame)
incidence = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/incidence.csv", DataFrame)
incidence_dc = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/incidence_dc.csv", DataFrame)
susceptance = CSV.read("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/susceptance.csv", DataFrame)

xf_renew = XLSX.readxlsx("_ModelsSchonheitAdjusted_LessTimeSteps_OBZ/data_renew_2015.xlsx")
df_pv = DataFrame(xf_renew["pv"][:][2:end,:], :auto)
rename!(df_pv, Dict(names(df_pv)[i] => Symbol.(xf_renew["pv"][:][1,:])[i] for i = 1:ncol(df_pv)))
df_wind = DataFrame(xf_renew["onshore"][:][2:end,:], :auto)
rename!(df_wind, Dict(names(df_wind)[i] => Symbol.(xf_renew["onshore"][:][1,:])[i] for i = 1:ncol(df_wind)))
df_wind_off = DataFrame(xf_renew["offshore"][:][2:end,:], :auto)
rename!(df_wind_off, Dict(names(df_wind_off)[i] => Symbol.(xf_renew["offshore"][:][1,:])[i] for i = 1:ncol(df_wind_off)))

a.ext[:sets] = Dict()
N = a.ext[:sets][:N] = df_node[:,:BusID] # Nodes
G = a.ext[:sets][:G] = [z for z = 1:maximum(df_gen[:,:GenID])] # Generators
Z = a.ext[:sets][:Z] = [z for z = 1:maximum(df_node[:,:Zone])] # Zones
L = a.ext[:sets][:L] = [l for l = 1:length(df_line[:,:BranchID])] # Lines
L_DC = a.ext[:sets][:L_DC] = [l for l = 1:length(df_DC[:,:BranchID])] # DC lines
T = a.ext[:sets][:T] = [t for t = 1:size(df_load)[1]] # Time steps
R = ["PV","Wind","Wind Offshore"]


df_node.ZoneRes = df_node.Zone
function create_res_table()
	res_temp = zeros(Float64, length(T), length(N), length(R))
	for n in N, r in R
		zone_temp = df_node.ZoneRes[df_node[:,:BusID].==n][1]
		cap_temp = sum(df_gen.Pmax[(df_gen[:,:Type].==r) .&
		                              (df_gen[:,:OnBus].==n)])
		if r == "PV"
			share_temp = df_pv[:,Symbol.(zone_temp)]
		elseif r == "Wind"
			share_temp = df_wind[:,Symbol.(zone_temp)]
		else
			share_temp = df_wind_off[:,Symbol.(zone_temp)]
		end
		res_temp[:, findfirst(N .== n), findfirst(R .== r)] = cap_temp*share_temp
	end

    res = zeros(Float64, length(T), length(N))
     for n in N
         for t in T
             res_tot = sum(res_temp[t,findfirst(N .== n),r] for r in [1,2,3])
             res[t,findfirst(N .== n)] = res_tot[1]
         end
     end
 	return res
 end

 res = create_res_table()

 function get_renew(t,n)
 	return sum(res[findfirst(T .== t), findfirst(N .== n), findfirst(R .== r)] for r in R)
 end

 
 nNodes = length(N) # amount of nodes
 nGenerators = length(G) #amount of generators
 nLines = length(L) # number of lines
 nLines_DC = length(L_DC) # number of DC lines
 nZones = maximum(df_node[:,:Zone]) # number of zones

 a.ext[:parameters] = Dict()
 GENCAP = a.ext[:parameters][:GENCAP] = df_gen[:,:Pmax]
 RES = a.ext[:parameters][:RES] = res
 LF_PV = a.ext[:parameters][:LF_PV] = Matrix(df_pv)
 LF_WINDON = a.ext[:parameters][:WINDON] = Matrix(df_wind)
 LF_WINDOFF = a.ext[:parameters][:WINDOFF] = Matrix(df_wind_off)
 MC = a.ext[:parameters][:MC] = df_gen[:,:Costs]
 DEM = a.ext[:parameters][:DEM] = Matrix(df_load)
 n_in_z = a.ext[:parameters][:n_in_z] = Dict(map(z -> z => [n for n in N if df_node[df_node[:,:BusID].==n, :Zone][1] == z], Z))
 g_in_n = a.ext[:parameters][:g_in_n] = Dict(map(n -> n => [g for g in G if df_gen[df_gen[:,:GenID].==g, :OnBus][1] == n && (df_gen[df_gen[:,:GenID].==g, :Type][1] == "Nuclear" || df_gen[df_gen[:,:GenID].==g, :Type][1] == "Gas/CCGT" || df_gen[df_gen[:,:GenID].==g, :Type][1] == "Lignite" || df_gen[df_gen[:,:GenID].==g, :Type][1] == "Hard Coal")], N))
 all_g_in_n = a.ext[:parameters][:all_g_in_n] = Dict(map(n -> n => [g for g in G if df_gen[df_gen[:,:GenID].==g, :OnBus][1] == n], N))
 l_in_z = a.ext[:parameters][:l_in_z] = Dict(map(z -> z => [l for l in df_line[:,:BranchID] if df_node[findfirst(N .== df_line[findfirst(df_line[:,:BranchID] .== l), :FromBus][1]), :Zone][1] == z || df_node[findfirst(N .== df_line[findfirst(df_line[:,:BranchID] .== l), :ToBus][1]), :Zone][1] == z], Z))
 TC = a.ext[:parameters][:TC] = 100000*df_line[:,:Pmax]
 TC_DC = a.ext[:parameters][:TC_DC] = 100000*df_DC[:,:Pmax]
 α = a.ext[:parameters][:α] = 0
 CC = a.ext[:parameters][:CC] = 0
 MinRAM = 0.7
 NTC_12, NTC_23, NTC_34, NTC_15, NTC_25, NTC_35, NTC_45 = a.ext[:parameters][:NTC_12], a.ext[:parameters][:NTC_23], a.ext[:parameters][:NTC_34], a.ext[:parameters][:NTC_15], a.ext[:parameters][:NTC_25], a.ext[:parameters][:NTC_35], a.ext[:parameters][:NTC_45] = 0.5*TC_DC[11], 0.5*TC_DC[12], 0.5*TC_DC[13], 0.5*TC_DC[2], 0.5*TC_DC[1]+0.5*TC_DC[6], 0.5*TC_DC[4]+0.5*TC_DC[5], 0.5*TC_DC[3]

#  a.ext[:parameters][:n_in_windz] = Dict()
#  a.ext[:parameters][:n_in_windz][1] = filter!(e->e∉[119, 120, 121, 122, 123, 124],a.ext[:parameters][:n_in_z][1])
#  a.ext[:parameters][:n_in_windz][2] = filter!(e->e∉[119, 120, 121, 122, 123, 124],a.ext[:parameters][:n_in_z][2])
#  a.ext[:parameters][:n_in_windz][3] = filter!(e->e∉[119, 120, 121, 122, 123, 124],a.ext[:parameters][:n_in_z][3])
#  a.ext[:parameters][:n_in_windz][4] = filter!(e->e∉[119, 120, 121, 122, 123, 124],a.ext[:parameters][:n_in_z][4])
#  a.ext[:parameters][:n_in_windz][5] = [119, 120, 121, 122, 123, 124]

t0 = 1
t1 = Int(floor(length(T)/4*1))
t2 = Int(floor(length(T)/4*2))
t3 = Int(floor(length(T)/4*3))
t4 = Int(floor(length(T)))

up_temp = collect(0:1:t1)./t1
down_temp = collect(t1:-1:0)./t1
down_temp_small = collect(t1:-1:1)./t1

RES[t0:t1+1,findfirst(N.== 119)], RES[t0:t1+1,findfirst(N.== 121)], RES[t0:t1+1,findfirst(N.== 122)] = up_temp.*GENCAP[findfirst(G .== all_g_in_n[119])], up_temp.*GENCAP[findfirst(G .== all_g_in_n[121])], up_temp.*GENCAP[findfirst(G .== all_g_in_n[122])] 
RES[t0:t1+1,findfirst(N.== 120)], RES[t0:t1+1,findfirst(N.== 123)], RES[t0:t1+1,findfirst(N.== 124)] = zeros(length(up_temp)), zeros(length(up_temp)), zeros(length(up_temp))
RES[t1+1:t2+1,findfirst(N.== 119)], RES[t1+1:t2+1,findfirst(N.== 121)], RES[t1+1:t2+1,findfirst(N.== 122)] = ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[119])], ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[121])], ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[122])]
RES[t1+1:t2+1,findfirst(N.== 120)], RES[t1+1:t2+1,findfirst(N.== 123)], RES[t1+1:t2+1,findfirst(N.== 124)] = up_temp.*GENCAP[findfirst(G .== all_g_in_n[120])], up_temp.*GENCAP[findfirst(G .== all_g_in_n[123])], up_temp.*GENCAP[findfirst(G .== all_g_in_n[124])]
RES[t2+1:t3+1,findfirst(N.== 119)], RES[t2+1:t3+1,findfirst(N.== 121)], RES[t2+1:t3+1,findfirst(N.== 122)] = down_temp.*GENCAP[findfirst(G .== all_g_in_n[119])], down_temp.*GENCAP[findfirst(G .== all_g_in_n[121])], down_temp.*GENCAP[findfirst(G .== all_g_in_n[122])]
RES[t2+1:t3+1,findfirst(N.== 120)], RES[t2+1:t3+1,findfirst(N.== 123)], RES[t2+1:t3+1,findfirst(N.== 124)] = ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[120])], ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[123])], ones(length(up_temp)).*GENCAP[findfirst(G .== all_g_in_n[124])]
RES[t3+1:t4,findfirst(N.== 119)], RES[t3+1:t4,findfirst(N.== 121)], RES[t3+1:t4,findfirst(N.== 122)] = zeros(length(up_temp)-1), zeros(length(up_temp)-1), zeros(length(up_temp)-1)
RES[t3+1:t4,findfirst(N.== 120)], RES[t3+1:t4,findfirst(N.== 123)], RES[t3+1:t4,findfirst(N.== 124)] = down_temp_small.*GENCAP[findfirst(G .== all_g_in_n[120])], down_temp_small.*GENCAP[findfirst(G .== all_g_in_n[123])], down_temp_small.*GENCAP[findfirst(G .== all_g_in_n[124])]

 # Merge zone 1 and 4 because they are one island
 z14 = vcat(l_in_z[1], l_in_z[4])
 z14 = unique(z14)
 delete!(a.ext[:parameters][:l_in_z], 1)
 delete!(a.ext[:parameters][:l_in_z], 4)
 delete!(a.ext[:parameters][:l_in_z], 5)
 a.ext[:parameters][:l_in_z][1] = z14
 l_in_z = a.ext[:parameters][:l_in_z]

  # Merge zone 1 and 4 because they are one island
 z14 = vcat(n_in_z[1], n_in_z[4])
 temp1 = n_in_z[2].*1
 temp2 = n_in_z[3].*1
 filter!(e->e∉[119, 120, 121, 122, 123, 124],temp1)
 filter!(e->e∉[119, 120, 121, 122, 123, 124],temp2)
 filter!(e->e∉[119, 120, 121, 122, 123, 124],z14)
 a.ext[:parameters][:n_in_z_island] = Dict()
 n_in_z_island = a.ext[:parameters][:n_in_z_island]
 n_in_z_island[1] = z14
 n_in_z_island[2] = temp1
 n_in_z_island[3] = temp2
 a.ext[:parameters][:n_in_z_island] = n_in_z_island

# Create incidence_zonal matrix, which connects lines with zones
incidence_zonal = zeros(length(L), length(Z))
for l in L

    FromBus = df_line[:,:FromBus][l]
    ToBus = df_line[:,:ToBus][l]

    FromZone = df_node[:,:Zone][findfirst(N .== FromBus)]
    ToZone = df_node[:,:Zone][findfirst(N .== ToBus)]

    if FromZone != ToZone
        incidence_zonal[l,FromZone] =  1
        incidence_zonal[l,ToZone] = -1
    end
    
end
incidence_zonal_DC = zeros(length(L_DC), length(Z))
for l_dc in L_DC

    FromBus = df_DC[:,:FromBus][l_dc]
    ToBus = df_DC[:,:ToBus][l_dc]

    FromZone = df_node[:,:Zone][findfirst(N .== FromBus)]
    ToZone = df_node[:,:Zone][findfirst(N .== ToBus)]

    if FromZone != ToZone
        incidence_zonal_DC[l_dc,FromZone] =  1
        incidence_zonal_DC[l_dc,ToZone] = -1
    end
    
end

for l in L
    if sum(abs(incidence_zonal[l,z]) for z in Z) == 2
        @info "AC-line ",l, "is a cross-border line"
    end 
end
for l_dc in L_DC
    if sum(abs(incidence_zonal_DC[l_dc,z]) for z in Z) == 2
        @info "DC-line ",l_dc, "is a cross-border line"
    end 
end

# Build Nodal PTDF
MWBase = 380^2
slack_node1 = 1
slack_position1 = findfirst(n_in_z_island[1] .== slack_node1)
slack_node2 = 39
slack_position2 = findfirst(n_in_z_island[2] .== slack_node2)
slack_node3 = 24
slack_position3 = findfirst(n_in_z_island[3] .== slack_node3)

inc = a.ext[:parameters][:inc] = Matrix(incidence)
inc_dc = a.ext[:parameters][:inc_dc] = Matrix(incidence_dc)
i = a.ext[:parameters][:i] = Matrix(incidence_zonal)
i_DC = a.ext[:parameters][:i_DC] = Matrix(incidence_zonal_DC)

# line_sus_mat = Matrix(susceptance)./MWBase*Matrix(incidence)
# node_sus_mat = transpose(Matrix(incidence))*Matrix(susceptance)./MWBase*Matrix(incidence)
# line_sus_mat_ = line_sus_mat[:, 1:end .!= slack_position]
# node_sus_mat_ = node_sus_mat[1:end .!= slack_position, 1:end .!= slack_position]
# PTDF = line_sus_mat_*inv(node_sus_mat_) #PTDF matrix without slack node
# zero_column = zeros(Float64, length(L), 1)
# PTDF = hcat(PTDF[:,1:(slack_position-1)], zero_column, PTDF[:,slack_position:end]) #PTDF matrix with slack node (zero column)


incidence1, susceptance1 = zeros(length(l_in_z[1]),length(n_in_z_island[1])), zeros(length(l_in_z[1]),length(l_in_z[1]))
incidence2, susceptance2 = zeros(length(l_in_z[2]),length(n_in_z_island[2])), zeros(length(l_in_z[2]),length(l_in_z[2]))
incidence3, susceptance3 = zeros(length(l_in_z[3]),length(n_in_z_island[3])), zeros(length(l_in_z[3]),length(l_in_z[3]))
for l in df_line[:,:BranchID]
    if l in l_in_z[1]
        susceptance1[findfirst(l_in_z[1] .== l),findfirst(l_in_z[1] .== l)] = Matrix(susceptance)[findfirst(df_line[:,:BranchID] .== l),findfirst(df_line[:,:BranchID] .== l)]
    end    
    for n in N
        if l in l_in_z[1] && n in n_in_z_island[1]
            incidence1[findfirst(l_in_z[1] .== l),findfirst(n_in_z_island[1] .== n)] = inc[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)]
        end  
    end  
end    
for l in df_line[:,:BranchID]
    if l in l_in_z[2]
        susceptance2[findfirst(l_in_z[2] .== l),findfirst(l_in_z[2] .== l)] = Matrix(susceptance)[findfirst(df_line[:,:BranchID] .== l),findfirst(df_line[:,:BranchID] .== l)]
    end    
    for n in N
        if l in l_in_z[2] && n in n_in_z_island[2]
            incidence2[findfirst(l_in_z[2] .== l),findfirst(n_in_z_island[2] .== n)] = inc[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)]
        end  
    end  
end   
for l in df_line[:,:BranchID]
    if l in l_in_z[3]
        susceptance3[findfirst(l_in_z[3] .== l),findfirst(l_in_z[3] .== l)] = Matrix(susceptance)[findfirst(df_line[:,:BranchID] .== l),findfirst(df_line[:,:BranchID] .== l)]
    end    
    for n in N
        if l in l_in_z[3] && n in n_in_z_island[3]
            incidence3[findfirst(l_in_z[3] .== l),findfirst(n_in_z_island[3] .== n)] = inc[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)]
        end  
    end  
end     

line_sus_mat1 = Matrix(susceptance1)./MWBase*Matrix(incidence1)
node_sus_mat1 = transpose(Matrix(incidence1))*Matrix(susceptance1)./MWBase*Matrix(incidence1)
line_sus_mat_1 = line_sus_mat1[:, 1:end .!= slack_position1]
node_sus_mat_1 = node_sus_mat1[1:end .!= slack_position1, 1:end .!= slack_position1]
PTDF1 = line_sus_mat_1*inv(node_sus_mat_1) #PTDF matrix without slack node
zero_column1 = zeros(Float64, length(l_in_z[1]), 1)
PTDF1 = hcat(PTDF1[:,1:(slack_position1-1)], zero_column1, PTDF1[:,slack_position1:end]) #PTDF matrix with slack node (zero column)

line_sus_mat2 = Matrix(susceptance2)./MWBase*Matrix(incidence2)
node_sus_mat2 = transpose(Matrix(incidence2))*Matrix(susceptance2)./MWBase*Matrix(incidence2)
line_sus_mat_2 = line_sus_mat2[:, 1:end .!= slack_position2]
node_sus_mat_2 = node_sus_mat2[1:end .!= slack_position2, 1:end .!= slack_position2]
PTDF2 = line_sus_mat_2*inv(node_sus_mat_2) #PTDF matrix without slack node
zero_column2 = zeros(Float64, length(l_in_z[2]), 1)
PTDF2 = hcat(PTDF2[:,1:(slack_position2-1)], zero_column2, PTDF2[:,slack_position2:end]) #PTDF matrix with slack node (zero column)

line_sus_mat3 = Matrix(susceptance3)./MWBase*Matrix(incidence3)
node_sus_mat3 = transpose(Matrix(incidence3))*Matrix(susceptance3)./MWBase*Matrix(incidence3)
line_sus_mat_3 = line_sus_mat3[:, 1:end .!= slack_position3]
node_sus_mat_3 = node_sus_mat3[1:end .!= slack_position3, 1:end .!= slack_position3]
PTDF3 = line_sus_mat_3*inv(node_sus_mat_3) #PTDF matrix without slack node
zero_column3 = zeros(Float64, length(l_in_z[3]), 1)
PTDF3 = hcat(PTDF3[:,1:(slack_position3-1)], zero_column3, PTDF3[:,slack_position3:end]) #PTDF matrix with slack node (zero column)

NPTDF = zeros(length(L),length(N))
for l in l_in_z[1]
    for n in n_in_z_island[1]
        NPTDF[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)] = PTDF1[findfirst(l_in_z[1] .== l),findfirst(n_in_z_island[1] .== n)]
    end
end
for l in l_in_z[2]
    for n in n_in_z_island[2]
        NPTDF[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)] = PTDF2[findfirst(l_in_z[2] .== l),findfirst(n_in_z_island[2] .== n)]
    end
end
for l in l_in_z[3]
    for n in n_in_z_island[3]
        NPTDF[findfirst(df_line[:,:BranchID] .== l),findfirst(N .== n)] = PTDF3[findfirst(l_in_z[3] .== l),findfirst(n_in_z_island[3] .== n)]
    end
end

NPTDF = a.ext[:parameters][:NPTDF] = transpose(NPTDF)

# # Array that indicates whether DC line is linking the AC area and the DC area
# global border_ACDC = zeros(length(L_DC))
# for l_dc in L_DC
#     if abs(df_DC[l_dc,:ToDirection]) == 1 | abs(df_DC[l_dc,:FromDirection]) == 1
#         global border_ACDC[l_dc] = 1
#     end
# end
# a.ext[:parameters][:border_ACDC] = border_ACDC


# # Add columns for visualisation of grid
# totalline = zeros(length(L)+length(L_DC),3)
# totalline[1:length(L),1] = df_line[:,:BranchID]
# totalline[1:length(L),2] = df_line[:,:FromBus]
# totalline[1:length(L),3] = df_line[:,:ToBus]
# totalline[length(L)+1:length(L)+length(L_DC),1] = df_DC[:,:BranchID]
# totalline[length(L)+1:length(L)+length(L_DC),2] = df_DC[:,:FromBus]
# totalline[length(L)+1:length(L)+length(L_DC),3] = df_DC[:,:ToBus]
# x_df = zeros(2)
# y_df = zeros(2)
# bus_df = zeros(2)
# for l in 1:length(L)+length(L_DC)
#     x_df = zeros(2)
#     y_df = zeros(2)
#     bus_df = zeros(2)
#     bus_df[1] = totalline[l,2]
#     bus_df[2] = totalline[l,3]
#     for i in 1:2     
#         x_df[i] = df_node[:,:Latitude][findfirst(N .== Int(bus_df[i]))]
#         y_df[i] = df_node[:,:Longitude][findfirst(N .== Int(bus_df[i]))]
#     end
#     df_network = DataFrame(x = x_df, y = y_df)
#     CSV.write("C://Users//KENISM//OneDrive - VITO//Documents//_Research_2022b_OffShoreBiddingZone//NetworkVisualisation//df_network$l.csv", df_network, delim=';', decimal='.')
# end


@info "Data processed."