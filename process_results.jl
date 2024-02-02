# Optimal values for day-ahead market outcome
GC_DA = a.ext[:parameters][:GC_DA]
RDC_perzone_opt = value.(b.ext[:expressions][:RDC_per_zone])
GC_opt = a.ext[:parameters][:GC_opt] = zeros(length(a.ext[:sets][:T]))
RDC_opt = a.ext[:parameters][:RDC_opt] = zeros(length(a.ext[:sets][:T]))
for t in T
    GC_opt[t] = a.ext[:parameters][:GC_opt][t] = sum(GC_DA[z,t] for z in Z)
    RDC_opt[t] = a.ext[:parameters][:RDC_opt][t] = sum(RDC_perzone_opt[z,t] for z in Z)
end
v_opt = value.(a.ext[:variables][:v])
p_opt = value.(a.ext[:variables][:p])
curt_opt = a.ext[:parameters][:curt_DA] 
F_opt = a.ext[:parameters][:F_DA] 
λ_opt = JuMP.dual.(a.ext[:constraints][:con1])
λ_opt2 = JuMP.dual.(a.ext[:constraints][:con2])
F_DC_opt = a.ext[:parameters][:F_DC_DA]

# Optimal values for redispatch outcome 
global redispatch_needed = zeros(length(a.ext[:sets][:T])) 
global overload = zeros(length(a.ext[:sets][:T]),length(L))
global overload_DC = zeros(length(a.ext[:sets][:T]),length(L_DC))
for l in 1:length(a.ext[:sets][:L])
    for t in 1:length(a.ext[:sets][:T])

    if abs(a.ext[:parameters][:F_DA][l,t]) > a.ext[:parameters][:TC][l] 
        global redispatch_needed[t] = 1
        global overload[t,l] = 1
        @info "There is congestion on an AC-line l.", l
    end 

    end
end 
for l_dc in 1:length(a.ext[:sets][:L_DC])
    for t in 1:length(a.ext[:sets][:T])

    if abs(a.ext[:parameters][:F_DC_DA][l_dc,t]) > a.ext[:parameters][:TC_DC][l_dc] 
        global redispatch_needed[t] = 1
        global overload_DC[t,l_dc] = 1
        @info "There is congestion on a DC-line l_dc.", l_dc
    end 

    end
end 

UP_opt = value.(b.ext[:variables][:UP])
DOWN_opt = value.(b.ext[:variables][:DOWN])
F_beforeRD_opt = value.(b.ext[:expressions][:F_beforeRD])
F_afterRD_opt = value.(b.ext[:expressions][:F_afterRD])
curt_delta_opt = value.(b.ext[:variables][:curt_delta])
F_DC_delta_opt = value.(b.ext[:variables][:F_DC_delta])
TOT_opt = zeros(length(T))

for t in 1:length(a.ext[:sets][:T])
    if redispatch_needed[t] == 1
    else
        UP_opt[:,t] = zeros(length(G))
        DOWN_opt[:,t] = zeros(length(G))
        F_afterRD_opt[:,t] = F_beforeRD_opt[:,t]
        curt_delta_opt[:,t] = zeros(length(N))
        F_DC_delta_opt[:,t] = zeros(length(L_DC))
        RDC_opt[t] = 0.0
    end 

TOT_opt[t] = GC_opt[t] + RDC_opt[t]
end

# global CR_opt = zeros(length(L),length(a.ext[:sets][:T])) 
# global CR_DC_opt = zeros(length(L_DC),length(a.ext[:sets][:T]))
# for l in 1:length(a.ext[:sets][:L])
#     for t in 1:length(a.ext[:sets][:T])

#         if abs(a.ext[:parameters][:F_DA][l,t]) >= a.ext[:parameters][:TC][l] 
#         global CR_opt[l,t] = abs(λ_opt[df_node[findfirst(N .== df_line[l,:FromBus]),:Zone],t] - λ_opt[df_node[findfirst(N .== df_line[l,:ToBus]),:Zone],t]) * abs(a.ext[:parameters][:F_DA][l,t])
#         @info "There is congestion rent on AC-line.", l
#     end 

#     end
# end 
# for l_dc in 1:length(a.ext[:sets][:L_DC])
#     for t in 1:length(a.ext[:sets][:T])

#     if abs(a.ext[:parameters][:F_DC_DA][l_dc,t]) >= a.ext[:parameters][:TC_DC][l_dc] 
#         global CR_DC_opt[l_dc,t] = abs(λ_opt[df_node[findfirst(N .== df_DC[l_dc,:FromBus]),:Zone],t] - λ_opt[df_node[findfirst(N .== df_DC[l_dc,:ToBus]),:Zone],t]) * abs(a.ext[:parameters][:F_DC_DA][l_dc,t])
#         @info "There is congestion rent on DC-line.", l_dc
#     end 

#     end
# end

global CR_opt = zeros(length(L),length(a.ext[:sets][:T])) 
global CR_DC_opt = zeros(length(L_DC),length(a.ext[:sets][:T]))
# for l in 1:length(a.ext[:sets][:L])
#     for t in 1:length(a.ext[:sets][:T])

#         if abs(a.ext[:parameters][:F_DA][l,t]) >= a.ext[:parameters][:TC][l]             
#         global CR_opt[l,t] = abs(λ_opt[df_node[:,:Zone][N .== df_line[l,:FromBus]][1],t] - λ_opt[df_node[:,:Zone][N .== df_line[l,:ToBus]][1],t]) * abs(a.ext[:parameters][:F_DA][l,t])
#         @info "There is congestion rent on AC-line.", l
#     end 

#     end
# end 
# for l_dc in 1:length(a.ext[:sets][:L_DC])
#     for t in 1:length(a.ext[:sets][:T])

#     if abs(a.ext[:parameters][:F_DC_DA][l_dc,t]) >= a.ext[:parameters][:TC_DC][l_dc] 
#         global CR_DC_opt[l_dc,t] = abs(λ_opt[df_node[:,:Zone][N .== df_DC[l_dc,:FromBus]][1],t] - λ_opt[df_node[:,:Zone][N .== df_DC[l_dc,:ToBus]][1],t]) * abs(a.ext[:parameters][:F_DC_DA][l_dc,t])
#         @info "There is congestion rent on DC-line.", l_dc
#     end 

#     end
# end

cd("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results") 


# Write results to CSV-files
df1 = DataFrame(cost_generation = GC_opt, cost_redispatch = RDC_opt)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Cost_total.csv", df1)
writedlm("Cost_total.txt", ["cost_generation" "cost_redispatch" ; GC_opt RDC_opt])

df2 = DataFrame(GC_DA.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Cost_generation.csv", df2)
writedlm("Cost_generation.txt", GC_DA)

df3 = DataFrame(RDC_perzone_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Cost_redispatch.csv", df3)
writedlm("Cost_redispatch.txt", RDC_perzone_opt)

df4 = DataFrame(λ_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Price.csv", df4)
writedlm("Price.txt", λ_opt)

df5 = DataFrame(v_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Dispatch.csv", df5)
writedlm("Dispatch.txt", v_opt)

df6 = DataFrame(p_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\NetPosition.csv", df6)
writedlm("NetPosition.txt", p_opt)

df7 = DataFrame(UP_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\UpwardRedispatch.csv", df7)
writedlm("UpwardRedispatch.txt", UP_opt)

df8 = DataFrame(DOWN_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\DownwardRedispatch.csv", df8)
writedlm("DownwardRedispatch.txt", DOWN_opt)

df9 = DataFrame(F_beforeRD_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\FlowBeforeRD.csv", df9)
writedlm("FlowBeforeRD.txt", F_beforeRD_opt)

df10 = DataFrame(F_afterRD_opt.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\FlowAfterRD.csv", df10)
writedlm("FlowAfterRD.txt", F_afterRD_opt)

df11 = DataFrame(curt_DA.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Curtailment.csv", df11)
writedlm("Curtailment.txt", curt_DA)

df12 = DataFrame(F_DC_DA.data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\DCflow.csv", df12)
writedlm("DCflow.txt", F_DC_DA)

df13 = DataFrame(b.ext[:variables][:curt_delta].data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Redispatch_curtailment.csv", df13)
writedlm("Redispatch_curtailment.txt", b.ext[:variables][:curt_delta])

df14 = DataFrame(b.ext[:variables][:F_DC_delta].data,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Redispatch_DCflow.csv", df14)
writedlm("Redispatch_DCflow.txt", b.ext[:variables][:F_DC_delta])

df15 = DataFrame(overload,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Overload.csv", df15)
writedlm("Overload.txt", overload)

df16 = DataFrame(overload_DC,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\Overload_DC.csv", df16)
writedlm("Overload_DC.txt", overload_DC)

df17 = DataFrame(i,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\i.csv", df17)
writedlm("i.txt", i)

df18 = DataFrame(i_DC,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\i_DC.csv", df18)
writedlm("i_DC.txt", i_DC)

df19 = DataFrame(CR_opt,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\CongestionRent.csv", df19)
writedlm("CongestionRent.txt", CR_opt)

df20 = DataFrame(CR_DC_opt,:auto)
CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022a_OffShoreBiddingZone\\_Models\\_Results\\CongestionRentDC.csv", df20)
writedlm("CongestionRent_DC.txt", CR_DC_opt)


# for l in 1:length(a.ext[:sets][:L])

#     if abs(a.ext[:parameters][:F_DA][l,2]) >= a.ext[:parameters][:TC][l] 
#         @info "There is congestion on an AC-line l.", l
#     end 

# end 
# for l_dc in 1:length(a.ext[:sets][:L_DC])

#     if abs(a.ext[:parameters][:F_DC_DA][l_dc,2]) >= a.ext[:parameters][:TC_DC][l_dc] 
#         @info "There is congestion on a DC-line l_dc.", l_dc
#     end 

# end 

# RES_wind = zeros(length(T),6)
# for t in T
#     RES_wind[t,1] = RES[t,findfirst(N .== 119)]
#     RES_wind[t,2] = RES[t,findfirst(N .== 120)] 
#     RES_wind[t,3] = RES[t,findfirst(N .== 121)] 
#     RES_wind[t,4] = RES[t,findfirst(N .== 122)] 
#     RES_wind[t,5] = RES[t,findfirst(N .== 123)] 
#     RES_wind[t,6] = RES[t,findfirst(N .== 124)] 
# end
# df21 = DataFrame(RES_wind,:auto)
# CSV.write("C:\\Users\\KENISM\\OneDrive - VITO\\Documents\\_Research_2022b_OffShoreBiddingZone\\_Models\\_Results\\RES_wind.csv", df21)
