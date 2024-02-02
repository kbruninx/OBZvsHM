m = Model(Gurobi.Optimizer)

function market_clearing!(a::Model)
    a.ext[:variables] = Dict()
    a.ext[:expressions] = Dict()
    a.ext[:constraints] = Dict()
    a.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    L_DC = a.ext[:sets][:L_DC]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]
    MC = a.ext[:parameters][:MC]
    DEM = a.ext[:parameters][:DEM]
    NPTDF = a.ext[:parameters][:NPTDF]
    n_in_z = a.ext[:parameters][:n_in_z]
    g_in_n = a.ext[:parameters][:g_in_n]
    TC = a.ext[:parameters][:TC]
    TC_DC = a.ext[:parameters][:TC_DC]
    inc = a.ext[:parameters][:inc] 
    inc_dc = a.ext[:parameters][:inc_dc]
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]

    # Create variables
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = a.ext[:variables][:p] = @variable(a, [n in N, t in T], base_name = "position")
    F = a.ext[:variables][:F] = @variable(a, [l in L, t in T], base_name = "flows")
    F_dc = a.ext[:variables][:F_dc] = @variable(a, [l_dc in L_DC, t in T], base_name = "DC flows")

    # Create expressions
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
        
    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for t in T for n in N))

    # Constraints 
    a.ext[:constraints][:con1] = @constraint(a, con1[n in N, t in T], sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] - p[n,t] == 0 )
    a.ext[:constraints][:con2] = @constraint(a, con2[n in N, t in T], sum(F[l,t]*inc[l,findfirst(N .== n)] for l in L) + sum(F_dc[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC) - p[n,t] == 0 )
    a.ext[:constraints][:con3] = @constraint(a, con3[l in L, t in T], F[l,t] == sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) -DEM[t,findfirst(N .== n)] - sum(F_dc[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    a.ext[:constraints][:con4] = @constraint(a, con4[l in L, t in T], -TC[l] <= F[l,t] <= TC[l] )
    a.ext[:constraints][:con5] = @constraint(a, con5[l_dc in L_DC, t in T], -TC_DC[l_dc] <= F_dc[l_dc,t] <= TC_DC[l_dc] )
    a.ext[:constraints][:con6] = @constraint(a, con6[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )

    #enter laws of Kirchoff
    a.ext[:constraints][:con9] = @constraint(a, con9[t in T], -F_dc[1,t]-(RES[t,findfirst(N .== 122)]-curt[122,t]) == 0)
    a.ext[:constraints][:con10] = @constraint(a, con10[t in T], -F_dc[2,t]-(RES[t,findfirst(N .== 119)]-curt[119,t]) +F_dc[8,t] +F_dc[7,t] == 0)
    a.ext[:constraints][:con11] = @constraint(a, con11[t in T], -F_dc[5,t]-(RES[t,findfirst(N .== 120)]-curt[120,t]) -F_dc[7,t] +F_dc[9,t] == 0)
    a.ext[:constraints][:con12] = @constraint(a, con12[t in T], -F_dc[6,t]-(RES[t,findfirst(N .== 124)]-curt[124,t]) -F_dc[8,t] -F_dc[9,t] == 0)
    a.ext[:constraints][:con13] = @constraint(a, con13[t in T], -F_dc[3,t]-(RES[t,findfirst(N .== 121)]-curt[121,t]) +F_dc[10,t] == 0)
    a.ext[:constraints][:con14] = @constraint(a, con14[t in T], -F_dc[4,t]-(RES[t,findfirst(N .== 123)]-curt[123,t]) -F_dc[10,t] == 0)
    
    return a
end

# Build your model
market_clearing!(a)
status = optimize!(a)

GC_DA = a.ext[:parameters][:GC_DA] = value.(a.ext[:expressions][:GC_per_zone])
v_DA = a.ext[:parameters][:v_DA] = value.(a.ext[:variables][:v])
curt_DA = a.ext[:parameters][:curt_DA] = value.(a.ext[:variables][:curt])
p_DA = a.ext[:parameters][:p_DA] = value.(a.ext[:variables][:p])
F_DA = a.ext[:parameters][:F_DA] = value.(a.ext[:variables][:F])
F_DC_DA = a.ext[:parameters][:F_DC_DA] = value.(a.ext[:variables][:F_dc])