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
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]

    # Create variables
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    vbar = a.ext[:variables][:vbar] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatchbar")
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    curtbar = a.ext[:variables][:curtbar] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailmentbar")
    p = a.ext[:variables][:p] = @variable(a, [z in Z, t in T], base_name = "position")
    pbar = a.ext[:variables][:pbar] = @variable(a, [z in Z, t in T], base_name = "positionbar")

    # Create expressions
    F = a.ext[:expressions][:F] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) -DEM[t,findfirst(N .== n)]) for n in N ))
    Fbar = a.ext[:expressions][:Fbar] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(vbar[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curtbar[n,t]) -DEM[t,findfirst(N .== n)]) for n in N ))
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )

    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N))

    # Constraints 
    a.ext[:constraints][:con1] = @constraint(a, con1[z in Z, t in T], -sum(- sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) - (RES[t,findfirst(N .== n)]-curt[n,t]) + DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )
    a.ext[:constraints][:con2] = @constraint(a, con2[t in T], sum(p[z,t] for z in Z) == 0 )
    a.ext[:constraints][:con3] = @constraint(a, con3[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )

    a.ext[:constraints][:con4] = @constraint(a, con4[z in Z, t in T], sum(sum(GENCAP[g] * vbar[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curtbar[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )
    a.ext[:constraints][:con5] = @constraint(a, con5[t in T], sum(pbar[z,t] for z in Z) == 0 )
    a.ext[:constraints][:con6] = @constraint(a, con6[l in L, t in T], Fbar[l,t] + TC[l] >= 0 )
    a.ext[:constraints][:con7] = @constraint(a, con7[l in L, t in T], TC[l] - Fbar[l,t] >= 0 )
    a.ext[:constraints][:con8] = @constraint(a, con8[n in N, t in T], curtbar[n,t] <= RES[t,findfirst(N .== n)] )

    return a
end

# Build your model
market_clearing!(a)
status = optimize!(a)

GC_DA = a.ext[:parameters][:GC_DA] = value.(a.ext[:expressions][:GC_per_zone])
v_DA = a.ext[:parameters][:v_DA] = value.(a.ext[:variables][:v])
curt_DA = a.ext[:parameters][:curt_DA] = value.(a.ext[:variables][:curt])
p_DA = a.ext[:parameters][:p_DA] = value.(a.ext[:variables][:p])
F_DA = a.ext[:parameters][:F_DA] = value.(a.ext[:expressions][:F])