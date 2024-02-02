b = Model(Gurobi.Optimizer)

function redispatch!(a::Model)

    b.ext[:variables] = Dict()
    b.ext[:expressions] = Dict()
    b.ext[:constraints] = Dict()
    b.ext[:objective] = Dict()

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
    α = a.ext[:parameters][:α]
    v_DA = a.ext[:parameters][:v_DA] 
    curt_DA = a.ext[:parameters][:curt_DA]
    p_DA = a.ext[:parameters][:p_DA] 
    F_DA = a.ext[:parameters][:F_DA] 
    F_DC_DA = a.ext[:parameters][:F_DC_DA] 
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]
    GC_DA = a.ext[:parameters][:GC_DA]

    # Create variables
    UP = b.ext[:variables][:UP] = @variable(b, [g in G, t in T], lower_bound = 0, base_name = "upward redispatch")
    DOWN = b.ext[:variables][:DOWN] = @variable(b, [g in G, t in T], lower_bound = 0, base_name = "downward redispatch")
    curt_delta = b.ext[:variables][:curt_delta] = @variable(b, [n in N, t in T], lower_bound = 0, base_name = "curtailment") 
    F_DC_delta = b.ext[:variables][:F_DC_delta] = @variable(b, [l_dc in L_DC, t in T], base_name = "DC flow adjustment")

    # Create expressions
    F_beforeRD = b.ext[:expressions][:F_beforeRD] = @expression(b, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v_DA[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt_DA[n,t]) -DEM[t,findfirst(N .== n)] - sum(F_DC_DA[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    F_afterRD = b.ext[:expressions][:F_afterRD] = @expression(b, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v_DA[g,t]+UP[g,t]-DOWN[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt_DA[n,t]-curt_delta[n,t]) -DEM[t,findfirst(N .== n)] - sum((F_DC_DA[l_dc,t]+F_DC_delta[l_dc,t])*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    RDC_per_zone = b.ext[:expressions][:RDC_per_zone] = @expression(b, [z in Z, t in T], sum(sum((1+α)*UP[g,t]*GENCAP[g]*MC[g]-(1-α)*DOWN[g,t]*GENCAP[g]*MC[g] for g in g_in_n[n]) + CC*curt_delta[n,t] for n in n_in_z[z] )) 

    # Objective
    RDC = b.ext[:objective][:RDC] = @objective(b, Min, sum( sum((1+α)*UP[g,t]*GENCAP[g]*MC[g]-(1-α)*DOWN[g,t]*GENCAP[g]*MC[g] for g in g_in_n[n]) + CC*curt_delta[n,t] for t in T for n in N))

    # Constraints 
    b.ext[:constraints][:con1] = @constraint(b, con1[t in T], sum(sum(UP[g,t]*GENCAP[g] for g in g_in_n[n]) for n in N) == sum(sum(DOWN[g,t]*GENCAP[g] for g in g_in_n[n]) + curt_delta[n,t] for n in N))
    b.ext[:constraints][:con2] = @constraint(b, con2[l in L, t in T], F_afterRD[l,t] + TC[l] >= 0 )
    b.ext[:constraints][:con3] = @constraint(b, con3[l in L, t in T], TC[l] - F_afterRD[l,t] >= 0 )
    b.ext[:constraints][:con4] = @constraint(b, con4[g in G, t in T], (1-v_DA[g,t]) - UP[g,t] >= 0 )
    b.ext[:constraints][:con5] = @constraint(b, con5[g in G, t in T], v_DA[g,t] - DOWN[g,t] >= 0 )
    b.ext[:constraints][:con6] = @constraint(b, con6[l_dc in L_DC, t in T], TC_DC[l_dc] >= (F_DC_DA[l_dc,t]+F_DC_delta[l_dc,t]) >= -TC_DC[l_dc] )
    b.ext[:constraints][:con7] = @constraint(b, con7[l_dc in L_DC, t in T], - F_DC_DA[l_dc,t] <= F_DC_delta[l_dc,t] <= TC_DC[l_dc] - F_DC_DA[l_dc,t])
    b.ext[:constraints][:con8] = @constraint(b, con8[n in N, t in T], -curt_DA[n,t] <= curt_delta[n,t] <= RES[t,findfirst(N .== n)]-curt_DA[n,t] )

    # Laws of Kirchoff in DC-area
    b.ext[:constraints][:con9] = @constraint(b, con9[t in T], -F_DC_DA[1,t]-F_DC_delta[1,t]-(RES[t,findfirst(N .== 122)]-curt_DA[122,t]-curt_delta[122,t]) == 0)
    b.ext[:constraints][:con10] = @constraint(b, con10[t in T], -F_DC_DA[2,t]-F_DC_delta[2,t]-(RES[t,findfirst(N .== 119)]-curt_DA[119,t]-curt_delta[119,t]) +F_DC_DA[8,t] +F_DC_DA[7,t] +F_DC_delta[8,t] +F_DC_delta[7,t] == 0)
    b.ext[:constraints][:con11] = @constraint(b, con11[t in T], -F_DC_DA[5,t]-F_DC_delta[5,t]-(RES[t,findfirst(N .== 120)]-curt_DA[120,t]-curt_delta[120,t]) -F_DC_DA[7,t] +F_DC_DA[9,t] -F_DC_delta[7,t] +F_DC_delta[9,t] == 0)
    b.ext[:constraints][:con12] = @constraint(b, con12[t in T], -F_DC_DA[6,t]-F_DC_delta[6,t]-(RES[t,findfirst(N .== 124)]-curt_DA[124,t]-curt_delta[124,t]) -F_DC_DA[8,t] -F_DC_DA[9,t] -F_DC_delta[8,t] -F_DC_delta[9,t] == 0)
    b.ext[:constraints][:con13] = @constraint(b, con13[t in T], -F_DC_DA[3,t]-F_DC_delta[3,t]-(RES[t,findfirst(N .== 121)]-curt_DA[121,t]-curt_delta[121,t]) +F_DC_DA[10,t] +F_DC_delta[10,t] == 0)
    b.ext[:constraints][:con14] = @constraint(b, con14[t in T], -F_DC_DA[4,t]-F_DC_delta[4,t]-(RES[t,findfirst(N .== 123)]-curt_DA[123,t]-curt_delta[123,t]) -F_DC_DA[10,t] -F_DC_delta[10,t] == 0)


    return status
end

# Build your model
redispatch!(a)
@info "Redispatch model built."
status = optimize!(b)
@info "Redispatch model optimised."
