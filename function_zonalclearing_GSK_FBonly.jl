m = Model(Gurobi.Optimizer)

function market_clearing!(a::Model)
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()
    m.ext[:objective] = Dict()

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
    v = m.ext[:variables][:v] = @variable(m, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = m.ext[:variables][:curt] = @variable(m, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = m.ext[:variables][:p] = @variable(m, [n in N, t in T], base_name = "position")
    F = m.ext[:variables][:F] = @variable(m, [l in L, t in T], base_name = "flows")
    F_dc = m.ext[:variables][:F_dc] = @variable(m, [l_dc in L_DC, t in T], base_name = "DC flows")
    
    # Create expressions
    GC_per_time_and_zone = m.ext[:expressions][:GC_per_zone] = @expression(m, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
           
    # Objective
    GC = m.ext[:objective][:GC] = @objective(m, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N))

    # Constraints 
    m.ext[:constraints][:con1] = @constraint(m, con1[n in N, t in T], sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] - p[n,t] == 0 )
    m.ext[:constraints][:con2] = @constraint(m, con2[n in N, t in T], sum(F[l,t]*inc[l,findfirst(N .== n)] for l in L) + sum(F_dc[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC) - p[n,t] == 0 )
    m.ext[:constraints][:con3] = @constraint(m, con3[l in L, t in T], F[l,t] == sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t])-DEM[t,findfirst(N .== n)] - sum(F_dc[l_dc,t]*inc_dc[l_dc,findfirst(N .== n)] for l_dc in L_DC)) for n in N ))
    m.ext[:constraints][:con4] = @constraint(m, con4[l in L, t in T], -TC[l] <= F[l,t] <= TC[l] )
    m.ext[:constraints][:con5] = @constraint(m, con5[l_dc in L_DC, t in T], -TC_DC[l_dc] <= F[l_dc,t] <= TC_DC[l_dc] )
    m.ext[:constraints][:con6] = @constraint(m, con6[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )

    return m
end

# Build your model
market_clearing!(a)
status = optimize!(m)

v_bc = a.ext[:parameters][:v_bc] = value.(m.ext[:variables][:v])
curt_bc = a.ext[:parameters][:curt_bc] = value.(m.ext[:variables][:curt])
p_bc = a.ext[:parameters][:p_bc] = value.(m.ext[:variables][:p])
F_bc = a.ext[:parameters][:F_bc] = value.(m.ext[:variables][:F])

function get_GSK!(a::Model)
    a.ext[:variables] = Dict()
    a.ext[:expressions] = Dict()
    a.ext[:constraints] = Dict()
    a.ext[:objective] = Dict()

    # Extract sets
    N = a.ext[:sets][:N]
    G = a.ext[:sets][:G]
    Z = a.ext[:sets][:Z]
    L = a.ext[:sets][:L]
    T = a.ext[:sets][:T]

    # Extract parameters
    GENCAP = a.ext[:parameters][:GENCAP]
    n_in_z = a.ext[:parameters][:n_in_z]
    g_in_n = a.ext[:parameters][:g_in_n]
    all_g_in_n = a.ext[:parameters][:all_g_in_n]

    cap, gsk = zeros(length(Z)), zeros(length(N),length(Z))
    for z in Z
        cap[z] = sum(sum(GENCAP[g] for g in all_g_in_n[n]) for n in n_in_z[z])
        for n in n_in_z[z]
            gsk[findfirst(N .== n),z] = sum(GENCAP[g] for g in all_g_in_n[n])/cap[z]
        end
    end 
    
    a.ext[:parameters][:GSK] = gsk

    return a 
end

# Build your model
get_GSK!(a)

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
    v_bc = a.ext[:parameters][:v_bc]
    p_bc = a.ext[:parameters][:p_bc]
    F_bc = a.ext[:parameters][:F_bc]
    GSK = a.ext[:parameters][:GSK]
    i = a.ext[:parameters][:i]
    RES = a.ext[:parameters][:RES]
    CC = a.ext[:parameters][:CC]

    # Create variables
    v = a.ext[:variables][:v] = @variable(a, [g in G, t in T], lower_bound = 0, upper_bound = 1, base_name = "dispatch")
    curt = a.ext[:variables][:curt] = @variable(a, [n in N, t in T], lower_bound = 0, base_name = "curtailment")
    p = a.ext[:variables][:p] = @variable(a, [z in Z, t in T], base_name = "position")
    F_FBMC = a.ext[:variables][:F_FBMC] = @variable(a, [l in L, t in T], base_name = "commercial flow")

    # Create expressions
    F = a.ext[:expressions][:F] = @expression(a, [l in L, t in T], sum(NPTDF[findfirst(N .== n),l]*(sum(GENCAP[g]*(v[g,t]) for g in g_in_n[n])+ (RES[t,findfirst(N .== n)]-curt[n,t])-DEM[t,findfirst(N .== n)]) for n in N ))
    ZPTDF = a.ext[:expressions][:ZPTDF] = @expression(a, [z in Z, l in L], sum(NPTDF[findfirst(N .== n),l]*GSK[findfirst(N .== n),z] for n in N )    )
    RAM_pos = a.ext[:expressions][:RAM_pos] = @expression(a, [l in L, t in T], TC[l] - (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )) )
    RAM_neg = a.ext[:expressions][:RAM_neg] = @expression(a, [l in L, t in T], TC[l] + (F_bc[l,t] - sum(ZPTDF[z,l] * p_bc[z,t] for z in Z )) )
    GC_per_time_and_zone = a.ext[:expressions][:GC_per_zone] = @expression(a, [z in Z, t in T], sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) + CC*curt[n,t] for n in n_in_z[z]) )
    
    # Objective
    GC = a.ext[:objective][:GC] = @objective(a, Min, sum(sum(MC[g]*GENCAP[g]*v[g,t] for g in g_in_n[n]) +CC*curt[n,t] for t in T for n in N))

    # Constraints 
    a.ext[:constraints][:con1] = @constraint(a, con1[z in Z, t in T], sum(sum(GENCAP[g] * v[g,t] for g in g_in_n[n]) + (RES[t,findfirst(N .== n)]-curt[n,t]) - DEM[t,findfirst(N .== n)] for n in n_in_z[z]) - p[z,t] == 0 )
    a.ext[:constraints][:con2] = @constraint(a, con2[z in Z, t in T], sum(F_FBMC[l,t]*i[l,z] for l in L) - p[z,t] == 0 )
    a.ext[:constraints][:con3] = @constraint(a, con3[l in L, t in T], F_FBMC[l,t] + RAM_neg[l,t] >= 0 )
    a.ext[:constraints][:con4] = @constraint(a, con4[l in L, t in T], RAM_pos[l,t] - F_FBMC[l,t] >= 0 )
    a.ext[:constraints][:con5] = @constraint(a, con5[l in L, t in T], F_FBMC[l,t] == sum(ZPTDF[z,l]*p[z,t] for z in Z))
    a.ext[:constraints][:con6] = @constraint(a, con6[n in N, t in T], curt[n,t] <= RES[t,findfirst(N .== n)] )

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