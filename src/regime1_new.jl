
module regime1
using Markdown, OrdinaryDiffEq, DiffEqBase, Plots, StochasticDiffEq, DifferentialEquations, ModelingToolkit, NonlinearSolve
@variables NH4_aq Na_aq H_aq Cl_aq NO3_aq SO4_aq HNO3_aq NH3_aq HCl_aq HSO4_aq Ca_aq K_aq Mg_aq
@variables P_NH3 
# unit: atm

@variables NH4HSO4_s NaHSO4_s CaSO4_s KHSO4_s  
# unit: mol/m3(air)

@variables I W_
@parameters T Na_i H2SO4_i NH3_i HNO3_i HCl_i Ca_i K_i Mg_i aw # m_all

begin
    K1 = 6.067*10^5
    K2 = 7.974*10^11
    K3 = 4.319*10^-5
    K4 = 1.569*10^-2
    K5 = 24.016
    K6 = 0.872
    K7 = 8.680
    K8 = 1.079*10^5
    K9 = 2.507*10^15
    K10 = 9.557*10^21
    K11 = 1.015*10^-2
    K12 = 5.764*10^1
    K13 = 1.805*10^-5
    K14 = 2.511*10^6
    K15 = 2.1*10^5
    K16 = 1.971*10^6
    K17 = 2.5*10^3
    K18 = 1.01.*10^-14
    K19 = 4.799*10^-1
    K20 = 1.817
    K21 = 1.086*10^-16
    K22 = 1.197*10
    K23 = 3.766*10
    K24 = 2.413*10^4
    K25 = 4.199*10^-17
    K26 = 1.383
    K27 = 2.972*10
    
    H4 = -9.589
    H5 = -8.423
    H6 = 14.075
    H7 = -6.167
    H8 = 36.798
    H9 = -8.754
    H10 = -1.347
    H11 = 8.85
    H12 = 13.79
    H13 = -1.5
    H14 = 29.17
    H15 = 29.17
    H16 = 30.20
    H17 = 30.20
    H18 = -22.52
    H19 = 0.98
    H20 = -2.65
    H21 = -71
    H22 = -8.22
    H23 = -1.56
    H24 = 0.79
    H25 = -74.735
    H26 = -2.87
    H27 = -5.19
    
    C4 = 45.807
    C5 = 17.964
    C6 = 19.388
    C7 = 19.953
    C11 = 25.14
    C12 = -5.39
    C13 = 26.92
    C14 = 16.83
    C15 = 16.83
    C16 = 19.91
    C17 = 19.91
    C18 = 26.92
    C19 = 39.75
    C20 = 38.57
    C21 = 2.40
    C22 = 16.01
    C23 = 16.90
    C24 = 14.75
    C25 = 6.025
    C26 = 15.83
    C27 = 54.40
    end

    struct ion
        z::Int
        m::Num
    end
    K = ion(1,K_aq)
    NH4 = ion(1,NH4_aq)
    Na = ion(1,Na_aq)
    H =ion(1,H_aq)
    Ca =ion(2,Ca_aq)
    Mg = ion(1,Mg_aq)
    NO3 = ion(1,NO3_aq)
    SO4 = ion(2,SO4_aq)
    HSO4 = ion(1,HSO4_aq)
    Cl = ion(1,Cl_aq)


    struct salt
    cation::ion
    anion::ion
    name::String
    q::Float64
    vec_same_anion::Vector{}
    vec_same_cation::Vector{}
    end
    KCl = salt(K,Cl,"KCl",0.92,[NH4,Na,H,Ca,Mg],[NO3,SO4]) #HSO4 
    KNO3 = salt(K,NO3,"KNO3",-2.33,[NH4,Na,H,Ca,Mg],[Cl,SO4]) #HSO4
    K2SO4 = salt(K,SO4,"K2SO4",-0.25,[NH4,Na,H,Mg],[NO3,Cl]) #Ca
    NH4Cl = salt(NH4,Cl,"NH4Cl",0.82,[K,Na,H,Ca,Mg],[NO3,SO4]) #NH4HSO4
    NaCl = salt(Na,Cl,"NaCl",2.23,[K,NH4,H,Ca,Mg],[NO3,SO4]) #NaHSO4
    HCl = salt(H,Cl,"HCl",6.00,[K,NH4,Na,Ca,Mg],[NO3,SO4,HSO4])
    CaCl2 = salt(Ca,Cl,"CaCl2",2.4,[K,NH4,Na,H,Mg],[NO3]) #CaSO4
    MgCl2 = salt(Mg,Cl,"MgCl2",2.90,[K,NH4,Na,H,Ca],[NO3,SO4])
    NH4NO3 = salt(NH4,NO3,"NH4NO3",-1.15,[K,Na,H,Ca,Mg],[Cl,SO4]) #NH4NO3
    NaNO3 = salt(Na,NO3,"NaNO3",-0.39,[K,NH4,H,Ca,Mg],[Cl,SO4])
    HNO3 = salt(H,NO3,"HNO3",2.60,[K,NH4,Na,Ca,Mg],[Cl,SO4,HSO4])
    CaNO32 = salt(Ca,NO3,"CaNO32",0.93,[K,NH4,Na,H,Mg],[Cl]) #CaSO4
    MgNO32 = salt(Mg,NO3,"MgNO32",2.32,[K,NH4,Na,H,Ca],[Cl,SO4])
    NH42SO4 = salt(NH4,SO4,"NH42SO4",-0.25,[K,Na,H,Mg],[NO3,Cl]) #NH4HSO4
    Na2SO4 = salt(Na,SO4,"Na2SO4",-0.19,[K,NH4,H,Mg],[NO3,Cl])
    H2SO4 = salt(H,SO4,"H2SO4",-0.1,[K,NH4,Na,Mg],[NO3,Cl,HSO4])
    MgSO4 = salt(Mg,SO4,"MgSO4",0.15,[K,NH4,Na,H],[NO3,Cl])
    HHSO4 = salt(H,HSO4,"HHSO4",8.00,[],[NO3,Cl,SO4])

    q_table_anion = Dict(
        "KCl"=>[NH4Cl.q,NaCl.q,HCl.q,CaCl2.q,MgCl2.q],
        "KNO3"=>[NH4NO3.q,NaNO3.q,HNO3.q,CaNO32.q,MgNO32.q],
        "K2SO4"=>[NH42SO4.q,Na2SO4.q,H2SO4.q,MgSO4.q],
        "NH4Cl"=>[KCl.q,NaCl.q,HCl.q,CaCl2.q,MgCl2.q],
        "NaCl"=>[KCl.q,NH4Cl.q,HCl.q,CaCl2.q,MgCl2.q],
        "HCl"=>[KCl.q,NH4Cl.q,NaCl.q,CaCl2.q,MgCl2.q],
        "CaCl2"=>[KCl.q,NH4Cl.q,NaCl.q,HCl.q,MgCl2.q],
        "MgCl2"=>[KCl.q,NH4Cl.q,NaCl.q,HCl.q,CaCl2.q],
        "NH4NO3"=>[KNO3.q,NaNO3.q,HNO3.q,CaNO32.q,MgNO32.q],
        "NaNO3"=>[KNO3.q,NH4NO3.q,HNO3.q,CaNO32.q,MgNO32.q],
        "HNO3"=>[KNO3.q,NH4NO3.q,NaNO3.q,CaNO32.q,MgNO32.q],
        "CaNO32"=>[KNO3.q,NH4NO3.q,NaNO3.q,HNO3.q,MgNO32.q],
        "MgNO32"=>[KNO3.q,NH4NO3.q,NaNO3.q,HNO3.q,CaNO32.q],
        "NH42SO4"=>[K2SO4.q,Na2SO4.q,H2SO4.q,MgSO4.q],
        "Na2SO4"=>[K2SO4.q,NH42SO4.q,H2SO4.q,MgSO4.q],
        "H2SO4"=>[K2SO4.q,NH42SO4.q,Na2SO4.q,MgSO4.q],
        "MgSO4"=>[K2SO4.q,NH42SO4.q,Na2SO4.q,H2SO4.q],
        "HHSO4"=>[0]
    )

    q_table_cation = Dict(
        "KCl"=>[KNO3.q,K2SO4.q],
        "KNO3"=>[KCl.q,K2SO4.q],
        "K2SO4"=>[KNO3.q,KCl.q],
        "NH4Cl"=>[NH4NO3.q,NH42SO4.q],
        "NaCl"=>[NaNO3.q,Na2SO4.q],
        "HCl"=>[HNO3.q,H2SO4.q,HHSO4.q],
        "CaCl2"=>[CaNO32.q],
        "MgCl2"=>[MgNO32.q,MgSO4.q],
        "NH4NO3"=>[NH4Cl.q,NH42SO4.q],
        "NaNO3"=>[NaCl.q,Na2SO4.q],
        "HNO3"=>[HCl.q,H2SO4.q,HHSO4.q],
        "CaNO32"=>[CaCl2.q],
        "MgNO32"=>[MgCl2.q,MgSO4.q],
        "NH42SO4"=>[NH4NO3.q,NH4Cl.q],
        "Na2SO4"=>[NaNO3.q,NaCl.q],
        "H2SO4"=>[HNO3.q,HCl.q,HHSO4.q],
        "MgSO4"=>[MgNO32.q,MgCl2.q],
        "HHSO4"=>[HNO3.q,HCl.q,H2SO4.q]
    )


    function get_cation_vec(s::salt)
        r1 = []
        r2 = []
        k = s.vec_same_cation
        for i in 1:length(k)
            push!(r1,k[i].z)
            push!(r2,k[i].m)
        end
        return r1,r2
    end
    function get_anion_vec(s::salt)
        r1 = []
        r2 = []
        k = s.vec_same_anion
        for i in 1:length(k)
            push!(r1,k[i].z)
            push!(r2,k[i].m)
        end
        return r1,r2
    end

    get_cation_vec(KCl)
    get_anion_vec(KCl)   
    
    function Equilibriumconst(K₀,T,H_group,c_group) # to calculate equilibrium constant at given temperature
        # input variables：
        # K₀ --> equilibrium constant at 298.15K(Tₒ), T --> given temperature,K, H_group --> ΔH₀/(R*T₀), c_group --> Δcₒ/R 
    
        # Equation: 
        # K = K₀*exp(-ΔH₀/(R*T₀)*(T₀/T-1)-Δcₒ/R*(1+log(T₀/T)-T₀/T))
    
        # K₀: equilibrium constant at 298.15K(Tₒ)
        # T: given temperature, K
        # R: universal gas constant 
        # ΔH₀:  enthalpy change of the reaction at 298.15K(Tₒ)
        # Δc₀: change of molar heat capacity of products minus reactants
    
        K = K₀*exp(-H_group*(298.15/T-1)-c_group*(1+log(298.15/T)-298.15/T))
    end
    
    function cal_I(z_all, m_all) # to compute the ionic strength of the multicomponent solution
        # input variables: 
        # mi --> molalities of ions, zi --> valence of ions
        # I += 1 / 2 * (m_all[i] * z_all[i]^2)
        return 1/2*sum(m_all.*z_all.^2)
    end

    function log_gamma_0(z1, z2, q, I)
        C = 1 + 0.055 * q * exp(-0.023 * I^3)
        Γ_ = exp(-0.5107 * I^0.5 / (1 + C * I^0.5))
        B = 0.75 - 0.065 * q
        Γ0 = (1 + B * (1 + 0.1 * I)^q - B) * Γ_
        log_gamma_12_0 = z1 * z2 * log(Γ0)
        return log_gamma_12_0
    end

    function F(s::salt, I)
        z1_ = get_cation_vec(s)[1]
        m1_ = get_cation_vec(s)[2]
        q1_ = q_table_cation[s.name]
        z2_ = get_anion_vec(s)[1]
        m2_ = get_anion_vec(s)[2]
        q2_ = q_table_anion[s.name]
        z1 = s.cation.z
        z2 = s.anion.z
        m1 = s.cation.m
        m2 = s.anion.m
        q = s.q
    
        Aγ = 0.511 # kg^0.5 mol^-0.5
        n = length(z2_)
        Y_21 = ((z1 + z2) / 2)^2 * m2 / I
        log_γ_12 = log_gamma_0(z1, z2, q, I)
        A_group = Aγ * I^0.5 / (1 + I^0.5)
        F1 = Y_21 * log_γ_12 + A_group * z1 * z2 * Y_21
        for i in 1:n
            Y_i1 = ((z1 + z2_[i]) / 2)^2 * m2_[i] / I
            log_γ_1i = log_gamma_0(z1, z2_[i], q2_[i], I)
            F1 += Y_i1 * log_γ_1i + A_group * z1 * z2_[i] * Y_i1
        end

        m = length(z1_)
        X_12 = ((z1 + z2) / 2)^2 * m1 / I
        F2 = X_12 * log_γ_12 + A_group * z1 * z2 * X_12
        for i in 1:m
            X_i2 = ((z1_[i] + z2) / 2)^2 * m1_[i] / I
            log_γ_i2 = log_gamma_0(z1_[i], z2, q1_[i], I)
            F2 += X_i2 * log_γ_i2 + A_group * z1_[i] * z2 * X_i2
        end
        return F1, F2
    end
    
    F(KCl,2)

    function ActivityCoefficient(s::salt, T, I)
        Aγ = 0.511 # kg^0.5 mol^-0.5
        z1 = s.cation.z
        z2 = s.anion.z
        F1,F2 = F(s,I)
    
        log_γ₁₂_T0 = exp(-Aγ * (z1 * z2 * I^0.5 / (1 + I^0.5)) + z1 * z2 / (z1 + z2) * (F1 / z1 + F2 / z2))
        A = -0.41 * I^0.5 / (1 + I^0.5) + 0.039*I^0.92
        γ₁₂ = exp(log_γ₁₂_T0 * (1.125 - 0.005 * (T - 273.15)) - (0.125 - 0.005 * (T - 273.15)) * A)
        return γ₁₂
    end

    AC_special = Dict(
        "CaSO4"=>0,
        "KHSO4"=>(ActivityCoefficient(HHSO4,T,I)*ActivityCoefficient(KCl,T,I)/ActivityCoefficient(HCl,T,I))^0.5,
        "NH4HSO4"=>(ActivityCoefficient(HHSO4,T,I)*ActivityCoefficient(NH4Cl,T,I)/ActivityCoefficient(HCl,T,I))^0.5,
        "NaHSO4"=>(ActivityCoefficient(HHSO4,T,I)*ActivityCoefficient(NaCl,T,I)/ActivityCoefficient(HCl,T,I))^0.5,
        "NH43HSO42"=>(((ActivityCoefficient(HHSO4,T,I)*ActivityCoefficient(NH4Cl,T,I)/ActivityCoefficient(HCl,T,I))^0.5)*ActivityCoefficient(NH42SO4,T,I)^3)^0.2
    )

    polynomial_fit = [
        [36.356, -165.66, 447.46, -673.55, 510.91, -155.56, 0],
        [20.847, -97.599, 273.220, -422.120, 331.160, -105.450, 0],
        [1.061, -0.101, 1.579 * 10^-2, -1.950 * 10^-3, 9.515 * 10^-5, -1.547 * 10^-6, 0],
        [1061.51, -4748.97, 8096.16, -6166.16, 1757.47, 0, 0],
        [1.2141 * 10^4, -5.1173 * 10^4, 8.1252 * 10^4, -5.7527 * 10^4, 1.5305 * 10^4, 0, 0],
        [179.721, -721.266, 1161.03, -841.479, 221 / 943, 0, 0],
        [-0.778, 177.74, -719.79, 1174.6, -863.44, 232.31, 0],
        [12.166, -16.154, 0, 10.886, 0, -6.815, 0],
        [11.505, -26.518, 34.937, -19.829, 0, 0, 0],
        [0.9988, -2.6947 * 10^-2, 1.9610 * 10^-4, 2.8154 * 10^-5, -6.1359 * 10^-7, 0, 0],
        [1.0614, -0.1014, 1.5796 * 10^-2, -1.9501 * 10^-3, 9.5147 * 10^-5, -1.5473 * 10^-6, 0],
        [1.0084, -4.9390 * 10^-2, 8.888 * 10^-3, -2.1570 * 10^-3, 1.6170 * 10^-4, -1.99 * 10^-6, -1.142 * 10^-7],
        [1.0052, -6.4840 * 10^-2, 3.519 * 10^-2, -1.319 * 10^-2, 1.9250 * 10^-3, -1.224 * 10^-4, 2.87 * 10^-6],
        [0.9968, -2.969 * 10^-2, 1.735 * 10^-5, -3.253 * 10^-4, 3.571 * 10^-5, -9.787 * 10^-7, 0],
        [0.9968, -2.611 * 10^-2, -1.599 * 10^-3, 1.355 * 10^-4, -2.3170 * 10^-6, -1.113 * 10^-8, 0],
        [1.0053, -2.4991 * 10^-2, 4.4688 * 10^-4, 1.6453 * 10^-5, -3.8940 * 10^-7, -4.7668 * 10^-8, 1.3753 * 10^-9],
        [1.0261, -4.9766 * 10^-2, 3.2757 * 10^-3, -2.4477 * 10^-4, 1.0766 * 10^-5, -1.8329 * 10^-7, 0],
        [1.0088, -5.3730 * 10^-2, 1.4201 * 10^-3, -9.2484 * 10^-4, 2.2796 * 10^-4, -1.5445 * 10^-5, 0],
    ]

    
    function Wateruptake(Mi::Vector{Num},a_w) 
        r = [] #unit: mol/kg
        a = [0,1,2,3,4,5,6]
        b = a_w .^a

        for i in 1:18
            push!(r,sum(polynomial_fit[i] .* b))
        end
        return sum(Mi ./ r)
    end

    Wateruptake([0.0, 0, KHSO4_s, 0, 0, 0, 0, 0, 0, 0, NaHSO4_s, 0, 0, 0, 0, 0, NH4HSO4_s, 0],aw)

    z_all = [1, 1, 1, 1, 1, 2, 2, 1, 1]
    Mg_aq = 0
    
    eqs = [
        # Ionic strength
        I ~ cal_I(z_all,[NH4_aq, Na_aq, H_aq, Cl_aq, NO3_aq, SO4_aq, Ca_aq, K_aq, HSO4_aq])
        
        # Equilibrium Equations
        0 ~ Ca_aq*SO4_aq
        Equilibriumconst(K5,T,H5,C5) ~ K_aq*HSO4_aq*AC_special["KHSO4"]^2
        Equilibriumconst(K11,T,H11,C11) ~ (H_aq*SO4_aq/HSO4_aq)*(ActivityCoefficient(H2SO4,T,I)^3)/(ActivityCoefficient(HHSO4,T,I)^2)
        Equilibriumconst(K24,T,H24,C24) ~ Na_aq*HSO4_aq*AC_special["NaHSO4"]^2
        Equilibriumconst(K26,T,H26,C26) ~ NH4_aq*HSO4_aq*AC_special["NH4HSO4"]^2

        # Charge Conservation
        0 ~ Na_aq + H_aq + NH4_aq + K_aq + 2*Ca_aq - NO3_aq - Cl_aq - HSO4_aq - 2*SO4_aq 
        
        # Water uptake kg/m3
        W_ ~ Wateruptake([0.0, 0, KHSO4_s, 0, 0, 0, 0, 0, 0, 0, NaHSO4_s, 0, 0, 0, 0, 0, NH4HSO4_s, 0],0.8)
        #runnable if give aw a value
        
        # Mass conservation
        Na_i ~ Na_aq*W_ + NaHSO4_s # unit: mol/m3
        Ca_i ~ Ca_aq*W_ + CaSO4_s
        K_i ~ K_aq*W_ + KHSO4_s 
        H2SO4_i ~ SO4_aq*W_ + HSO4_aq*W_ + CaSO4_s + KHSO4_s + NaHSO4_s + NH4HSO4_s
        HCl_i ~ Cl_aq*W_ 
        HNO3_i ~ NO3_aq*W_ 
        NH3_i ~ NH4_aq*W_ + NH4HSO4_s
    ]

    @named ns = NonlinearSystem(eqs,[NH4_aq, Na_aq, H_aq, Cl_aq, NO3_aq, SO4_aq, HSO4_aq, Ca_aq, K_aq, NH4HSO4_s, NaHSO4_s, CaSO4_s, KHSO4_s, I, W_],[T, Na_i, H2SO4_i, NH3_i, HNO3_i, HCl_i, Ca_i, K_i, Mg_i])

    nlsys_func = generate_function(ns)[2]
    f = eval(nlsys_func)

    # [CaNO32_s, CaCl2_s, KHSO4_s, K2SO4_s, KNO3_s, KCl_s, MgSO4_s, MgNO32_s, MgCl2_s, NaNO3_s, NaHSO4_s, NaCl_s, Na2SO4_s, NH42SO4_s, NH4Cl_s, NH4NO3_s, NH4HSO4_s, NH43HSO42_s]
    # M_mass(g/mol) = [164.09, 110.98, 136.169, 174.259, 101.1032, 74.5513, 120.366, 148.3, 95.211, 84.9947, 120.06, 58.44, 142.04, 132.14, 53.491, 80.043, 115.11, 247.25] 

    M_mass_input = (1, 22.99*10^6, 98.079*10^6, 17.031*10^6, 63.01*10^6, 36.458*10^6, 40.08*10^6, 39.1*10^6, 24.31*10^6)
    params = (298.15, 0.0, 10.0, 3.4, 2.0, 0.0, 0.4, 0.33, 0.0)./M_mass_input
    # T Na_i H2SO4_i NH3_i HNO3_i HCl_i Ca_i K_i Mg_i aw

    j_func = generate_jacobian(ns)[2] # second is in-place
    #	takes 20 minutes to run!
    j! = eval(j_func)

    using NLsolve
    guess = Complex.([params[4],params[2],10^(-7),params[6],params[5],params[3],params[3],params[7],params[8],0,0,0,0,1,0])
    #NH4_aq, Na_aq, H_aq, Cl_aq, NO3_aq, SO4_aq, HSO4_aq, Ca_aq, K_aq, NH4HSO4_s, NaHSO4_s, CaSO4_s, KHSO4_s, I, W_

    nlsolve((out, x) -> f(out, x, params), (out, x) -> j!(out, x, params), guess,show_trace=true, extended_trace=true, iterations = 1000, method = :trust_region)

end

