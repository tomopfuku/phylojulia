using Distributions

function calc_total_dif(seq1,seq2)
    dif = 0
    for i in [1:length(seq1)]
        if seq1[i] != seq2[i]
            dif = dif + 1
        end
    end
    return dif
end

function dist_prior (dist::Float64)
    d = Normal(0.3,0.3)
    return log(pdf(d,dist))
end

function jc_loglike(theta, dif::Float64, sites::Float64)
    dist = theta[1]  
    if dist < 0
        return 1000000
    end
    first = dif * log( 0.75 - (0.75 * (e^(-4.*(dist/3.))) ))
    second = (sites-dif) * log( 0.25 + (0.75 * (e^(-4.*(dist/3.))) ))
    return (first + second)
end

function calc_SV_diff(seq1,seq2)
    dif = 0
    S=0
    V=0
    seqlen=length(seq1)
    purines = ['A','G']
    pyrimidines = ['C','T']
    for i in [1:length(seq1)]
        if seq1[i] != seq2[i]
            if seq1[i] in purines && seq2[i] in purines || seq1[i] in pyrimidines && seq2[i] in pyrimidines
                S = S + 1
            else 
                V = V + 1
            end
        end
    end
    return [float(S/seqlen),float(V/seqlen)] 
end            

function kappa_prior(kappa::Float64)
    d=Normal(3,2)
    return log(pdf(d,kappa))
end

function k80_loglike(theta,S,V,sites)
    d = theta[1]
    kappa = theta[2]
    if kappa < 0 || d < 0 
        return 1000000000
    end
    p0 = (1/4.)+((1/4.)*exp(-4*d/(kappa+2.)))+((1/2.)*exp(-2*d*(kappa+1.)/(kappa+2.)))
    p1 = (1/4.)+((1/4.)*exp(-4*d/(kappa+2.)))-((1/2.)*exp(-2*d*(kappa+1.)/(kappa+2.)))
    p2 = (1/4.)-((1/4.)*exp(-4*d/(kappa+2.)))
    return -((sites-S-V)*log(p0/4.)+S*log(p1/4.)+V*log(p2/4.))
end

function multiplier_prop(theta::Float64,epsilon=0.2)
    u = rand()
    c = e^((u-0.5)*epsilon)
    return theta * c
end

function sliding_window(theta, win::Float64 = 0.1)
    f = theta-win/2.
    b = theta+win/2.
    return f+rand()*(b-f)
end

function rj_proposal(theta)
    if length(theta) == 1
        d = theta[1]
        gamma = Gamma(5,5)
        kappa = rand(gamma)
        theta_star = [d, kappa]
    elseif length(theta) == 2
        theta_star = [theta[1]]
    end
    return theta_star
end

function rj_multiplier(theta)
    if length(theta) == 1
        d = theta[1]
        u = rand()
        kappa = d*(e^(5*(u-0.5))) 
        theta_star = [d, kappa]
    elseif length(theta) == 2
        d = theta[1]
        kappa = theta[2]
        if d/(e^(5/2)) < kappa && kappa < d*(e^(5/2))
            theta_star = [theta[1]]
        else
            theta_star = theta
        end
    end
    return theta_star
end

function mcmc ()
    srand(12345)
    seq1="ACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTAACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTA"
    seq2="ACACTCCGCCTACGTAACGTAGTGCCGTAACGCTCCACCTACGTAACGTAGTGCCGTAACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTA"
    printfreq = 1000
    outfl = open("outfile.test","w")
    dif = float(calc_total_dif(seq1,seq2))
    SV = calc_SV_diff(seq1,seq2)
    dif = 10.
    SV = [5.,5.]
    sites = float(length(seq1))
    theta = [0.3]
    ll = jc_loglike(theta[1],dif,sites)
    lp = dist_prior(theta[1])
    gen = 1000000
    for i in [1:gen]
        ##first update parameters using standard MH 
        if length(theta) == 1
            theta_star = multiplier_prop(theta[1],0.1)
            ll_star = jc_loglike(theta[1],dif,sites)
            lp_star = dist_prior(theta[1])
            alpha = exp(ll_star+lp_star-ll-lp)
            if alpha > rand()
                theta = theta_star
                ll = ll_star
                lp = lp_star
        end
        elseif length(theta) == 2
            if rand() > 0.5
                theta_star = [multiplier_prop(theta[1],0.1),theta[2]]
                ll_star = k80_loglike(theta_star,SV[1],SV[2],sites)
                lp_star = dist_prior(theta_star[1])
                alpha = exp(ll_star+lp_star-ll-lp)
            else
                theta_star = [theta[1],multiplier_prop(theta[2],0.4)]
                ll_star = k80_loglike(theta_star,SV[1],SV[2],sites)
                lp_star = kappa_prior(theta_star[2])
                alpha = exp(ll_star+lp_star-ll-lp)
            end
            if alpha > rand()
                theta = theta_star
                ll = ll_star
                lp = lp_star
            end
        end
        theta_star = rj_proposal(theta)
        if length(theta) == 1
            ll_star = k80_loglike(theta_star,SV[1],SV[2],sites)
            #alpha = 1/theta_star[2]
            lp_star = log(0.5)
            alpha = exp(ll_star+lp_star-ll-lp)*theta_star[2]   #(2/((2*1)-2-1))
            #println(alpha)
        else
            ll_star = jc_loglike(theta_star,dif,sites)
            #alpha = theta[2]
            lp_star = log(0.5)
            #alpha = exp(ll_star+lp_star-ll-lp)*   #(((2*2)-2-2+1)/2)
            alpha = exp(ll_star+lp_star-ll-lp)*1/theta[2] #(((2*2)-2-2+1)/2)
            #println(alpha)
        end
        if rand() < alpha
            theta = theta_star
            ll = ll_star
            lp = lp_star
        end
        if i%printfreq == 0
            model = length(theta)
            dist = theta[1]
            println(outfl,"$(i)\t$(model)\t$(theta)\t$(ll)\t$(lp)")
            #println("$(i)\t$(model)\t$(theta)\t$(ll)\t$(lp)")
        end
    end
    close(outfl)
end

@time mcmc()









