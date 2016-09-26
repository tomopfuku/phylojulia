using Distributions

function calc_dif(seq1,seq2)
    dif = 0
    for i in [1:length(seq1)]
        if seq1[i] != seq2[i]
            dif = dif + 1
        end
    end
    return dif
end

function prior (theta::Float64)
    d = Normal(0.3,0.1)
    return log(pdf(d,theta))
end

function likelihood(theta::Float64,x::Float64,n::Float64)
    return log((((3/4.)-((e^((-4.*theta)/3))*(3./4)))^(x))*(((1/4.)+((e^((-4.*theta)/3))*(3./4)))^(n-x)))
end

function loglikelihood(theta::Float64, x::Float64,n::Float64)
  if theta < 0
    return 1000000
  end
  first = x * log( 0.75 - (0.75 * (e^(-4.*(theta/3.))) ))
  second = (n-x) * log( 0.25 + (0.75 * (e^(-4.*(theta/3.))) ))
  return (first + second)
end

function multiplier_prop(theta::Float64)
    u = rand()
    epsilon = 0.2
    c = e^((u-0.5)*epsilon)
    return theta * c
end

function rand_prop(theta::Float64)
    d = Normal(theta,0.005)
    jump = abs(rand(d))
    return jump
end

function sliding_window(theta::Float64,win::Float64)
    f = theta-win/2.
    b = theta+win/2.
    return f+rand()*(b-f)
end

function mcmc ()
    srand(12345)
    seq1="ACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTAACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTA"
    seq2="ACGCTCCACCTACGTAACGTAGTGCCGTAACGCTCCACCTACGTAACGTAGTGCCGTAACGTTGCACCTGCGTGACGAAGTGACGTAACGCTCCACCTACGTAACGTAGTGCCGTA"
    printfreq = 1000
    outfl = open("outfile.test","w")
    x = float(calc_dif(seq1,seq2))
    n = float(length(seq1))
    theta = 0.5
    ll = loglikelihood(theta,x,n)
    lp = prior(theta)
    gen = 1000000
    for i in [1:gen]
        theta_star = multiplier_prop(theta)
        ll_star = loglikelihood(theta_star,x,n)
        lp_star = prior(theta_star)
        alpha = exp(ll_star+lp_star-ll-lp)
        if rand() < alpha
            theta = theta_star
            ll = ll_star
            lp = lp_star
        end
        if i%printfreq == 0
            println(outfl,"$(i)\t$(theta)\t$(ll)\t$(lp)")
            println("$(i)\t$(theta)\t$(ll)\t$(lp)")
        end
    end
    close(outfl)
end

@time mcmc()









