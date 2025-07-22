# Julia module for Q distribution calculations
# See the driver programs q20dist.jl and qsolve.jl

module qdist

using Cubature, NLsolve, QuadGK, Distributed

function chi(A)
    3 * 1.2^2 * A^(5/3) / sqrt(5*pi)
end

function g0(x)
    # The function exp(-|x|) sinh(x) / x
    if x == 0
        return 1.
    else
        return -sign(x) *expm1(-2*abs(x))/(2*x)
    end
end


function g01(x)
    # The function exp(-|x|) sinh(x) / x
    if x == 0
        return 1.
    else
        return -sign(x) * (exp(-2*abs.(x))-1)/(2*x)
    end
end


function g1(x)
    # The function exp(-|x|) (sinh(x) - x cosh(x)) / x^2"""
    # The combined piecewise function is accurate to ~14 digits.
    coeff = [1/3, 1/3, 1/5, 4/45, 2/63, 1/105, 1/405, 8/14175, 2/17325]
    res = 0.
    if abs(x) < 0.0625
        # Series expansion for small x
        x1 = -abs(x)
        for i=1:length(coeff)
            res = res + x1 * coeff[i]
            x1 = -x1 * abs(x)
        end
        res = res * sign(x)
    else
        res = (-sign(x) * expm1(-2*abs(x)) - x*(exp(-2*abs(x)) + 1.)) /
        (2. * x^2) 
    end
    return res
end


function g2(x)
    # The function exp(-|x|) ((2 + x^2) sinh(x) - 2x cosh(x))/x^3
    # The combined piecewise function is accurate to ~14 digits.
    coeff = [1/3, 1/3, 4/15, 7/45, 22/315, 8/315, 22/2835, 29/14175,
             74/155925, 46/467775, 16/868725]
    res = 0.
    if abs(x) < 0.125
        # Series expansion for small x
        x1 = 1.
        for i=1:length(coeff)
            res = res + x1 * coeff[i]
            x1 = x1 * (-abs(x))
        end
    else
        res = ((2. + x^2) * (-sign(x))*expm1(-2*abs(x)) -
               2*x*(exp(-2*abs(x)) + 1.)) / (2*x^3)
    end
    return res
end


function bscale(a,b,c)
    # Approximate scale of β beyond which exp(-a β^2 - b β^3 - c β^4) decays.

    @assert c > 0 || (a > 0 && b == 0 && c == 0)

    s = 0.
    tau = a*c/b^2
    if tau < 9/32
        # Find the largest root of df/dβ
        β1 = 3*b/(8*c) * (1 + sqrt(1 - 32*tau/9))
        β2 = 3*b/(8*c) * (1 - sqrt(1 - 32*tau/9))
        s = max(abs(β1),abs(β2))
    end

    # Add the width of the exponential decay (ignoring β^3)
    # (Seems to work ok -- should look at some plots to confirm.)
    if a > 0 && c > 0
        s += max(1/sqrt(2*a), 0.58/c^(1/4))
    elseif a<=0 && c > 0
        s += 0.58/c^(1/4)
    elseif c == 0
        s += 1/sqrt(2*a)
    end
    
    return s
end


function pint(a,b,c,n,m)
    # Integral of β^n cos(3γ)^m with respect to the distribution
    #
    #  p(β,γ) β^4 dβ dγ dΩ
    #       = exp(-a β^2 - b β^3 cos(3γ) - c β^4) sin(3γ) β^4 dβ dγ dΩ
    #
    # So that <β^n cos(3γ)^m> = pint(a,b,c,n,m)/pint(a,b,c,0,0).
    #
    # Input:
    #  a,b,c: Parameters of the distribution (as above).
    #  n:     Power of β (n >= 0)
    #  m:     Power of cos(3γ) (must be m = 0, 1, or 2)
    #
    # Output:
    #  The integral \int p(β,γ) β^(n+4) cos(3γ)^m dβ dγ dΩ from β=0 to inf.,
    #  γ = 0 to π/3 and over the Euler angles Ω.
    
    # write(STDOUT, "a = $a, b = $b, c = $c\n");
    # flush(STDOUT)
    if !(c > 0 || (a > 0 && b == 0 && c == 0))
        error("ASSERTION FAILED ON LINE 109 of qdist.jl:\n a = $a, b = $b, c = $c")
    end

    # Scale s such that the integrand decays for beta*s ≳ 1.
       s = bscale(a,b,c)
    #s = 1
    function h(β)
        if m == 0
            return β^(n+4) * exp(-a*(β*s)^2 + abs(b)*(β*s)^3 - c*(β*s)^4) *
            g0(b*(β*s)^3)
        elseif m == 1
            return β^(n+4) * exp(-a*(β*s)^2 + abs(b)*(β*s)^3 - c*(β*s)^4) *
            g1(b*(β*s)^3)
        elseif m == 2
            return β^(n+4) * exp(-a*(β*s)^2 + abs(b)*(β*s)^3 - c*(β*s)^4) *
            g2(b*(β*s)^3)
        else
            error("Invalid value of m for in pint()")
            return 0
        end
    end


    beta1(t) = t / (1. - t)
    dbeta1(t) = 1. / (1. - t)^2

    g(t) = h(beta1(t))*dbeta1(t)

    res, err = quadgk(t->g(t),0,1,rtol=1E-13)
    res = res * s^(n+5) * 8*pi^2 / 3
    err = err * s^(n+5) * 8*pi^2 / 3
    return res
end


function f!(fx,x,q2,q3,q4,chi)
    # Vector-valued function f to minimize in order to determine a,b,c.
    #
    # Input:
    #   x:         Parameters [a,b,c] of the intrinsic quadrupole distribution.
    #   q2,q3,q4:  Desired second, third, and fourth moments of Q_{2,0}.
    #   chi:       Constant connecting the intrinsic deformation β^2 to Q^2.
    #
    # Output:
    #   [f1, f2, f3], where
    #
    #      f1 = χ^2 <β^2> /5 - <Q_{2,0}^2>
    #      f2 = χ^3 <β^3 cos(3γ)> * 2/35 - <Q_{2,0}^3>
    #      f3 = χ^4 <β^4> * 3/35 - <Q_{2,0}^4>

    a, b, c = x[1], x[2], x[3]

    I00 = pint(a,b,c,0,0)

    if I00 > 10000
        N = 0
    else
        N = 1/I00
    end
        
  
        
    
    E20 = pint(a,b,c,2,0)*N
    E31 = pint(a,b,c,3,1)*N
    E40 = pint(a,b,c,4,0)*N

    fx[1] = chi^2 * E20 / 5 - q2
    fx[2] = chi^3 * E31 * 2/35 - q3
    fx[3] = chi^4 * E40 * 3/35 - q4


# print("a = ",a,"\n")
## print("b = ",b,"\n")
# print("c = ",c,"\n")
# print("Norm  = ",I00,"\n")
# print("Res. of 2nd moment = ",fx[1],"\n")
# print("Res. of 3rd moment = ",fx[2],"\n")
# print("Res. of 4th moment = ",fx[3],"\n")
end


function J!(Jx,x,q2,q3,q4,chi)
    # Jacobian of f.
    #
    # Input:
    #   x:         Parameters [a,b,c] of the intrinsic quadrupole distribution.
    #   q2,q3,q4:  Desired second, third, and fourth moments of Q_{2,0}.
    #   chi:       Constant connecting the intrinsic deformation β^2 to Q^2.
    #
    # Output:
    #   The Jacobian matrix of f.

    a, b, c = x[1], x[2], x[3]
    I00 = pint(a,b,c,0,0)
    E20 = pint(a,b,c,2,0)/I00
    E31 = pint(a,b,c,3,1)/I00
    E40 = pint(a,b,c,4,0)/I00
    E51 = pint(a,b,c,5,1)/I00
    E60 = pint(a,b,c,6,0)/I00
    E62 = pint(a,b,c,6,2)/I00
    E71 = pint(a,b,c,7,1)/I00
    E80 = pint(a,b,c,8,0)/I00

    # Jacobian
    Jx[1,1] = -chi^2 * (E40 - E20^2)
    Jx[1,2] = -chi^2 * (E51 - E20*E31)
    Jx[1,3] = -chi^2 * (E60 - E20*E40)
    Jx[2,1] = -chi^3 * (E51 - E31*E20)
    Jx[2,2] = -chi^3 * (E62 - E31^2)
    Jx[2,3] = -chi^3 * (E71 - E31*E40)
    Jx[3,1] = -chi^4 * (E60 - E40*E20)
    Jx[3,2] = -chi^4 * (E71 - E40*E31)
    Jx[3,3] = -chi^4 * (E80 - E40^2)
end


function fradial(q20,r1,r2,theta1,a,b,c,chi)
    # The integrand for the marginal distribution,
    #   8π/5 exp(-a q^2/χ^2 + √(7/2) b (q x q).q/χ^3 - c(q.q)^2/χ^4)
    # using radial coordinates to represent q_{2,mu≠0}.
    q2 = q20^2 + 2*(r1^2 + r2^2)
    q3 = sqrt(2/7)*(-q20^3 - 3*r1^2*q20 + 6*r2^2*q20) -
         6*sqrt(3/7)*r1^2*r2*cos(theta1)
    q4 = q2^2
    return exp(-a*q2/chi^2 + sqrt(7/2)*b*q3/chi^3 - c*q4/chi^4)*r1*r2*8*pi
end



function pq20int_radial(q20,a,b,c,chi)
    # The integral for the marginal distribution
    #    P(q_{2,0}) =  8π/5 ∫ P(q_{2,mu≠0}) dq_{2,mu≠0},
    # where
    #    P(q_{2,mu}) = exp(-a q^2/χ^2 + √(7/2) b (q x q).q/χ^3 - c(q.q)^2/χ^4)
    # using radial coordinates to represent each q_{2,mu}.
    #
    # Input:
    #   a,b,c:  Parameters of the distribution
    #   chi:    Coefficient which converts between β and q
    #
    # See also pint_radial below.
    g(r1,r2,theta1) = fradial(q20,r1,r2,theta1,a,b,c,chi)

    # Change of variables to integrate from 0 to +inf
  #  x1(t) = t / (1. - t)
   # dx1(t) = 1. / (1. - t)^2
    
    #h(t) = g(x1(t[1]), x1(t[2]), t[3]) * dx1(t[1]) * dx1(t[2])

    r1(t) = t / (1. - t)
    dr1(t) = 1. / (1. - t)^2
    
    h(t) = g(r1(t[1]), r1(t[2]), t[3]) * (dr1(t[1]) * dr1(t[2]))

   
   # res, err = hcubature(h, [0.,0.,0.], [Inf,Inf,2*pi], reltol=1E-4)
    res, err = hcubature(h,[0.,0.,0.],[1.,1.,2*pi],reltol=1E-5)
    return res
end


function pq20(qvals,a,b,c,chi)
    # Compute the marginal distribution P(q_{2,0}) from a,b,c.
    #
    # Input:
    #   qvals:   q_{2,0} values (equally spaced)
    #   a,b,c:   Parameters of the distribution
    #   chi:     coefficient which converts between β and q
    #
    # Output:
    #   returns: Array of P(q_{2,0}) values
    
    # Normalization
    I00 = pint(a,b,c,0,0)*chi^5
 #   print((I00/pint_radial(a,b,c,chi)))
 #   I00 = pint_radial(a,b,c,chi)
    # Compute the distribution (in parallel)
    return pmap(x -> pq20int_radial(x,a,b,c,chi)/I00, qvals)
end





function fradial2(q2p2,q2m2,q20,r1,theta1,a,b,c,chi)
    # The integrand for the marginal distribution,
    #   8π exp(-a q^2/χ^2 + √(7/2) b (q x q).q/χ^3 - c(q.q)^2/χ^4)
    # using radial coordinates to represent q_{2,±1}.
    q2 = q20^2 + 2*r1^2 + 1/2*(q2p2^2 + q2m2^2)
    q3 = sqrt(2/7)*(-q20^3 - 3*r1^2*q20 + 3/2*q20*(q2p2^2 + q2m2^2)) -
         3*sqrt(3/7)*r1^2* (q2p2*cos(2*theta1) - q2m2*sin(2*theta1))
    q4 = q2^2
    return exp(-a*q2/chi^2 + sqrt(7/2)*b*q3/chi^3 - c*q4/chi^4)*r1
end


function pq22int_radial(q22,a,b,c,chi,pm)
    # The integral for the marginal distribution
    #    P(q_{2,0}) =  8π/5 ∫ P(q_{2,mu≠0}) dq_{2,mu≠0},
    # where
    #    P(q_{2,mu}) = exp(-a q^2/χ^2 + √(7/2) b (q x q).q/χ^3 - c(q.q)^2/χ^4)
    # using radial coordinates to represent each q_{2,mu}.
    #
    # Input:
    #   a,b,c:  Parameters of the distribution
    #   chi:    Coefficient which converts between β and q
    #
    # See also pint_radial below.


             # Change of variables to integrate from 0 to +inf
    x1(t) = t / (1. - t)
    dx1(t) = 1. / (1. - t)^2

    # Change of variables to integrate from -inf to +inf
    x(t) = t / (1. - t^2)
    dx(t) = (1. + t^2) / (1. - t^2)^2

    
    if pm == "p"
        g1(q2,q20,r1,theta1) = fradial2(q22,q2,q20,r1,theta1,a,b,c,chi)

        h1(t) = g1(x(t[1]),x(t[2]),x1(t[3]),t[4]) * dx(t[1]) * dx(t[2]) * dx1(t[3])
   # res, err = hcubature(h, [0.,0.,0.], [Inf,Inf,2*pi], reltol=1E-4)
        res, err = hcubature(h1,[-1.,-1.,0.,0.],[1.,1.,1.,2*pi],reltol=1E-3)

        
    elseif pm == "m"
        g2(q2,q20,r1,theta1) = fradial2(q2,q22,q20,r1,theta1,a,b,c,chi)
         h2(t) = g2(x(t[1]),x(t[2]),x1(t[3]),t[4]) * dx(t[1]) * dx(t[2]) * dx1(t[3])
   # res, err = hcubature(h, [0.,0.,0.], [Inf,Inf,2*pi], reltol=1E-4)
        res, err = hcubature(h2,[-1.,-1.,0.,0.],[1.,1.,1.,2*pi],reltol=1E-3)
    end

         
    return res
end



function pq22(qvals,a,b,c,chi,pm)
    # Compute the marginal distribution P(q_{2,0}) from a,b,c.
    #
    # Input:
    #   qvals:   q_{2,0} values (equally spaced)
    #   a,b,c:   Parameters of the distribution
    #   chi:     coefficient which converts between β and q
    #
    # Output:
    #   returns: Array of P(q_{2,0}) values
    
    # Normalization
    I00 = pint(a,b,c,0,0)*chi^5
 #   print((I00/pint_radial(a,b,c,chi)))
 #   I00 = pint_radial(a,b,c,chi)
    # Compute the distribution (in parallel)
    return pmap(x -> pq22int_radial(x,a,b,c,chi,pm)/I00, qvals)
end





################################################################################
# Some test routines                                                           #
################################################################################

## Testing the distribution P(q_{2,0})

function pint_radial(a,b,c,chi)
    # The integral of the entire unnormalized distribution,
    #    8π/5 ∫ P(q_{2,mu}) dq_{2,mu},
    # where
    #    P(q_{2,mu}) = exp(-a q^2/χ^2 + √(7/2) b (q x q).q/χ^3 - c(q.q)^2/χ^4)
    # using radial coordinates to represent each q_{2,mu}.
    #
    # Input:
    #   a,b,c:  Parameters of the distribution
    #   chi:    Coefficient which converts between β and q
    #
    # This is a test routine, compare to pint(a,b,c,0,0) * chi^5 (seems to
    # work).
    g(q20,r1,r2,theta1) = fradial(q20,r1,r2,theta1,a,b,c,chi)

    # Change of variables to integrate from 0 to +inf
    x1(t) = t / (1. - t)
    dx1(t) = 1. / (1. - t)^2
    # Change of variables to integrate from -inf to +inf
    x(t) = t / (1. - t^2)
    dx(t) = (1. + t^2) / (1. - t^2)^2
    
    h(t) = g(x(t[1]), x1(t[2]), x1(t[3]), t[4]) *
            dx(t[1]) * dx1(t[2]) * dx1(t[3])

    res, err = hcubature(h, [-1.,0.,0.,0.], [1.,1.,1.,2*pi], reltol=1E-4)

    return res
end


# Testing pint(a,b,c,n,m)

function testint()
    # Compare some integrals against Mathematica results.
    ints = Float64[]
    for a in [400.,1000.]
        for b in [-100.,-10.,10.,100.]
            for c in [10.,50.]
                for m in [0,1,2]
                    x = NaN
                    try
                        x = pint(a,b,c,0,m)
                    catch
                    end
                    append!(ints,[x])
                end
            end
        end
    end

    # Load Mathematica data from file
    mints = Float64[]
    open("$(ENV["HOME"])/Projects/qdist/misc/ints_MMA.dat") do file
        for line in readlines(file)
            append!(mints, [float(line)])
        end
    end

    # Compare my integrals to Mathematica
    maxerr = 0.
    for i=1:min(length(ints),length(mints))
        err = abs((mints[i]-ints[i])/mints[i])
        maxerr = max(err, maxerr)
    end

    return maxerr
end


end
