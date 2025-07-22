push!(LOAD_PATH, ".")
using qdist, NLsolve, QuadGK, Printf

function solve(filename,A)
    # Solve for the Q distribution parameters a,b,c given data on the moments
    # of Q_{2,0}.
    #
    # Input:
    #   filename:   File with Q_{2,0} moment data. (For format see below.)
    #   A:          Number of nucleons for this nucleus.
    #
    # Output:
    #   Writes a table of beta and a,b,c to stdout."""

    @printf("# Q distribution coefficients\n")
    @printf("# log(P) = log N - a β^2 - b β^3 cos(3γ) - c β^4\n")
    @printf("#%9s  %12s  %12s  %12s %12s\n", "beta",
        "a", "b", "c","tau")

    open(filename) do infile
        for line in readlines(infile)
            # Preserve empty lines of input in the output
            # (the jackknife blocks are separated by empty lines
            if length(strip(line)) == 0
                @printf("\n")
                continue
            end
            # Ignore comment lines
            if strip(line)[1] == '#'
                continue
            end
            # <beta> <q1>  <q2>  ... <q5> 
            vals = map(x->parse(Float64,x),split(line))
			# map(float,split(line))
            beta = vals[1]
         	#q2 = 4*vals[3]
         	#q3 = 8*vals[4]
         	#q4 = 16*vals[5]
			q2 = vals[3]
			q3 = vals[4]
			q4 = vals[5]
		  
		
            
            a0 = qdist.chi(A)^2/(2*q2)
            b0 = -5.
            c0 = 5.
	    x0 = [a0; b0; c0]
            a = b = c = 0.

           
            f1!(fx,x) = qdist.f!(fx,x,q2,q3,q4,qdist.chi(A))
            J1!(Jx,x) = qdist.J!(Jx,x,q2,q3,q4,qdist.chi(A))
            r = nlsolve(f1!, J1!,x0, xtol=1E-1,iterations=1000)
            a,b,c = r.zero
            tau = a*c/b^2
			
#	    print(converged(r))
			
			

			
            @printf("%10.6f  %12.10f  %12.10f  %12.10f  %12.10f\n",
                beta, a, b, c, tau)
            flush(stdout)
        end
    end
end

solve(ARGS[1],parse(Int,ARGS[2]))
