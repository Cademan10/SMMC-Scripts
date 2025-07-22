# Compute q_{2,2}+q_{2,-2} and i(q_{2,2} - q_{2,-2}) distribution from a,b,c values.
# Usage: julia -p <np> q2dist.jl <data file> <Qmax> <M> <A> < 0 or plus(p) or minus(m)>
# where
#    <np> is the number of processes
#    <data file> is the output of qsolve.jl
#    <Qmax> is the desired maxium Q_{2,0} value
#    <M> is s.t. 2 M+1 = number of desired q points
#    <A> is the mass number.
# 	< 0 or plus(p) or minus(m)> determines which observable to calculate the distribution of (0 for Q_20, p for Q_22 + Q_2-2 and m for i * (Q_22 - Q2-2))
push!(LOAD_PATH, ".")
import qdist
using Printf

function distributions(filename,qmax,M,A,pm0)
    # Compute P(q_{2,0}) given a file with the parameters a,b,c for different
    # temperatures.
    #
    # Input:
    #   filename:   File with parameters a,b,c
    #   qmax:       Max value of |q_{2,0}|
    #   M:          2 M + 1 = number of q points
    #   A:          Number of nucleons for this nucleus.
    #  pm0:         Type of operator to project onto (0 = q20, p = q_{2,2} + q_{2,-2} and m = i * (q_{2,2} - q_{2,-2})
    #
    # Output:
    #   Writes a file "q2dist.dat" with all the data to disk.

    dq = 2*qmax/(2*M+1)
    qvals = range(-qmax-dq/2,qmax+dq/2,2*M+1)

    outfile = stdout
    if pm0 == "0"
        @printf(outfile,"# q_{2,0} distribution\n")
    elseif pm0 == "p"
        
        @printf(outfile,"# q_{2,2} + q_{2,-2} distribution\n")
    elseif pm0 == "m"
        @printf(outfile,"# i * (q_{2,2} - q_{2,-2}) distribution\n")
    end
    @printf(outfile,"# input file: %s\n", filename)
    @printf(outfile,"# M = %d, qmax = %f\n", M, qmax)
    
    index = 0
    open(filename) do infile
        
        for line in readlines(infile)
            if strip(line)[1] == '#' || length(strip(line)) == 0
                continue
            end
            # <beta> <a> <b> <c>
            vals = map(x->parse(Float64,x),split(line))
            beta, a, b, c = vals[1:4]
            @printf(outfile,"# index %d\n", index)
            @printf(outfile,"# beta = %f\n", beta)
            @printf(outfile,"#%14s %12s\n", "q_{2,0}", "P(q)")

            if pm0=="0"

                pvals = qdist.pq20(qvals,a,b,c,qdist.chi(A))
            elseif pm0 in ["p","m"]
                pvals = qdist.pq22(qvals,a,b,c,qdist.chi(A),pm0)
             end  
                    
            I00 = qdist.pint(a,b,c,0,0)*qdist.chi(A)^5
            for iq=1:length(qvals)
                @printf(outfile,"%15.8E %12.4E\n", qvals[iq], pvals[iq])
            end
            @printf(outfile,"\n\n")
            @printf(stderr,"Done with beta = %f, I00 = %f\n", beta, I00)
            flush(outfile)
            flush(stdout)
            index += 1
        end
    end
    nothing
end


#distributions("qdist.dat", 1000., 250, 152)
distributions(ARGS[1], parse(Float64,ARGS[2]), parse(Int,ARGS[3]), parse(Int,ARGS[4]), ARGS[5])
