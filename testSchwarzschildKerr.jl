# use modules to scope the namespaes
module Schwarzschild
    include("Schwarzschild.jl")
end

module Kerr
    include("Kerr.jl")
end

Schwarzschild.g0
Kerr.g0

gammaKerr = Kerr.compute_christoffel_analytical(100.0, deg2rad(90))
gammaSchwarz = Schwarzschild.compute_christoffel_analytical(100.0, deg2rad(90))

gammaKerr[1, 1, 2]
gammaSchwarz[1, 1, 2]

gammaKerr[2, 1, 1]
gammaSchwarz[2, 1, 1]

gammaKerr[2, 2, 2]
gammaSchwarz[2, 2, 2]

gammaKerr[2, 3, 3]
gammaSchwarz[2, 3, 3]

gammaKerr[2, 4, 4]
gammaSchwarz[2, 4, 4]

gammaKerr[3, 2, 3]
gammaSchwarz[3,2,3]

gammaKerr[3, 4, 4]
gammaSchwarz[3, 4, 4]

gammaKerr[4, 2, 4]
gammaSchwarz[4, 2, 4]

gammaKerr[4, 3, 4]
gammaSchwarz[4, 3, 4]
