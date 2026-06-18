# install.jl has to be in the same folder as integrate.jl, math.jl, etc. in order for package environment to successfully work!

using Pkg

Pkg.activate(".")

"""Checks if the package is installed. If not, installs the package."""
function ensure_package_installed(pkg_name::String)
    installed_packages = Pkg.dependencies()
    if !any(p -> p.name == pkg_name, values(installed_packages))
        println("Package '$pkg_name' not found. Installing...")
        Pkg.add(pkg_name)
        println("Package '$pkg_name' installed successfully.")
    end
end

ensure_package_installed("ClusterManagers")
ensure_package_installed("QuadGK")
ensure_package_installed("Interpolations")
ensure_package_installed("Memoization")
ensure_package_installed("LinearAlgebra")
ensure_package_installed("Distributed")
ensure_package_installed("OrdinaryDiffEq")
ensure_package_installed("FLoops")
ensure_package_installed("DistributedArrays")
ensure_package_installed("ForwardDiff")
ensure_package_installed("CSV")
ensure_package_installed("Dates")
ensure_package_installed("NFFT")
ensure_package_installed("FFTW")
ensure_package_installed("DataFrames")

#Pkg.resolve()
#Pkg.instantiate()
#Pkg.precompile()