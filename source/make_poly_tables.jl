using Pkg
function ensure_package_installed(pkg_name::String)
    # Check if the package is installed
    installed_packages = Pkg.dependencies()
    if !any(p -> p.name == pkg_name, values(installed_packages))
        println("Package '$package_name' not found. Installing...")
        Pkg.add(package_name)
        println("Package '$package_name' installed successfully.")
    end
end

Pkg.activate(".")
ensure_package_installed("DifferentialEquations")
ensure_package_installed("QuadGK")
ensure_package_installed("Interpolations")
ensure_package_installed("CSV")

using CSV

include("polytrope.jl")