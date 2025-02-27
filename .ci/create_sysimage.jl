using PackageCompiler, Libdl

function create_sysimage_linear_algebra()
    sysimage = tempname() * "." * Libdl.dlext

    if haskey(ENV, "BUILDKITE")
        ncores = Sys.CPU_THREADS
    else
        ncores = ceil(Int, Sys.CPU_THREADS / 2)
    end

    withenv("JULIA_IMAGE_THREADS" => ncores) do
        create_sysimage(["LinearAlgebra", "Test", "Distributed", "Dates", "Printf", "Random"]; sysimage_path=sysimage, incremental=false, filter_stdlibs=true)
    end
    return sysimage, ncores
end
