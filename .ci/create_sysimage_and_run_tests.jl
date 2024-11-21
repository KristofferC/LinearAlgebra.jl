using PackageCompiler

sysimage = tempname()

if haskey(ENV, "BUILDKITE")
    ncores = Sys.CPU_THREADS
else
    ncores = ceil(Int, Sys.CPU_THREADS / 2)
end

withenv("JULIA_IMAGE_THREADS" => ncores) do
    create_sysimage(["LinearAlgebra", "Test", "Distributed", "Dates", "REPL", "Printf", "Random"]; sysimage_path=sysimage, incremental=false, filter_stdlibs=true)
end

current_dir = @__DIR__
cmd = """Base.runtests(["LinearAlgebra"]; propagate_project=true, ncores=$ncores)"""
withenv("JULIA_NUM_THREADS" => 1) do
    run(`$(Base.julia_cmd()) --sysimage=$sysimage --project=$current_dir -e $cmd`)
end
