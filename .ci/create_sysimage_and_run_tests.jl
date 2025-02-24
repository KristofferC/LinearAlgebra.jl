include("create_sysimage.jl")
sysimage, ncores = create_sysimage_linear_algebra()
current_dir = @__DIR__
cmd = """Base.runtests(["LinearAlgebra"]; propagate_project=true, ncores=$ncores)"""
withenv("JULIA_NUM_THREADS" => 1) do
    run(`$(Base.julia_cmd()) --sysimage=$sysimage --project=$current_dir -e $cmd`)
end
