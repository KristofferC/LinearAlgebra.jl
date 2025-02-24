include("create_sysimage.jl")
sysimage, ncores = create_sysimage_linear_algebra()
doc_file = joinpath(@__DIR__, "..", "docs", "make.jl")
withenv("JULIA_NUM_THREADS" => 1) do
    run(`$(Base.julia_cmd()) --sysimage=$sysimage --project=$(dirname(doc_file)) $doc_file`)
end
