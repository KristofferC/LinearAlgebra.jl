using LinearAlgebra
using Documenter: DocMeta, makedocs, deploydocs, HTML

DocMeta.setdocmeta!(LinearAlgebra, :DocTestSetup, :(using LinearAlgebra); recursive=true)

makedocs(
    modules = [LinearAlgebra],
    sitename = "LinearAlgebra",
    pages = Any[
        "LinearAlgebra" => "index.md",
        ];
    warnonly = [:missing_docs, :cross_references],
    format = HTML(size_threshold = 600 * 2^10 #=600 KiB=#),
    )

deploydocs(repo = "github.com/KristofferC/LinearAlgebra.jl.git", push_preview = true)
