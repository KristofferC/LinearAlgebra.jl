using PackageCompiler

sysimage = tempname()
create_sysimage(["LinearAlgebra", "Test", "Distributed", "Dates", "REPL", "Printf", "Random"]; sysimage_path=sysimage, incremental=false, filter_stdlibs=true)

current_dir = @__DIR__
run(`$(Base.julia_cmd()) --sysimage=$sysimage --project=$current_dir -e 'Base.runtests(["LinearAlgebra"]; propagate_project=true)'`)
