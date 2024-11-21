using PackageCompiler

sysimage = tempname()
create_sysimage(["LinearAlgebra", "Test", "Distributed", "Dates", "REPL", "Printf", "Random"]; sysimage_path=sysimage, incremental=false, filter_stdlibs=true)

current_dir = @__DIR__
ncores = Sys.isapple() ? Sys.CPU_THREADS : ceil(Int, Sys.CPU_THREADS / 2)
cmd = """Base.runtests(["LinearAlgebra"]; propagate_project=true, ncores=$ncores)"""
run(`$(Base.julia_cmd()) --sysimage=$sysimage --project=$current_dir -e $cmd`)
