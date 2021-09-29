using Documenter, ExternalBallistics
push!(LOAD_PATH,"../src/")
makedocs(sitename="My Documentation")

makedocs(
         sitename = "ExternalBallistics.jl",
         modules  = [ExternalBallistics],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/YayeIrene/ExternalBallistics.jl.git",
)
