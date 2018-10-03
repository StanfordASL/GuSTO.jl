isinstalled(pkg::AbstractString) =
  pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)

asl_pkgs = ["BulletCollision", "AstrobeeRobot", "PandaRobot"]
for pkg in asl_pkgs
  if isinstalled(pkg)
    println("$pkg already installed")
    continue
  else
    println("Cloning out $pkg")
  end
  Pkg.clone("git://github.com/StanfordASL/$pkg.jl.git")
  Pkg.build(pkg)
end
