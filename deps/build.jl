isinstalled(pkg::AbstractString) =
  pkg != "METADATA" && pkg != "REQUIRE" && pkg[1] != '.' && Pkg.cd(isdir, pkg)

asl_pkgs = ["BulletCollision", "AstrobeeRobot", "PandaRobot"]
for pkg in asl_pkgs
  isinstalled(pkg) && continue
  Pkg.clone("git://github.com/StanfordASL/$pkg.jl.git")
  Pkg.build(pkg)
end
