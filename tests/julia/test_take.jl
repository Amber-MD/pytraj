using PyCall

@pyimport pytraj as pt
@pyimport pytraj.sandbox as sb

traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
println(traj)

for frame in traj
    println(frame)
end

println(sb.take(traj, (1:5, "@CA")))
println(sb.get_top(traj))
println(pt.rmsd(traj))
println(pt.rmsd(traj, dtype="dict"))
