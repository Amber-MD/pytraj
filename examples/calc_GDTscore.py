import unittest
from pytraj import io as mdio
from pytraj.common_actions import calc_score

def main():
    # set frame index
    IDX0 = 9
    IDX1 = 8

    # load traj
    traj = mdio.load("../tests/data/md1_prod.Tc5b.x", "../tests/data/Tc5b.top")

    # write pdb files for TMalign program so we can compare our results to TMalign
    # ./TMalign -A test_gdt_0.pdb -B test_gdt_1.pdb
    # turn off those if you don't like
    mdio.writetraj(filename="./output/test_gdt_0.pdb", traj=traj[IDX0], top=traj.top, overwrite=True)
    mdio.writetraj(filename="./output/test_gdt_1.pdb", traj=traj[IDX1], top=traj.top, overwrite=True)

    # do our calculation
    # option:
    #     score = 'gdtscore', 'tmscore' or 'maxsub'
    #     "@CA" : follow AMBER mask, change it if needed
    score = 'tmscore'
    print('%s = %s ' % (score, calc_score(traj[IDX0], traj[IDX1], "@CA", traj.top, score=score)/1000.))
    print(calc_score(traj[IDX1], traj[IDX0], "@CA", traj.top, score=score)/1000.)

if __name__ == "__main__":
    main()
