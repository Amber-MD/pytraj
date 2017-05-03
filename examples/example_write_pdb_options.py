import pytraj as pt
from pytraj.utils import tempfolder


def example_write_pdb():
    traj = pt.iterload("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")
    # multiple pdb in multiple files, using `save` method in traj
    with tempfolder():
        basename = "test_pdb_files.pdb"
        traj.save(basename, overwrite=True, options="multi")

    # multiple pdb in multiple files, using `pt.io.write_traj`
    with tempfolder():
        basename = "test_pdb_files_pt.io_write_traj.pdb"
        pt.io.write_traj(basename, traj, overwrite=True, options="multi")

    # multiple pdb in SINGLE file
    with tempfolder():
        basename = "test_pdb_files.pdb"
        traj.save(basename, overwrite=True)

    # multiple pdb in SINGLE file with `optionsl` keyword
    # write to output so we can manually check
    with tempfolder():
        basename = "test_pdb_files_model.pdb"
        traj.save(basename, overwrite=True, options='model')


if __name__ == "__main__":
    example_write_pdb()
