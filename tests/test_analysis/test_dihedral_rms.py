import pytraj as pt
from pytraj.testing import aa_eq, tempfolder
from utils import tz2_ortho_trajin, tz2_ortho_top


def test_dihedral_rms():
    cm = """
    parm {}
    trajin {}
    dihrms :1-10
    """.format(tz2_ortho_top, tz2_ortho_trajin)

    with tempfolder():
        # Run cpptraj state
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cpptraj_results = state.data[-1]

        # Run pytraj's dihedral_rms
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)
        pytraj_results = pt.dihedral_rms(
            traj,
            mask=":1-10",
        )

        # Compare results
        aa_eq(pytraj_results, cpptraj_results)


def test_dihedral_rms_to_first():
    with tempfolder():
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)

        # Dihedral RMSD to the first frame
        pytraj_results = pt.dihedral_rms(traj, mask="phi psi", ref=0)
        state = pt.datafiles.load_cpptraj_state("""
        parm {}
        trajin {}
        dihrms ToFirst out dihrms.dat noheader phi psi
        """.format(tz2_ortho_top, tz2_ortho_trajin)).run()
        cpptraj_results = state.data[-1]

        aa_eq(pytraj_results, cpptraj_results)


def test_dihedral_rms_to_reference():
    with tempfolder():
        traj = pt.iterload(tz2_ortho_trajin, tz2_ortho_top)

        # Dihedral RMSD to a reference frame
        pytraj_results = pt.dihedral_rms(traj, mask="phi psi", ref=5)
        state = pt.datafiles.load_cpptraj_state("""
        parm {}
        trajin {}
        reference {} 6 [MyRef]
        dihrms ToRef ref [MyRef] out toref.dat noheader phi psi
        """.format(tz2_ortho_top, tz2_ortho_trajin, tz2_ortho_trajin)).run()
        cpptraj_results = state.data[-1]

        aa_eq(pytraj_results, cpptraj_results)