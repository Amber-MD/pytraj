import pytraj as pt
from pytraj.testing import aa_eq, tempfolder, cpptraj_test_dir, assert_equal_dict


def _get_dict_from_data(cpptraj_result, name):
    out_dict = {}
    for d in cpptraj_result:
        if d.name == name:
            out_dict[d.key] = d.values
    return out_dict

def test_multipucker():
    # Test 1: Nucleic pucker analysis
    parm_file = f"{cpptraj_test_dir}/adh026.3.pdb"
    cm_nucleic = f"""
    parm {parm_file}
    trajin {parm_file}
    multipucker ADHas resrange 1-3 altona out nucleic.dat
    multipucker ADHcp resrange 1-3 cremer out nucleic.dat
    """
    with tempfolder():
        state_nucleic = pt.datafiles.load_cpptraj_state(cm_nucleic).run()
        cpptraj_results_nucleic = state_nucleic.data[:]

        traj = pt.iterload(parm_file, parm_file)
        pytraj_results_nucleic = pt.multipucker(
            traj,
            resrange="1-3",
            method="altona",
            out="nucleic.dat",
        )
        assert_equal_dict(pytraj_results_nucleic, _get_dict_from_data(cpptraj_results_nucleic, name="ADHas"))
        pytraj_results_nucleic_cp = pt.multipucker(
            traj,
            resrange="1-3",
            method="cremer",
            out="nucleic.dat",
        )
        assert_equal_dict(pytraj_results_nucleic_cp, _get_dict_from_data(cpptraj_results_nucleic, name="ADHcp"))

    # Test 2: Furanoid pucker analysis
    parm_file_furanoid = f"{cpptraj_test_dir}/Test_Pucker/Furanoid.mol2"
    cm_furanoid = f"""
    parm {parm_file_furanoid}
    trajin {parm_file_furanoid}
    multipucker Furanoid puckertype furanoid:C2:C3:C4:C5:O2 cremer \
      out furanoid.dat amplitude ampout furanoid.dat range360
    """
    with tempfolder():
        state_furanoid = pt.datafiles.load_cpptraj_state(cm_furanoid).run()
        cpptraj_results_furanoid = state_furanoid.data[:]

        traj_furanoid = pt.iterload(parm_file_furanoid, parm_file_furanoid)
        pytraj_results_furanoid = pt.multipucker(
            traj_furanoid,
            puckertype="furanoid:C2:C3:C4:C5:O2",
            method="cremer",
            out="furanoid.dat",
            amplitude=True,
            ampout="furanoid.dat",
            range360=True,
        )
        assert_equal_dict(pytraj_results_furanoid, _get_dict_from_data(cpptraj_results_furanoid, name="Furanoid"))

    # Test 3: Pyranoid pucker analysis
    parm_file_pyranoid = f"{cpptraj_test_dir}/Test_Pucker/Pyranoid.mol2"
    cm_pyranoid = f"""
    parm {parm_file_pyranoid}
    trajin {parm_file_pyranoid}
    multipucker Pyranoid puckertype pyranoid:C1:C2:C3:C4:C5:O5 cremer \
      out pyranoid.type.dat \
      amplitude ampout pyranoid.type.dat \
      theta thetaout pyranoid.type.dat range360
    multipucker Pyr pyranose cremer \
      out pyranoid.auto.dat \
      amplitude ampout pyranoid.auto.dat \
      theta thetaout pyranoid.auto.dat range360
    """
    with tempfolder():
        state_pyranoid = pt.datafiles.load_cpptraj_state(cm_pyranoid).run()
        cpptraj_results_pyranoid = state_pyranoid.data[:]

        traj_pyranoid = pt.iterload(parm_file_pyranoid, parm_file_pyranoid)
        pytraj_results_pyranoid = pt.multipucker(
            traj_pyranoid,
            puckertype="pyranoid:C1:C2:C3:C4:C5:O5",
            method="cremer",
            out="pyranoid.type.dat",
            amplitude=True,
            ampout="pyranoid.type.dat",
            theta=True,
            thetaout="pyranoid.type.dat",
            range360=True,
        )
        assert_equal_dict(pytraj_results_pyranoid, _get_dict_from_data(cpptraj_results_pyranoid, name="Pyranoid"))