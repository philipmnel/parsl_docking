from parsl.app.app import python_app, bash_app
from parsl.data_provider.files import File
from functools import lru_cache

@python_app
def dock_ligand(ligand, receptor, center, outputs=[]):
    """
    Dock ligand in Vina.
    Saves .pdbqt file of docked ligand.
    Returns list of Vina affinities of docked poses.
    """
    from vina import Vina
    from pathlib import Path

    v = Vina(sf_name="vina",
             cpu=0)
    docked_file = str(outputs[0])

    v.set_receptor(receptor)
    #if Path(docked_file).is_file():
    #    print("Found docked results.")
    #    with open(docked_file) as f:
    #        ret = f.read()
    #        ret = ret.split('\n')[1] 
    #        ret = ret.split(':')[1]
    #        ret = ret.strip().split(' ')[0]
    #        energies = [[ float(ret) ]] 
    #    print(energies)
    #else:
    v.set_ligand_from_file(ligand)
    v.compute_vina_maps(center=center, box_size=[10, 10, 10])

    # Dock the ligand
    v.dock(exhaustiveness=16, n_poses=10)
    #v.write_poses(docked_file, n_poses=5, overwrite=True)
    energies = v.energies(n_poses=5)
    return energies

@python_app
def process_poses(inputs=[]):
    """
    Take docked ligand file. 
    Return list of docked poses in xyz coordinates w/ charges.
    """
    import subprocess
    from rdkit import Chem
    
    docked_file = str(inputs[0])
    #ligand_file = f"{ligand_folder}/{ligand_file}"
    docked_sdf = f"{docked_file.split('.')[0]}.sdf"
    ret = subprocess.run([f"mk_export.py {docked_file} -o {docked_sdf}"], shell=True, capture_output=True, text=True)

    pose_list = []
    charge_list = []

    ret = subprocess.run([f"obabel -isdf {docked_sdf} -oxyz"], shell=True, capture_output=True, text=True)
    pose_str = ret.stdout.split('\n')
    natom = int(pose_str[0])
    #TODO: maybe grab multiple poses; not necessary :)
    pose_str = pose_str[2:natom+2]
    pose_str = "\n".join(pose_str)
    pose_list.append(pose_str)

    ret = subprocess.run([f"obabel -isdf {docked_sdf} -osmi"], shell=True, capture_output=True, text=True)
    smi_str = ret.stdout.split('\n')[0]
    rdk_mol = Chem.MolFromSmiles(smi_str)
    charge = Chem.GetFormalCharge(rdk_mol)
    charge_list.append(charge)

    return pose_list, charge_list

@lru_cache()
def load_model():
    import apnet
    model = apnet.PairModel.from_file("/theoryfs2/ds/pmnelson/chem/apnet/apnet/pair_models/pair0")
    print("initialized ap-net model :)")
    return model

@python_app
def rescore_ligand(protein_chg, protein_str, pose, charge):
    """
    Rescore docked pose with AP-Net.
    Returns AP-Net interaction energy.
    """
    import qcelemental as qcel

    dimer = qcel.models.Molecule.from_data(f"""
    {protein_chg} 1
    {protein_str}
    --
    {charge} 1
    {pose}
    """)

    model = load_model()    
    pred = model.predict(dimers=[dimer])

    return pred[0]

@python_app
def rescore_vina(receptor, center, inputs=[]):
    """
    Dock ligand in Vina.
    Saves .pdbqt file of docked ligand.
    Returns list of Vina affinities of docked poses.
    """
    from vina import Vina
    from pathlib import Path

    v = Vina(sf_name="vina",
             cpu=0)
    docked_file = str(inputs[0])

    with open(docked_file) as f:
        s = f.read()
        s = s.split("MODEL 1")[1]
        s = s.split("ENDMDL")[0]

    v.set_receptor(receptor)
    #if Path(docked_file).is_file():
    v.set_ligand_from_string(s)
    v.compute_vina_maps(center=center, box_size=[20, 20, 20])

    # Dock the ligand
    #v.dock(exhaustiveness=16, n_poses=10)
    #v.write_poses(docked_file, n_poses=5, overwrite=True)
    #energies = v.energies(n_poses=5)
    
    try:
        energies = v.score()
    except:
        energies = [0.0]

    return energies
