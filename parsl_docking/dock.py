from parsl.app.app import python_app, bash_app
from parsl.data_provider.files import File
from functools import lru_cache

from typing import List, Tuple

@python_app
def dock_ligand(ligand: str, receptor: str, center: List[float], outputs: List[str]=[]) -> List[float]:
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
    if Path(docked_file).is_file():
        print("Found docked results.")
        with open(docked_file) as f:
            ret = f.read()
            ret = ret.split('\n')[1] 
            ret = ret.split(':')[1]
            ret = ret.strip().split(' ')[0]
            energies = [[ float(ret) ]] 
        print(energies)
    else:
        v.set_ligand_from_file(ligand)
        v.compute_vina_maps(center=center, box_size=[40, 40, 40])

        # Dock the ligand
        v.dock(exhaustiveness=16, n_poses=10)
        v.write_poses(docked_file, n_poses=5, overwrite=True)
        energies = v.energies(n_poses=5)
    return energies

@python_app
def process_poses(inputs=[]) :
    """
    Take docked ligand file. 
    Return list of docked poses in xyz coordinates w/ charges.
    """
    import subprocess
    from rdkit import Chem
    from pathlib import Path
    
    docked_file = str(inputs[0])
    #ligand_file = f"{ligand_folder}/{ligand_file}"
    docked_sdf = f"{docked_file.split('.')[0]}.sdf"
    if not Path(docked_sdf).is_file():
        ret = subprocess.run([f"mk_export.py {docked_file} -o {docked_sdf}"], shell=True, capture_output=True, text=True)
        if ret.returncode != 0:
            return "", 0
    #ret = subprocess.run([f"mk_export.py {docked_file} -o {docked_sdf}"], shell=True, capture_output=True, text=True)

    ret = subprocess.run([f"obabel -isdf {docked_sdf} -oxyz"], shell=True, capture_output=True, text=True)
    if ret.returncode != 0:
        return "", 0
    else:
        pose_str = ret.stdout.split('\n')
        natom = int(pose_str[0])
        #TODO: maybe grab multiple poses; not necessary :)
        pose_str = pose_str[2:natom+2]
        pose_str = "\n".join(pose_str)

    ret = subprocess.run([f"obabel -isdf {docked_sdf} -osmi"], shell=True, capture_output=True, text=True)
    if ret.returncode != 0:
        return "", 0
    else:
        smi_str = ret.stdout.split('\n')[0]
        rdk_mol = Chem.MolFromSmiles(smi_str)
        charge = Chem.GetFormalCharge(rdk_mol)

    return pose_str, charge

@lru_cache()
def load_model(model_num: int) -> object:
    import apnet
    #model = apnet.PairModel.from_file("/theoryfs2/ds/pmnelson/chem/apnet/apnet/pair_models/pair0")
    model = apnet.PairModel.pretrained(model_num)
    print("initialized ap-net model :)")
    return model

@python_app
def rescore_ligand(protein_chg: int, protein_str: str, charge: int, pose: str, model_path: str) -> float:
    """
    Rescore docked pose with AP-Net.
    Returns AP-Net interaction energy.
    """
    import qcelemental as qcel
    import numpy as np

    try:
        dimer = qcel.models.Molecule.from_data(f"""
        {protein_chg} 1
        {protein_str}
        --
        {charge} 1
        {pose}
        """)
        pred_list = []
        for model_num in range(5):
            model = load_model(model_num)
            pred = model.predict(dimers=[dimer])
            pred_list.append(sum(pred[0]))
        
        ret = np.mean(pred_list)
        std = np.std(pred_list)

    except:
        ret = 0.0
        std = 0.0

    return (ret, std)

