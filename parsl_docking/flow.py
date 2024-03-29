import pandas as pd
from parsl.data_provider.files import File
import subprocess
from pathlib import Path
from typing import List

#local imports
from .dock import dock_ligand, process_poses, rescore_ligand


def prepare_ligands(input_df, data_dir):
    """
    Convert SDF files to PDBQT files
    """

    sdf_file_list = input_df["SDF File"]
    pdbqt_file_list = sdf_file_list.apply(lambda x: x.replace("sdf", "pdbqt"))
    pdbqt_file_list2 = []

    for sdf_file, pdbqt_file in zip(sdf_file_list, pdbqt_file_list):
        path = Path(f"{data_dir}/{pdbqt_file}")
        if path.is_file():
            pdbqt_file_list2.append(pdbqt_file)
            continue
        else:
            ret = subprocess.run([f"mk_prepare_ligand.py -i {data_dir}/{sdf_file} -o {data_dir}/{pdbqt_file}"], shell=True, capture_output=True, text=True)
        if path.is_file():
            pdbqt_file_list2.append(pdbqt_file)
        else:
            pdbqt_file_list2.append("Failed")

    return pdbqt_file_list2



def flow(input_df: pd.DataFrame, receptor: str, center: list, protein_file: str, protein_chg: int, model_path: str):
    """
    Main function to dock ligands and rescore with AP-Net

    input_df: dataframe with ligand information. Required columns: "PDBQT File", with path to ligand file
    receptor: path to pdbqt file of receptor
    center: center of docking box
    protein_file: path to xyz file of receptor
    protein_chg: charge of receptor
    model_path: path to AP-Net model   

    returns: dataframe with docking and APNet rescored energies
    """


    #receptor = "data/5htr/target/receptor.pdbqt"
    #center = [110.735, 129.173, 119.183]
    result_df = input_df.copy()

    ligand_list = list(result_df["PDBQT File"])
    #name_list = list(map(lambda f: f.split('/')[-1].split('.')[0], ligand_list))
    name_list = [name.split('/')[-1].split('.')[0] for name in ligand_list]

    docked_file_list = []
    docked_file_path = ligand_list[0].split('/')[:-2]
    docked_file_path.append('docked')
    docked_file_path = '/'.join(docked_file_path)

    for name in name_list:
        docked_file = f'{docked_file_path}/{name}-docked.pdbqt'
        docked_file_list.append(docked_file)

    result_df["Docked File"] = docked_file_list

    #dock ligands
    task_list = []
    for ligand, docked_file in zip(ligand_list, docked_file_list):
        task_list.append(dock_ligand(ligand, receptor, center, outputs=[File(docked_file)]))

    results = [task.result() for task in task_list]
    energies = [e[0][0] for e in results]
    result_df["Vina Energy"] = energies

    #process results
    task_list = []
    for docked_file in result_df["Docked File"]:
        task_list.append(process_poses(inputs=[File(docked_file)]))

    results = [task.result() for task in task_list]
    poses = [p for p, _ in results]
    charges = [c for _, c in results]

    #rescore docked ligands
    #protein_file = open("data/5htr/target/rec.xyz")
    with open(protein_file, 'r') as f:
        protein_str = f.read()
    #protein_chg = 10

    task_list = []
    for pose, charge in zip(poses,charges):
        task_list.append(rescore_ligand(protein_chg, protein_str, charge, pose, model_path))

    results = [task.result() for task in task_list]
    print(results)
    result_df["AP Net Energy"] = [r[0] for r in results]
    result_df["AP Net STD"] = [r[1] for r in results]

    #print(result_df)
    
    return result_df
