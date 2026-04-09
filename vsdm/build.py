"""
vsdm.build
==========
Prepara o sistema para simulação de dinâmica molecular usando OpenMM.

Fluxo:
  1. Gera parâmetros do ligante (acpype/GAFF2 ou OpenFF)
  2. Monta sistema solúvel (proteína + ligante + água + íons)
  3. Ou: monta sistema de membrana (proteína + ligante + bicamada + água + íons)
  4. Serializa sistema OpenMM (XML) e coordenadas (PDB)

Requer: openmm, acpype, packmol-memgen (para membrana)
"""

import os
import subprocess
import shutil
import tempfile
from pathlib import Path


# ---------------------------------------------------------------------------
# Parâmetros de ligante via acpype (GAFF2)
# ---------------------------------------------------------------------------

def _prepare_ligand_file(ligand_file: str, outdir: Path) -> Path:
    """
    Prepara o arquivo do ligante para o acpype.
    Converte PDBQT -> SDF via obabel se necessário.
    Remove sais e pega apenas o maior fragmento.
    """
    ligand_path = Path(ligand_file)
    suffix = ligand_path.suffix.lower()

    # Formatos aceitos diretamente pelo acpype
    if suffix in [".mol2", ".sdf", ".mol"]:
        return ligand_path

    # Converte PDBQT -> SDF via obabel
    if suffix == ".pdbqt":
        print(f"[build] Convertendo {ligand_path.name} de PDBQT para SDF...")

        # Extrai primeiro modelo do PDBQT
        AD_TO_ELEMENT = {
            "C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N",
            "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
            "H": "H", "HD": "H", "HS": "H", "P": "P", "F": "F",
            "Cl": "Cl","CL": "Cl","Br": "Br","BR": "Br","I": "I",
        }

        lines = []
        in_model = False
        for line in ligand_path.read_text().splitlines():
            if line.startswith("MODEL") and not in_model:
                in_model = True; continue
            if line.startswith("ENDMDL") and in_model:
                break
            if in_model and line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                ad_type = parts[-1] if parts else "C"
                element = AD_TO_ELEMENT.get(ad_type, "C")
                new_line = "HETATM" + line[6:17] + "LIG" + line[20:76].ljust(56) + f"{element:>2}"
                lines.append(new_line)

        # PDB -> SDF via obabel com geração de H
        sdf_path = outdir / f"{ligand_path.stem}.sdf"
        cmd = ["obabel", str(ligand_path), "-O", str(sdf_path), "-h", "-f", "1", "-l", "1"]
        result = subprocess.run(cmd, capture_output=True)

        if sdf_path.exists() and sdf_path.stat().st_size > 10:
            return sdf_path

    # Converte PDB -> SDF
    if suffix == ".pdb":
        sdf_path = outdir / f"{ligand_path.stem}.sdf"
        cmd = ["obabel", str(ligand_path), "-O", str(sdf_path), "-h"]
        subprocess.run(cmd, capture_output=True)
        if sdf_path.exists():
            return sdf_path

    raise ValueError(f"Formato de ligante não suportado: {suffix}. Use SDF, MOL2 ou PDBQT.")


def generate_ligand_params(
    ligand_file: str,
    outdir: Path,
    charge_method: str = "bcc",
    net_charge: int = 0,
) -> dict:
    """
    Gera parâmetros GAFF2 para o ligante usando acpype.

    Parâmetros
    ----------
    ligand_file   : arquivo do ligante (PDB, SDF, MOL2)
    outdir        : pasta de saída
    charge_method : método de carga ("bcc" = AM1-BCC, "gas" = Gasteiger)
    net_charge    : carga total da molécula

    Retorna
    -------
    dict com paths: xml (OpenMM), mol2, frcmod
    """
    if not shutil.which("acpype"):
        raise EnvironmentError(
            "acpype não encontrado. Instale com: conda install -c conda-forge acpype"
        )

    acpype_dir = outdir / "acpype_ligand"
    acpype_dir.mkdir(parents=True, exist_ok=True)

    # Converte formato se necessário (PDBQT -> SDF, etc.)
    ligand_ready = _prepare_ligand_file(ligand_file, acpype_dir)

    print(f"[build] Gerando parâmetros GAFF2 para {ligand_ready.name}...")
    print(f"  Método de carga : {charge_method}")
    print(f"  Carga total     : {net_charge}")

    cmd = [
        "acpype",
        "-i", str(ligand_ready),
        "-c", charge_method,
        "-n", str(net_charge),
        "-a", "gaff2",
        "-f",
    ]

    # Usa paths absolutos para evitar duplicação de caminhos
    ligand_abs = ligand_ready.resolve()
    run_dir    = ligand_abs.parent

    cmd = [
        "acpype",
        "-i", str(ligand_abs),
        "-c", charge_method,
        "-n", str(net_charge),
        "-a", "gaff2",
        "-f",
    ]

    # Roda no diretório do SDF para que a pasta .acpype seja criada lá
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(run_dir)
    )

    # acpype gera pasta com nome baseado no arquivo de entrada
    stem = ligand_abs.stem
    acpype_output = run_dir / f"{stem}.acpype"

    if not acpype_output.exists():
        candidates = list(run_dir.glob("*.acpype"))
        if candidates:
            acpype_output = candidates[0]
        else:
            print(f"[build] ERRO acpype stdout:\n{result.stdout[-500:]}")
            print(f"[build] ERRO acpype stderr:\n{result.stderr[-500:]}")
            raise RuntimeError("acpype falhou — pasta .acpype não gerada.")

    # Localiza arquivos gerados
    mol2_files  = list(acpype_output.glob("*.mol2"))
    frcmod_files = list(acpype_output.glob("*.frcmod"))

    # Converte para formato OpenMM usando parmed
    xml_path = _acpype_to_openmm_xml(acpype_output, outdir)

    print(f"[build] Parâmetros GAFF2 gerados:")
    print(f"  MOL2   : {mol2_files[0] if mol2_files else 'N/A'}")
    print(f"  XML    : {xml_path}")

    return {
        "xml":    xml_path,
        "mol2":   mol2_files[0] if mol2_files else None,
        "frcmod": frcmod_files[0] if frcmod_files else None,
        "acpype_dir": acpype_output,
    }


def _acpype_to_openmm_xml(acpype_dir: Path, outdir: Path) -> Path:
    """Converte parâmetros acpype para XML do OpenMM via parmed."""
    try:
        import parmed as pmd

        # Encontra arquivos AMBER (prmtop + inpcrd)
        prmtop_files = list(acpype_dir.glob("*.prmtop"))
        inpcrd_files = list(acpype_dir.glob("*.inpcrd"))

        if not prmtop_files or not inpcrd_files:
            raise FileNotFoundError("Arquivos AMBER (.prmtop, .inpcrd) não encontrados.")

        prmtop = prmtop_files[0]
        inpcrd = inpcrd_files[0]

        # Carrega com parmed e exporta para OpenMM XML
        struct = pmd.load_file(str(prmtop), str(inpcrd))
        xml_path = outdir / "ligand.xml"

        # Serializa force field
        from openmm.app import AmberPrmtopFile
        from openmm import XmlSerializer
        import openmm.app as app

        prmtop_omm = app.AmberPrmtopFile(str(prmtop))
        system = prmtop_omm.createSystem()

        with open(xml_path, "w") as f:
            f.write(XmlSerializer.serialize(system))

        return xml_path

    except Exception as e:
        print(f"[build] Aviso: conversão para XML falhou ({e}). Usando AMBER direto.")
        # Retorna path dos arquivos AMBER como fallback
        prmtop_files = list(acpype_dir.glob("*.prmtop"))
        return prmtop_files[0] if prmtop_files else None


# ---------------------------------------------------------------------------
# Sistema solúvel (proteína + ligante + água + íons)
# ---------------------------------------------------------------------------

def build_soluble_system(
    protein_file: str,
    ligand_xml: str = None,
    ligand_mol2: str = None,
    water_model: str = "tip3p",
    box_padding: float = 1.0,
    forcefield: str = "amber14-all",
    ions: str = "NaCl",
    ion_conc: float = 0.15,
    out: str = "system",
) -> dict:
    """
    Monta sistema solúvel com OpenMM.

    Parâmetros
    ----------
    protein_file : PDB da proteína (receptor_clean.pdb do vsdock)
    ligand_xml   : XML do ligante gerado pelo acpype
    ligand_mol2  : MOL2 do ligante (para posicionamento)
    water_model  : modelo de água (tip3p, spce, tip4pew)
    box_padding  : distância proteína-borda em nm
    forcefield   : campo de força da proteína
    ions         : sal para neutralização
    ion_conc     : concentração iônica em M
    out          : prefixo dos arquivos de saída

    Retorna
    -------
    dict com paths: system_xml, coords_pdb, topology
    """
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    print(f"[build] Montando sistema solúvel...")
    print(f"  Proteína        : {protein_file}")
    print(f"  Campo de força  : {forcefield}")
    print(f"  Modelo de água  : {water_model}")
    print(f"  Padding         : {box_padding} nm")
    print(f"  Íons            : {ions} {ion_conc} M")

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Carrega proteína
    pdb = app.PDBFile(protein_file)

    # Configura force field
    ff_files = [f"{forcefield}.xml", f"amber14/{water_model}.xml"]

    # Se tiver ligante com XML, adiciona
    if ligand_xml and Path(ligand_xml).exists():
        ff_files.append(str(ligand_xml))

    try:
        forcefield_obj = app.ForceField(*ff_files)
    except Exception as e:
        # Fallback para amber14 padrão
        print(f"[build] Aviso FF: {e}. Usando amber14-all + tip3p.")
        forcefield_obj = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    # Modeller para adicionar H e solvatar
    modeller = app.Modeller(pdb.topology, pdb.positions)

    print("[build] Adicionando hidrogênios...")
    modeller.addHydrogens(forcefield_obj)

    print("[build] Solvatando sistema...")
    modeller.addSolvent(
        forcefield_obj,
        model=water_model,
        padding=box_padding * unit.nanometers,
        ionicStrength=ion_conc * unit.molar,
    )

    print(f"[build] Sistema: {modeller.topology.getNumAtoms()} átomos")

    # Cria sistema OpenMM
    print("[build] Criando sistema OpenMM...")
    system = forcefield_obj.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )

    # Serializa sistema
    system_xml = f"{out}.xml"
    coords_pdb = f"{out}.pdb"

    print(f"[build] Salvando sistema...")
    with open(system_xml, "w") as f:
        f.write(mm.XmlSerializer.serialize(system))

    with open(coords_pdb, "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    print(f"[build] Sistema salvo:")
    print(f"  XML  : {system_xml}")
    print(f"  PDB  : {coords_pdb}")
    print(f"[build] Pronto para minimização: vsdm minimize --system {system_xml} --coords {coords_pdb}")

    return {
        "system_xml": system_xml,
        "coords_pdb": coords_pdb,
        "topology":   modeller.topology,
        "n_atoms":    modeller.topology.getNumAtoms(),
    }


# ---------------------------------------------------------------------------
# Sistema de membrana (proteína + bicamada + água + íons)
# ---------------------------------------------------------------------------

def build_membrane_system(
    protein_file: str,
    membrane: str = "POPC:100",
    insertion_region: str = None,
    ligand_mol2: str = None,
    water_model: str = "tip3p",
    box_padding: float = 1.0,
    forcefield: str = "amber14-all",
    ions: str = "NaCl",
    ion_conc: float = 0.15,
    out: str = "system_membrane",
) -> dict:
    """
    Monta sistema de membrana usando PACKMOL-Memgen.

    Parâmetros
    ----------
    protein_file      : PDB da proteína transmembrana
    membrane          : composição lipídica (ex: "POPC:100,DOPC:50")
    insertion_region  : resíduos transmembrana para alinhamento (ex: "1-20")
    ligand_mol2       : ligante a incluir no sistema
    water_model       : modelo de água
    box_padding       : padding em nm
    forcefield        : campo de força
    ions              : sal
    ion_conc          : concentração iônica em M
    out               : prefixo de saída

    Retorna
    -------
    dict com paths gerados pelo PACKMOL-Memgen
    """
    if not shutil.which("packmol-memgen"):
        raise EnvironmentError(
            "packmol-memgen não encontrado. "
            "Instale com: conda install -c conda-forge packmol-memgen"
        )

    print(f"[build] Montando sistema de membrana com PACKMOL-Memgen...")
    print(f"  Proteína    : {protein_file}")
    print(f"  Membrana    : {membrane}")
    print(f"  Padding     : {box_padding} nm")

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Monta lipid_ratio para packmol-memgen
    # Formato: "POPC:100,DOPC:50" -> "--lipid_ratio POPC:100 DOPC:50"
    lipid_parts = []
    for part in membrane.split(","):
        part = part.strip()
        if ":" in part:
            lipid, count = part.split(":")
            lipid_parts.extend([lipid.strip(), count.strip()])

    cmd = [
        "packmol-memgen",
        "--pdb", protein_file,
        "--lipids", membrane.replace(",", " ").replace(":", " "),
        "--ratio", " ".join(part.split(":")[1] for part in membrane.split(",") if ":" in part),
        "--dist", str(int(box_padding * 10)),   # Å
        "--salt", "0.15",
        "--salt_c", "Na+",
        "--salt_a", "Cl-",
        "--nottrim",
        "--overwrite",
        "--output", f"{out}_membrane.pdb",
    ]

    # Adiciona região de inserção se fornecida
    if insertion_region:
        cmd += ["--zmem", insertion_region]

    print(f"[build] Rodando PACKMOL-Memgen...")
    result = subprocess.run(cmd, capture_output=True, text=True)

    membrane_pdb = f"{out}_membrane.pdb"

    if not Path(membrane_pdb).exists():
        print(f"[build] ERRO PACKMOL-Memgen:\n{result.stderr}")
        raise RuntimeError("PACKMOL-Memgen falhou.")

    print(f"[build] Sistema de membrana gerado: {membrane_pdb}")
    print(f"[build] Próximo passo: vsdm minimize --system {out}.xml --coords {membrane_pdb}")

    # Agora cria o sistema OpenMM a partir do PDB de membrana
    # usando CHARMM36 que tem parâmetros de lipídeos
    try:
        result_omm = _membrane_to_openmm(membrane_pdb, forcefield, water_model, out)
        return result_omm
    except Exception as e:
        print(f"[build] Aviso: conversão OpenMM falhou ({e}). PDB de membrana disponível em {membrane_pdb}.")
        return {"membrane_pdb": membrane_pdb}


def _membrane_to_openmm(membrane_pdb: str, forcefield: str, water_model: str, out: str) -> dict:
    """Converte PDB de membrana para sistema OpenMM com CHARMM36."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    pdb = app.PDBFile(membrane_pdb)

    # CHARMM36 tem parâmetros para lipídeos comuns
    ff = app.ForceField("charmm36.xml", "charmm36/water.xml")

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.2 * unit.nanometers,
        constraints=app.HBonds,
        switchDistance=1.0 * unit.nanometers,
    )

    system_xml = f"{out}.xml"
    coords_pdb = f"{out}_final.pdb"

    with open(system_xml, "w") as f:
        f.write(mm.XmlSerializer.serialize(system))

    with open(coords_pdb, "w") as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

    print(f"[build] Sistema OpenMM salvo: {system_xml}")

    return {
        "system_xml":    system_xml,
        "coords_pdb":    coords_pdb,
        "membrane_pdb":  membrane_pdb,
        "n_atoms":       modeller.topology.getNumAtoms(),
    }


# ---------------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------------

def build_system(
    protein: str,
    ligand: str = None,
    ligand_ff: str = None,
    membrane: str = None,
    insertion_region: str = None,
    water: str = "tip3p",
    box_type: str = "cubic",
    box_padding: float = 1.0,
    ions: str = "NaCl",
    forcefield: str = "amber14-all",
    out: str = "system",
) -> dict:
    """
    Entry point principal do módulo build.
    Detecta automaticamente o tipo de sistema e prepara para simulação.
    """
    out_dir = Path(out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    ligand_params = None

    # 1. Gera parâmetros do ligante se necessário
    if ligand and not ligand_ff:
        print("[build] Gerando parâmetros de ligante automaticamente (GAFF2/acpype)...")
        ligand_params = generate_ligand_params(
            ligand_file=ligand,
            outdir=out_dir / "ligand_params",
        )
        ligand_ff = str(ligand_params["xml"]) if ligand_params.get("xml") else None
    elif ligand_ff:
        print(f"[build] Usando parâmetros de ligante fornecidos: {ligand_ff}")

    # 2. Monta sistema
    if membrane:
        return build_membrane_system(
            protein_file=protein,
            membrane=membrane,
            insertion_region=insertion_region,
            ligand_mol2=str(ligand_params["mol2"]) if ligand_params else ligand,
            water_model=water,
            box_padding=box_padding,
            forcefield=forcefield,
            ions=ions,
            out=out,
        )
    else:
        return build_soluble_system(
            protein_file=protein,
            ligand_xml=ligand_ff,
            ligand_mol2=str(ligand_params["mol2"]) if ligand_params else ligand,
            water_model=water,
            box_padding=box_padding,
            forcefield=forcefield,
            ions=ions,
            out=out,
        )
