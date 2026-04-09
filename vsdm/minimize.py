"""
vsdm.minimize
=============
Minimização de energia do sistema usando OpenMM.

Remove clashes estéricos gerados durante a solvatação e adição de íons,
preparando o sistema para equilibração.
"""

from pathlib import Path


def minimize_system(
    system: str,
    coords: str,
    tolerance: float = 10.0,
    max_steps: int = 10000,
    out: str = "minimized",
) -> dict:
    """
    Minimiza a energia do sistema.

    Parâmetros
    ----------
    system    : arquivo XML do sistema OpenMM
    coords    : arquivo PDB com coordenadas iniciais
    tolerance : tolerância de convergência em kJ/mol/nm. Default: 10.0
    max_steps : número máximo de passos. Default: 10000
    out       : prefixo dos arquivos de saída

    Retorna
    -------
    dict com paths: coords_pdb, energy_initial, energy_final
    """
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    print(f"[minimize] Iniciando minimização de energia...")
    print(f"  Sistema     : {system}")
    print(f"  Coordenadas : {coords}")
    print(f"  Tolerância  : {tolerance} kJ/mol/nm")
    print(f"  Max passos  : {max_steps}")

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Carrega sistema e coordenadas
    with open(system) as f:
        sys_obj = mm.XmlSerializer.deserialize(f.read())

    pdb = app.PDBFile(coords)

    # Configura integrador (necessário para criar contexto)
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # Detecta plataforma disponível
    platform = _get_best_platform()
    print(f"  Plataforma  : {platform.getName()}")

    simulation = app.Simulation(
        pdb.topology, sys_obj, integrator, platform
    )
    simulation.context.setPositions(pdb.positions)

    # Energia inicial
    state = simulation.context.getState(getEnergy=True)
    e_initial = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"\n[minimize] Energia inicial : {e_initial:.2f} kJ/mol")

    # Minimização
    print(f"[minimize] Minimizando...")
    simulation.minimizeEnergy(
        tolerance=tolerance * unit.kilojoules_per_mole / unit.nanometers,
        maxIterations=max_steps,
    )

    # Energia final
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    e_final = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    print(f"[minimize] Energia final   : {e_final:.2f} kJ/mol")
    print(f"[minimize] Redução         : {e_initial - e_final:.2f} kJ/mol ({(1 - e_final/e_initial)*100:.1f}%)")

    # Salva coordenadas minimizadas
    coords_out = f"{out}.pdb"
    positions  = state.getPositions()

    with open(coords_out, "w") as f:
        app.PDBFile.writeFile(simulation.topology, positions, f)

    print(f"\n[minimize] Coordenadas salvas em: {coords_out}")
    print(f"[minimize] Pronto para equilibração:")
    print(f"  vsdm equilibrate --system {system} --coords {coords_out}")

    return {
        "coords_pdb":    coords_out,
        "energy_initial": e_initial,
        "energy_final":   e_final,
    }


def _get_best_platform():
    """Seleciona a melhor plataforma disponível (CUDA > OpenCL > CPU)."""
    import openmm as mm

    for platform_name in ["CUDA", "OpenCL", "CPU"]:
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            return platform
        except Exception:
            continue

    return mm.Platform.getPlatformByName("CPU")
