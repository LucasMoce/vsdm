"""
vsdm CLI — entry point principal
Molecular Dynamics pipeline for vsdock hits.
"""

import argparse
from vsdm import __version__


def cmd_build(args):
    from vsdm.build import build_system
    build_system(
        protein=args.protein,
        ligand=args.ligand,
        ligand_ff=args.ligand_ff,
        membrane=args.membrane,
        insertion_region=args.insertion_region,
        water=args.water,
        box_type=args.box_type,
        box_padding=args.box_padding,
        ions=args.ions,
        forcefield=args.ff,
        out=args.out,
    )


def cmd_minimize(args):
    from vsdm.minimize import minimize_system
    minimize_system(
        system=args.system,
        coords=args.coords,
        tolerance=args.tolerance,
        max_steps=args.max_steps,
        out=args.out,
    )


def cmd_equilibrate(args):
    from vsdm.equilibrate import equilibrate_system
    equilibrate_system(
        system=args.system,
        coords=args.coords,
        temperature=args.temp,
        pressure=args.pressure,
        nvt_time=args.nvt_time,
        npt_time=args.npt_time,
        restraints=args.restraints,
        out=args.out,
    )


def cmd_simulate(args):
    from vsdm.simulate import run_simulation
    run_simulation(
        system=args.system,
        coords=args.coords,
        checkpoint=args.checkpoint,
        time_ns=args.time,
        dt=args.dt,
        out=args.out,
    )


def cmd_analyze(args):
    from vsdm.analyze import analyze_trajectory
    analyze_trajectory(
        topology=args.topology,
        trajectory=args.trajectory,
        module=args.module,
        analysis_type=args.type,
        selection1=args.selection1,
        selection2=args.selection2,
        out=args.out,
    )


def cmd_free_energy(args):
    from vsdm.free_energy import calc_free_energy
    calc_free_energy(
        topology=args.topology,
        trajectory=args.trajectory,
        method=args.method,
        receptor=args.receptor,
        ligand=args.ligand,
        salt_conc=args.salt_conc,
        frames=args.frames,
        out=args.out,
    )


def cmd_alchemical(args):
    from vsdm.alchemical import run_alchemical
    run_alchemical(
        action=args.action,
        system=args.system,
        ligand_a=args.ligand_a,
        ligand_b=args.ligand_b,
        lambda_vdw=args.lambda_vdw,
        lambda_elec=args.lambda_elec,
        estimator=args.estimator,
        window_dir=args.window_dir,
        out=args.out,
    )


def cmd_permeation(args):
    from vsdm.permeation import run_permeation
    run_permeation(
        system=args.system,
        coords=args.coords,
        pull_group=args.pull_group,
        ref_group=args.ref_group,
        axis=args.axis,
        rate=args.rate,
        force_const=args.force_const,
        windows=args.windows,
        wham=args.wham,
        out=args.out,
    )


def main():
    parser = argparse.ArgumentParser(
        prog="vsdm",
        description="Molecular Dynamics pipeline for vsdock hits",
    )
    parser.add_argument("--version", action="version", version=f"vsdm {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    # build
    p = subparsers.add_parser("build", help="Prepara sistema para simulação")
    p.add_argument("--protein", required=True, help="Arquivo PDB da proteína")
    p.add_argument("--ligand", default=None, help="Arquivo do ligante (PDB, SDF, MOL2)")
    p.add_argument("--ligand-ff", default=None, dest="ligand_ff",
                   help="Parâmetros do ligante (XML OpenMM). Se omitido, gera automaticamente via acpype/OpenFF")
    p.add_argument("--membrane", default=None,
                   help="Composição lipídica (ex: 'POPC:100,DOPC:50')")
    p.add_argument("--insertion-region", default=None, dest="insertion_region",
                   help="Região transmembrana para alinhamento (ex: '1-20')")
    p.add_argument("--water", default="tip3p",
                   help="Modelo de água (tip3p, spce, tip4pew). Default: tip3p")
    p.add_argument("--box-type", default="cubic", dest="box_type",
                   help="Geometria da caixa (cubic, dodecahedron). Default: cubic")
    p.add_argument("--box-padding", type=float, default=1.0, dest="box_padding",
                   help="Distância proteína-borda em nm. Default: 1.0")
    p.add_argument("--ions", default="NaCl",
                   help="Sal para neutralização (ex: 'NaCl:0.15'). Default: NaCl")
    p.add_argument("--ff", default="amber14-all",
                   help="Campo de força da proteína. Default: amber14-all")
    p.add_argument("--out", default="system",
                   help="Prefixo dos arquivos de saída. Default: system")
    p.set_defaults(func=cmd_build)

    # minimize
    p = subparsers.add_parser("minimize", help="Minimização de energia")
    p.add_argument("--system", required=True, help="Arquivo XML do sistema")
    p.add_argument("--coords", required=True, help="Coordenadas iniciais (PDB)")
    p.add_argument("--tolerance", type=float, default=10.0,
                   help="Tolerância de convergência kJ/mol/nm. Default: 10.0")
    p.add_argument("--max-steps", type=int, default=10000, dest="max_steps",
                   help="Máximo de passos. Default: 10000")
    p.add_argument("--out", default="minimized", help="Prefixo de saída. Default: minimized")
    p.set_defaults(func=cmd_minimize)

    # equilibrate
    p = subparsers.add_parser("equilibrate", help="Equilibração NVT + NPT")
    p.add_argument("--system", required=True)
    p.add_argument("--coords", required=True)
    p.add_argument("--temp", type=float, default=300.0, help="Temperatura K. Default: 300")
    p.add_argument("--pressure", type=float, default=1.0, help="Pressão bar. Default: 1.0")
    p.add_argument("--nvt-time", type=float, default=0.1, dest="nvt_time",
                   help="Duração NVT em ns. Default: 0.1")
    p.add_argument("--npt-time", type=float, default=1.0, dest="npt_time",
                   help="Duração NPT em ns. Default: 1.0")
    p.add_argument("--restraints", type=float, default=1000.0,
                   help="Constante de restrição kJ/mol/nm². Default: 1000")
    p.add_argument("--out", default="equilibrated")
    p.set_defaults(func=cmd_equilibrate)

    # simulate
    p = subparsers.add_parser("simulate", help="Dinâmica molecular de produção")
    p.add_argument("--system", required=True)
    p.add_argument("--coords", required=True)
    p.add_argument("--checkpoint", default=None, help="Checkpoint da equilibração")
    p.add_argument("--time", type=float, required=True, help="Duração em ns")
    p.add_argument("--dt", type=float, default=2.0, help="Passo de integração em fs. Default: 2.0")
    p.add_argument("--out", default="production")
    p.set_defaults(func=cmd_simulate)

    # analyze
    p = subparsers.add_parser("analyze", help="Análise de trajetória")
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--module", required=True,
                   choices=["protein", "membrane", "ligand", "complex"],
                   help="Categoria de análise")
    p.add_argument("--type", required=True,
                   help="Métrica de análise (rmsd, rmsf, sasa, pca, dssp, rg, dccm, "
                        "apl, thickness, scd, tilt, lipid-diffusion, annular-lipids, "
                        "deformation, symmetry-rmsd, contacts, water-bridges, pharmacophore)")
    p.add_argument("--selection1", default="protein", dest="selection1")
    p.add_argument("--selection2", default=None, dest="selection2")
    p.add_argument("--out", default="analysis")
    p.set_defaults(func=cmd_analyze)

    # free-energy
    p = subparsers.add_parser("free-energy", help="Energia livre de ligação (MM-PB/GBSA)")
    p.add_argument("--topology", required=True)
    p.add_argument("--trajectory", required=True)
    p.add_argument("--method", required=True, choices=["pbsa", "gbsa"])
    p.add_argument("--receptor", default="protein")
    p.add_argument("--ligand", default="resname LIG")
    p.add_argument("--salt-conc", type=float, default=0.15, dest="salt_conc")
    p.add_argument("--igb", type=int, default=2)
    p.add_argument("--frames", type=int, default=100)
    p.add_argument("--out", default="free_energy")
    p.set_defaults(func=cmd_free_energy)

    # alchemical
    p = subparsers.add_parser("alchemical", help="Energia livre alquímica")
    p.add_argument("--action", required=True, choices=["setup", "run", "analyze"])
    p.add_argument("--system", default=None)
    p.add_argument("--ligand-a", default=None, dest="ligand_a")
    p.add_argument("--ligand-b", default=None, dest="ligand_b")
    p.add_argument("--lambda-vdw", type=int, default=11, dest="lambda_vdw")
    p.add_argument("--lambda-elec", type=int, default=11, dest="lambda_elec")
    p.add_argument("--estimator", default="mbar", choices=["mbar", "ti"])
    p.add_argument("--window-dir", default="lambda_windows", dest="window_dir")
    p.add_argument("--out", default="alchemical")
    p.set_defaults(func=cmd_alchemical)

    # permeation
    p = subparsers.add_parser("permeation", help="PMF por SMD + Umbrella Sampling")
    p.add_argument("--system", required=True)
    p.add_argument("--coords", required=True)
    p.add_argument("--pull-group", required=True, dest="pull_group")
    p.add_argument("--ref-group", required=True, dest="ref_group")
    p.add_argument("--axis", default="z", choices=["x", "y", "z"])
    p.add_argument("--rate", type=float, required=True, help="Taxa de pulling nm/ps")
    p.add_argument("--force-const", type=float, required=True, dest="force_const",
                   help="Constante de mola kJ/mol/nm²")
    p.add_argument("--windows", type=int, required=True)
    p.add_argument("--wham", action="store_true")
    p.add_argument("--out", default="permeation")
    p.set_defaults(func=cmd_permeation)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
