# vsdm

Molecular Dynamics pipeline para hits do vsdock.
Integra OpenMM, MDTraj, MDAnalysis e PACKMOL-Memgen.

## Instalação

```bash
# Dependências via conda (recomendado)
conda install -c conda-forge openmm mdtraj mdanalysis packmol-memgen acpype

# Instala o vsdm
pip install git+https://github.com/LucasMoce/vsdm.git
```

## Fluxo de uso

```bash
# 1. Prepara sistema (proteína solúvel + ligante)
vsdm build --protein receptor_clean.pdb \
           --ligand docking/pdbqt/CHEMBL1443.pdbqt \
           --ff amber14-all --out system

# 2. Prepara sistema de membrana
vsdm build --protein receptor_clean.pdb \
           --ligand docking/pdbqt/CHEMBL1443.pdbqt \
           --membrane "POPC:100" \
           --insertion-region "1-20" \
           --out system_membrane

# 3. Minimização
vsdm minimize --system system.xml --coords system.pdb --out minimized

# 4. Equilibração
vsdm equilibrate --system system.xml --coords minimized.pdb --out equilibrated

# 5. Produção
vsdm simulate --system system.xml --coords equilibrated.pdb \
              --checkpoint equilibrated.chk --time 100 --out production

# 6. Análise
vsdm analyze --topology system.xml --trajectory production.dcd \
             --module complex --type rmsd --out analysis

# 7. Energia livre (MM-GBSA)
vsdm free-energy --topology system.xml --trajectory production.dcd \
                 --method gbsa --out binding_energy
```

## Integração com vsdock

O vsdm lê diretamente os outputs do vsdock:

```
vsdock_project/
├── docking/
│   ├── docking_results.csv   ← ranking dos hits
│   └── pdbqt/                ← estruturas 3D
├── receptor_clean.pdb        ← receptor preparado
└── vsdock_state.yaml         ← parâmetros do projeto
```

## Subcomandos

| Comando | Descrição |
|---------|-----------|
| `build` | Prepara sistema (topologia, solvação, membrana) |
| `minimize` | Minimização de energia |
| `equilibrate` | Equilibração NVT + NPT |
| `simulate` | Produção MD |
| `analyze` | Análise de trajetória (RMSD, RMSF, PCA, etc.) |
| `free-energy` | MM-PB/GBSA |
| `alchemical` | Energia livre alquímica |
| `permeation` | PMF por SMD + Umbrella Sampling |
