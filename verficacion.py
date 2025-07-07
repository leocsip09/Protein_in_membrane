import re


def verificar_atomtypes_faltantes(path_top_or_itp, path_ffnonbonded):
    with open(path_top_or_itp, "r") as f:
        contenido_top = f.read()

    with open(path_ffnonbonded, "r") as f:
        contenido_ff = f.read()

    usados = set()
    in_atoms_section = False
    for line in contenido_top.splitlines():
        line = line.strip()
        if line.startswith("[ atoms ]"):
            in_atoms_section = True
            continue
        if in_atoms_section:
            if line.startswith("[") and not line.startswith("[ atoms ]"):
                break
            if line and not line.startswith(";"):
                parts = line.split()
                if len(parts) >= 2:
                    usados.add(parts[1])  # type is the second column

    definidos = set()
    in_atomtypes = False
    for line in contenido_ff.splitlines():
        if line.strip().startswith("[ atomtypes ]"):
            in_atomtypes = True
            continue
        if in_atomtypes:
            if line.strip().startswith("["):
                break
            if line.strip() and not line.strip().startswith(";"):
                definidos.add(line.split()[0])

    faltantes = sorted(usados - definidos)
    if faltantes:
        print("Atomtypes faltantes en ffnonbonded.itp:")
        for at in faltantes:
            print(f" - {at}")
    else:
        print("✅ Todos los atomtypes usados están definidos.")

    return faltantes


verificar_atomtypes_faltantes("/Users/leonardoponzebellido/Documents/Protein_in_membrane/membranas/dppc128/lipid.itp", "/Users/leonardoponzebellido/Documents/Simulacion/gromos53a6_lipid.ff/ffnonbonded.itp")