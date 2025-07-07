import os
import subprocess
import shutil
import re

def run(cmd, msg=None):
    print(f"\n> Ejecutando: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    if msg:
        print(msg)

def elegir_membrana():
    opciones = {
        "1": "dppc128.pdb",
    }
    print("\nSelecciona la membrana (Lista disponible por el momento):")
    for k, v in opciones.items():
        print(f"{k}. {v}")
    while True:
        eleccion = input("Opción: ").strip()
        if eleccion in opciones:
            return opciones[eleccion]
        print("Opción no válida")

def generar_topologia(prot, pdb, pdb_dir):
    print("\n====== Paso 1: Generar topología de la proteína ======\n")
    print("Para esta ocasión vamos a escoger GROMOS96 53A6 para la topología de la proteína (Opción 13)")
    print("Para los grupos acetilo y amida, selecciona 'None' (Opción 2 para ambas) \n¿Estás de acuerdo? (s/n): ", end="")
    if input().strip().lower() != 's':
        print("Operación cancelada.")
        return

    output_gro = os.path.join(pdb_dir, f"{prot}_processed.gro")
    run(f"gmx pdb2gmx -f \"{pdb}\" -o \"{output_gro}\" -ignh -ter -water spc",
        "Topología generada, ¡yeiiii!")
    
    for fname in ["topol.top", "posre.itp"]:
        if os.path.exists(fname):
            shutil.move(fname, os.path.join(pdb_dir, fname))
    
    for f in os.listdir("."):
        if f.endswith(".ff") and os.path.isdir(f):
            shutil.move(f, os.path.join(pdb_dir, f))

def insertar_en_seccion(filepath, seccion, nuevas_lineas, limpiar=False):
    """
    Inserta nuevas líneas en la sección indicada de un archivo .itp.
    Si limpiar=True, elimina ciertas líneas innecesarias como comentarios obsoletos o duplicados.
    """
    with open(filepath, 'r') as f:
        lineas = f.readlines()

    dentro = False
    nueva_seccion = []
    for i, line in enumerate(lineas):
        if line.strip().lower().startswith(f"[ {seccion.lower()} ]"):
            dentro = True
            nueva_seccion.append(i)
        elif dentro and line.strip().startswith('['):
            nueva_seccion.append(i)
            break

    if not nueva_seccion:
        print(f"No se encontró la sección [{seccion}] en {filepath}")
        return

    inicio = nueva_seccion[0] + 1
    fin = nueva_seccion[1] if len(nueva_seccion) > 1 else len(lineas)

    contenido = lineas[:inicio]

    originales = lineas[inicio:fin]
    if limpiar and seccion.lower() == "nonbond_params":
        originales = [
            l for l in originales
            if ";; parameters for lipid-GROMOS interactions" not in l
            and "HW" not in l
        ]

    contenido += originales

    for l in nuevas_lineas:
        if l.strip().startswith("[") or l.strip().startswith(";") or l.strip() == "":
            continue
        contenido.append(l)

    contenido += lineas[fin:]

    with open(filepath, "w") as f:
        f.writelines(contenido)


def modificar_topologia(prot, pdb_dir):
    print("\n====== Paso 2: Modificar la topología de la proteína ======\n")
    print("Vamos a modificar el campo de fuerza para incluir parámetros de lípidos Berger.")
    print("¿Deseas continuar? (s/n): ", end="")
    if input().strip().lower() != 's':
        print("Operación cancelada.")
        return
    origen_ff = "/usr/local/gromacs/share/gromacs/top/gromos53a6.ff"
    destino_ff = os.path.join(pdb_dir, "gromos53a6_lipid.ff")
    lipid_itp = "membranas/dppc128/lipid.itp"
    print("\nModificando campo de fuerza (forcefield.doc) para incluir parámetros de lípidos Berger...")

 
    run(f"cp -r \"{origen_ff}\" \"{destino_ff}\"")

    with open(os.path.join(destino_ff, "forcefield.doc"), "w") as f:
        f.write("GROMOS96 53a6 force field, extended to include Berger lipid parameters\n")

    print("\nforcefield.doc modificado correctamente ^^\n")

    if not os.path.exists(lipid_itp):
        print("Error: No se encontró lipid.itp")
        return

    with open(lipid_itp, "r") as f:
        lipid_lines = f.readlines()

    def extraer_bloque(tag):
        start, end = None, None
        for i, line in enumerate(lipid_lines):
            if line.strip().lower().startswith(f"[ {tag.lower()} ]"):
                start = i
                continue
            if start and line.strip().startswith("[") and not line.strip().lower().startswith(f"[ {tag.lower()} ]"):
                end = i
                break
        return lipid_lines[start:end] if start else []

    atomtypes = extraer_bloque("atomtypes")
    nonbond_params = extraer_bloque("nonbond_params")
    pairtypes = extraer_bloque("pairtypes")
    dihedraltypes = extraer_bloque("dihedraltypes")

    tabla_atomtypes = {
        "LO":   "   LO    8    15.9994      0.000     A  2.36400e-03 1.59000e-06 ;carbonyl O, OPLS\n",
        "LOM":  "  LOM    8    15.9994      0.000     A  2.36400e-03 1.59000e-06 ;carboxyl O, OPLS\n",
        "LNL":  "  LNL    7    14.0067      0.000     A  3.35300e-03 3.95100e-06 ;Nitrogen, OPLS\n",
        "LC":   "   LC    6    12.0110      0.000     A  4.88800e-03 1.35900e-05 ;Carbonyl C, OPLS\n",
        "LH1":  "  LH1    6    13.0190      0.000     A  4.03100e-03 1.21400e-05 ;CH1, OPLS\n",
        "LH2":  "  LH2    6    14.0270      0.000     A  7.00200e-03 2.48300e-05 ;CH2, OPLS\n",
        "LP":   "   LP   15    30.9738      0.000     A  9.16000e-03 2.50700e-05 ;phosphor, OPLS\n",
        "LOS":  "  LOS    8    15.9994      0.000     A  2.56300e-03 1.86800e-06 ;ester oxygen, OPLS\n",
        "LP2":  "  LP2    6    14.0270      0.000     A  5.87400e-03 2.26500e-05 ;RB CH2, Bergers LJ\n",
        "LP3":  "  LP3    6    15.0350      0.000     A  8.77700e-03 3.38500e-05 ;RB CH3, Bergers LJ\n",
        "LC3":  "  LC3    6    15.0350      0.000     A  9.35700e-03 3.60900e-05 ;CH3, OPLS\n",
        "LC2":  "  LC2    6    14.0270      0.000     A  5.94700e-03 1.79000e-05 ;CH2, OPLS\n",
    }

    atomtypes_mod = []
    for line in atomtypes:
        if line.strip().startswith("[") or line.strip().startswith(";") or line.strip() == "":
            continue
        key = line.split()[0]
        atomtypes_mod.append(tabla_atomtypes.get(key, line))

    ffnb = os.path.join(destino_ff, "ffnonbonded.itp")
    ffbonded = os.path.join(destino_ff, "ffbonded.itp")

    insertar_en_seccion(ffnb, "atomtypes", atomtypes_mod)
    print("Sección [ atomtypes ] modificada correctamente.")
    insertar_en_seccion(ffnb, "nonbond_params", nonbond_params, limpiar=True)
    print("Sección [ nonbond_params ] modificada correctamente.")
    insertar_en_seccion(ffnb, "pairtypes", pairtypes)
    print("Sección [ pairtypes ] modificada correctamente.")
    insertar_en_seccion(ffbonded, "dihedraltypes", dihedraltypes)
    print("Sección [ dihedraltypes ] modificada correctamente.")

    topol = os.path.join(pdb_dir, "topol.top")
    if os.path.exists(topol):
        with open(topol, "r") as f:
            contenido = f.read()
        contenido = contenido.replace(
            'gromos53a6.ff/forcefield.itp',
            'gromos53a6_lipid.ff/forcefield.itp'
        ).replace(
            'gromos53a6.ff/spc.itp',
            'gromos53a6_lipid.ff/spc.itp'
        ).replace(
            'gromos53a6.ff/ions.itp',
            'gromos53a6_lipid.ff/ions.itp'
        )
        with open(topol, "w") as f:
            f.write(contenido)

        with open(topol, "r") as f:
            lines = f.readlines()

        new_lines = []
        inserted_dppc = False
        for i, line in enumerate(lines):
            new_lines.append(line)

            if not inserted_dppc and line.strip() == "#endif":
                if i >= 2 and lines[i - 2].strip() == "#ifdef POSRES":
                    new_lines.append("\n; Include DPPC chain topology\n")
                    new_lines.append('#include "dppc.itp"\n')
                    inserted_dppc = True

        with open(topol, "w") as f:
            f.writelines(new_lines)

        print("\ntopol.top actualizado con include de DPPC, SPC e IONS corregidos.")



    print("\nCampo de fuerza y topología modificados correctamente ^v^.")


def main():
    print("\n====== Proteina en membrana con GROMACS automatizada c; ======\n")

    prot = input("¿Cuál es el nombre de la proteína? (sin extensión, el nombre que va antes de .pdb): ").strip()
    if not prot:
        print("Error: nombre de la proteína no puede estar vacío.")
        return
    
    print("\nFriendly reminder: toda la simulación se hará en el directorio donde se encuentre el .pdb de la proteína.")
    pdb = input("Introduce la ruta al archivo .pdb de la proteína: ").strip()
    if not os.path.isfile(pdb):
        print("Error: archivo no encontrado.")
        return
    
    pdb_dir = os.path.dirname(os.path.abspath(pdb))
    
    generar_topologia(prot, pdb, pdb_dir)

    modificar_topologia(prot, pdb_dir)

    # membrana = elegir_membrana()

    



if __name__ == "__main__":
    main()
