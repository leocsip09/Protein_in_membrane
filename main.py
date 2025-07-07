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

def modificar_topologia(prot, pdb_dir):
    print("\n====== Paso 2: Modificar la topología de la proteína ======\n")
    destino_ff = os.path.join(pdb_dir, "gromos53a6_lipid.ff")
    print("Modificando campo de fuerza (forcefield.doc) para incluir parámetros de lípidos Berger...")

    run(f"cp -r /usr/local/gromacs/share/gromacs/top/gromos53a6.ff \"{destino_ff}\"")

    forcefield_doc_path = os.path.join(destino_ff, "forcefield.doc")
    with open(forcefield_doc_path, "w") as f:
        f.write("GROMOS96 53a6 force field, extended to include Berger lipid parameters\n")

    print("forcefield.doc modificado correctamente ^v^.")


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
