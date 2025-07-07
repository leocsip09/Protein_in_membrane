import os
import subprocess
import shutil

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

    # membrana = elegir_membrana()

    



if __name__ == "__main__":
    main()
