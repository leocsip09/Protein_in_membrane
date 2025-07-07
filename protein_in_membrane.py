import os
import subprocess

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

def main():
    print("====== Proteina en membrana con GROMACS automatizada c; ======\n")

    prot = input("¿Cuál es el nombre de la proteína? (sin extensión, el nombre que va antes de .pdb): ").strip()

    pdb = input("Ruta al archivo .pdb de la proteína: ").strip()
    if not os.path.isfile(pdb):
        print("Error: archivo no encontrado.")
        return


    # membrana = elegir_membrana()

    # Paso 1: Generar topología
    print("\n====== Paso 1: Generar topología de la proteína ======")
    print("Para esta ocasión vamos a escoger GROMOS96 53A6 para la topología de la proteína (Opción 13)")
    print("Para los grupos acetilo y amida, selecciona 'None' (Opción 2 para ambas) \n ¿Estás de acuerdo? (s/n): ", end="")
    if input().strip().lower() != 's':
        print("Operación cancelada.")
        return
    run(f"gmx pdb2gmx -f {pdb} -o {prot}_processed.gro -ignh -ter -water spc",
        "Topología generada, yeiiii!")



if __name__ == "__main__":
    main()
