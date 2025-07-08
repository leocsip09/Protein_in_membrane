import os
import subprocess

def run(cmd, msg=None):
    print(f"\n> Ejecutando: {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    if msg:
        print(msg)


def analizar_sistema_membrana(pdb_dir):
    print("\n====== Paso 8: Análisis del sistema de membrana ======\n")
    print("Vamos a realizar algunos análisis del sistema de membrana.")
    os.chdir(pdb_dir)

    #Este es temporal, luego usaré solo run()
    def ejecutar_comando(comando):
        print(f"\n> Ejecutando: {' '.join(comando)}\n")
        subprocess.run(comando)

    def visualizar_grafica():
        if input("¿Deseas visualizar la gráfica con xmgrace? (s/n): ").lower() == 's':
            archivo = input("Nombre del archivo .xvg a visualizar: ")
            run(f"xmgrace {archivo}")

    def visualizar_sistema():
        vis = input("¿Deseas visualizar el sistema? (s/n): ").lower()
        if vis == 's':
            visor = input("¿Con qué visor? (1) VMD o (2) ChimeraX: ")
            archivo = input("Archivo a visualizar (ej. md_0_1.gro): ")
            if visor == '1':
                ejecutar_comando(['vmd', archivo])
            elif visor == '2':
                run(f"gmx edit conf -f {archivo} -o sistema.pdb")
                print("Sistema convertido a .pdb (Ejecutar desde el directorio) :3")
            else:
                print("Opción inválida")

    def deuterium_order():
        ejecutar_comando(['gmx', 'make_ndx', '-f', 'md_0_1.tpr', '-o', 'sn1.ndx'])
        ejecutar_comando(['gmx', 'order', '-s', 'md_0_1.tpr', '-f', 'md_0_1.xtc', '-n', 'sn1.ndx', '-d', 'z', '-od', 'deuter_sn1.xvg'])
        visualizar_grafica()

    def densidad_membrana():
        ejecutar_comando(['gmx', 'make_ndx', '-f', 'md_0_1.tpr', '-o', 'density_groups.ndx'])
        grupos = [('Headgroups', 'dens_headgroups.xvg'), ('Glycerol_Ester', 'dens_glyc.xvg'), ('Acyl_Chains', 'dens_acyl.xvg')]
        for nombre, archivo in grupos:
            print(f"\n→ Analizando densidad del grupo: {nombre}")
            ejecutar_comando(['gmx', 'density', '-s', 'md_0_1.tpr', '-f', 'md_0_1.xtc', '-n', 'density_groups.ndx', '-o', archivo, '-d', 'Z'])
            visualizar_grafica()

    def diffusion_lateral():
        ejecutar_comando(['gmx', 'make_ndx', '-f', 'md_0_1.tpr', '-o', 'p8.ndx'])
        ejecutar_comando(['gmx', 'msd', '-s', 'md_0_1.tpr', '-f', 'md_0_1.xtc', '-n', 'p8.ndx', '-lateral', 'z'])
        visualizar_grafica()

    def menu():
        while True:
            print("\nAnálisis disponibles:")
            print("1. Ver sistema en VMD o ChimeraX")
            print("2. Deuterium Order Parameters")
            print("3. Densidad de la Membrana")
            print("4. Difusión lateral de lípidos")
            print("5. Salir")
            opcion = input("Selecciona una opción: ")

            if opcion == '1':
                visualizar_sistema()
            elif opcion == '2':
                deuterium_order()
            elif opcion == '3':
                densidad_membrana()
            elif opcion == '4':
                diffusion_lateral()
            elif opcion == '5':
                print("Saliendo del análisis.")
                break
            else:
                print("Opción inválida")

    menu()

pdb = input("Introduce la ruta al archivo .pdb de la proteína: ").strip()
if not os.path.isfile(pdb):
    print("Error: archivo no encontrado.")
pdb_dir = os.path.dirname(os.path.abspath(pdb))
analizar_sistema_membrana(pdb_dir)