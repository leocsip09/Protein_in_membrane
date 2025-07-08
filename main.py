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
    Si limpiar=True y la sección es 'nonbond_params', reemplaza 'HW' por 'H' en las líneas.
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
        originales = [l.replace('HW', 'H') for l in originales]

    contenido += originales

    for l in nuevas_lineas:
        if l.strip().startswith("[") or l.strip().startswith(";") or l.strip() == "":
            continue
        contenido.append(l)

    contenido += lineas[fin:]

    with open(filepath, "w") as f:
        f.writelines(contenido)

def modificar_topologia(pdb_dir):
    print("\n====== Paso 2: Modificar la topología de la proteína ======\n")
    print("Vamos a modificar el campo de fuerza para incluir parámetros de lípidos Berger.")
    print("¿Deseas continuar? (s/n): ", end="")
    if input().strip().lower() != 's':
        print("Operación cancelada.")
        return

    origen_ff = "/usr/local/gromacs/share/gromacs/top/gromos53a6.ff"
    destino_ff = os.path.join(pdb_dir, "gromos53a6_lipid.ff")
    lipid_itp = "membranas/dppc128/lipid.itp"

    print("\nModificando campo de fuerza (forcefield.doc) para incluir parámetros de lípidos Berger...\n")
    run(f"cp -r \"{origen_ff}\" \"{destino_ff}\"")

    with open(os.path.join(destino_ff, "forcefield.doc"), "w") as f:
        f.write("GROMOS96 53a6 force field, extended to include Berger lipid parameters\n")
    print("forcefield.doc modificado correctamente ^^\n")

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
            if start is not None and line.strip().startswith("["):
                end = i
                break
        return lipid_lines[start:end] if start is not None else []

    def limpiar_nonbond_params(lines):
        limpio = []
        borrar = False
        for line in lines:
            if ";; parameters for lipid-GROMOS interactions" in line:
                borrar = True
                continue
            if ";; lipid-SPC/SPCE interactions" in line:
                borrar = False
                limpio.append(line)
                continue
            if not borrar:
                limpio.append(line)
        return limpio

    atomtypes = extraer_bloque("atomtypes")
    nonbond_params = limpiar_nonbond_params(extraer_bloque("nonbond_params"))
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
    tipos_en_lipid = set()
    for line in atomtypes:
        if line.strip().startswith("[") or line.strip().startswith(";") or line.strip() == "":
            continue
        key = line.split()[0]
        tipos_en_lipid.add(key)
        atomtypes_mod.append(tabla_atomtypes.get(key, line))
    for key, linea in tabla_atomtypes.items():
        if key not in tipos_en_lipid:
            atomtypes_mod.append(linea)

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

    with open(ffnb, "r") as f:
        ff_lines = f.readlines()

    ff_lines = [line.replace(" HW", " H") for line in ff_lines]

    if len(ff_lines) >= 2123:
        del ff_lines[1583:2123]

    with open(ffnb, "w") as f:
        f.writelines(ff_lines)

    topol = os.path.join(pdb_dir, "topol.top")
    if os.path.exists(topol):
        with open(topol, "r") as f:
            lines = f.readlines()

        new_lines = []
        inserted_dppc = False
        for i, line in enumerate(lines):
            line = line.replace('gromos53a6.ff/forcefield.itp', 'gromos53a6_lipid.ff/forcefield.itp')
            line = line.replace('gromos53a6.ff/spc.itp', 'gromos53a6_lipid.ff/spc.itp')
            line = line.replace('gromos53a6.ff/ions.itp', 'gromos53a6_lipid.ff/ions.itp')
            new_lines.append(line)

            if not inserted_dppc and line.strip() == "#endif":
                if i >= 2 and lines[i - 2].strip() == "#ifdef POSRES":
                    new_lines.append("\n\n; Strong position restraints for InflateGRO\n")
                    new_lines.append("#ifdef STRONG_POSRES\n")
                    new_lines.append('#include "strong_posre.itp"\n')
                    new_lines.append("#endif\n")
                    new_lines.append("\n; Include DPPC chain topology\n")
                    new_lines.append('#include "dppc.itp"\n')
                    inserted_dppc = True

        with open(topol, "w") as f:
            f.writelines(new_lines)

        print("\ntopol.top actualizado con include de DPPC, SPC e IONS corregidos.")

    print("\nCampo de fuerza y topología modificados correctamente ^v^.")

def actualizar_topologia_con_dppc(topol_path, gro_path_system, gro_path_protein, n_atoms_por_dppc=50):
    with open(gro_path_system, 'r') as f:
        n_atoms_total = len(f.readlines()) - 3

    with open(gro_path_protein, 'r') as f:
        n_atoms_protein = len(f.readlines()) - 3

    print(f"\nSe toman: {n_atoms_total} átomos del sistema, de los cuales {n_atoms_protein} son de la proteína.")
    n_dppc = (n_atoms_total - n_atoms_protein) // n_atoms_por_dppc
    print(f"\nSe asume que cada DPPC tiene {n_atoms_por_dppc} átomos, por lo que se añadirán {n_dppc} moléculas de DPPC.")

    with open(topol_path, 'r') as f:
        lines = f.readlines()

    mol_section_start = None
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            mol_section_start = i
            break
    if mol_section_start is None:
        return

    protein_index = None
    for i in range(mol_section_start + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped.startswith(";") or stripped == "":
            continue
        if "protein" in stripped.lower():
            protein_index = i
            break
    if protein_index is None:
        return

    if any("dppc" in line.lower() for line in lines[protein_index+1:]):
        return

    lines.insert(protein_index + 1, f"DPPC\t{n_dppc}\n")

    with open(topol_path, "w") as f:
        f.writelines(lines)

def eliminar_agua_y_actualizar_topologia(pdb_dir, topol_path):
    print("\n¿Desea eliminar el excedente de moléculas de agua? (Es necesario) (s/n)\n")
    if input().strip().lower() == 's':
        cmd = [
            "perl", "files/water_deletor.pl",
            "-in", os.path.join(pdb_dir, "system_solv.gro"),
            "-out", os.path.join(pdb_dir, "system_solv_fix.gro"),
            "-ref", "O33",
            "-middle", "C50",
            "-nwater", "3"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)  # Mostrar mensaje del script Perl

        match = re.search(r"(\d+)\s+water molecules remain", result.stdout)
        if not match:
            raise ValueError("No se pudo encontrar el número de moléculas de agua restantes en la salida.")

        n_sol = int(match.group(1))

        with open(topol_path, "r") as f:
            lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].strip().startswith("SOL"):
                lines[i] = f"SOL    {n_sol}\n"
                break
        else:
            lines.append(f"SOL    {n_sol}\n")

        with open(topol_path, "w") as f:
            f.writelines(lines)


def caja_y_solvatar(prot, pdb_dir):
    print("\n====== Paso 3: Caja de solvatación ======\n")

    print("Vamos a preparar la caja con la membrana DPPC y solvatación.")
    if input("¿Deseas continuar? (s/n): ").strip().lower() != 's':
        print("Operación cancelada.")
        return
    
    origen = os.path.join("membranas", "dppc128", "dppc.itp")
    destino = os.path.join(pdb_dir, "dppc.itp")
    shutil.copyfile(origen, destino)

    minim_path = os.path.join(os.getcwd(), "files", "minim.mdp")
    dppc_path = os.path.join(os.getcwd(), "membranas", "dppc128", "dppc128.pdb")
    top_path = os.path.join(os.getcwd(), "files", "topol_dppc.top")
    tpr_out = os.path.join(pdb_dir, "dppc.tpr")
    whole_out = os.path.join(pdb_dir, "dppc128_whole.gro")

    print("\nOrientando la proteina y la membrana...")

    topol_dppc_in_pdb = os.path.join(pdb_dir, "topol_dppc.top")
    shutil.copyfile(top_path, topol_dppc_in_pdb)

    run(f"gmx grompp -f {minim_path} -c {dppc_path} -p {topol_dppc_in_pdb} -o {tpr_out} -maxwarn 3")

    print(f"\nSeleccionar grupo 0, 'System' para la siguiente pregunta.")
    input("¿Entendido? Presiona Enter para continuar...")
    run(f"gmx trjconv -s {tpr_out} -f {dppc_path} -o {whole_out} -pbc mol -ur compact")

    print("\nInsertando la proteína en la caja...")

    usar_defecto = input("¿Deseas usar las dimensiones por defecto de la caja de DPPC128? (s/n): ").strip().lower()
    if usar_defecto == 's':
        box_coords = "6.41840 6.44350 6.59650"
    else:
        while True:
            coords_input = input("Introduce las dimensiones de la caja (x y z en nm, separadas por espacios): ").strip()
            partes = coords_input.split()
            if len(partes) == 3:
                try:
                    float(partes[0]); float(partes[1]); float(partes[2])
                    box_coords = coords_input
                    break
                except ValueError:
                    print("Error: Las coordenadas deben ser números.")
            else:
                print("Error: Debes ingresar 3 valores separados por espacio.")

    prot_gro = os.path.join(pdb_dir, f"{prot}_processed.gro")
    newbox_gro = os.path.join(pdb_dir, f"{prot}_newbox.gro")

    run(f"gmx editconf -f {prot_gro} -o {newbox_gro} -c -box {box_coords}")

    print(f"\nCaja creada con dimensiones: {box_coords} y proteína centrada.\n")
    print("Empaquetando los lípidos alrededor de la proteína...\n")

    origenInf = os.path.join("files", "inflategro.pl")
    destinoInf = os.path.join(pdb_dir, "inflategro.pl")
    shutil.copyfile(origenInf, destinoInf)

    run(f"cat {pdb_dir}/{prot}_newbox.gro {pdb_dir}/dppc128_whole.gro > {pdb_dir}/system.gro")
    print("\nGenerando nuevo archivo de restricción de posición usando genrestr...\n")
    print("En la siguiente opción, selecciona el Grupo 1 (1), que corresponde a restringir la posición de la proteína nvn/.\n")
    input("\n¿Entendido? Presiona Enter para continuar...")
    run(f"gmx genrestr -f {pdb_dir}/{prot}_newbox.gro -o {pdb_dir}/strong_posre.itp -fc 100000 100000 100000")

    with open(minim_path, "r") as f:
        contenido = f.readlines()

    if not contenido or "define = -DSTRONG_POSRES" not in contenido[0]:
        contenido.insert(0, "define = -DSTRONG_POSRES\n")

    with open(minim_path, "w") as f:
        f.writelines(contenido)
    

    origen_minimInf = os.path.join("files", "minim_inflategro.mdp")
    destino_minimInf = os.path.join(pdb_dir, "minim_inflategro.mdp")
    shutil.copyfile(origen_minimInf, destino_minimInf)

    print("\nInflando el sistema...\n")

    run(f"perl {pdb_dir}/inflategro.pl {pdb_dir}/system.gro 4 DPPC 14 {pdb_dir}/system_inflated.gro 5 area.dat")
    actualizar_topologia_con_dppc(
        topol_path=os.path.join(pdb_dir, "topol.top"),
        gro_path_system=os.path.join(pdb_dir, "system_inflated.gro"),
        gro_path_protein=os.path.join(pdb_dir, f"{prot}_processed.gro")
    )

    origen_area = os.path.join("area.dat")
    destino_area = os.path.join(pdb_dir, "area.dat")
    shutil.move(origen_area, destino_area)

    run(f"gmx grompp -f {pdb_dir}/minim_inflategro.mdp -c {pdb_dir}/system_inflated.gro -p {pdb_dir}/topol.top -r {pdb_dir}/system_inflated.gro -o {pdb_dir}/system_inflated_em.tpr -maxwarn 3")
    run(f"gmx mdrun -deffnm {pdb_dir}/system_inflated_em")

    print("\nVamos a corregir el PBC (Periodic Boundary Conditions) del sistema inflado...\n")
    print("En la siguiente opción, selecciona el Grupo 0 (0), que corresponde a todo el sistema.\n")
    input("\n¿Entendido? Presiona Enter para continuar...")

    run(f"gmx trjconv -s {pdb_dir}/system_inflated_em.tpr -f {pdb_dir}/system_inflated_em.gro -o {pdb_dir}/tmp.gro -pbc mol")
    run(f"mv {pdb_dir}/tmp.gro {pdb_dir}/system_inflated_em.gro")
    run(f"perl {pdb_dir}/inflategro.pl {pdb_dir}/system_inflated_em.gro 0.95 DPPC 0 {pdb_dir}/system_shrink1.gro 5 {pdb_dir}/area_shrink1.dat")

    print("\nReduciendo el tamaño de la caja...\nVamos a hacer lo mismo 26 veces :p")
    input("¿Entendido? Presiona Enter para continuar...")

    run("gcc -o shrink_loop files/shrink_loop.c")

    run(f"./shrink_loop {pdb_dir}")

    print("\nEmpezando la solvatación con agua...\n")
    run(f"gmx solvate -cp {pdb_dir}/system_shrink26_em.gro -cs spc216.gro -o {pdb_dir}/system_solv.gro -p {pdb_dir}/topol.top")
    print("\n¿Deseas visualizar el sistema con VMD? (s/n): ", end="")
    if input().strip().lower() == 's':
        if shutil.which("vmd"):
            run(f"vmd {pdb_dir}/system_solv.gro")
        else:
            print("VMD no está instalado o no está en el PATH.")

    eliminar_agua_y_actualizar_topologia(pdb_dir, os.path.join(pdb_dir, "topol.top"))

    print("\n¿Deseas visualizar cómo se ve el sistema actualmente con VMD? (s/n): ", end="")
    if input().strip().lower() == 's':
        if shutil.which("vmd"):
            run(f"vmd {pdb_dir}/system_solv_fix.gro")
        else:
            print("VMD no está instalado o no está en el PATH.")


    origen_mdout = os.path.join("mdout.mdp")
    destino_mdout = os.path.join(pdb_dir, "mdout.mdp")
    shutil.move(origen_mdout, destino_mdout)




def add_ions(pdb_dir):
    print("\n====== Paso 4: Añadir iones al sistema ======\n")
    print("\nVamos a añadir iones al sistema para neutralizar las cargas.")
    input("¿Entendido? Presiona Enter para continuar...")

    run(f"gmx grompp -f files/ions.mdp -c {pdb_dir}/system_solv_fix.gro -p {pdb_dir}/topol.top -o {pdb_dir}/ions.tpr -maxwarn 3")
    print("\nSelecciona el grupo 15 (SOL) para la siguiente pregunta.")
    run(f"gmx genion -s {pdb_dir}/ions.tpr -o {pdb_dir}/system_solv_ions.gro -p {pdb_dir}/topol.top -pname NA -nname CL -neutral")
    print("\nIones añadidos al sistema y topología actualizada correctamente.\n")

def verificar_minimizacion(pdb_dir):
    print("\nVerificando la minimización energética (Epot y Fmax)...")

    xvg_path = os.path.join(pdb_dir, "em_ener.xvg")
    edr_path = os.path.join(pdb_dir, "em.edr")

    subprocess.run(f"echo 10 11 | gmx energy -f {edr_path} -o {xvg_path}", shell=True, check=True)

    with open(xvg_path, "r") as f:
        lines = [line.strip() for line in f if not line.startswith(("#", "@"))]
        if not lines:
            print("No se encontraron datos en em_ener.xvg")
            return
        last_line = lines[-1].split()
        time = float(last_line[0])
        epot = float(last_line[1])
        fmax = float(last_line[2])
        print(f"\nÚltimo paso de minimización:\nTiempo = {time} ps\nEpot = {epot:.2f} kJ/mol\nFmax = {fmax:.2f} kJ/mol/nm")
    print("\n ¿Deseas visualizar el gráfico de energía potencial y fuerza máxima? (s/n): ", end="")
    if input().strip().lower() == 's':
        print("\nAbriendo gráfico con xmgrace...\n")
        subprocess.run(f"xmgrace {xvg_path}", shell=True)

def minimizacion_energia(pdb_dir):
    print("\n====== Paso 5: Minimización energética ======\n")
    print("\nVamos a minimizar la energía del sistema para ser simulable c:.")
    input("¿Entendido? Presiona Enter para continuar...")
    run(f"gmx grompp -f files/minim_en.mdp -c {pdb_dir}/system_solv_ions.gro -p {pdb_dir}/topol.top -o {pdb_dir}/em.tpr -maxwarn 3")
    run(f"gmx mdrun -deffnm {pdb_dir}/em")
    print("Minimización completada, verificando que todo esté en orden...")
    verificar_minimizacion(pdb_dir)
    print("\nMinimización energética completada correctamente. El sistema está listo para la simulación.\n")



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

    modificar_topologia(pdb_dir)

    caja_y_solvatar(prot, pdb_dir)

    add_ions(pdb_dir)

    minimizacion_energia(pdb_dir)

    
if __name__ == "__main__":
    main()
