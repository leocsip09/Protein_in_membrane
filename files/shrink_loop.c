#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// A quien esté leyendo esto:
// ¿Esperabas un script de bash? 
// Se entiende, solo quería usar C por los jajas. ^^

int main(int argc, char *argv[]) {
    char command[2048];
    int i, prev;

    if (argc != 2) {
        fprintf(stderr, "Uso: %s <ruta_a_pdb_dir>\n", argv[0]);
        return 1;
    }

    const char *pdb_dir = argv[1];
    setenv("GMX_MAXBACKUP", "-1", 1);

    sprintf(command, "gmx grompp -f %s/minim_inflategro.mdp -c %s/system_inflated.gro -p %s/topol.top -r %s/system_inflated.gro -o %s/system_inflated_em.tpr -maxwarn 3", 
            pdb_dir, pdb_dir, pdb_dir, pdb_dir, pdb_dir);
    system(command);

    sprintf(command, "gmx mdrun -deffnm %s/system_inflated_em", pdb_dir);
    system(command);

    sprintf(command, "echo 0 | gmx trjconv -s %s/system_inflated_em.tpr -f %s/system_inflated_em.gro -o %s/tmp.gro -pbc mol", 
            pdb_dir, pdb_dir, pdb_dir);
    system(command);

    sprintf(command, "mv %s/tmp.gro %s/system_inflated_em.gro", pdb_dir, pdb_dir);
    system(command);

    for (i = 1; i <= 26; i++) {
        printf("########################################\n");
        printf("#\n# EJECUCIÓN DE LA %d° ITERACIÓN...\n#\n", i);
        printf("########################################\n");

        prev = i - 1;

        if (i == 1) {
            char path[1024];
            sprintf(path, "%s/system_inflated_em.gro", pdb_dir);
            if (access(path, F_OK) != 0) {
                fprintf(stderr, "%s no existe. Saliendo.\n", path);
                return 1;
            }
            sprintf(command, "perl %s/inflategro.pl %s/system_inflated_em.gro 0.95 DPPC 0 %s/system_shrink%d.gro 5 %s/area_shrink%d.dat", 
                    pdb_dir, pdb_dir, pdb_dir, i, pdb_dir, i);
            system(command);
        } else {
            char path[1024];
            sprintf(path, "%s/system_shrink%d_em.gro", pdb_dir, prev);
            if (access(path, F_OK) != 0) {
                fprintf(stderr, "%s no existe. Saliendo.\n", path);
                return 1;
            }
            sprintf(command, "perl %s/inflategro.pl %s/system_shrink%d_em.gro 0.95 DPPC 0 %s/system_shrink%d.gro 5 %s/area_shrink%d.dat", 
                    pdb_dir, pdb_dir, prev, pdb_dir, i, pdb_dir, i);
            system(command);
        }

        sprintf(command, "gmx grompp -f %s/minim_inflategro.mdp -c %s/system_shrink%d.gro -r %s/system_shrink%d.gro -p %s/topol.top -o %s/system_shrink%d_em.tpr -maxwarn 3", 
                pdb_dir, pdb_dir, i, pdb_dir, i, pdb_dir, pdb_dir, i);
        system(command);

        sprintf(command, "gmx mdrun -deffnm %s/system_shrink%d_em", pdb_dir, i);
        system(command);

        sprintf(command, "gmx trjconv -s %s/system_shrink%d_em.tpr -f %s/system_shrink%d_em.gro -o %s/tmp.gro -pbc mol", 
                pdb_dir, i, pdb_dir, i, pdb_dir);
        system(command);

        sprintf(command, "mv %s/tmp.gro %s/system_shrink%d_em.gro", pdb_dir, pdb_dir, i);
        system(command);
    }

    return 0;
}
