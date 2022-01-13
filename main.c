#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/random.h>

#include <picosat.h>

#define STR_(x) #x
#define STR(x) STR_(x)
#define MAX_SOLUTIONS 1000

#define MAX_ATOMS 512
#define MAX_BONDS 512
#define MAX_CUT_EDGES (64*1024*1024)

#define min(a_, b_) ((a_) < (b_) ? (a_) : (b_))
#define max(a_, b_) ((a_) > (b_) ? (a_) : (b_))
#define array_length(arr_) (sizeof(arr_) / sizeof((arr_)[0]))
#define DEBUG 0

typedef enum {
    ATOM_UNSPECIFIED,
    ATOM_H,
    ATOM_O,
    ATOM_N,
    ATOM_C,
    ATOM_MAX
} AtomKind;

typedef struct {
    int x;
    int y;
    AtomKind kind;
    int mark;
} Atom;

Atom atoms[512] = {};
int  atoms_count = 0;
int unspecified_atom_ids[512] = {};
int  unspecified_atoms_count = 0;


typedef enum {
    BOND_MINUS,
    BOND_SLASH,
    BOND_BACKSLASH
} BondKind;

typedef struct {
    int x;
    int y;
    BondKind kind;
    int atom_id1;
    int atom_id2;
    int solution;
} Bond;

Bond bonds[512];
int bonds_count = 0;

typedef struct {
    int count;
    int bond_ids[3];
    int atom_ids[3];
} Link;

void add_atom(int x, int y, AtomKind kind) {
    atoms[atoms_count] = (Atom){
        .x = x,
        .y = y,
        .kind = kind
    };
    if(kind == ATOM_UNSPECIFIED) unspecified_atom_ids[unspecified_atoms_count++] = atoms_count;
    atoms_count++;
}

int get_atom(int x, int y) {
    for(int i = 0; i < atoms_count; i++) {
        if(atoms[i].x == x && atoms[i].y == y) {
            return i;
        }
    }
    return -1;
}

void add_bond(int x, int y, BondKind kind) {
    bonds[bonds_count] = (Bond){
        .x = x,
        .y = y,
        .kind = kind
    };
    bonds_count++;
}

int cut_edges[MAX_CUT_EDGES];
int cut_edges_count = 0;
void add_to_cut_set(int bond) {
    assert(cut_edges_count < MAX_CUT_EDGES);
    cut_edges[cut_edges_count++] = bond;
}

Link get_atom_link(int atom_id) {
    Link result = {};
    for(int i = 0; i < bonds_count; i++) {
        if(bonds[i].atom_id1 == atom_id || bonds[i].atom_id2 == atom_id) {
            result.bond_ids[result.count] = i;
            result.atom_ids[result.count] = bonds[i].atom_id1 == atom_id ? bonds[i].atom_id2 : bonds[i].atom_id1;
            result.count++;
        }
    }
    return result;
}


void parse(char* puzzle, int len) {
    int x = 0;
    int y = 0;

    //printf("%s\n", puzzle);
    
    for(int i = 0; i < len; i++) {
        char c = puzzle[i];
        switch(c) {
        case ' ' : {                                   x++; break; }
        case '\n': { x = 0;                            y++; break; }
        case 'X' : { add_atom(x, y, ATOM_UNSPECIFIED); x++; break; }
        case 'H' : { add_atom(x, y, ATOM_H);           x++; break; }
        case 'O' : { add_atom(x, y, ATOM_O);           x++; break; }
        case 'N' : { add_atom(x, y, ATOM_N);           x++; break; }
        case 'C' : { add_atom(x, y, ATOM_C);           x++; break; }
        case '-' : { add_bond(x, y, BOND_MINUS);       x++; break; }
        case '/' : { add_bond(x, y, BOND_SLASH);       x++; break; }
        case '\\': { add_bond(x, y, BOND_BACKSLASH);   x++; break; }
        default  : { printf("\"%c\"=%02x\n", c, c);    assert(0); break; }
        }
    }

    for(int i = 0; i < bonds_count; i++) {
        Bond* b = &bonds[i];
        switch(b->kind) {
        case BOND_MINUS: {
            b->atom_id1 = get_atom(b->x - 1, b->y);
            b->atom_id2 = get_atom(b->x + 1, b->y);
            break;
        }
        case BOND_SLASH: {
            b->atom_id1 = get_atom(b->x - 1, b->y + 1);
            b->atom_id2 = get_atom(b->x + 1, b->y - 1);
            break;
        }
        case BOND_BACKSLASH: {
            b->atom_id1 = get_atom(b->x - 1, b->y - 1);
            b->atom_id2 = get_atom(b->x + 1, b->y + 1);
            break;
        }
        }
    }
}


void print(int show_solution) {
    char* bond_names_simple[3] = { "-", "/", "\\" };
    char* bond_names[4][4] = {
        { " ", "-", "=", "≡" },
        { " ", "/", "⫽", "⫻" },
        { " ", "\\", "⑊", "3" },
        };
    char* atom_names[] = { "X", "H", "O", "N", "C" };
    for(int i = 0; i < atoms_count; i++) {
        printf("\033[%i;%iH%s", atoms[i].y+1, atoms[i].x+1, atom_names[atoms[i].kind]);
    }
    for(int i = 0; i < bonds_count; i++) {
        if(show_solution) {
            printf("\033[%i;%iH%s", bonds[i].y+1, bonds[i].x+1, bond_names[bonds[i].kind][bonds[i].solution]);
        } else {
            printf("\033[%i;%iH%s", bonds[i].y+1, bonds[i].x+1, bond_names_simple[bonds[i].kind]);
        }
    }
}

void print_file(FILE* fp) {
    char* bond_names[3] = { "-", "/", "\\" };
    char* atom_names[] = { "X", "H", "O", "N", "C" };

    // Find bounds:
    int min_x = 99999;
    int max_x = 0;
    int min_y = 99999;
    int max_y = 0;

    for(int i = 0; i < atoms_count; i++) {
        min_x = min(atoms[i].x, min_x);
        min_y = min(atoms[i].y, min_y);
        max_x = max(atoms[i].x, max_x);
        max_y = max(atoms[i].y, max_y);
    }
    
    for(int y = min_y; y <= max_y; y++) {
        for(int x = min_x; x <= max_x; x++) {
            int space = 1;
            for(int i = 0; i < atoms_count; i++) {
                if(atoms[i].x == x && atoms[i].y == y) { fprintf(fp, "%s", atom_names[atoms[i].kind]); space = 0; }
            }
            for(int i = 0; i < bonds_count; i++) {
                if(bonds[i].x == x && bonds[i].y == y) { fprintf(fp, "%s", bond_names[bonds[i].kind]); space = 0; }
            }
            if(space) fprintf(fp, " ");
        }
        fprintf(fp, "\n");
    }
}

void print_latex(char* file_name, int show_solution) {
    FILE* fp = fopen(file_name, "wb");
    char a_names[] = { 'X', 'H', 'O', 'N', 'C' };
    fprintf(fp,
            "\\documentclass{amsart}\n"
            "\\usepackage[utf8]{inputenc}\n"
            "\\usepackage[ngerman]{babel}\n"
            "\\pagestyle{empty}\n"
            "\\usepackage[a4paper,top=1.2cm,centering]{geometry}\n"
            "\\usepackage{tikz}\n"
            "\\setlength{\\parindent}{0pt}\n"
            "\\newlength\\triplesep\n"
            "\\newlength\\triplelinewidth\n"
            "\\setlength\\triplesep{1pt}\n"
            "\\setlength\\triplelinewidth{0.4pt}\n"
            "\\tikzset{s0/.style={line width=\\triplelinewidth,white}}\n"
            "\\tikzset{s1/.style={line width=\\triplelinewidth,black}}\n"
            "\\tikzset{s2/.style={line width=\\triplesep,white, preaction={\n"
            "      draw,line width=\\triplesep+2\\triplelinewidth,black} } }\n"
            "\\tikzset{s3/.style={line width=\\triplelinewidth,black, preaction={\n"
            "      preaction={draw,line width=2\\triplesep+3\\triplelinewidth,black},\n"
            "      draw,line width=2\\triplesep+\\triplelinewidth,white} } }\n"
            "\\begin{document}\n"
            "\\begin{tikzpicture}[scale=0.5, yscale=0.866]\n"
        );
    for(int i = 0; i < atoms_count; i++) {
        int x = atoms[i].x;
        x = (x/2) + (x/4);
        int y = atoms[i].y;
        fprintf(fp, "\\node (a%i) at (%i, %i) {%c};\n", i, x, y, a_names[atoms[i].kind]);
    }

    if(show_solution) {
        for(int i = 0; i < bonds_count; i++) {
            fprintf(fp, "\\draw[s%i] (a%i) -- (a%i);\n", bonds[i].solution, bonds[i].atom_id1, bonds[i].atom_id2);
        }
    } else {
        for(int i = 0; i < bonds_count; i++) {
            fprintf(fp, "\\draw[dotted] (a%i) -- (a%i);\n", bonds[i].atom_id1, bonds[i].atom_id2);
        }
    }

    fprintf(fp, "\\end{tikzpicture}\n"
           "\\end{document}\n");
    fclose(fp);
}

int get_bond_literal(int bond_id, int bond_order) {
    return bond_id * 3 + bond_order;
}

void synthesize_sum_ieq(PicoSAT* ps, int* arr, int n, int r, int sign) {
    if(n == 2) {
        for(int i1 = 0; i1 <= 4; i1++) {
            for(int i2 = 0; i2 <= 4; i2++) {
                int b1 = sign*get_bond_literal(arr[0], i1);
                int b2 = sign*get_bond_literal(arr[1], i2);
                int lits[3] = { 0 };
                int k = 0;
                if(i1 > 0 && i1 < 4) lits[k++] = b1;
                if(i2 > 0 && i2 < 4) lits[k++] = b2;
                if(i1 + i2 == r + 1) {
                    if((sign < 0 && i1 < 4 && i2 < 4)  || (sign > 0 && i1 != 0 && i2 != 0)) {
                        if(DEBUG) printf("%2i: %i %i\n", sign*r, i1, i2);
                        picosat_add_lits(ps, lits);
                    }
                }
            }
        }
    } else {
        assert(n == 3);
        for(int i1 = 0; i1 <= 4; i1++) {
            for(int i2 = 0; i2 <= 4; i2++) {
                for(int i3 = 0; i3 <= 4; i3++) {
                    int b1 = sign*get_bond_literal(arr[0], i1);
                    int b2 = sign*get_bond_literal(arr[1], i2);
                    int b3 = sign*get_bond_literal(arr[2], i3);
                    int lits[4] = { 0 };
                    int k = 0;
                    if(i1 != 0 && i1 < 4) lits[k++] = b1;
                    if(i2 != 0 && i2 < 4) lits[k++] = b2;
                    if(i3 != 0 && i3 < 4) lits[k++] = b3;
                    if(i1 + i2 + i3 == r + 1 && sign < 0 && i1  < 4 && i2  < 4 && i3  < 4) {
                        if(DEBUG) printf("%2i: %i %i %i | %2i %2i %2i | %2i %2i %2i\n", sign*r, i1, i2, i3, b1, b2, b3, lits[0], lits[1], lits[2]);
                        picosat_add_lits(ps, lits);
                    }
                    if(i1 + i2 + i3 == r + 2 && sign > 0 && i1 != 0 && i2 != 0 && i3 != 0) {
                        if(DEBUG) printf("%2i: %i %i %i | %2i %2i %2i | %2i %2i %2i\n", sign*r, i1, i2, i3, b1, b2, b3, lits[0], lits[1], lits[2]);
                        picosat_add_lits(ps, lits);
                    }
                }
            }
        }
    }
}

#if DEBUG
void test_synthesize_sum_ieq(void) {
    PicoSAT* ps;
    int bonds[] = { 0, 1, 2 };
    for(int n = 2; n < 4; n++) {
        for(int r = 1; r <= 4; r++) {
            ps = picosat_init();
            picosat_adjust(ps, 3*n);
            for(int i = 0; i < n; i++) {
                int e1 = get_bond_literal(bonds[i], 1);
                int e2 = get_bond_literal(bonds[i], 2);
                int e3 = get_bond_literal(bonds[i], 3);
                picosat_add_arg(ps, e1, -e2, 0);        
                picosat_add_arg(ps, e2, -e3, 0);        
            }

            for(int sign = -1; sign < 2; sign += 2) {
                synthesize_sum_ieq(ps, bonds, n, r, sign);
            }
            //picosat_print(ps, stdout);
            int variable_count = picosat_variables(ps);
            while(1) {
                int result = picosat_sat(ps, -1);
                if(result != PICOSAT_SATISFIABLE) break;
                int negative_assignment[variable_count + 1];
                int set_variables = 0;
                for(int i = 1; i <= variable_count; i++) {
                    int x = picosat_deref(ps, i);
                    
                    if(x ==  1) negative_assignment[set_variables++] = -i;
                    if(x == -1) negative_assignment[set_variables++] =  i;
                    printf("%c ", x == 1 ? '+' : x == -1 ? '-' : '?');
                }
                negative_assignment[set_variables++] = 0;
                printf("\n");
                picosat_add_lits(ps, negative_assignment);
            }
            picosat_reset(ps);
        }
    }
}
#endif


typedef struct {
    int num_solutions;
    int num_decisions;
} SolveValue;

SolveValue solve(int max_solutions) {
    PicoSAT* ps = picosat_init();
    picosat_set_seed(ps, time(0));
    picosat_set_global_default_phase(ps, 3);

    for(int i = 0; i < cut_edges_count; i++) {
        picosat_add(ps, cut_edges[i]);
    }
    
    for(int i = 0; i < bonds_count; i++) {
        int e1 = get_bond_literal(i, 1);
        int e2 = get_bond_literal(i, 2);
        int e3 = get_bond_literal(i, 3);
        picosat_add_arg(ps, e1, -e2, 0);        
        picosat_add_arg(ps, e2, -e3, 0);        
    }
    for(int i = 0; i < atoms_count; i++) {
        Link link = get_atom_link(i);
        if(atoms[i].kind != ATOM_UNSPECIFIED) {
            synthesize_sum_ieq(ps, link.bond_ids, link.count, atoms[i].kind,  1);
            synthesize_sum_ieq(ps, link.bond_ids, link.count, atoms[i].kind, -1);
        } else {
            synthesize_sum_ieq(ps, link.bond_ids, link.count, 4, -1);
        }
    }

    int variable_count = picosat_variables(ps);
    int num_decisions = 0;
    int solutions_count = 0;
    while(solutions_count < max_solutions) {
        int result = picosat_sat(ps, -1);

        if(result != PICOSAT_SATISFIABLE) break;
        for(int i = 0; i < atoms_count; i++) {
            atoms[i].mark = 0;
        }
        int connected_atoms[MAX_ATOMS] = { 0 };
        int connected_atoms_count = 1;
        atoms[0].mark = 1;

        int new_cut_edges[MAX_BONDS + 1] = { };
        int new_cut_edges_count = 0;
        int visited_atoms = 0;
        while(visited_atoms < connected_atoms_count) {
            int atom = connected_atoms[visited_atoms];
            Link link = get_atom_link(atom);
            for(int i = 0; i < link.count; i++) {
                if(!atoms[link.atom_ids[i]].mark) {
                    if(picosat_deref(ps, get_bond_literal(link.bond_ids[i], 1)) == 1) {
                        atoms[link.atom_ids[i]].mark = 1;
                        connected_atoms[connected_atoms_count++] = link.atom_ids[i];
                    } else {
                        new_cut_edges[new_cut_edges_count++] = link.bond_ids[i];
                    }
                }
            }
            visited_atoms++;
        }
        for(int i = 0; i < new_cut_edges_count; ) {
            int a1 = bonds[new_cut_edges[i]].atom_id1;
            int a2 = bonds[new_cut_edges[i]].atom_id2;
            int sum = atoms[a1].mark + atoms[a2].mark;
            if(sum == 2) {
                new_cut_edges[i] = new_cut_edges[--new_cut_edges_count];
            } else {
                assert(sum == 1);
                new_cut_edges[i] = get_bond_literal(new_cut_edges[i], 1);
                i++;
            }
        }
        
        if(connected_atoms_count == atoms_count) {
            if(solutions_count == 0) num_decisions = picosat_decisions(ps);
            for(int i = 0; i < bonds_count; i++) {
                bonds[i].solution = 0;
                for(int j = 1; j <= 3; j++) {
                    if(picosat_deref(ps, get_bond_literal(i, j)) == 1) {
                        bonds[i].solution = j;
                    }
                }
            }
            int negative_assignment[variable_count + 1];
            int set_variables = 0;
            solutions_count++;
            for(int i = 1; i <= variable_count; i++) {
                int x = picosat_deref(ps, i);
                
                if(x ==  1) negative_assignment[set_variables++] = -i;
                if(x == -1) negative_assignment[set_variables++] =  i;
            }
            negative_assignment[set_variables++] = 0;
            picosat_add_lits(ps, negative_assignment);
        } else {
            for(int i = 0; i < new_cut_edges_count; i++) {
                add_to_cut_set(new_cut_edges[i]);
                picosat_add(ps, new_cut_edges[i]);
            }
            if(new_cut_edges_count > 0) {
                add_to_cut_set(0);
                picosat_add(ps, 0);
            }
        }
        
    }
    picosat_reset(ps);
    return (SolveValue) {
        .num_solutions = solutions_count,
        .num_decisions = num_decisions,
    };
}

int sample_distribution(int* distribution, int length, int total) {
    int r = rand() % total;
    for(int i = 0; i < length; i++) {
        if(r < distribution[i]) {
            return i;
        } else {
            r -= distribution[i];
        }
    }
    assert(0);
    return 0;
}

void usage(const char* argv0) {
    printf("Usage: %s size [#H #O #N #C]\n"
           "where size is the name of a template file (for example 'medium')\n"
           "and #H, #O, #N, #C are integers indicating the probabilities"
           "that the respective atoms are chosen.\n\n"
           "The program generates a puzzle of the respective size"
           "and writes it into 'puzzle.txt'.\n", argv0);
}

int main(int argc, const char** argv) {
#if DEBUG
    test_synthesize_sum_ieq();
    return 0;
#endif    

    int len = 0;
    char* puzzle = NULL;

    if(argc >= 2) {
        FILE* fp = fopen(argv[1], "r");
        if(fp) {
            fseek(fp, 0, SEEK_END);
            len = ftell(fp);
            fseek(fp, 0, SEEK_SET);
            puzzle = malloc(len);
            fread(puzzle, len, 1, fp);
            fclose(fp);
        }
    }
    if(puzzle == NULL) {
        usage(argv[0]);
        return 1;
    }
    
    system("clear");
    parse(puzzle, len);
    unsigned int seed;
    getrandom(&seed, sizeof(seed), 0);
    srand(seed);

    int atom_scale = 1; (void) atom_scale;
    
    int iterations = 0;
    int old_num_solutions = atom_scale*atoms_count; (void) old_num_solutions;


    int distribution[] = {
        0, 1, 5, 8, 3
    };
    int distribution_length = array_length(distribution);

    if(argc >= 6) {
        distribution[1] = atoi(argv[2]);
        distribution[2] = atoi(argv[3]);
        distribution[3] = atoi(argv[4]);
        distribution[4] = atoi(argv[5]);
    }
    
    int distribution_total = 0;
    for(int i = 0; i < distribution_length; i++) {
        distribution_total += distribution[i];
    }
#define NUM_CHOICES 2
    while(1) {
        AtomKind old_kinds[NUM_CHOICES];
        int indices[NUM_CHOICES];
        for(int i = 0; i < NUM_CHOICES; i++) {
            if(unspecified_atoms_count > 0) {
                indices[i] = unspecified_atom_ids[unspecified_atoms_count-1];
                unspecified_atoms_count--;
            } else {
                indices[i] = rand() % atoms_count;
            }
            old_kinds[i] = atoms[indices[i]].kind;
        }


        for(int i = 0; i < NUM_CHOICES; i++) {
            while(old_kinds[i] == atoms[indices[i]].kind) {
                atoms[indices[i]].kind = sample_distribution(distribution, distribution_length, distribution_total);
            }
        }
        SolveValue solve_value = solve(2);
        int new_num_solutions = solve_value.num_solutions;

        if(new_num_solutions == 1) {
            break;
        } else if(1 <= new_num_solutions) {
            print(1);
            old_num_solutions = new_num_solutions;
        } else {
            for(int i = 0; i < NUM_CHOICES; i++) {
                atoms[indices[i]].kind = old_kinds[i];
                if(old_kinds[i] == ATOM_UNSPECIFIED) unspecified_atoms_count++;
            }
        }
        int atom_kinds_count[5] = { };
        for(int i = 0; i < atoms_count; i++) {
            atom_kinds_count[atoms[i].kind]++;
        }
        printf("\033[34;1H\033[K+ %i %i/%i %i (%i %i %i %i %i)\n", iterations, new_num_solutions, old_num_solutions, cut_edges_count,
               atom_kinds_count[0], atom_kinds_count[1], atom_kinds_count[2], atom_kinds_count[3], atom_kinds_count[4]);
        iterations++;
    }

    system("clear");
    print(1);
    printf("\n");
    FILE* fp = fopen("puzzle.txt", "wb");
    print_file(fp);
    fclose(fp);
    return 0;
}
