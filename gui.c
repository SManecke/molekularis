#include <raylib.h>
#include <rlgl.h>

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>

#include "picosat/picosat.h"

#define ArrayLength(a) (sizeof(a) / sizeof((a)[0]))

#define MAX_ATOMS 512
#define MAX_BONDS 512
#define MAX_CUT_EDGES (64*1024*1024)

typedef float f32;

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

    Vector2 p;
} Atom;

Atom atoms[512] = {};
int  atoms_count = 0;
int  unspecified_atom_ids[512] = {};
int  unspecified_atoms_count = 0;
char* atom_names[] = { "X", "H", "O", "N", "C" };
int  atom_name_widths[] = { 8, 8, 8, 8, 8 };
typedef enum {
    BOND_0   = 0,
    BOND_1   = 1,
    BOND_2   = 2,
    BOND_3   = 3,
    BOND_FIX = 1 << 2,
} Solution;

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
    Solution solution;
} Bond;

Bond bonds[512];
int bonds_count = 0;

f32 zoom = 31;

typedef struct {
    int count;
    int bond_ids[3];
    int atom_ids[3];
} Link;

Rectangle bounding_rect = { };
float width;
float height;

Rectangle RectangleEnlargeToContain(Rectangle rect, Vector2 pos) {
    float minx = fmin(rect.x, pos.x);
    float miny = fmin(rect.y, pos.y);
    float maxx = fmax(rect.x + rect.width, pos.x);
    float maxy = fmax(rect.y + rect.height, pos.y);
    return (Rectangle) { .x = minx, .y = miny, .width = maxx - minx, .height = maxy - miny };
}

void add_atom(int x, int y, AtomKind kind) {
    atoms[atoms_count] = (Atom){
        .x = x,
        .y = y,
        .kind = kind,
        .p.x  = zoom * ((x/2) + (x/4)),
        .p.y  = zoom * 0.866 * y,
    };
    bounding_rect = RectangleEnlargeToContain(bounding_rect, atoms[atoms_count].p);
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

void parse(char* puzzle, int len) {
    int x = 0;
    int y = 0;
    
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

int rewindState = 0;
Solution rewindSolution[64*1024];

Solution ToggleSolution(Solution solution, int inc) {
    int count = solution & 0x3;
    int fixed = solution & 0x4;
    return ((count+inc) & 0x3) | fixed;
}

Solution ToggleFixed(Solution solution) {
    int count = solution & 0x3;
    int fixed = solution & 0x4;
    return count | ((~fixed) & 0x4);
}


int GetAtomAt(Vector2 p) {
    Rectangle rect = bounding_rect;
    float scale = fmax(1, fmax(rect.width / (width - 10), rect.height / (height - 40)));
    p = Vector2Scale(p, scale);
    for(int i = 0; i < atoms_count; i++) {
        if(CheckCollisionPointCircle(p, atoms[i].p, 10)) {
            return i;
        }
    }
    return -1;
}

f32 f32_clamp(f32 a, f32 t, f32 b) {
    if(t < a) return a;
    if(t > b) return b;
    return t;
}

int GetBondAt(Vector2 p) {
    Rectangle rect = bounding_rect;
    float scale = fmax(1, fmax(rect.width / (width - 10), rect.height / (height - 40)));
    p = Vector2Scale(p, scale);
    for(int i = 0; i < bonds_count; i++) {
        Vector2 a = atoms[bonds[i].atom_id1].p;
        Vector2 b = atoms[bonds[i].atom_id2].p;
        Vector2 pa = Vector2Subtract(p, a);
        Vector2 ba = Vector2Subtract(b, a);
        float h = f32_clamp(0, Vector2DotProduct(pa,ba) / Vector2DotProduct(ba,ba), 1);
        if(Vector2Length(Vector2Subtract(pa, Vector2Scale(ba, h))) < 7) {
            return i;
        }
    }
    return -1;
}

Vector2 CenteredPosition(Vector2 pos) {
    Rectangle rect = bounding_rect;
    float scale = fmin(1, fmin((width - 10) / rect.width, (height - 40) / rect.height));
    return (Vector2) {
        .x = width  / 2 + scale*( -rect.width  / 2 + pos.x),
        .y = height / 2 + scale*( -rect.height / 2 + pos.y),
    };  
}

void DrawEdge(Vector2 p1, Vector2 p2, Solution s) {
    Vector2 dir = Vector2Subtract(p2, p1);
    Vector2 n = Vector2Normalize(dir);
    Vector2 o = { .x = -n.y, .y = n.x };

    f32 length = Vector2Length(dir);
    int count = s & 0x3;
    int fix   = s & 0x4;
    if(count == 0 && fix == 0) {
        for(int i = 3; 5*i + 10 < length; i++) {
            Vector2 p = Vector2Add(p1, Vector2Scale(n, 5*i));
            p.x = lroundf(p.x + 0.5) - 0.5;
            p.y = lroundf(p.y + 0.5) - 0.5;

            Vector2 strip[4] = {
                { p.x - 0.5, p.y - 0.5 },
                { p.x - 0.5, p.y + 0.5 },
                { p.x + 0.5, p.y - 0.5 },
                { p.x + 0.5, p.y + 0.5 }
            };

            DrawTriangleStrip(strip, 4, BLACK);
        }
    } else {
        for(int i = 0; i < count; i++) {
            Vector2 q1 = Vector2Add(Vector2Add(p1, Vector2Scale(n,  13)), Vector2Scale(o, 4*i - 2*(count - 1)));
            Vector2 q2 = Vector2Add(Vector2Add(p2, Vector2Scale(n, -13)), Vector2Scale(o, 4*i - 2*(count - 1)));
            if(q1.y == q2.y) {
                q1.y = q2.y = lroundf(q1.y + 0.5) - 0.5;
            }
            DrawLineEx(q1, q2, 0.8, BLACK);
        }
    }
    if(count > 0 && fix) {
        Vector2 center = Vector2Scale(Vector2Add(p1, p2), 0.5);
        Vector2 q1 = Vector2Add(center, Vector2Scale(o,  10));
        Vector2 q2 = Vector2Add(center, Vector2Scale(o, -10));

        if(q1.x == q2.x) {
            q1.x = q2.x = lroundf(q1.x + 0.5) - 0.5;
        }
        DrawLineEx(q1, q2, 0.8, BLACK);
    }
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


int main() {
    FILE* fp = fopen("puzzle.txt", "rb");
    fseek(fp, 0, SEEK_END);
    int length = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char puzzle[4096];
    fread(puzzle, length, 1, fp);
    fclose(fp);
    parse(puzzle, length);

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetTraceLogLevel(LOG_WARNING);
    InitWindow(0, 0, "Molecularis");
    rlDisableBackfaceCulling();
    SetTargetFPS(60);
    SetWindowState(FLAG_WINDOW_RESIZABLE);
    ToggleFullscreen();
    ToggleFullscreen();
    int keepRunning = 1;
    int frameCount = 0;

    Font font = LoadFontEx("/usr/share/fonts/truetype/lyx/cmr10.ttf", 20, NULL, 95);
    
    while(keepRunning) {
        if(WindowShouldClose()) keepRunning = 0;
        width = GetScreenWidth();
        height = GetScreenHeight();
        
        Vector2 mouse = GetMousePosition();

        if(IsKeyPressed(KEY_Q)) keepRunning = 0;
        if(IsKeyPressed(KEY_ESCAPE)) keepRunning = 0;
        if(IsKeyPressed(KEY_F10)) ToggleFullscreen();
        if(IsKeyPressed(KEY_F12)) TakeScreenshot("screenshot.png");
        
        Vector2 offset = CenteredPosition((Vector2){});

        int currentAtom = GetAtomAt(Vector2Subtract(mouse, offset));
        int currentBond = currentAtom == -1 ? GetBondAt(Vector2Subtract(mouse, offset)) : -1;

        int change = 0;
        if(currentBond != -1) {
            if(IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                bonds[currentBond].solution = ToggleSolution(bonds[currentBond].solution, 1);
                change = 1;
            }
            
            if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
                bonds[currentBond].solution = ToggleSolution(bonds[currentBond].solution, -1);
                change = 1;
            }
        
            if(IsMouseButtonPressed(MOUSE_MIDDLE_BUTTON)) {
                bonds[currentBond].solution = ToggleFixed(bonds[currentBond].solution);
                change = 1;
            }
        }

        if(change) {
            rewindState++;
            for(int i = 0; i < bonds_count; i++) {
                rewindSolution[rewindState * bonds_count + i] = bonds[i].solution;
            }
        }

        if(IsKeyPressed(KEY_Z) && rewindState > 0) {
            rewindState--;
            for(int i = 0; i < bonds_count; i++) {
                bonds[i].solution = rewindSolution[rewindState * bonds_count + i];
            }
        }

        
        BeginDrawing();

        int isSolved = 0;
        Color backgroundColor = isSolved ? GREEN : YELLOW;
        ClearBackground(backgroundColor);
        
        for(int i = 0; i < atoms_count; i++) {
            Vector2 center = CenteredPosition(atoms[i].p);
            int kind = atoms[i].kind;
            char* name = atom_names[kind];
            Vector2 pos = Vector2Subtract(center, (Vector2) { atom_name_widths[kind], 7 });
            Color textColor = BLACK;
            //Link link = get_atom_link(atoms(
            DrawTextEx(font, name, pos, 20, 30, textColor);
        }

        for(int i = 0; i < bonds_count; i++) {
            Vector2 p1 = CenteredPosition(atoms[bonds[i].atom_id1].p);
            Vector2 p2 = CenteredPosition(atoms[bonds[i].atom_id2].p);
            Solution s = bonds[i].solution;
            DrawEdge(p1, p2, s);
        }

        EndDrawing();

        frameCount++;
    }
}
