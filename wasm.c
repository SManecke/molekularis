typedef unsigned char u8; 
typedef int           i32; 
typedef unsigned int  u32; 
typedef float         f32; 

enum {
    MOUSE_INVALID = -1,
    MOUSE_LEFT = 0,
    MOUSE_MIDDLE = 1,
    MOUSE_RIGHT = 2,
};

void clear(i32 r, i32 g, i32 b);
void draw_line(f32 x0, f32 y0, f32 x1, f32 y1, f32 thickness);
void draw_dot(f32 x, f32 y);
void draw_char(f32 x, f32 y, i32 c, f32 h);
void debug_break(void);

#define ArrayLength(a) (sizeof(a) / sizeof((a)[0]))
#define assert(cond) do { if(!(cond)) { debug_break; } } while(0)


void *memcpy(void *restrict dest, const void *restrict src, u32 n) {
    for(int i = 0; i < n; i++) {
        ((char*)dest)[i] = ((char*)src)[i];
    }
    return dest;
}
#define MAX_ATOMS 512
#define MAX_BONDS 512
#define MAX_CUT_EDGES (64*1024*1024)

f32 f32_min(f32 a, f32 b) {
    if(a < b) return a;
    return b;
}

f32 f32_max(f32 a, f32 b) {
    if(a > b) return a;
    return b;
}

static const f32 tiny = 1.0e-30;

// From musl:
/* Get a 32 bit int from a float.  */
#define GET_FLOAT_WORD(w,d)                                                     \
    do {                                                                        \
        union {f32 f; u32 i;} __u;                                              \
        __u.f = (d);                                                            \
        (w) = __u.i;                                                            \
    } while (0)

                    /* Set a float from a 32 bit int.  */
#define SET_FLOAT_WORD(d,w)                                                     \
    do {                                                                        \
        union {f32 f; u32 i;} __u;                                              \
        __u.i = (w);                                                            \
        (d) = __u.f;                                                            \
    } while (0)
                                    

f32 f32_sqrt(f32 x) {
    f32 z;
    i32 sign = (int)0x80000000;
    i32 ix,s,q,m,t,i;
    u32 r;

    GET_FLOAT_WORD(ix, x);

    /* take care of Inf and NaN */
    if ((ix&0x7f800000) == 0x7f800000)
        return x*x + x; /* sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN */

    /* take care of zero */
    if (ix <= 0) {
        if ((ix&~sign) == 0)
            return x;  /* sqrt(+-0) = +-0 */
        if (ix < 0)
            return (x-x)/(x-x);  /* sqrt(-ve) = sNaN */
    }
    /* normalize x */
    m = ix>>23;
    if (m == 0) {  /* subnormal x */
        for (i = 0; (ix&0x00800000) == 0; i++)
            ix<<=1;
        m -= i - 1;
    }
    m -= 127;  /* unbias exponent */
    ix = (ix&0x007fffff)|0x00800000;
    if (m&1)  /* odd m, double x to make it even */
        ix += ix;
    m >>= 1;  /* m = [m/2] */

    /* generate sqrt(x) bit by bit */
    ix += ix;
    q = s = 0;       /* q = sqrt(x) */
    r = 0x01000000;  /* r = moving bit from right to left */

    while (r != 0) {
        t = s + r;
        if (t <= ix) {
            s = t+r;
            ix -= t;
            q += r;
        }
        ix += ix;
        r >>= 1;
    }

    /* use floating add to find out rounding direction */
    if (ix != 0) {
        z = 1.0f - tiny; /* raise inexact flag */
        if (z >= 1.0f) {
            z = 1.0f + tiny;
            if (z > 1.0f)
                q += 2;
            else
                q += q & 1;
        }
    }
    ix = (q>>1) + 0x3f000000;
    ix += m << 23;
    SET_FLOAT_WORD(z, ix);
    return z;
}

typedef enum {
    ATOM_UNSPECIFIED,
    ATOM_H,
    ATOM_O,
    ATOM_N,
    ATOM_C,
    ATOM_MAX
} AtomKind;

typedef struct {
    f32 x, y;
} Vector2;

Vector2 vector2_add(Vector2 v1, Vector2 v2) {
    return (Vector2) { v1.x + v2.x, v1.y + v2.y };
}

Vector2 vector2_sub(Vector2 v1, Vector2 v2) {
    return (Vector2) { v1.x - v2.x, v1.y - v2.y };
}

Vector2 vector2_scale(Vector2 v, f32 scale) {
    return (Vector2) { v.x*scale, v.y*scale };
}

f32 vector2_dot(Vector2 v1, Vector2 v2) {
    return v1.x*v2.x + v1.y*v2.y;
}

f32 vector2_length(Vector2 v) {
    return f32_sqrt((v.x*v.x) + (v.y*v.y));
}

Vector2 vector2_noz(Vector2 v) {
    Vector2 result = { 0 };
    f32 length = f32_sqrt((v.x*v.x) + (v.y*v.y));
    if(length > 0.0001f) {
        result.x = v.x*1.0f/length;
        result.y = v.y*1.0f/length;
    }
    return result;
}

i32 check_circle_collision(Vector2 center1, f32 radius1, Vector2 center2, f32 radius2) {
    f32 dx = center2.x - center1.x;      // X distance between centers
    f32 dy = center2.y - center1.y;      // Y distance between centers

    f32 distance = f32_sqrt(dx*dx + dy*dy); // Distance between centers
    return distance <= radius1 + radius2;
}


typedef struct {
    f32 x, y, width, height;
} Rectangle;

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
char atom_names[] = { 'X', 'H', 'O', 'N', 'C' };
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
} Bond;

#define MAX_BONDS 512
Bond bonds[MAX_BONDS];

u8 bond_solutions[MAX_BONDS];
int bonds_count = 0;

f32 zoom = 31;

typedef struct {
    int count;
    int bond_ids[3];
    int atom_ids[3];
} Link;


Rectangle bounding_rect = { };
f32 width;
f32 height;

Rectangle rectangle_enlarge_to_contain(Rectangle rect, Vector2 pos) {
    f32 minx = f32_min(rect.x, pos.x);
    f32 miny = f32_min(rect.y, pos.y);
    f32 maxx = f32_max(rect.x + rect.width, pos.x);
    f32 maxy = f32_max(rect.y + rect.height, pos.y);
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
    bounding_rect = rectangle_enlarge_to_contain(bounding_rect, atoms[atoms_count].p);
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

Solution toggle_solution(Solution solution, int inc) {
    int count = solution & 0x3;
    int fixed = solution & 0x4;
    return ((count+inc) & 0x3) | fixed;
}

Solution toggle_fixed(Solution solution) {
    int count = solution & 0x3;
    int fixed = solution & 0x4;
    return count | ((~fixed) & 0x4);
}

int get_atom_at(Vector2 p) {
    Rectangle rect = bounding_rect;
    f32 scale = f32_max(1, f32_max(rect.width / (width - 10), rect.height / (height - 40)));
    p = vector2_scale(p, scale);
    for(int i = 0; i < atoms_count; i++) {
        if(check_circle_collision(p, 0, atoms[i].p, 10)) {
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

int get_bond_at(Vector2 p) {
    Rectangle rect = bounding_rect;
    f32 scale = f32_max(1, f32_max(rect.width / (width - 10), rect.height / (height - 40)));
    p = vector2_scale(p, scale);
    for(int i = 0; i < bonds_count; i++) {
        Vector2 a = atoms[bonds[i].atom_id1].p;
        Vector2 b = atoms[bonds[i].atom_id2].p;
        Vector2 pa = vector2_sub(p, a);
        Vector2 ba = vector2_sub(b, a);
        f32 h = f32_clamp(0, vector2_dot(pa,ba) / vector2_dot(ba,ba), 1);
        if(vector2_length(vector2_sub(pa, vector2_scale(ba, h))) < 7) {
            return i;
        }
    }
    return -1;
}

f32 get_scale(void) {
    Rectangle rect = bounding_rect;
    return f32_min(1, f32_min((width - 20) / rect.width, (height - 40) / rect.height));
}

Vector2 centered_position(Vector2 pos) {
    Rectangle rect = bounding_rect;
    f32 scale = get_scale();
    return (Vector2) {
        .x = width  / 2 + scale*( -rect.width  / 2 + pos.x),
        .y = height / 2 + scale*( -rect.height / 2 + pos.y),
    };  
}

void draw_edge(Vector2 p1, Vector2 p2, Solution s) {
    Vector2 dir = vector2_sub(p2, p1);
    Vector2 n = vector2_noz(dir);
    Vector2 o = { .x = -n.y, .y = n.x };

    f32 scale  = get_scale();
    f32 length = vector2_length(dir);
    int count = s & 0x3;
    int fix   = s & 0x4;
    if(count == 0 && fix == 0) {
        for(int i = 3; 5*i + 10 < length / scale; i++) {
            Vector2 p = vector2_add(p1, vector2_scale(n, 5*i*scale));
            draw_dot(p.x, p.y);
        }
    } else {
        for(int i = 0; i < count; i++) {
            Vector2 q1 = vector2_add(vector2_add(p1, vector2_scale(n,  13*scale)), vector2_scale(o, scale*(4*i - 2*(count - 1))));
            Vector2 q2 = vector2_add(vector2_add(p2, vector2_scale(n, -13*scale)), vector2_scale(o, scale*(4*i - 2*(count - 1))));
            draw_line(q1.x, q1.y, q2.x, q2.y, 0.8);
        }
    }
    if(count > 0 && fix) {
        Vector2 center = vector2_scale(vector2_add(p1, p2), 0.5);
        Vector2 q1 = vector2_add(center, vector2_scale(o,  scale*10));
        Vector2 q2 = vector2_add(center, vector2_scale(o, -scale*10));
        draw_line(q1.x, q1.y, q2.x, q2.y, 0.8);
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

char puzzle_text[4096];
i32  puzzle_text_length;
void init() {
    parse(puzzle_text, puzzle_text_length);
}

typedef struct {
    short bond_id;
    u8 diff;
} Rewind;

Rewind rewind_stack[4 * 1024] = { };
int rewind_stack_size = 0;
int rewind_stack_size_max = 0;

void undo() {
    if(rewind_stack_size > 0) {
        rewind_stack_size--;
        Rewind rw = rewind_stack[rewind_stack_size];
        bond_solutions[rw.bond_id] ^= rw.diff;
    }
}

void redo() {
    if(rewind_stack_size < rewind_stack_size_max) {
        Rewind rw = rewind_stack[rewind_stack_size];
        bond_solutions[rw.bond_id] ^= rw.diff;
        rewind_stack_size++;
    }
}

enum {
    VALUE_KIND_PUZZLE_TEXT = 0,
    VALUE_KIND_PUZZLE_TEXT_LENGTH = 1,
    VALUE_KIND_REWIND_STACK = 2,
    VALUE_KIND_REWIND_STACK_SIZE = 3,
    VALUE_KIND_REWIND_STACK_SIZE_MAX = 4,
    VALUE_KIND_BOND_SOLUTIONS = 5,
    VALUE_KIND_BONDS_COUNT = 6,
};

i32 get_value(i32 kind) {
    if(kind == VALUE_KIND_PUZZLE_TEXT)           return (i32)puzzle_text;
    if(kind == VALUE_KIND_PUZZLE_TEXT_LENGTH)    return (i32)puzzle_text_length;
    if(kind == VALUE_KIND_REWIND_STACK)          return (i32)rewind_stack;
    if(kind == VALUE_KIND_REWIND_STACK_SIZE)     return (i32)rewind_stack_size;
    if(kind == VALUE_KIND_REWIND_STACK_SIZE_MAX) return (i32)rewind_stack_size_max;
    if(kind == VALUE_KIND_BOND_SOLUTIONS)        return (i32)bond_solutions;
    if(kind == VALUE_KIND_BONDS_COUNT)           return (i32)bonds_count;
    assert(0);
    return -1;
}


void set_value(i32 kind, i32 value) {
    if(kind == VALUE_KIND_PUZZLE_TEXT_LENGTH)    puzzle_text_length    = value;
    else if(kind == VALUE_KIND_REWIND_STACK_SIZE)     rewind_stack_size     = value;
    else if(kind == VALUE_KIND_REWIND_STACK_SIZE_MAX) rewind_stack_size_max = value;
    else if(kind == VALUE_KIND_BONDS_COUNT)           bonds_count           = value;
    else assert(0);
}


void update(i32 w, i32 h, i32 x, i32 y, i32 button) {
    width = w;
    height = h;
    Vector2 offset = centered_position((Vector2){});
    Vector2 mouse = {x, y};
    
    int currentAtom = get_atom_at(vector2_sub(mouse, offset));
    int currentBond = currentAtom == -1 ? get_bond_at(vector2_sub(mouse, offset)) : -1;

    if(currentBond != -1) {
        Solution diff = 0; (void)diff;
        Solution current = bond_solutions[currentBond];
        if(button == MOUSE_LEFT) {
            Solution new = toggle_solution(current, 1);
            diff = (new ^ bond_solutions[currentBond]);
            bond_solutions[currentBond] = new;
        }
        if(button == MOUSE_RIGHT) {
            Solution new = toggle_solution(current, -1);
            diff = (new ^ bond_solutions[currentBond]);
            bond_solutions[currentBond] = new;
        }

        if(button == MOUSE_MIDDLE) {
            Solution new = toggle_fixed(current);
            diff = (new ^ bond_solutions[currentBond]);
            bond_solutions[currentBond] = new;
        }

        if(diff) {
            if(rewind_stack_size < ArrayLength(rewind_stack)) {
                rewind_stack[rewind_stack_size].bond_id = currentBond;
                rewind_stack[rewind_stack_size].diff    = diff;
                rewind_stack_size++;
                rewind_stack_size_max = rewind_stack_size;
            } else {
                for(int i = 1; i < ArrayLength(rewind_stack); i++) {
                    rewind_stack[i-1] = rewind_stack[i];
                }
                rewind_stack[rewind_stack_size - 1].bond_id = currentBond;
                rewind_stack[rewind_stack_size - 1].diff    = diff;
            }
        }
    }
    int isSolved = 0;
    if(isSolved) {
        clear(  0, 228,  48);
    } else {
        clear(180, 180, 180);
    }

    f32 scale = get_scale();
    for(int i = 0; i < atoms_count; i++) {
        Vector2 center = centered_position(atoms[i].p);
        int kind = atoms[i].kind;
        char name = atom_names[kind];
        draw_char(center.x, center.y, name, scale);
    }

    for(int i = 0; i < bonds_count; i++) {
        Vector2 p1 = centered_position(atoms[bonds[i].atom_id1].p);
        Vector2 p2 = centered_position(atoms[bonds[i].atom_id2].p);
        Solution s = bond_solutions[i];
        draw_edge(p1, p2, s);
    }
}