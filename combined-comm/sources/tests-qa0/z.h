/* This file contains what qa0->h should do to z.qa0 */
/* comments here are for clarification only */
struct UpLink {
    int u; /* built-in type int converts to int */
    int f;
};

struct DownLink {
    int f;
};

typedef int Links[4]; /* array is converted to typedef, *dim* is computed */

struct Neighbor {
    int mask;
    int up_gauge;
    Links up_fermion; /* links has Links as C-name */
    Links down_gauge;
    Links down_fermion;
};

/* macro over d and p/m is expanded, notice procedure name construction
from the stem */
int Up_face0plus(void);
int Down_face0_plus(int size, struct F* r, struct F* s);
int Up_face0minus(void);
int Down_face0_minus(int size, struct F* r, struct F* s);
int Up_face1plus(void);
int Down_face1_plus(int size, struct F* r, struct F* s);
int Up_face1minus(void);
int Down_face1_minus(int size, struct F* r, struct F* s);
int Up_face2plus(void);
int Down_face2_plus(int size, struct F* r, struct F* s);
int Up_face2minus(void);
int Down_face2_minus(int size, struct F* r, struct F* s);
int Up_face3plus(void);
int Down_face3_plus(int size, struct F* r, struct F* s);
int Up_face3minus(void);
int Down_face3_minus(int size, struct F* r, struct F* s);

/* notice that zzz is of type struct bar *, not Neighbor! */
void funny_name(struct bar * zzz);

