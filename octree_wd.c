/* Taken from https://gist.github.com/astarasikov/c037846009a5a6ec4b93 */
/* Added field for data by https://github.com/alexdany657 */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Octree;
typedef struct Octree Octree;

enum {
    OCTREE_NODE_COUNT = 8,
    OCTREE_COORD_COUNT = 3,

    OCTREE_POINT_LIMIT = 10,

    COORD_X = 0,
    COORD_Y = 1,
    COORD_Z = 2,
};

typedef enum Retcode {
    RETCODE_OK = 1,
    RETCODE_FAIL = 2,
} Retcode;

#define CHECK_OR_FAIL(cond, code) do {\
    if (!(cond)) {\
        fprintf(stderr, "%s: failed to check " #cond "\n", __func__);\
        return code ;\
    }\
} while (0)

#if 0
#define OCTREE_DEBUGF(fmt, args...) do { fprintf(stderr, fmt, ##args); } while(0)
#else
#define OCTREE_DEBUGF(fmt, args...) do {} while (0)
#endif

struct Octree {
    Octree *childNodes[OCTREE_NODE_COUNT];
    Octree *parentNode;
    bool isLeaf;
    double center[OCTREE_COORD_COUNT];
    double size[OCTREE_COORD_COUNT];

    double points[OCTREE_POINT_LIMIT * OCTREE_COORD_COUNT];
    double vecData[OCTREE_POINT_LIMIT * OCTREE_COORD_COUNT];
    size_t pointCount;
};

#define CHECK_OCTREE_COORD(coord) do { \
    double dist = fabs(position[COORD_ ##coord] - tree->center[COORD_ ##coord]); \
    OCTREE_DEBUGF("%s: dist=%f radius=%f\n", __func__, dist, tree->size[COORD_ ##coord]); \
    OCTREE_DEBUGF("\t position=[%f %f %f] center=[%f %f %f]\n",\
            position[0], position[1], position[2],\
            tree->center[0], tree->center[1], tree->center[2]);\
    CHECK_OR_FAIL(dist <= tree->size[COORD_ ##coord], RETCODE_FAIL); \
} while (0)

#define OCTREE_CMP_COORD(coord) (position[COORD_ ##coord] > tree->center[COORD_ ##coord])

static Octree *makeTree(double *at, double *size)
{
    Octree *tree = (Octree*)malloc(sizeof(Octree));
    CHECK_OR_FAIL(tree != NULL, NULL);
    memset(tree, 0, sizeof(Octree));
    memcpy(tree->center, at, OCTREE_COORD_COUNT * sizeof(double));
    memcpy(tree->size, size, OCTREE_COORD_COUNT * sizeof(double));
    tree->isLeaf = true;
    return tree;
}

static inline Retcode calculateChildIndex(Octree *tree, double *position, unsigned *out)
{
    CHECK_OR_FAIL(tree != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(position != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(out != NULL, RETCODE_FAIL);
    
    CHECK_OCTREE_COORD(X);
    CHECK_OCTREE_COORD(Y);
    CHECK_OCTREE_COORD(Z);
    
    bool cx = OCTREE_CMP_COORD(X);
    bool cy = OCTREE_CMP_COORD(Y);
    bool cz = OCTREE_CMP_COORD(Z);
    
    unsigned mask = (cx << COORD_X) | (cy << COORD_Y) | (cz << COORD_Z);
    *out = mask;
    return RETCODE_OK;
}

static Retcode lookup(Octree* tree, double *position, Octree **out)
{
    CHECK_OR_FAIL(tree != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(position != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(out != NULL, RETCODE_FAIL);
    *out = NULL;

    unsigned mask = 0;
    CHECK_OR_FAIL(RETCODE_FAIL != calculateChildIndex(tree, position, &mask), RETCODE_FAIL);

    if (tree->isLeaf) {
        *out = tree;
        return RETCODE_OK;
    }

    return lookup(tree->childNodes[mask], position, out);
};

static inline Retcode subdivide(Octree *tree)
{
    const double half = 0.5;
    const double nhalf = -0.5;

    for (size_t i = 0; i < OCTREE_NODE_COUNT; i++) {
        double sign_x = (i & (1 << COORD_X)) ? (half) : (nhalf);
        double sign_y = (i & (1 << COORD_Y)) ? (half) : (nhalf);
        double sign_z = (i & (1 << COORD_Z)) ? (half) : (nhalf);

        double x = tree->center[COORD_X] + (sign_x * tree->size[COORD_X]);
        double y = tree->center[COORD_Y] + (sign_y * tree->size[COORD_Y]);
        double z = tree->center[COORD_Z] + (sign_z * tree->size[COORD_Z]);

        double origin[OCTREE_COORD_COUNT] = {x, y, z};
        double size[OCTREE_COORD_COUNT] = { 
            half * tree->size[COORD_X],
            half * tree->size[COORD_Y],
            half * tree->size[COORD_Z]
        };
        Octree *t = makeTree(origin, size);
        CHECK_OR_FAIL(t != NULL, RETCODE_FAIL);
        tree->childNodes[i] = t;
        t->parentNode = tree;
    }
    return RETCODE_OK;
}

static Retcode insert(Octree *tree, double *position, double *vec)
{
    CHECK_OR_FAIL(tree != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(position != NULL, RETCODE_FAIL);
    
    unsigned mask = 0;
    CHECK_OR_FAIL(RETCODE_FAIL != calculateChildIndex(tree, position, &mask), RETCODE_FAIL);

    if (tree->isLeaf && (tree->pointCount + 1 < OCTREE_POINT_LIMIT)) {
        //insert locally
        memcpy(tree->points + OCTREE_COORD_COUNT * tree->pointCount,
            position, OCTREE_COORD_COUNT * sizeof(double));
        memcpy(tree->vecData + OCTREE_COORD_COUNT * tree->pointCount,
            vec, OCTREE_COORD_COUNT * sizeof(double));
        ++tree->pointCount;
        return RETCODE_OK;
    }

    if (tree->isLeaf) {
        //initialize child nodes
        CHECK_OR_FAIL(RETCODE_FAIL != subdivide(tree), RETCODE_FAIL);

        //mark self as non-leaf and recurse to insert child nodes
        tree->isLeaf = false;

        //insert own points to child nodes
        for (size_t i = 0; i < tree->pointCount; i++) {
            Retcode rc = insert(tree,
                    tree->points + OCTREE_COORD_COUNT * i, tree->vecData + OCTREE_COORD_COUNT * i);
            CHECK_OR_FAIL(RETCODE_FAIL != rc, RETCODE_FAIL);
        }

        memset(tree->points, 0,
                OCTREE_POINT_LIMIT * OCTREE_COORD_COUNT * sizeof(double));
        memset(tree->vecData, 0,
                OCTREE_POINT_LIMIT * OCTREE_COORD_COUNT * sizeof(double));
        tree->pointCount = 0;
    }

    Octree *child = tree->childNodes[mask];
    CHECK_OR_FAIL(tree != NULL, RETCODE_FAIL);
    CHECK_OR_FAIL(RETCODE_FAIL != insert(child, position, vec), RETCODE_FAIL);

    return RETCODE_OK;
}

static Retcode remove_octree(Octree *tree) {
    for (size_t i = 0; i < OCTREE_NODE_COUNT; ++i) {
        if (tree->childNodes[i]) {
            remove_octree(tree->childNodes[i]);
        }
    }
    free(tree);
    return RETCODE_OK;
}

static Retcode dfs(Octree *tree) {
    for (size_t i = 0; i < OCTREE_NODE_COUNT; ++i) {
        if (tree->childNodes[i]) {
            dfs(tree->childNodes[i]);
        }
    }
    printf("node: %p with %i children:\n", tree, tree->pointCount);
    for (size_t i = 0; i < OCTREE_NODE_COUNT; ++i) {
        printf("\t%p\n", tree->childNodes[i]);
    }
    return RETCODE_OK;
}

/*
int main() {
    const double tsize = 400;
    double origin[OCTREE_COORD_COUNT] = {
        tsize / 2, tsize / 2, tsize / 2
    };

    const size_t points_to_insert = 40;
    const size_t seed = 0xdead;
    srand(seed);

    Octree *root = makeTree(origin, origin);

    for (size_t c = 0; c < points_to_insert; c++) {
        double x = (rand() % (int)tsize);
        double y = (rand() % (int)tsize);
        double z = (rand() % (int)tsize);
        double position[OCTREE_COORD_COUNT] = { x, y, z };

        CHECK_OR_FAIL(insert(root, position, position) != RETCODE_FAIL, RETCODE_FAIL);

    }
    fprintf(stderr, "inserted points\n");
    
    for (size_t c = 0; c < points_to_insert; c++) {
        double x = (rand() % (int)tsize);
        double y = (rand() % (int)tsize);
        double z = (rand() % (int)tsize);
        double position[OCTREE_COORD_COUNT] = { x, y, z };

        Octree *t = NULL;
        Retcode rc = lookup(root, position, &t);
        if ((rc == RETCODE_OK) && (t != NULL)) {
            printf("found tree for [%f %f %f]\n",
                    position[0], position[1], position[2]);
            for (size_t i = 0; i < t->pointCount; i++) {
                printf("\t neighbour [%f %f %f]\n",
                        t->points[OCTREE_COORD_COUNT * i],
                        t->points[OCTREE_COORD_COUNT * i + 1],
                        t->points[OCTREE_COORD_COUNT * i + 3]);
            }
        }
        else {
            printf("not found tree for [%f %f %f]\n",
                    position[0], position[1], position[2]);
        }
    }

    remove_octree(root);

    return 0;
}
*/
