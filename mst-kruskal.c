// ====================================================================
// disjoint set data structure, with path compression and union by rank
// ====================================================================

typedef struct disjoint_set {
  int *P; // parents
  int *R; // ranks
} disjoint_set;


int find(int x, int *P) {
  int y = x;
  int t = P[x];

  // find the root of the tree
  while (y != P[y])
    y = P[y];

  // follow the chain of ancestors of x
  // set the root as parent of every one of those
  while (t != y) {
    P[x] = y;
    x = P[t];
    t = P[x];
  }

  return y;
}


// I rEaLly want to call my function `union`
// we are assuming x and y are already roots
void union_(int x, int y, disjoint_set *F) {
  // same set
  if (x == y)
    return;

  int rx = F->R[x];
  int ry = F->R[y];

  if (rx < ry)
    F->P[x] = y;
  else if (ry < rx)
    F->P[y] = x;
  else {
    F->P[y] = x;
    F->R[x]++;
  }
}


// ====================================================================
// kruskal implementation, used for both sequential & parallel versions
// ====================================================================

// holds maximum spanning forest information
// that needs to be sent between processors
typedef struct msf {
  int size;
  edge *edges;
} msf;


void kruskal(edge *edges, int nb_edge, int count, msf *forest) {
  disjoint_set sets;
  int roota, rootb;
  edge *e;

  forest->edges = malloc(count * sizeof(edge));
  int nb = 0;

  sets.P = malloc(count * sizeof(int));
  sets.R = calloc(count, sizeof(int));

  // initialize data structure with singletons
  for (int i = count; i--;)
    sets.P[i] = i;

  for (int k = 0; k < nb_edge && nb < count; k++) {
    // find an edge joining 2 disjoint sets
    e = &edges[k];
    roota = find(e->a, sets.P);
    rootb = find(e->b, sets.P);

    if (roota == rootb)
      continue;

    // add it to the forest
    forest->edges[nb++] = *e;

    // and merge the two sets
    union_(roota, rootb, &sets);
  }

  forest->size = nb;

  free(sets.P);
  free(sets.R);
}


void merge(msf *a, msf *b, int count) {
  disjoint_set sets;
  sets.P = malloc(count * sizeof(int));
  sets.R = calloc(count, sizeof(int));

  for (int i = count; i--;)
    sets.P[i] = i;

  int size = a->size + b->size;
  int nb = 0;
  edge *edges = malloc(size * sizeof(edge));
  edge *e;

  int roota, rootb;
  int ia = 0;
  int ib = 0;

  while (size--) {
    if (ia >= a->size)
      e = &b->edges[ib++];
    else if (ib >= b->size)
      e = &a->edges[ia++];
    else {
      int c = cmp_edges(&a->edges[ia], &b->edges[ib]);
      // we cannot have c == 0 since every processor consider different edges
      if (c < 0)
        e = &a->edges[ia++];
      else
        e = &b->edges[ib++];
    }

    roota = find(e->a, sets.P);
    rootb = find(e->b, sets.P);

    if (roota == rootb)
      continue;

    // add it to the forest
    edges[nb++] = *e;

    // and merge the two sets
    union_(roota, rootb, &sets);
  }

  free(sets.P);
  free(sets.R);
  free(a->edges);

  a->edges = edges;
  a->size = nb;
}


