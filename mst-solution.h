typedef struct edge {
  int a;
  int b;
  int w;
} edge;

typedef struct p_edge {
  int v;
  int t;
  edge e;
} p_edge;


// comparison function b/w 2 edges, + lexicographic priority
int cmp_edges(const void *a, const void *b) {
  edge *u = (edge *)a;
  edge *v = (edge *)b;

  if (u->w == v->w)
    return (u->a < v->a || u->a == v->a && u->b < v->b) ? -1 : 1;
  else
    return u->w - v->w;
}
