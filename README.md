# Voronoi Diagrams Using Quad-Edge Data Structure #

This document presents the paper "Primitives for the manipulation of general subdivisions and the computation of Voronoi." ACM transactions on graphics (TOG) 4.2 (1985): 74-123.  I found that the paper is about 50-page long and there is only a scanned version available online.  I hope this document makes the process of implementing your own quad-edge and Voronoi Diagram algorithm easier.  The implementation is in rust but the pseudo-code in this document will be written in C++ for your convenience.

This is also the course project of CS274 Computational Geometry in UC Berkely.


## Quad-Edge Data Structure ##
First let's review the definition of an edge that may occur in a old traditional edge DCEL:
```c++
struct edge_t {
   vertex_t *org, *dst;
   face_t *left, *right;

   edge_t *prev, *next; // Onext in the paper
   edge_t *left_prev, *left_next; // Lnext in the paper
   edge_t *flip, sym;
   edge_t *dual;
   edge_t *rot;
} *e;
```
* `org` and `dst` represent the two end-points of this edge;
* `left` and `right` represents the two adjacent faces of edge in the direction `org` to `dst`;
* `sym` represents the reverse of the edge, i.e., swap of `org`/`dst` and `left`/`right`;
* `flip` is a bit magic here.  Imagine that the edge is a bug standing on a piece of paper.  Then `e->flip` represents the bug hanging upside down from the other side of the paper with the same `org` and `dst`, i.e., the "normal" of edge `e` is flipped;
* `prev` and `next` represent the immediate edge following `e` in counterclockwise order with respect to `e->org`;
* `left_prev` and `left_next` represent the next counterclockwise edge with the same left face.

For duality of an edge, we have the following definition (property):
```c++
assert(e->dual->dual == e);
assert(e->sym->dual == e->dual->sym);
assert(e->flip->dual == e->dual->flip->sym);
assert(e->left_next->dual == e->dual->org_prev);
```

Stare at the following duality example and you will understand why:

![subdivision dual](https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Duals_graphs.svg/313px-Duals_graphs.svg.png)

For convenience, we define the dual of flipped edge:
```c++
assert(e->rot == d->flip->dual);
```

`rot` is a very powerful operator, e.g., `e->rot->rot->rot->rot == e`.  Geometrically, `rot` can be understand as rotating the edge 90 degree counterclockwise (and flip its normal direction).  With the help of `e->rot`, we only need three edge reference of `edge_t` to traverse the edge list:
```c++
struct quad_edge_intro_t {
   vertex_t *org, *dst;
   face_t *left, *right;

   quad_edge_intro_t *rot, *flip, *next;
} *e;
```
For example, if we want to visit `e->sym` and `e->left_next` (the notation is abused here), we have
```c++
assert(e->sym == d->rot->rot);
assert(e->left_next == d->rot->rot->rot->next->rot);
```

We can make several improvements:
1. since we often refer `rot` many time, we may want to pre-process it and store it into an array;
2. we do not need to store `sym` and `dual` pointers;
3. we do not need to store pointer to vertexs as vertexs can be referred by any edge that starts with that vertex.  The same can be done for faces using its vertex dual.

If the manifold (plane) is orientable, i.e., the manifold is not a Möbius strip, which normally is the case, then the normal direction of the edge is not important because `e` and `e->flip` cannot co-exists in the same subdivision.  We do not even need to store the `flip` pointer.  As a result, you also don't need to store the flipped edge themselves.

The final elegant data structure with these optimization look like this
```c++
struct quad_edge_t {
   data_t *o;  // additional data if necessary
   quad_edge_t *rot[3], *next;

   quad_edge_t *sym() { return rot[1]; }                   // reverse edge
   quad_edge_t *Onext() { return next; }                   // next edge w.r.t. org
   quad_edge_t *Oprev() { return rot[0]->next->rot[0]; }   // prev edge w.r.t. org
   quad_edge_t *Dnext() { return rot[1]->next->rot[1]; }   // next edge w.r.t. dst
   quad_edge_t *Dprev() { return rot[2]->next->rot[2]; }   // prev edge w.r.t. dst
   quad_edge_t *Lnext() { return rot[2]->next->rot[0]; }   // next edge w.r.t. left face
   quad_edge_t *Rnext() { return rot[0]->next->rot[2]; }   // next edge w.r.t. right face
   quad_edge_t *Lprev() { return rot[1]->next; }           // prev edge w.r.t. left face
   quad_edge_t *Rprev() { return next->rot[1]; }           // prev edge w.r.t. right face
} *e;
```

## Topological Operations of Quad-Edge ##

### Create a New Edge ###
The first operator is to create a new non-loop edge on a empty sphere with its `sym` and `rot`.

![create edge](https://raw.githubusercontent.com/zhou13/voronoi/master/image/make_edge.png)

```c++
quad_edge_t *make_edge()
{
    quad_edge_t *e = new quad_edge_t[4];
    e[0]->rot[2] = e[1]->rot[1] = e[2]->rot[0] = e[3];
    e[1]->rot[2] = e[2]->rot[1] = e[0]->rot[0] = e[0];
    e[2]->rot[2] = e[3]->rot[1] = e[0]->rot[0] = e[1];
    e[3]->rot[2] = e[0]->rot[1] = e[1]->rot[0] = e[2];
    e[0].next = e[0];
    e[1].next = e[3]
    e[2].next = e[2];
    e[3].next = e[1]
    return e;
}
```

### Splice Edge ###

The second operator is to split/merge "vertex":

![splice](https://raw.githubusercontent.com/zhou13/voronoi/master/image/splice.png)

```c++
// a and b must be both primal or both dual
void splice(quad_edge_t *a, quad_edge_t *b)
{
    std::swap(a->next, b->next);
    auto alph = a->next->rot[0];
    auto beta = b->next->rot[0];
    std::swap(alph->next, beta->next);
}
```
The above figure explains this function's behavior:
* if edge a and b share the same org, the vertex will be split into two and two faces will be merged;
* if edge a and b do not share the same org, their org vertices will be merged and two faces will be split.

## Topological Operations For Voronoi Diagram ##

### Connect ###

The connect function connects destination of a to the origin of a using two splice operation.

![connect](https://cdn.rawgit.com/zhou13/voronoi/master/image/connect.svg)

```c++
// a and b must share the same left faces
void connect(quad_edge_t *a, quad_edge_t *b)
{
    auto e = make_edge();
    splice(e, a->Lnext());
    splice(e->sym(), b);
    update_meta(e);
}
```
Here `update_meta(e)` will fill the necessary information (such as the coordinate of the original of `e`).

### Delete ###

Edge deletion is just an inverse operation of `connect`, which removes an edge from a subdivision.

```c++
void delete(quad_edge_t *e)
{
    splice(e, e->Oprev());
    splice(e->sym, e->sym()->Oprev());
}
```

### Flip ###

Flip is the edge-flip operation for generating delaunay diagram.
![create edge](https://raw.githubusercontent.com/zhou13/voronoi/master/image/flip.png)

```c++
void flip(quad_edge_t *e)
{
    a = e->Oprev();
    b = e->sym()->Oprev();
    splice(e, a);
    splice(e->sym(), b);
    splice(e, a.Lnext());
    splice(e->sym(), b.Lnext());
    update_meta(e);
}
```

##  Geometric Predicate ##

![predicate](https://cdn.rawgit.com/zhou13/voronoi/master/image/preds.gif)

##  Delaunay Triangulation Property ##

1.  Adding new points to the plane will never introduce new edges between old points

##  Delaunay Triangulation Divide-and-Conquer ##

For split the point set S into L and R according to the lexicographical order.
