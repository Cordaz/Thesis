#ifndef DATA_STRUCTURES_H
	#define DATA_STRUCTURES_H

	//// PROTOTYPES
	struct node_s;
	struct edge_s;
	struct list_edge_s;


	typedef struct list_edge_s {
		struct edge_s * e;
		struct list_edge_s * next;
	} list_edge_t;

	typedef struct tree_edge_s {
		struct edge_s * e;
		struct tree_edge_s * p;
		struct tree_edge_s * left;
		struct tree_edge_s * right;
	} tree_edge_t;

	typedef struct node_s {
		int id;
		char * seq;
		list_edge_t * in;
		list_edge_t * out;
		tree_edge_t * in_kstep;
		tree_edge_t * out_kstep;
	} node_t;

	typedef struct edge_s {
		int id;
		int count;
		int input_count;
		node_t * from;
		node_t * to;
	} edge_t;

	typedef struct graph_s {
		node_t ** nodes; //hash table of nodes
		edge_t ** edges; //hash table of edges
	} graph_t;

	//// FUNCTIONS PROTOTYPES
	/*
	 * Creates (and allocates) a new node containing:
	 * - a sequence
	 * - corresponding id
	 * - NULL pointer for in/out edges
	 *
	 * Return NULL if can't create node.
	 *
	 */
	node_t * create_node(int, char *);

	/*
	 * Add edge to in-list.
	 *
	 */
	node_t * add_in_edges(node_t *, edge_t *);

	/*
	 * Add edge to out-list.
	 *
	 */
	 node_t * add_out_edges(node_t *, edge_t *);

	 /*
	  * Add edge to in-list.
	  *
	  */
	 node_t * add_in_kstep_edges(node_t *, edge_t *);

	 /*
	  * Add edge to out-list.
	  *
	  */
	  node_t * add_out_kstep_edges(node_t *, edge_t *);

	 /*
	  * If there is an edge out of the first node into the second (k-step)
	  * returns the edge, NULL otherwise.
	  *
	  */
	  edge_t * exist_edge(node_t *, node_t *);

	 /*
	  * Creates (and allocates) a new edge between two nodes
	  * - corresponding id
	  * - count to 1
	  *
	  * Return NULL if can't create edge.
	  *
	  */
	  edge_t * create_edge(node_t *, node_t *, int);

#endif
