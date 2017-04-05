#ifndef DATA_STRUCTURES_H
	#define DATA_STRUCTURES_H

	//// PROTOTYPES
	struct node_s;

	typedef struct node_s {
		int id;
		char * seq;
		struct edge_s * in[4];
		struct edge_s * out[4];
		struct edge_s ** in_kstep;
		struct edge_s ** out_kstep;
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
	node_t * create_node(int, char *, int);


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
