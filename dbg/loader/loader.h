#ifndef DBG_READER_H
	#define DBG_READER_H

	#include <stdio.h>
	#include <string.h>
	#include <malloc.h>
	#include <stdlib.h>
	#include <math.h>
	#include "../data_structures.h"

	/*
	 * Alloc the external structure of the graph.
	 * The two doubles are the dimension of nodes and edges
	 * It also init all content to NULL
	 *
	 */
	graph_t * alloc_graph(double, double);

	/*
	 * Loads nodes from the *.nodes file prompted
	 * MUST be invoked after alloc_graph().
	 *
	 */
	graph_t * load_nodes(graph_t *, char *);

	/*
	 * Load edges from the *.edges file prompted
	 * It also loads edges into nodes.
	 * MUST be invoked after load_nodes().
	 *
	 */
	graph_t * load_edges(graph_t *, char *);

#endif
