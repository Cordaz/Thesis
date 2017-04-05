#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "../utilities/data_structures.h"
#include "../utilities/my_lib.h"
#include "../utilities/FIFO.h"
#include "../utilities/argparse.h"

#ifdef SILENT
	#define printf(...)
#endif

#define BUFFER 256
#define MIN_ARGS 2

static const char* const usage[] = {
	"dsc [options] -p path_pattern -s search_size",
	NULL
};

graph_t * load_graph(const char *, int, int, int);
int get_base_index(char);

int main(int argc, const char * argv[]) {
	//Setup
	const char * pattern = NULL;
	int s = 0;
	int k = 10;
	struct argparse_option options[] = {
		OPT_HELP(),
		OPT_GROUP("Mandatory"),
		OPT_STRING('p', "pattern", &pattern, "pattern to .graph.*, .psm and .approx.cnt files"),
		OPT_INTEGER('s', "search", &s, "size of motifs to retrieve"),
		OPT_GROUP("Optional"),
		OPT_INTEGER('k', "kmer", &k, "kmer length of graph (default 10)"),
		OPT_END()
	};

	struct argparse argparse;
	argparse_init(&argparse, options, usage, 0);
	argparse_describe(&argparse, "\nSearch for approximate occurences of k-mer of length s in the adjDBG.", "\nNote that the grah is intended to be structured on k/2");
	if(argc < MIN_ARGS) {
		argparse_usage(&argparse);
		return 0;
	}
	argc = argparse_parse(&argparse, argc, argv);

	if(!s) {
		fprintf(stdout, "[ERROR] search size must be specified\n");
		argparse_usage(&argparse);
		return 0;
	}

	if(!pattern) {
		fprintf(stdout, "[ERROR] pattern must be specified\n");
		argparse_usage(&argparse);
		return 0;
	}

	int pid = getpid();
	time_t rawtime;
	struct tm * timeinfo;
	//End setup

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Starting with PID %d\n", pid);

	int nodes = (int)pow((double)4, k/2);
	int edges = (int)pow((double)4, k/2+1);
	graph_t * dbg = load_graph(pattern, k/2, nodes, edges);
	if(!dbg) {
		return 1;
	}
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf(stdout, "[%02d:%02d:%02d][%5d] ", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec, pid);
	fprintf(stdout, "Loaded graph in memory\n");


	


	return 0;
}


graph_t * load_graph(const char * pattern, int k, int nodes, int edges) {
	int i;
	char buf[BUFFER+1];
	char * token;
	const char sep[2] = "\t";
	char args[5][BUFFER+1];
	//Allocate graph
	graph_t * dbg;
	if( !(dbg = (graph_t*)malloc(sizeof(graph_t))) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}
	if( !(dbg->nodes = (node_t**)malloc(sizeof(node_t*) * nodes)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	if( !(dbg->edges = (edge_t**)malloc(sizeof(edge_t*) * edges)) ) {
		fprintf(stdout, "[ERROR] couldn't allocate memory\n");
		return NULL;
	}

	FILE * fp;
	char file_path[BUFFER+1];
	node_t * n;
	int node_id;
	strcpy(file_path, pattern);
	strcat(file_path, ".graph.nodes");
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Load nodes
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<2 ) { //Expecting two args
			strcpy(args[i], token);

			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 2 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0'; //Remove '\n'
		node_id = atoi(args[0]);
		if( !(n = create_node(node_id, args[1], nodes)) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		dbg->nodes[node_id] = n;

		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	edge_t * e;
	int id;
	int n0_hash;
	int n1_hash;

	//Load one-step out edges
	strcpy(file_path, pattern);
	strcat(file_path, ".graph.edges.out");
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<3 ) { //Expecting three args
			strcpy(args[i], token);
			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 3 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0';
		args[2][k] = '\0';
		id = atoi(args[0]);
		n0_hash = hash(args[1], k);
		n1_hash = hash(args[2], k);
		if( !( e = create_edge( dbg->nodes[n0_hash], dbg->nodes[n1_hash], id ) ) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		dbg->edges[id] = e;
		dbg->nodes[n0_hash]->out[get_base_index(args[2][k-1])] = e;
		dbg->nodes[n1_hash]->in[get_base_index(args[1][0])] = e;

		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	//Load k/2-step out edges
	strcpy(file_path, pattern);
	sprintf(buf, ".graph.edges.out.%dstep", k);
	strcat(file_path, buf);
	if( !(fp = fopen(file_path, "r")) ) {
		fprintf(stdout, "[ERROR] couldn't open %s\n", file_path);
		return NULL;
	}
	//Read and ignore header
	fgets(buf, BUFFER, fp);
	//Read first line
	fgets(buf, BUFFER, fp);
	while(!feof(fp)) {
		token = strtok(buf, sep);
		i=0;
		while( token && i<5 ) { //Expecting five args
			strcpy(args[i], token);
			token = strtok(NULL, sep);
			i++;
		}
		if(token) {
			fprintf(stdout, "[ERROR] expecting 5 arguments, more found. Exit.\n");
			return NULL;
		}
		args[1][k] = '\0';
		args[2][k] = '\0';
		id = atoi(args[0]);
		n0_hash = hash(args[1], k);
		n1_hash = hash(args[2], k);
		if( !( e = create_edge( dbg->nodes[n0_hash], dbg->nodes[n1_hash], id ) ) ) {
			fprintf(stdout, "[ERROR] couldn't allocate memory\n");
			return NULL;
		}
		e->count = atoi(args[3]);
		e->input_count = atoi(args[4]);
		dbg->nodes[n0_hash]->out_kstep[n1_hash] = e;
		dbg->nodes[n1_hash]->in_kstep[n0_hash] = e;
		fgets(buf, BUFFER, fp);
	}
	fclose(fp);

	return dbg;
}


int get_base_index(char b) {
	switch(b) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			return -1;
	}
}
