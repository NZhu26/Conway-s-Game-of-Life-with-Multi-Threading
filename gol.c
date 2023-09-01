/*
* Neil Zhu
*
* This file implements the Game of Life utilizing parallelization
*/

/*
 * To run:
 * ./gol file1.txt  0  # run with config file file1.txt, do not print board
 * ./gol file1.txt  1  # run with config file file1.txt, ascii animation
 * ./gol file1.txt  2  # run with config file file1.txt, ParaVis animation
 *
 */
#include <pthreadGridVisi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include "colors.h"

/****************** Definitions **********************/
/* Three possible modes in which the GOL simulation can run */
#define OUTPUT_NONE   (0)   // with no animation
#define OUTPUT_ASCII  (1)   // with ascii animation
#define OUTPUT_VISI   (2)   // with ParaVis animation

/* Used to slow down animation run modes: usleep(SLEEP_USECS);
 * Change this value to make the animation run faster or slower
 */
#define SLEEP_USECS    (200000)

/* A global variable to keep track of the number of live cells in the
 * world (this is the ONLY global variable you may use in your program)
 */
static int total_live = 0;

/* Mutex*/
static pthread_mutex_t my_mutex = PTHREAD_MUTEX_INITIALIZER;

/* Barrier */
static pthread_barrier_t my_barrier;

/* This struct represents all the data you need to keep track of your GOL
 * simulation.  Rather than passing individual arguments into each function,
 * we'll pass in everything in just one of these structs.
 * this is passed to play_gol, the main gol playing loop
 *
 * NOTE: You will need to use the provided fields here, but you'll also
 *       need to add additional fields. (note the nice field comments!)
 * NOTE: DO NOT CHANGE THE NAME OF THIS STRUCT!!!!
 */
struct gol_data {

    // NOTE: DO NOT CHANGE the names of these 4 fields (but USE them)
    int rows;  // the row dimension
    int cols;  // the column dimension
    int iters; // number of iterations to run the gol simulation
    int output_mode; // set to:  OUTPUT_NONE, OUTPUT_ASCII, or OUTPUT_VISI
    int num_threads; //number of threads

    //set to(0: row-wise grid cell allocation, 1: column-wise grid cell
    //allocation)
    int threads_partition;

    //set to (0: don't print config info, 1: print allocation info) 
    int print_partition; 

    int initalive; //initial amount of cells alive
    int *graph; //graph 1
    int *graph2; //graph 2
    int *currentGraph; //integer pointer to which graph currently available, 
                       //0 = graph, 1 = graph2
    int *work_per_thread; //amount of rows/cols per thread
    pthread_barrier_t *done;  // pointer to single shared barrier

    /* fields used by ParaVis library (when run in OUTPUT_VISI mode). */
    // NOTE: DO NOT CHANGE their definitions BUT USE these fields
    visi_handle handle;
    color3 *image_buff;
};

/* This struct contains the information each thread need to do work*/
struct thread_info{
    int tid;
    struct gol_data *data;
};

/****************** Function Prototypes **********************/

/* the main gol game playing loop (prototype must match this) */
void *play_gol(void *args);

/* init gol data from the input file and run mode cmdline args */
int init_game_data_from_args(struct gol_data *data, char **argv);

// A mostly implemented function, but a bit more for you to add.
/* print board to the terminal (for OUTPUT_ASCII mode) */
void print_board(struct gol_data *data, int round);

/*the play one round function plays one round of the game*/
void play_round(struct gol_data *data);

/* gets index of a 2d graph point for a 2d array*/
int coord_index(struct gol_data *data, int i, int j);

/*checks if coordinate is alive
* returns 1 if alive, returns 0 if dead
*/
int check_alive_dead(struct gol_data *data, int i, int j);

/*checks eight neighbors of a cell to see how many neighbors are alive
* returns number of neighbors of a cell that are alive
*/
int check_neighbors(struct gol_data *data, int i, int j);

/*decides whether a cell should live or die
* returns 1 if cell should live, returns 0 if cell should die
*/
int celllives_celldies(struct gol_data *data, int i, int j);

/* use updated data to set colors for visualization */
void update_colors(struct gol_data *data, int startRow, int endRow, int startCol, int endCol, int myid);

/* partitions the graphs of data for parallelization*/
int *partition_graph(struct gol_data *data);

/* prints partitioning information*/
void print_partitioning(struct gol_data *data);

// thread main function
//   arg: an int value specifying a logical tid value
void *thread_function(void *arg);


/************ Definitions for using ParVisi library ***********/
/* initialization for the ParaVisi library (DO NOT MODIFY) */
int setup_animation(struct gol_data* data);
/* register animation with ParaVisi library (DO NOT MODIFY) */
int connect_animation(void (*applfunc)(struct gol_data *data),
        struct gol_data* data);
/* name for visi (you may change the string value if you'd like) */
static char visi_name[] = "GOL!";





int main(int argc, char **argv) {
    int ret;
    struct gol_data data;
    struct timeval start_time;
    struct timeval stop_time;
    double secs;
    pthread_t *tids;
    int i, ntids;
    struct thread_info *tid_args;

    /* check number of command line arguments */
    if (argc < 6) {
        printf("usage: %s infile.txt output_mode[0,1,2] number of threads num_threads partition[0,1] print_partition[0,1]\n", argv[0]);
        printf("output mode - (0: no visualization, 1: ASCII, 2: ParaVisi)\n");
        printf("num_threads partition - (0: row-wise grid cell allocation, 1: column-wise grid cell allocation)\n");
        printf("print_partition - (0: do not print configuration information, 1: print allocation information)\n");
        exit(1);
    }

    /* Initialize game state (all fields in data) from information
     * read from input file */
    ret = init_game_data_from_args(&data, argv);
    if (ret != 0) {
        printf("Initialization error: file %s, mode %s\n", argv[1], argv[2]);
        exit(1);
    }

    /* initialize ParaVisi animation (if applicable) */
    if (data.output_mode == OUTPUT_VISI) {
        setup_animation(&data);
    }

    /* ASCII output: clear screen & print the initial board */
    if (data.output_mode == OUTPUT_ASCII) {
        if (system("clear")) { perror("clear"); exit(1); }
        ret = gettimeofday(&start_time, NULL);
        print_board(&data, 0);
        usleep(SLEEP_USECS);
    }

    /* Invoke play_gol in different ways based on the run mode */
    if (data.output_mode == OUTPUT_NONE) {  // run with no animation
        ret = gettimeofday(&start_time, NULL);
        play_gol(&data);
        ret = gettimeofday(&stop_time, NULL);
    }
    else if (data.output_mode == OUTPUT_ASCII) { // run with ascii animation
        play_gol(&data);

        // clear the previous print_board output from the terminal:
        // (NOTE: you can comment out this line while debugging)
        if (system("clear")) { perror("clear"); exit(1); }

        // NOTE: DO NOT modify this call to print_board at the end
        //       (it's to help us with grading your output)
        print_board(&data, data.iters);
        ret = gettimeofday(&stop_time, NULL);
    }
    else {  // OUTPUT_VISI: run with ParaVisi animation
        ntids = data.num_threads;

        tids = malloc(sizeof(pthread_t) * ntids);
        if (!tids) {
            perror("malloc: pthread_t array");
            exit(1);
        }

        tid_args = malloc(sizeof(struct thread_info) * ntids);
        if (!tid_args) {
            perror("malloc: tid_args array");
            exit(1);
        }

        for (i = 0; i < ntids; i++) {
            tid_args[i].tid = i;
            tid_args[i].data = &data;

            ret = pthread_create(&tids[i], 0, thread_function, &tid_args[i]);

            if (ret) {
                perror("Error pthread_create\n");
                exit(1);
            }
        }

        //start ParaVisi animation
        run_animation(data.handle, data.iters);

        // wait for all worker threads to exit
        for (i = 0; i < ntids; i++) {
            pthread_join(tids[i], 0);
        }

        /* Clean up memory. */
        printf("Main thread done\n");
        free(tid_args);
        free(data.done);
        tid_args = NULL;
        free(tids);
        tids = NULL;
    }

    // NOTE: you need to determine how and where to add timing code
    //       in your program to measure the total time to play the given
    //       number of rounds played.
    if (data.output_mode != OUTPUT_VISI) {
        secs = ((stop_time.tv_sec - start_time.tv_sec)*1000000.0 + stop_time.tv_usec - start_time.tv_usec)/1000000.0;

        /* Print the total runtime, in seconds. */
        // NOTE: do not modify these calls to fprintf
        fprintf(stdout, "Total time: %0.3f seconds\n", secs);
        fprintf(stdout, "Number of live cells after %d rounds: %d\n\n",
                data.iters, total_live);
    }

    /* Cleaning up memory. */
    free(data.work_per_thread);
    free(data.graph);
    free(data.graph2);
    free(data.done);
    
    /* Release the synchronization variables. */
    if (pthread_mutex_destroy(&my_mutex)) {
        printf("pthread_mutex_destroy error\n");
        exit(1);
    }

    if (pthread_barrier_destroy(&my_barrier)) {
        printf("pthread_barrier_destroy error\n");
        exit(1);
    }

    return 0;
}

/* Function passed into each thread to perform game of life calculations
*  
*  arg: thread_info struct
*
*  return: NULL
*/
void *thread_function(void *arg) {
    //gets logical thread id out of thread_info struct passed in as arg
    int myid = ((struct thread_info *)arg)->tid;

    //gets gol_data struct out of thread_info struct passed in as arg
    struct gol_data *data = ((struct thread_info *)arg)->data; 

    int startRow; //starting row of part of graph that thread works on
    int startCol; //starting col of part of graph that thread works on
    int endRow; //ending row of part of graph that thread works on
    int endCol; //ending col of part of graph that thread works on

    //Prints thread's partitioning information
    if(data->print_partition == 1){ //prints info only if user passes argument
        if(data->threads_partition == 0){ //prints row-wise partitioning info
            int offset = 0;
            for(int i=0; i<data->num_threads; i++){
                if(myid == i){
                    printf("tid  %d: ", i);
                    printf("rows:  %d:%d (%d) ", offset, data->work_per_thread[i]-1+offset, data->work_per_thread[i]);
                    printf("cols:  0:%d (%d)\n", data->cols-1, data->cols);
                    fflush(stdout);
                }
                offset = offset + data->work_per_thread[i];
            }
        }
        else{ //prints column-wise partitioning info
            int offset = 0;
            for(int i=0; i<data->num_threads; i++){
                if(myid == i){
                    printf("tid  %d: ", i);
                    printf("rows:  0:%d (%d) ", data->rows-1, data->rows);
                    printf("cols:  %d:%d (%d)\n", offset, data->work_per_thread[i]-1+offset, data->work_per_thread[i]);
                    fflush(stdout);
                }
                offset = offset + data->work_per_thread[i];
            }
        }
    }


    //Gets which rows or columns the thread needs to work on
    if(data->threads_partition == 0){ //Gets row-wise partitioning info
        int offset = 0;
        for(int i=0; i<data->num_threads; i++){
            if(myid == i){
                startRow = offset;
                endRow = data->work_per_thread[i]-1+offset;
                startCol = 0;
                endCol = data->cols-1;
            }
            offset = offset + data->work_per_thread[i];
        }
    }
    else{ //Gets column-wise partitioning info
        int offset = 0;
        for(int i=0; i<data->num_threads; i++){
            if(myid == i){
                startRow = 0;
                endRow = data->rows-1;
                startCol = offset;
                endCol = data->work_per_thread[i]-1+offset;
            }
            offset = offset + data->work_per_thread[i];
        }
    }


    if(data->output_mode == 0){ //No display mode
        for(int round = 1; round <= data->iters; round++){
            int live_per_round = 0;
            if(myid == 0){
                total_live = 0;
            }
            pthread_barrier_wait(&my_barrier);
            for (int i = startRow; i <= endRow; ++i) {
                for (int j = startCol; j <= endCol; ++j) {
                    if(celllives_celldies(data, i, j) == 1){
                        data->currentGraph[coord_index(data, i, j)] = 1;
                        live_per_round++;
                    }
                    else{
                        data->currentGraph[coord_index(data, i, j)] = 0;
                    }
                } 
            }
            pthread_barrier_wait(&my_barrier);

            pthread_mutex_lock(&my_mutex);
            total_live = total_live + live_per_round;
            pthread_mutex_unlock(&my_mutex);
            
            pthread_barrier_wait(&my_barrier);

            if(myid == 0){
                if (data->currentGraph == data->graph) {
                    data->currentGraph = data->graph2;
                }
                else {
                    data->currentGraph = data->graph;
                } 
            }
        }
    }
    else if(data->output_mode == 1){ //ASCII mode
        for(int round = 1; round <= data->iters; round++){
            int live_per_round = 0;
            if(myid == 0){
                total_live = 0;
            }
            pthread_barrier_wait(&my_barrier);
            if (system("clear")) { perror("clear"); exit(1); }
            for (int i = startRow; i <= endRow; ++i) {
                for (int j = startCol; j <= endCol; ++j) {
                    if(celllives_celldies(data, i, j) == 1){
                        data->currentGraph[coord_index(data, i, j)] = 1;
                        live_per_round++;
                    }
                    else{
                        data->currentGraph[coord_index(data, i, j)] = 0;
                    }
                } 
            }
            pthread_mutex_lock(&my_mutex);
            total_live = total_live + live_per_round;
            pthread_mutex_unlock(&my_mutex);

            pthread_barrier_wait(&my_barrier);
            if(myid == 0){
                print_board(data, round);
                usleep(SLEEP_USECS);
            }
            pthread_barrier_wait(&my_barrier);
            
            if(myid == 0){
                if (data->currentGraph == data->graph) {
                    data->currentGraph = data->graph2;
                }
                else {
                    data->currentGraph = data->graph;
                } 
            }
        }
    }
    else{ //Visi mode
        for(int round = 1; round <= data->iters; round++){
            for (int i = startRow; i <= endRow; ++i) {
                for (int j = startCol; j <= endCol; ++j) {
                    if(celllives_celldies(data, i, j) == 1){
                        data->currentGraph[coord_index(data, i, j)] = 1;
                    }
                    else{
                        data->currentGraph[coord_index(data, i, j)] = 0;
                    }
                }
            }
            pthread_barrier_wait(&my_barrier);
            update_colors(data, startRow, endRow, startCol, endCol, myid);
            draw_ready(data->handle);
            pthread_barrier_wait(&my_barrier);
            if(myid == 0){
                usleep(SLEEP_USECS);
            }
            pthread_barrier_wait(&my_barrier);
            if(myid == 0){
                if (data->currentGraph == data->graph){
                    data->currentGraph = data->graph2;
                }  
                else {
                    data->currentGraph = data->graph;
                }
            }
        }
    }

    return NULL;
}

/* Partitions graphs row-wise or column-wise for parallelization
*       data: pointer to gol_data struct 
*/
int *partition_graph(struct gol_data *data){
    if(data->threads_partition == 0){ //partition row-wise
        int num_extra_rows = data->rows%data->num_threads;
        printf("\n");
        int *num_rows = malloc(sizeof(int) * data->num_threads);
        for(int i=0; i<data->num_threads; i++){
            num_rows[i] = data->rows/data->num_threads;
        }
        for(int j=0; j<data->num_threads; j++){
            if(num_extra_rows > 0){
                num_rows[j] = num_rows[j] + 1;
                num_extra_rows--;
            }
        }

        return num_rows;
    }
    else{ //partition column-wise
        int num_extra_cols = data->cols%data->num_threads;
        printf("\n");
        int *num_cols = malloc(sizeof(int) * data->num_threads);
        for(int i=0; i<data->num_threads; i++){
            num_cols[i] = data->cols/data->num_threads;
        }
        for(int j=0; j<data->num_threads; j++){
            if(num_extra_cols > 0){
                num_cols[j] = num_cols[j] + 1;
                num_extra_cols--;
            }
        }

        return num_cols;
    }
}


/* initialize the gol game state from command line arguments
 *       argv[1]: name of file to read game config state from
 *       argv[2]: run mode value
 * data: pointer to gol_data struct to initialize
 * argv: command line args
 *       argv[1]: name of file to read game config state from
 *       argv[2]: run mode
 * returns: 0 on success, 1 on error
 */
int init_game_data_from_args(struct gol_data *data, char **argv) {
    if(atoi(argv[2]) == 0){
        data->output_mode = OUTPUT_NONE;
    }
    else if(atoi(argv[2]) == 1){
        data->output_mode = OUTPUT_ASCII;
    }
    else{
        data->output_mode = OUTPUT_VISI;
    }

    FILE *infile;
    infile = fopen(argv[1], "r");
    if (infile == NULL) {
        printf("Error: failed to open file: %s\n", argv[1]);
        exit(1);
    }
    int ret = fscanf(infile, "%d", &data->rows);
    if (ret == 0) {
        printf("Improper file format.\n");
        exit(1);
    }
    fscanf(infile, "%d", &data->cols);
    fscanf(infile, "%d", &data->iters);
    fscanf(infile, "%d", &data->initalive);

    //partitions work for each thread
    data->num_threads = atoi(argv[3]);
    data->threads_partition = atoi(argv[4]);
    data->print_partition = atoi(argv[5]);
    data->work_per_thread = partition_graph(data);
    
    /* Initialize the mutex. */
    if (pthread_mutex_init(&my_mutex, NULL)) {
        printf("pthread_mutex_init error\n");
        exit(1);
    }

    /* Initialize the barrier with the number of threads that will
     * synchronize on it: */
    if (pthread_barrier_init(&my_barrier, NULL, data->num_threads)) {
        printf("pthread_barrier_init error\n");
        exit(1);
    }

    //initializes graphs
    int *graph = (malloc((data->rows)*(data->cols)*sizeof(int)));
    int *graph2 = (malloc((data->rows)*(data->cols)*sizeof(int)));
    data->graph = graph;
    data->graph2 = graph2;
    
    for (int i = 0; i<(data->rows)*(data->cols); i++){
        graph[i] = 0;
        graph2[i] = 0;
    }
    for (int i = 0; i<data->initalive;i++){
        int row = 0;
        int col = 0;
        fscanf(infile,"%d%d", &row, &col);
        int index = coord_index(data, row, col);
        data->graph[index] = 1;
        data->graph2[index] = 1;
    }

    /* Barrier malloc*/
    data->done = malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(data->done, NULL, data->num_threads);
    data->currentGraph = data->graph;
    if (!data->graph) {
        printf("ERROR barrier malloc\n");
        return 1;
    }


    return 0;
}


/* gets index of a 2d graph point for a 2d array*/
int coord_index(struct gol_data *data, int i, int j) {
    return (i * data->cols) + j;
}

/* checks to see if the cell is still alive
data: our gol_data variable
i: y coordinate
j: x coordinate
returns 1 if it is alive and 0 if it is not
*/
int check_alive_dead(struct gol_data *data, int i, int j) {
    if (data->currentGraph == data->graph){
        if(data->graph2[coord_index(data, i, j)] == 1){
            return 1;
        }
    }
    else {
        if(data->graph[coord_index(data, i, j)] == 1){
            return 1;
        }
    }
    return 0;
}

/* looks at neighbors of a specified location.
    gol_data: struct of data about the game
    i: y coordinate
    j: x coordinate
    returns the number of cells still alive
*/
int check_neighbors(struct gol_data *data, int i, int j) {
    int num_alive = 0;

    int top_neighbors_outofbounds = 0;
    int bottom_neighbors_outofbounds = 0;
    int left_neighbors_outofbounds = 0;
    int right_neighbors_outofbounds = 0;

    int topleft_corner_case = 0;
    int topright_corner_case = 0;
    int bottomleft_corner_case = 0;
    int bottomright_corner_case = 0;

    //checks top three neighboring cells if they are in the bounds of the graph
    if(i == 0){
        top_neighbors_outofbounds = 1;
    }

    //checks bottom three neighboring cells if they are in the bounds of the graph
    if(i == (data->rows - 1)){
        bottom_neighbors_outofbounds = 1;
    }

    //checks left three neighboring cells if they are in the bounds of the graph
    if(j == 0){
        left_neighbors_outofbounds = 1;
    }

    //checks right three neighboring cells if they are in the bounds of the graph
    if(j == (data->cols - 1)){
        right_neighbors_outofbounds = 1;
    }



    //checks if the cell is in the top left corner of the graph
    if((top_neighbors_outofbounds == 1) && (left_neighbors_outofbounds == 1)){
        topleft_corner_case = 1;
    }

    //checks if the cell is in the top right corner of the graph
    if((top_neighbors_outofbounds == 1) && (right_neighbors_outofbounds == 1)){
        topright_corner_case = 1;
    }

    //checks if the cell is in the bottom left corner of the graph
    if((bottom_neighbors_outofbounds == 1) && (left_neighbors_outofbounds == 1)){
        bottomleft_corner_case = 1;
    }

    //checks if the cell is in the bottom right corner of the graph
    if((bottom_neighbors_outofbounds == 1) && (right_neighbors_outofbounds == 1)){
        bottomright_corner_case = 1;
    }

    //checks top left neighbor cell
    if((top_neighbors_outofbounds == 0) && (left_neighbors_outofbounds == 0)){
        if(check_alive_dead(data, i-1, j-1) == 1){
            num_alive++;
        }
    }
    else{
        if(topleft_corner_case == 1){
            if(check_alive_dead(data, data->rows-1, data->cols-1) == 1){
                num_alive++;
            }
        }
        else if(top_neighbors_outofbounds == 1){
            if(check_alive_dead(data, data->rows-1, j-1) == 1){
                num_alive++;
            }
        }
        else{
            if(check_alive_dead(data, i-1, data->cols-1) == 1){
                num_alive++;
            }
        }
    }

    //checks top middle neighbor cell
    if(top_neighbors_outofbounds == 0){
        if(check_alive_dead(data, i-1, j) == 1){
            num_alive++;
        }
    }
    else{
        if(check_alive_dead(data, data->rows-1, j) == 1){
            num_alive++;
        }
    }

    //checks top right neighbor cell
    if((top_neighbors_outofbounds == 0) && (right_neighbors_outofbounds == 0)){
        if(check_alive_dead(data, i-1, j+1) == 1){
            num_alive++;
        }
    }
    else{
        if(topright_corner_case == 1){
            if(check_alive_dead(data, data->rows-1, 0) == 1){
                num_alive++;
            }
        }
        else if(top_neighbors_outofbounds == 1){
            if(check_alive_dead(data, data->rows-1, j+1) == 1){
              num_alive++;
            }
        }
        else{
            if(check_alive_dead(data, i-1, 0) == 1){
                num_alive++;
            }
        }
    }

    //checks middle left neighbor cell
    if(left_neighbors_outofbounds == 0){
        if(check_alive_dead(data, i, j-1) == 1){
            num_alive++;
        }
    }
    else{
        if(check_alive_dead(data, i, data->cols-1) == 1){
            num_alive++;
        }
    }

    //checks middle right neighbor cell
    if(right_neighbors_outofbounds == 0){
        if(check_alive_dead(data, i, j+1) == 1){
            num_alive++;
        }
    }
    else{
        if(check_alive_dead(data, i, 0) == 1){
            num_alive++;
        }
    }

    //checks bottom left neighbor cell
    if((bottom_neighbors_outofbounds == 0) && (left_neighbors_outofbounds == 0)){
        if(check_alive_dead(data, i+1, j-1) == 1){
            num_alive++;
        }
    }
    else{
        if(bottomleft_corner_case == 1){
            if(check_alive_dead(data, 0, data->cols-1) == 1){
                num_alive++;
            }
        }
        else if(bottom_neighbors_outofbounds == 1){
            if(check_alive_dead(data, 0, j-1) == 1){
                num_alive++;
            }
        }
        else{
            if(check_alive_dead(data, i+1, data->cols-1) == 1){
                num_alive++;
            }
        }
    }

    //checks bottom middle neighbor cell
    if(bottom_neighbors_outofbounds == 0){
        if(check_alive_dead(data, i+1, j) == 1){
            num_alive++;
        }
    }
    else{
        if(check_alive_dead(data, 0, j) == 1){
            num_alive++;
        }
    }

    //checks bottom right neighbor cell
    if((bottom_neighbors_outofbounds == 0) && (right_neighbors_outofbounds == 0)){
        if(check_alive_dead(data, i+1, j+1) == 1){
            num_alive++;
        }
    }
    else{
        if(bottomright_corner_case == 1){
            if(check_alive_dead(data, 0, 0) == 1){
                num_alive++;
            }
        }
        else if(bottom_neighbors_outofbounds == 1){
            if(check_alive_dead(data, 0, j+1) == 1){
                num_alive++;
            }
        }
        else{
            if(check_alive_dead(data, i+1, 0) == 1){
                num_alive++;
            }    
        }
    }
    return num_alive;
}

/*
takes in the data and a coordinate, uses other functions to decide whether the
cell lives or dies
data: game data
i: y coordinate
j: x coordinate
*/
int celllives_celldies(struct gol_data *data, int i, int j){
    int num_neighbors_alive = check_neighbors(data, i, j);
    if (data->currentGraph == data->graph) {
        //dead cell lives if it has 3 alive neighbors, otherwise it stays dead
        if(data->graph2[coord_index(data, i, j)] == 0){
            if(num_neighbors_alive == 3){
                return 1;
            }
        }
        else {
        //alive cell stays alive if it has 2 or 3 alive neighbors, otherwise it 
        //dies    
           if((num_neighbors_alive == 2) || (num_neighbors_alive == 3)){
                return 1;
            }  
        }
    }
    else{
        if(data->graph[coord_index(data, i, j)] == 0){
        //dead cell lives if it has 3 alive neighbors, otherwise it stays dead
            if(num_neighbors_alive == 3){
                return 1;
            }
        }
        else {
        //alive cell stays alive if it has 2 or 3 alive neighbors, otherwise it 
        //dies
           if((num_neighbors_alive == 2) || (num_neighbors_alive == 3)){
                return 1;
            }  
        }
    }
    return 0;
}

/*
Used for option 2 when we show the animation, updates the colors of the
animation as the data is updated in the play gol function.
data: game data
i: y coordinate
j: x coordinate
startRow: starting row of part of graph to be worked on
endRow: ending row of part of graph to be worked on
startCol: starting col of part of graph to be worked on
endCol: ending col of part of graph to be worked on
*/
void update_colors(struct gol_data *data, int startRow, int endRow, int startCol, int endCol, int myid) {
    int i, j, r, c, index, buff_i;
    color3 *buff;

    buff = data->image_buff;  // just for readability
    r = data->rows;
    c = data->cols;

    for (i = startRow; i <= endRow; i++) {
        for (j = startCol; j <= endCol; j++) {
            index = i*c + j;
            // translate row index to y-coordinate value
            // in the image buffer, r,c=0,0 is assumed to be the _lower_ left
            // in the grid, r,c=0,0 is _upper_ left.
            buff_i = (r - (i+1))*c + j;

            // update animation buffer
            if (data->currentGraph[index] == 0) {
                buff[buff_i] = colors[((myid)%8)];
            } else {
                buff[buff_i] = c3_green;
            }
        }
    }

    // force threads to wait until all are done before each starts next
    pthread_barrier_wait(data->done);
}

/* the gol application main loop function:
 *  runs rounds of GOL,
 *    * updates program state for next round (world and total_live)
 *    * performs any animation step based on the output/run mode
 *
 *   data: pointer to a struct gol_data  initialized with
 *         all GOL game playing state
 */
void *play_gol(void *args) {
    struct gol_data *data;
    data = (struct gol_data *)args;
    pthread_t *tids;
    int ret, i, ntids;
    struct thread_info *tid_args;

    ntids = data->num_threads;

    tids = malloc(sizeof(pthread_t) * ntids);
    if (!tids) {
        perror("malloc: pthread_t array");
        exit(1);
    }

    tid_args = malloc(sizeof(struct thread_info) * ntids);
    if (!tid_args) {
        perror("malloc: tid_args array");
        exit(1);
    }

    for (i = 0; i < ntids; i++) {
        tid_args[i].tid = i;
        tid_args[i].data = data;

        ret = pthread_create(&tids[i], 0, thread_function, &tid_args[i]);

        if (ret) {
            perror("Error pthread_create\n");
            exit(1);
        }
    }

    // wait for all worker threads to exit
    for (i = 0; i < ntids; i++) {
        pthread_join(tids[i], 0);
    }

    /* Clean up memory. */
    free(tid_args);
    tid_args = NULL;
    free(tids);
    tids = NULL;

    return NULL;
}

/* Print the board to the terminal.
 *   data: gol game specific data
 *   round: the current round number
 *
 * NOTE: You may add extra printfs if you'd like, but please
 *       leave these fprintf calls exactly as they are to make
 *       grading easier!
 */
void print_board(struct gol_data *data, int round) {
    int i, j;

    /* Print the round number. */
    fprintf(stderr, "Round: %d\n", round);
    for (i = 0; i < data->rows; ++i) {
        for (j = 0; j < data->cols; ++j) {
            if(round == 0){
                if(data->currentGraph[coord_index(data, i, j)] == 1){
                    fprintf(stderr, " @"); //leave alone
                    pthread_mutex_lock(&my_mutex);
                    total_live++; 
                    pthread_mutex_unlock(&my_mutex);   
                }
                else{
                    fprintf(stderr, " ."); //leave alone
                }
            }
            else{
                if(data->currentGraph[coord_index(data, i, j)] == 1){
                    fprintf(stderr, " @"); //leave alone
                }
                else{
                    fprintf(stderr, " ."); //leave alone
                }
            }
        }
        fprintf(stderr, "\n");
    }

    /* Print the total number of live cells. */
    fprintf(stderr, "Live cells: %d\n\n", total_live);
}


/**********************************************************/
/***** START: DO NOT MODIFY THIS CODE *****/
/* initialize ParaVisi animation */
int setup_animation(struct gol_data* data) {
    /* connect handle to the animation */
    data->handle = init_pthread_animation(data->num_threads, data->rows,
            data->cols, visi_name);
    if (data->handle == NULL) {
        printf("ERROR init_pthread_animation\n");
        exit(1);
    }
    // get the animation buffer
    data->image_buff = get_animation_buffer(data->handle);
    if(data->image_buff == NULL) {
        printf("ERROR get_animation_buffer returned NULL\n");
        exit(1);
    }
    return 0;
}

/* sequential wrapper functions around ParaVis library functions */
void (*mainloop)(struct gol_data *data);

void* seq_do_something(void * args){
    mainloop((struct gol_data *)args);
    return 0;
}

int connect_animation(void (*applfunc)(struct gol_data *data),
        struct gol_data* data)
{
    pthread_t pid;

    mainloop = applfunc;
    if( pthread_create(&pid, NULL, seq_do_something, (void *)data) ) {
        printf("pthread_created failed\n");
        return 1;
    }
    return 0;
}
/***** END: DO NOT MODIFY THIS CODE *****/
/******************************************************/
