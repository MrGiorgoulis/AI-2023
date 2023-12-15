// -------------------------------------------------------------
// This program is an adaptation of the original program, which solved
// 3x3 puzzles using four algorithms. The changes made altered the type
// of problem the program solves as following:
//
// This program takes as input 2 integers and shows the steps
// it takes to make the first integer into the second integer
// with the use of 6 distinct allowed operations on the integer-register:
// - Increase the number by 1
// - Decrease the number by 1
// - Double the number
// - Half the number
// - Number to the power of 2
// - Square root of number
//
// The program solves the problems using four algorithms:
// - Depth first search
// - Breadth first search
// - Best first search
// - A*
//
// Author: Georgios Katzos, December 2023
// Original Author: Ioannis Refanidis, January 2008
//
// --------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define breadth 1		// Constants denoting the four algorithms
#define depth	2
#define best	3
#define astar	4

#define increment 1
#define decrease 2
#define doubling 3
#define halfing 4
#define square 5
#define sqroot 6

struct tree_node
{
    int reg;
    int h;				// the value of the heuristic function for this node
    int g;				// the cost of this node from the root of the search tree
    int f;              // The estimated cost to the destination node.f=0 or f=h or f=h+g, depending on the search algorithm used.
    int d;          //The depth of the tree where the node is at
    struct tree_node *parent;	// pointer to the parent node (NULL for the root).
    int action;			// The direction of the last move
};

// A node of the frontier. Frontier is kept as a double-linked list,
// for efficiency reasons for the breadth-first search algorithm.
struct frontier_node
{
    struct tree_node *n;			// pointer to a search-tree node
    struct frontier_node *previous;		// pointer to the previous frontier node
    struct frontier_node *next;		// pointer to the next frontier node
};

struct frontier_node *frontier_head=NULL;	// The one end of the frontier
struct frontier_node *frontier_tail=NULL;	// The other end of the frontier

clock_t t1;				// Start time of the search algorithm
clock_t t2;				// End time of the search algorithm
#define TIMEOUT		60	// Program terminates after TIMOUT secs

int solution_length;	// The lenght of the solution table.
int **solution;		// Pointer to a dynamic table with the moves of the solution.

int input_value;
int goal;


// Reading run-time parameters.
int get_method(char* s)
{
    if (strcmp(s,"breadth")==0)
        return  breadth;
    else if (strcmp(s,"depth")==0)
        return depth;
    else if (strcmp(s,"best")==0)
        return best;
    else if (strcmp(s,"astar")==0)
        return astar;
    else
        return -1;
}

// This function checks whether a node in the search tree
// holds exactly the same value with at least one of its
// predecessors. This function is used when creating the childs
// of an existing search tree node, in order to check for each one of the childs
// whether this appears in the path from the root to its parent.
// This is a moderate way to detect loops in the search.
// Inputs:
//		struct tree_node *new_node	: A search tree node (usually a new one)
// Output:
//		1 --> No coincidence with any predecessor
//		0 --> Loop detection
int check_with_parents(struct tree_node *new_node)
{
    struct tree_node *parent=new_node->parent;
    while (parent!=NULL)
    {
        if (new_node->reg == parent->reg)
            return 0;
        parent=parent->parent;
    }
    return 1;
}

// This function calculated the heuristic function value of the input value,
// given the goal, as the 1/4 of the absolute distance of the
// number/register and the goal number.
//
// Input:
//      int reg; An integer
// Output:
//      The value of the heuristic function for this node.
int heuristic(struct tree_node * node)
{
    int absolute_distance = abs(goal-node->reg);
    return ceil((double)absolute_distance / 4);
//    return log2(absolute_distance);
//        return log(absolute_distance)/log(2);
}

// This function adds a pointer to a new leaf search-tree node at the front of the frontier.
// This function is called by the depth-first search algorithm.
// Inputs:
//		struct tree_node *node	: A (leaf) search-tree node.
// Output:
//		0 --> The new frontier node has been added successfully.
//		-1 --> Memory problem when inserting the new frontier node .
int add_frontier_front(struct tree_node *node)
{
    // Creating the new frontier node
    struct frontier_node *new_frontier_node=(struct frontier_node*)
            malloc(sizeof(struct frontier_node));
    if (new_frontier_node==NULL)
        return -1;

    new_frontier_node->n = node;
    new_frontier_node->previous = NULL;
    new_frontier_node->next = frontier_head;

    if (frontier_head==NULL)
    {
        frontier_head=new_frontier_node;
        frontier_tail=new_frontier_node;
    }
    else
    {
        frontier_head->previous=new_frontier_node;
        frontier_head=new_frontier_node;
    }
    return 0;
}

// This function adds a pointer to a new leaf search-tree node at the back of the frontier.
// This function is called by the breadth-first search algorithm.
// Inputs:
//		struct tree_node *node	: A (leaf) search-tree node.
// Output:
//		0 --> The new frontier node has been added successfully.
//		-1 --> Memory problem when inserting the new frontier node.
int add_frontier_back(struct tree_node *node)
{
    // Creating the new frontier node
    struct frontier_node *new_frontier_node=(struct frontier_node*) malloc(sizeof(struct frontier_node));
    if (new_frontier_node==NULL)
        return -1;

    new_frontier_node->n=node;
    new_frontier_node->next=NULL;
    new_frontier_node->previous=frontier_tail;

    if (frontier_tail==NULL)
    {
        frontier_head=new_frontier_node;
        frontier_tail=new_frontier_node;
    }
    else
    {
        frontier_tail->next=new_frontier_node;
        frontier_tail=new_frontier_node;
    }
    return 0;
}

// This function adds a pointer to a new leaf search-tree node within the frontier.
// The frontier is always kept in increasing order wrt the f values of the corresponding
// search-tree nodes. The new frontier node is inserted in order.
// This function is called by the heuristic search algorithm.
// Inputs:
//		struct tree_node *node	: A (leaf) search-tree node.
// Output:
//		0 --> The new frontier node has been added successfully.
//		-1 --> Memory problem when inserting the new frontier node .
int add_frontier_in_order(struct tree_node *node)
{
    // Creating the new frontier node
    struct frontier_node *new_frontier_node=(struct frontier_node*)
            malloc(sizeof(struct frontier_node));
    if (new_frontier_node==NULL)
        return -1;

    new_frontier_node->n=node;
    new_frontier_node->previous=NULL;
    new_frontier_node->next=NULL;

    if (frontier_head==NULL)
    {
        frontier_head=new_frontier_node;
        frontier_tail=new_frontier_node;
    }
    else
    {
        struct frontier_node *pt;
        pt=frontier_head;

        // Search in the frontier for the first node that corresponds to either a larger f value
        // or to an equal f value but larger h value
        // Note that for the best first search algorithm, f and h values coincide.
        while (pt!=NULL && (pt->n->f<node->f || (pt->n->f==node->f && pt->n->h<node->h)))
            pt=pt->next;

        if (pt!=NULL)
        {
            // new_frontier_node is inserted before pt .
            if (pt->previous!=NULL)
            {
                pt->previous->next=new_frontier_node;
                new_frontier_node->next=pt;
                new_frontier_node->previous=pt->previous;
                pt->previous=new_frontier_node;
            }
            else
            {
                // In this case, new_frontier_node becomes the first node of the frontier.
                new_frontier_node->next=pt;
                pt->previous=new_frontier_node;
                frontier_head=new_frontier_node;
            }
        }
        else
        {
            // if pt==NULL, new_frontier_node is inserted at the back of the frontier
            frontier_tail->next=new_frontier_node;
            new_frontier_node->previous=frontier_tail;
            frontier_tail=new_frontier_node;
        }
    }
    return 0;
}

// This function checks whether an integer reg is a perfect square of sqrt(reg).
// Input:
//      int reg; An integer
// Output:
//      - 1 if the integer reg is a perfect square of sqrt(reg)
//      - 0 if the integer reg is not a perfect square of sqrt(reg)
int check_sqrt(int reg){
    double number = (double) reg;
    double squareRoot = sqrt(number);
    int squareRootInt = (int)squareRoot;

    if (squareRootInt * squareRootInt == (int)number)
        return 1;
    return 0;
}

// This function calculates and returns the cost of each operation.
// Inputs:
//      int reg;    An integer
//      int action; The operation performed
// Output:
//      The cost of the operation based on operation performed and the value of reg.
int calculate_cost(int reg, int action){

    if(action==1)   //Increase
        return 2;
    else if(action==2)  //Decrease
        return 2;
    else if(action==3)  //Double
        return ceil((double)reg/2) + 1;
    else if(action==4)  //Half
        return ceil((double)reg/4) + 1;
    else if(action==5)  //Square
        return ((pow(reg,2)-reg)/4 +1);
    else if(action==6)  //Square Root
        return (reg - sqrt(reg))/4 + 1;

    return -1;
}


// This function expands a leaf-node of the search tree.
// A leaf-node may have up to 4 childs. A table with 4 pointers
// to these childs is created, with NULLs for those childrens that do not exist.
// In case no child exists (due to loop-detections), the table is not created
// and a 'higher-level' NULL indicates this situation.
// Inputs:
//		struct tree_node *current_node	: A leaf-node of the search tree.
// Output:
//		The same leaf-node expanded with pointers to its children (if any).
int find_children(struct tree_node *current_node, int method)
{
    // Find the blank position in the current puzzle
    //find_blank(current_node->p,&i,&j);
    int reg = current_node->reg;

    //Increment by 1
    if(reg<pow(10,9)){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = increment;
        child->g = current_node->g + calculate_cost(current_node->reg, increment);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = reg+1;
        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);

            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }

    //Decrease by 1
    if(reg>0){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = decrease;
        child->g = current_node->g + calculate_cost(current_node->reg, decrease);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = reg-1;

        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);
            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }

    //Double
    if(reg>0 && 2*reg<=pow(10,9)){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = doubling;
        child->g = current_node->g + calculate_cost(current_node->reg, doubling);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = reg*2;

        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);
            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }

    //Half
    if(reg>0){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = halfing;
        child->g = current_node->g + calculate_cost(current_node->reg, halfing);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = reg/2;  //Integer division returns floor by definition

        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);
            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }

    //Square
    if(pow(reg,2)<=pow(10,9)){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = square;
        child->g = current_node->g + calculate_cost(current_node->reg, square);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = pow(reg, 2);

        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);
            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }

    //Sqrt
    if(reg>1 && check_sqrt(reg)){
        // Initializing the new child
        struct tree_node *child=(struct tree_node*) malloc(sizeof(struct tree_node));
        if (child==NULL) return -1;

        child->parent = current_node;
        child->action = sqroot;
        child->g = current_node->g + calculate_cost(current_node->reg, sqroot);
        child->d = current_node->d + 1;

        //Computing the value for the new child
        child->reg = (int)sqrt((double)reg);

        // Check for loops
        if (!check_with_parents(child))
            // In case of loop detection, the child is deleted
            free(child);
        else
        {
            // Computing the heuristic value
            child->h=heuristic(child);
            if (method==best)
                child->f = child->h;
            else if (method==astar)
                child->f = child->g + child->h;
            else
                child->f = 0;

            int err=0;
            if (method==depth)
                err=add_frontier_front(child);
            else if (method==breadth)
                err=add_frontier_back(child);
            else if (method==best || method==astar)
                err=add_frontier_in_order(child);
            if (err<0)
                return -1;
        }
    }
    return 1;
}

// Auxiliary function that displays a message in case of wrong input parameters.
void syntax_message()
{
    printf("register <method> <input-int> <output-int>\n\n");
    printf("where: ");
    printf("<method> = breadth|depth|best|astar\n");
    printf("<input-int> is an integer given to the program as input.\n");
    printf("<output-int> is the final output integer.\n");
}

// Giving a (solution) leaf-node of the search tree, this function computes
// the moves of the blank that have to be done, starting from the root puzzle,
// in order to go to the leaf node's puzzle.
// Inputs:
//		struct tree_node *solution_node	: A leaf-node
// Output:
//		The sequence of blank's moves that have to be done, starting from the root puzzle,
void extract_solution(struct tree_node *solution_node)
{
    int i;

    struct tree_node *temp_node=solution_node;
    solution_length = solution_node->d;
    solution= (int**) malloc(solution_length*sizeof(int*));
    for(i=0;i<solution_length;i++){
        solution[i] = (int*)malloc(3*sizeof(int));
    }
    temp_node=solution_node;
    i=solution_length;
    while (temp_node->parent!=NULL)
    {

        i--;
        solution[i][0] = temp_node->action;
        solution[i][1] = temp_node->g;
        solution[i][2] = temp_node->reg;
        temp_node=temp_node->parent;
    }
}

// This function writes the solution into a file
// Inputs:
//		char* filename	: The name of the file where the solution will be written.
// Outputs:
//		Nothing (apart from the new file)
void write_solution_to_file(char* filename, int solution_length, int**solution)
{
    int i;
    FILE *fout;
    fout=fopen(filename,"w");
    if (fout==NULL)
    {
        printf("Cannot open output file to write solution.\n");
        printf("Now exiting...");
        return;
    }

    printf("Solution Length: %d", solution_length);
    int total_cost = solution[solution_length-1][1];
    printf("\nTotal cost: %d", total_cost);

    fprintf(fout,"%d %d\n",solution_length, total_cost);
    for (i=0;i<solution_length;i++){

        if(i==0){
            int temp = solution[i][0];
            if(temp == increment)
                fprintf(fout,"increase   %d %d\n",input_value,solution[i][1]);
            else if(temp == decrease)
                fprintf(fout,"decrease   %d %d\n",input_value,solution[i][1]);
            else if(temp==doubling)
                fprintf(fout,"double     %d %d\n",input_value,solution[i][1]);
            else if(temp==halfing)
                fprintf(fout,"half       %d %d\n",input_value,solution[i][1]);
            else if(temp==square)
                fprintf(fout,"square     %d %d\n",input_value,solution[i][1]);
            else if(temp==sqroot)
                fprintf(fout,"sqrt       %d %d\n",input_value,solution[i][1]);
        }
        else{
            int temp = solution[i][0];
            if(temp == increment)
                fprintf(fout,"increase   %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
            else if(temp == decrease)
                fprintf(fout,"decrease   %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
            else if(temp==doubling)
                fprintf(fout,"double     %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
            else if(temp==halfing)
                fprintf(fout,"half       %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
            else if(temp==square)
                fprintf(fout,"square     %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
            else if(temp==sqroot)
                fprintf(fout,"sqrt       %d %d\n",solution[i-1][2],solution[i][1]-solution[i-1][1]);
        }

    }
    fclose(fout);
}

// This function initializes the search, i.e. it creates the root node of the search tree
// and the first node of the frontier.
void initialize_search(int reg, int method)
{
    struct tree_node *root=NULL;	// the root of the search tree.

    // Initialize search tree
    root=(struct tree_node*) malloc(sizeof(struct tree_node));
    root->parent=NULL;
    root->action=-1;
    root->reg = reg;
    root->g=0;
    root->d=0;
    root->h= heuristic(root);

    if (method==best)
        root->f=root->h;
    else if (method==astar)
        root->f=root->g+root->h;
    else
        root->f=0;

    // Initialize frontier
    add_frontier_front(root);
}

// This function implements at the higest level the search algorithms.
// The various search algorithms differ only in the way the insert
// new nodes into the frontier, so most of the code is commmon for all algorithms.
// Inputs:
//		Nothing, except for the global variables root, frontier_head and frontier_tail.
// Output:
//		NULL --> The problem cannot be solved
//		struct tree_node*	: A pointer to a search-tree leaf node that corresponds to a solution.
struct tree_node *search(int method)
{
    clock_t t;
    struct frontier_node *temp_frontier_node;
    struct tree_node *current_node;

    while (frontier_head!=NULL)
    {
        t=clock();
        if (t-t1 > CLOCKS_PER_SEC*TIMEOUT)
        {
            printf("Timeout\n");
            return NULL;
        }

        // Extract the first node from the frontier
        current_node = frontier_head->n;

        if (current_node->reg==goal)
            return current_node;

        // Delete the first node of the frontier
        temp_frontier_node=frontier_head;
        frontier_head = frontier_head->next;
        free(temp_frontier_node);
        if (frontier_head==NULL)
            frontier_tail=NULL;
        else
            frontier_head->previous=NULL;

        // Find the children of the extracted node
        int err=find_children(current_node, method);

        if (err<0)
        {
            printf("Memory exhausted while creating new frontier node. Search is terminated...\n");
            return NULL;
        }
    }
    return NULL;
}

int main(int argc, char** argv)
{
    int err;
    struct tree_node *solution_node;
    int reg;		// The initial puzzle read from a file
    int method;				// The search algorithm that will be used to solve the puzzle.

    if (argc!=5)
    {
        printf("Wrong number of arguments. Use correct syntax:\n");
        syntax_message();
        return -1;
    }
    printf("Argv[1]: %s\n", argv[1]);
    printf("Argv[2]: %s\n", argv[2]);
    printf("Argv[3]: %s\n", argv[3]);
    printf("Argv[4]: %s\n", argv[4]);

    method=get_method(argv[1]);
    if (method<0)
    {
        printf("Wrong method. Use correct syntax:\n");
        syntax_message();
        return -1;
    }

    reg = atoi(argv[2]);
    input_value = reg;
    goal = atoi(argv[3]);

    t1=clock();
    initialize_search(reg, method);

    solution_node = search(method);			// The main call

    t2=clock();

    if (solution_node!=NULL)
        extract_solution(solution_node);
    else
        printf("No solution found.\n");

    if (solution_node!=NULL)
    {
        printf("Solution found! (%d steps)\n",solution_length);
        printf("Time spent: %f secs\n",((float) t2-t1)/CLOCKS_PER_SEC);
        write_solution_to_file(argv[4], solution_length, solution);
    }

    return 0;
}
