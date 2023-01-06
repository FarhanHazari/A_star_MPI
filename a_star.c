#include "tools.h"
#include "heap.h"
#include "string.h"
#include <mpi.h>

#define MAX_NEIGHBORS 8
#define INIT_HEAP_CAPACITY 4

// A heuristic function is a function h() that returns a (double) distance
// between a start and end position of the grid. The function could also
// depend on the grid (e.g. the number of walls encountered by the start-finish segment),
// but the latter parameter need not be used. You can define your own heuristic.
typedef double (*heuristic)(position, position, grid *);

// "Bird's eye view" heuristic for A*.
double hvo(position s, position t, grid *G)
{
    double x = t.x - s.x;
    double y = t.y - s.y;
    return sqrt(x * x + y * y);
}

// "Alpha x Bird's eye view" heuristic for A*.
static double alpha = 0; // 0 = Djikstra, 1 = A*, 2 = Approximation ...
double halpha(position s, position t, grid *G)
{
    return alpha * hvo(s, t, G);
}

double weight[] = {
    1.0,   // V_FREE
    -99.9, // V_WALL
    3.0,   // V_SAND
    9.0,   // V_WATER
    2.3,   // V_MUD
    1.5,   // V_GRASS
    0.1,   // V_TUNNEL
};

// Function to compare the score of 2 nodes
int fcmp_nodescore(const void *u, const void *v)
{
    const double a = ((struct node *)u)->score;
    const double b = ((struct node *)v)->score;
    return (a < b) ? -1 : (a > b);
}

node createNode(grid G, position p, node parent, heuristic h)
{
    double diagonal_len = (p.x != parent->pos.x && p.y != parent->pos.y) ? 0.01 : 0.0;
    node n = malloc(sizeof(struct node));
    n->pos = p;
    n->parent = parent;
    n->cost = parent->cost + weight[G.value[p.x][p.y]];
    n->score = n->cost + h(n->pos, G.end, &G) + diagonal_len;
    return n;
}

mpi_node CreateMpiNode(grid G, position p, mpi_node *parent, int parent_rank, int parent_win_i, heuristic h)
{
    mpi_node n;
    double diagonal_len = (p.x != parent->pos.x && p.y != parent->pos.y) ? 0.1 : 0.0;
    n.pos = p;
    n.cost = parent->cost + weight[G.value[p.x][p.y]];
    n.score = n.cost + h(p, G.end, &G) + diagonal_len;
    n.parent_rank = parent_rank;
    n.parent_win_i = parent_win_i;
    return n;
}

mpi_node *mpiNodeToPtr(mpi_node n)
{
    mpi_node *n_ptr = malloc(sizeof(mpi_node));
    n_ptr->pos = n.pos;
    n_ptr->cost = n.cost;
    n_ptr->score = n.score;
    n_ptr->parent_rank = n.parent_rank;
    n_ptr->parent_win_i = n.parent_win_i;
    return n_ptr;
}

// Return the rank number of the core that is going to process the node with position p
int hda(position p, int world_size)
{
    return (p.x + p.y) % world_size;
}

void CreateMpiPositionDataType(MPI_Datatype *position_dt)
{
    int lengths[2] = {1, 1};

    // Calculate displacements
    // In C, by default padding can be inserted between fields. MPI_Get_address will allow
    // to get the address of each struct field and calculate the corresponding displacement
    // relative to that struct base address. The displacements thus calculated will therefore
    // include padding if any.
    MPI_Aint displacements[2];
    position dummy_position;
    MPI_Aint base_address;
    MPI_Get_address(&dummy_position, &base_address);
    MPI_Get_address(&dummy_position.x, &displacements[0]);
    MPI_Get_address(&dummy_position.y, &displacements[1]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, lengths, displacements, types, position_dt);
    MPI_Type_commit(position_dt);
}

void CreateMpiNodeDataType(MPI_Datatype *mpi_node_dt, MPI_Datatype mpi_position_dt)
{
    int lengths[5] = {1, 1, 1, 1, 1};
    MPI_Aint displacements[5];
    mpi_node dummy_node;
    MPI_Aint base_address;
    MPI_Get_address(&dummy_node, &base_address);
    MPI_Get_address(&dummy_node.pos, &displacements[0]);
    MPI_Get_address(&dummy_node.cost, &displacements[1]);
    MPI_Get_address(&dummy_node.score, &displacements[2]);
    MPI_Get_address(&dummy_node.parent_rank, &displacements[3]);
    MPI_Get_address(&dummy_node.parent_win_i, &displacements[4]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    displacements[3] = MPI_Aint_diff(displacements[3], base_address);
    displacements[4] = MPI_Aint_diff(displacements[4], base_address);
    MPI_Datatype types[5] = {mpi_position_dt, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT};
    MPI_Type_create_struct(5, lengths, displacements, types, mpi_node_dt);
    MPI_Type_commit(mpi_node_dt);
}

double A_star_mpi(grid G, heuristic h)
{
    int rank, world_size;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Create the datatypes
    MPI_Datatype mpi_position_dt, mpi_node_dt;
    CreateMpiPositionDataType(&mpi_position_dt);
    CreateMpiNodeDataType(&mpi_node_dt, mpi_position_dt);

    int dim = G.X * G.Y;

    const int ARRAY_SIZE = dim; // TODO: Improve
    mpi_node window_buffer[ARRAY_SIZE];
    int cur_win_i = 0;

    // Create a heap with a capacity of the dimension of the graph
    heap Q = heap_create(INIT_HEAP_CAPACITY, fcmp_nodescore);

    // Verify if destination is a wall
    if (G.value[G.end.x][G.end.y] == V_WALL)
    {
        fprintf(stderr, "DESTINATION ON WALL\n");
        heap_destroy(Q);
        return -1;
    }

    int starting_process_rank = hda(G.start, world_size);
    int ending_process_rank = hda(G.end, world_size);

    // Set tags
    int destination_reached_tag = 1;
    int node_tag = 2;
    int path_construction_tag = 3;
    int path_done_tag = 4;

    // Get the node that will process the origin node
    if (rank == starting_process_rank)
    {
        // Init origin node
        mpi_node *s = malloc(sizeof(mpi_node));
        s->pos = G.start;
        s->cost = 0;
        s->score = s->cost + h(s->pos, G.end, &G);
        s->parent_rank = -1;
        s->parent_win_i = -1;
        if (heap_add(Q, s)) // add s to heap Q
        {
            fprintf(stderr, "Heap cannot expand anymore\n");
            heap_destroy(Q);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        G.mark[s->pos.x][s->pos.y] = M_FRONT;
    }

    while (true)
    {
        do
        {
            int flag_found;

            // Check if destination has been reached
            MPI_Iprobe(ending_process_rank, destination_reached_tag, MPI_COMM_WORLD, &flag_found, MPI_STATUS_IGNORE);
            if (flag_found)
            {
                MPI_Recv(NULL, 0, MPI_INT, ending_process_rank, destination_reached_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Loop until path has been fully constructed
                while (true)
                {
                    // Check path is done
                    int flag_path_done;
                    MPI_Iprobe(ending_process_rank, path_done_tag, MPI_COMM_WORLD, &flag_path_done, MPI_STATUS_IGNORE);
                    if (flag_path_done)
                    {
                        MPI_Recv(NULL, 0, MPI_INT, ending_process_rank, path_done_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        heap_destroy(Q);
                        return 1;
                    }

                    // Check if ending process require info about a node
                    int flag_path;
                    MPI_Iprobe(ending_process_rank, path_construction_tag, MPI_COMM_WORLD, &flag_path, MPI_STATUS_IGNORE);
                    if (flag_path)
                    {
                        int win_i;
                        MPI_Recv(&win_i, 1, MPI_INT, ending_process_rank, path_construction_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Send(&window_buffer[win_i], 1, mpi_node_dt, ending_process_rank, path_construction_tag, MPI_COMM_WORLD);
                    }
                }
            }

            MPI_Status status_node;
            int flag_node;

            // Check if we received a node
            MPI_Iprobe(MPI_ANY_SOURCE, node_tag, MPI_COMM_WORLD, &flag_node, &status_node);

            // We may have received multiple ones from different sources
            while (flag_node)
            {
                // Get number of nodes received
                int number_nodes_receiving;
                MPI_Get_count(&status_node, mpi_node_dt, &number_nodes_receiving);

                mpi_node nodes[number_nodes_receiving];
                MPI_Recv(nodes, number_nodes_receiving, mpi_node_dt, status_node.MPI_SOURCE, status_node.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Add nodes to heap
                for (int i = 0; i < number_nodes_receiving; i++)
                {
                    if (heap_add(Q, mpiNodeToPtr(nodes[i])))
                    {
                        fprintf(stderr, "Heap cannot expand anymore\n");
                        heap_destroy(Q);
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    }
                    G.mark[nodes[i].pos.x][nodes[i].pos.y] = M_FRONT;
                }

                // Check if we received a node again
                MPI_Iprobe(MPI_ANY_SOURCE, node_tag, MPI_COMM_WORLD, &flag_node, &status_node);
            }
        } while (heap_empty(Q));

        mpi_node *u = heap_pop(Q); // extract the node with minimum score

        // Node already visited ?
        // TODO: Check if cost lower than actual cost
        if (G.mark[u->pos.x][u->pos.y] == M_USED)
        {
            continue;
        }

        // Check if we are on the destination position
        if (rank == ending_process_rank && u->pos.x == G.end.x && u->pos.y == G.end.y)
        {
            // Broadcast that destination has been reached
            MPI_Request req[world_size - 1]; // Minus itself
            int cur_req = 0;
            for (int dst = 0; dst < world_size; dst++)
            {
                if (dst != rank)
                    MPI_Isend(NULL, 0, MPI_INT, dst, destination_reached_tag, MPI_COMM_WORLD, &req[cur_req++]);
            }
            MPI_Waitall(cur_req, req, MPI_STATUSES_IGNORE);

            // Construct the path
            mpi_node path = *u;
            while (path.pos.x != G.start.x && path.pos.y != G.start.y)
            {
                // Draw the path
                G.mark[path.pos.x][path.pos.y] = M_PATH;

                // Get info about parent node
                mpi_node parent;
                if (path.parent_rank != rank)
                {
                    int win_i = path.parent_win_i;
                    MPI_Send(&win_i, 1, MPI_INT, path.parent_rank, path_construction_tag, MPI_COMM_WORLD);
                    MPI_Recv(&parent, 1, mpi_node_dt, path.parent_rank, path_construction_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                else
                {
                    parent = window_buffer[path.parent_win_i];
                }

                path = parent;
            }

            // Broadcast that path has been constructed
            cur_req = 0;
            for (int dst = 0; dst < world_size; dst++)
            {
                if (dst != rank)
                    MPI_Isend(NULL, 0, MPI_INT, dst, path_done_tag, MPI_COMM_WORLD, &req[cur_req++]);
            }
            MPI_Waitall(cur_req, req, MPI_STATUSES_IGNORE);

            heap_destroy(Q);
            return u->cost;
        }

        // Add node to P
        G.mark[u->pos.x][u->pos.y] = M_USED;

        window_buffer[cur_win_i] = *u;

        // Create a 2D array where the nodes are going to be stored before getting sent
        mpi_node node_storage[world_size][MAX_NEIGHBORS];

        // Create an array that counts the nb of nodes to send to each process and initialize it with 0
        int nb_nodes_per_process[world_size];
        for (int i = 0; i < world_size; i++)
        {
            nb_nodes_per_process[i] = 0;
        }

        // For every neighbor of u
        for (int y = -1; y <= 1; y++)
        {
            for (int x = -1; x <= 1; x++)
            {
                // Calculate position of the node
                position p;
                p.x = u->pos.x + x;
                p.y = u->pos.y + y;

                // If node has not been visited yet or is not a wall
                if (G.mark[p.x][p.y] == M_NULL && G.value[p.x][p.y] != V_WALL)
                {
                    // Create and add node to tsend it to its destination process
                    int dst_process = hda(p, world_size);
                    mpi_node n = CreateMpiNode(G, p, u, rank, cur_win_i, h);

                    if (dst_process != rank)
                    {
                        node_storage[dst_process][nb_nodes_per_process[dst_process]++] = n;
                    }
                    else
                    {
                        if (heap_add(Q, mpiNodeToPtr(n)))
                        {
                            fprintf(stderr, "Heap overloaded\n");
                            heap_destroy(Q);
                            MPI_Abort(MPI_COMM_WORLD, 1);
                        }
                    }

                    G.mark[n.pos.x][n.pos.y] = M_FRONT; // -> Broadcast
                }
            }
        }

        MPI_Request req[MAX_NEIGHBORS]; // Maximum request = Max neighbors = 8
        int cur_req = 0;
        for (int dst = 0; dst < world_size; dst++)
        {
            int nb_nodes = nb_nodes_per_process[dst];
            if (nb_nodes > 0)
            {
                MPI_Isend(&(node_storage[dst][0]), nb_nodes, mpi_node_dt, dst, node_tag, MPI_COMM_WORLD, &req[cur_req++]);
            }
        }
        MPI_Waitall(cur_req, req, MPI_STATUSES_IGNORE);

        cur_win_i++;
    }

    // If path not found
    heap_destroy(Q);
    return -1;
}

double A_star_sequential(grid G, heuristic h)
{
    heap Q = heap_create(INIT_HEAP_CAPACITY, fcmp_nodescore);

    // Init destination node
    node t = malloc(sizeof(struct node));
    t->pos = G.end;

    // Verify if t is a wall
    if (G.value[t->pos.x][t->pos.y] == V_WALL)
    {
        fprintf(stderr, "DESTINATION ON WALL\n");
        heap_destroy(Q);
        return -1;
    }

    // Init origin node
    node s = malloc(sizeof(struct node));
    s->pos = G.start;
    s->parent = NULL;
    s->cost = 0;
    s->score = s->cost + h(s->pos, t->pos, &G);

    if (heap_add(Q, s)) // add s to heap Q
    {
        fprintf(stderr, "Heap cannot expand anymore\n");
        heap_destroy(Q);
        return -1;
    }
    G.mark[s->pos.x][s->pos.y] = M_FRONT;

    while (!heap_empty(Q))
    {                         // As long as there are nodes in Q
        node u = heap_pop(Q); // extract the node with minimum score

        // Noeud already visited ?
        if (G.mark[u->pos.x][u->pos.y] == M_USED)
            continue;

        // Check if we are on the destination position
        if (u->pos.x == t->pos.x && u->pos.y == t->pos.y)
        {
            t = u;
            node path = t;
            while (path != s)
            {
                // Draw the path
                G.mark[path->pos.x][path->pos.y] = M_PATH;
                path = path->parent;
            }
            heap_destroy(Q);
            return t->cost;
        }

        // Add node to P
        G.mark[u->pos.x][u->pos.y] = M_USED;

        // For every neighbor of u
        for (int y = -1; y <= 1; y++)
        {
            for (int x = -1; x <= 1; x++)
            {
                // Calculate position of the node
                position p;
                p.x = u->pos.x + x;
                p.y = u->pos.y + y;

                if (G.mark[p.x][p.y] == M_NULL && G.value[p.x][p.y] != V_WALL)
                {
                    // Create and add node to the heap Q
                    node v = createNode(G, p, u, h);
                    if (heap_add(Q, v))
                    {
                        printf("Heap cannot expand anymore\n");
                        heap_destroy(Q);
                        return -1;
                    };
                    G.mark[v->pos.x][v->pos.y] = M_FRONT; // -> Broadcast
                }
            }
        }
    }

    heap_destroy(Q);
    return -1;
}

int main(int argc, char *argv[])
{

    if (argc < 2 || argc > 6)
    {
        fprintf(stderr, "Number of arguments should be 5\n"
                        "Usage: ./a_star <seed> <grid width> <grid height> <grid type [empty|walls|maze])> "
                        "<algorithm [0 (Djikstra)|1 (AStar)|2 (Approx)]>");
        return 1;
    }

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set seed, if seed equal 0 then generate random one
    unsigned seed = atoi(argv[1]);
    seed = (seed == 0) ? time(NULL) % 1000 : seed;
    srandom(seed);

    // Set width and height
    int width = atoi(argv[2]);
    int height = atoi(argv[3]);

    // Set type
    char *type = argv[4];

    // Sel algorithm (Djikstra, A*, Approx)
    alpha = atoi(argv[5]);

    // Set Grid according to type provided
    grid G;
    if (strcmp(type, "empty") == 0)
    {
        G = initGridPoints(width, height, V_FREE, 1);
    }
    else if (strcmp(type, "walls") == 0)
    {
        G = initGridPoints(width, height, V_WALL, 0.2);
    }
    else if (strcmp(type, "maze") == 0)
    {
        int cw = 3;
        G = initGridLaby(width / (cw + 1), height / (cw + 1), cw);
    }
    else
    {
        fprintf(stderr, "Unknown type provided: %s\nTypes allowed: empty, walls, maze", type);
        return 1;
    }

    double (*f)(grid, heuristic);
    if (world_size > 1)
        f = A_star_mpi;
    else
        f = A_star_sequential;

    double d, start, delta;
    start = MPI_Wtime();
    d = f(G, halpha);
    delta = MPI_Wtime() - start;

    // path found or not?
    if (d < 0)
    {
        printf("path not found!\n");
        freeGrid(G);
        return 1;
    }

    int dst_process = hda(G.end, world_size);

    if (rank == dst_process)
    {
        // TODO: MPI_Gather on Grid to update it with all other processes Grid

        // counts the number of vertices explored to compare the heuristics
        // int m = 0;
        // for (int i = 0; i < G.X; i++)
        //     for (int j = 0; j < G.Y; j++)
        //         m += (G.mark[i][j] != M_NULL);
        // printf("#nodes explored: %i\n", m);

        printf("Nb_cores: %d\nDimensions: %d\nBingo! Path found.. Cost: %g\tPerf: %lgs\n", world_size, width, d, delta);
    }

    freeGrid(G);
    MPI_Finalize();
    return 0;
}