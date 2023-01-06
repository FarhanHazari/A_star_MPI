#include "tools.h"

// Returns true if (i,j) is on the border of grid G.
static inline int onBorder(grid *G, int i, int j)
{
    return (i == 0) || (j == 0) || (i == G->X - 1) || (j == G->Y - 1);
}

// Allocates a grid with dimensions x,y as well as its image.
// We force x,y>=3 to have at least one point that is not on the border.
static grid allocGrid(int x, int y)
{
    grid G;
    position p = {-1, -1};
    G.start = G.end = p;
    if (x < 3)
        x = 3;
    if (y < 3)
        y = 3;
    G.X = x;
    G.Y = y;
    G.value = malloc(x * sizeof(*(G.value)));
    G.mark = malloc(x * sizeof(*(G.mark)));

    for (int i = 0; i < x; i++)
    {
        G.value[i] = malloc(y * sizeof(*(G.value[i])));
        G.mark[i] = malloc(y * sizeof(*(G.mark[i])));
        for (int j = 0; j < y; j++)
            G.mark[i][j] = M_NULL; // initialise
    }

    return G;
}

// Returns a random position on the grid that is uniform among all values of the grid of type t (excluding the borders of the grid).
// If no type t cases are found, the position {-1,-1} is returned.
position randomPosition(grid G, int t)
{
    int i, j, c;
    int n;                      // n: number of type t cases off the border
    int r = -1;                 // r: random number in [0,n[
    position p = {-1, -1};      // p: default position
    const int stop = G.X * G.Y; // stop: to exit the loops
    const int x1 = G.X - 1;
    const int y1 = G.Y - 1;

    // We do two passes: a 1st to count the number n of type t cases,
    // and a 2nd to randomly select the position among them.
    // At the end of the first pass, we know the number n of type t cases.
    // We then randomly draw a number r in [0,n[.
    // Then we start counting again (n=0) of type t cases and stop as soon as we reach case number r.

    c = 0;
    do
    {
        n = 0;
        for (i = 1; i < x1; i++)
            for (j = 1; j < y1; j++)
                if (G.value[i][j] == t)
                {
                    if (n == r)
                    {
                        p = (position){i, j};
                        i = j = stop; // always false on the first pass
                    }
                    n++;
                }
        c = 1 - c;
        if (c)
            r = random() % n;
    } while (c); // True 1st time, False the 2nd

    return p;
}

// Frees the pointers allocated by allocGrid().
void freeGrid(grid G)
{
    for (int i = 0; i < G.X; i++)
    {
        free(G.value[i]);
        free(G.mark[i]);
    }
    free(G.value);
    free(G.mark);
}

// Returns a grid of dimensions x,y initialized with random values.
grid initGridPoints(int x, int y, int type, double density)
{
    grid G = allocGrid(x, y); // Allocates the grid and its image

    // Verify correct type, default: M_NULL
    if ((type < 0))
        type = M_NULL;

    // Put the borders and fills the inside
    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++)
            G.value[i][j] =
                onBorder(&G, i, j) ? V_WALL : ((RAND01 <= density) ? type : V_FREE);

    // Random position start/end
    // G.start = randomPosition(G, V_FREE);
    // G.end = randomPosition(G, V_FREE);

    // Default position
    G.start = (position){.x = G.X - 2, .y = G.Y - 2};
    G.end = (position){.x = 1, .y = 1};

    return G;
}

// Returns a random grid of dimensions x,y (at least 3) corresponding to a random uniform spanning tree.
// We fix the start point = bottom right and end = top left. The width of the corridors is given by w>0.
// This is the Wilson algorithm by "random walks with loop erasure" (see https://bl.ocks.org/mbostock/11357811).
grid initGridLaby(int x, int y, int w)
{

    // Verify parameters
    if (x < 3)
        x = 3;
    if (y < 3)
        y = 3;
    if (w <= 0)
        w = 1;

    // Allocates the grid and its image
    int *value = malloc(x * y * sizeof(*value));

    // Allocates the grid and its image
    grid Gw = allocGrid(x * (w + 1) + 1, y * (w + 1) + 1);

    // Default position
    Gw.start = (position){.x = Gw.X - 2, .y = Gw.Y - 2};
    Gw.end = (position){.x = 1, .y = 1};

    // Initially walls only on the borders
    for (int i = 0; i < Gw.X; i++)
    {
        for (int j = 0; j < Gw.Y; j++)
        {
            Gw.value[i][j] =
                ((i % (w + 1) == 0) || (j % (w + 1) == 0)) ? V_WALL : V_FREE;
        }
    }

    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++)
            value[i * y + j] = -1;

    int count = 1;
    value[0] = 0;
    while (count < x * y)
    {
        int i0 = 0;
        while (i0 < x * y && value[i0] != -1)
            i0++;
        value[i0] = i0 + 1;
        while (i0 < x * y)
        {
            int x0 = i0 / y;
            int y0 = i0 % y;
            while (true)
            {
                int dir = random() & 3; // same as random()%4
                switch (dir)
                {
                case 0:
                    if (x0 <= 0)
                        continue;
                    x0--;
                    break;
                case 1:
                    if (y0 <= 0)
                        continue;
                    y0--;
                    break;
                case 2:
                    if (x0 >= x - 1)
                        continue;
                    x0++;
                    break;
                case 3:
                    if (y0 >= y - 1)
                        continue;
                    y0++;
                    break;
                }
                break;
            }
            if (value[x0 * y + y0] == -1)
            {
                value[x0 * y + y0] = i0 + 1;
                i0 = x0 * y + y0;
            }
            else
            {
                if (value[x0 * y + y0] > 0)
                {
                    while (i0 != x0 * y + y0 && i0 > 0)
                    {
                        int i1 = value[i0] - 1;
                        value[i0] = -1;
                        i0 = i1;
                    }
                }
                else
                {
                    int i1 = i0;
                    i0 = x0 * y + y0;
                    do
                    {
                        int x0 = i0 / y;
                        int y0 = i0 % y;
                        int x1 = i1 / y;
                        int y1 = i1 % y;
                        if (x0 < x1)
                            for (int i = 0; i < w; ++i)
                                Gw.value[x1 * (w + 1)][y0 * (w + 1) + i + 1] = V_FREE;
                        if (x0 > x1)
                            for (int i = 0; i < w; ++i)
                                Gw.value[x0 * (w + 1)][y0 * (w + 1) + i + 1] = V_FREE;
                        if (y0 < y1)
                            for (int i = 0; i < w; ++i)
                                Gw.value[x1 * (w + 1) + i + 1][y1 * (w + 1)] = V_FREE;
                        if (y0 > y1)
                            for (int i = 0; i < w; ++i)
                                Gw.value[x1 * (w + 1) + i + 1][y0 * (w + 1)] = V_FREE;
                        i0 = i1;
                        i1 = value[i0] - 1;
                        value[i0] = 0;
                        count++;
                    } while (value[i1] != 0);
                    break;
                }
            }
        }
    }

    free(value);

    return Gw;
}

void saveGridValueFile(grid G, char *filename)
{
    // MPI_File f;
    // if (MPI_File_open(MPI_COMM_SELF, "test-val.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f) != MPI_SUCCESS)
    // {
    //     printf("Failure in opening the file.\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // };
    char c;
    char new_line = '\n';
    fprintf(stderr, "#GRID_VALUE\n\n");
    for (int y = 0; y < G.Y; y++)
    {
        for (int x = 0; x < G.X; x++)
        {
            int v = G.value[x][y];
            switch (v)
            {
            case V_FREE:
                c = ' ';
                break;
            case V_WALL:
                c = '#';
                break;
            case V_SAND:
                c = ';';
                break;
            case V_WATER:
                c = '~';
                break;
            case V_MUD:
                c = ',';
                break;
            case V_GRASS:
                c = '.';
                break;
            case V_TUNNEL:
                c = '+';
                break;
            default:
                c = ' ';
            }
            if (x == G.start.x && y == G.start.y)
                c = 's';
            else if (x == G.end.x && y == G.end.y)
                c = 't';
            fprintf(stderr, "%c", c);
            // MPI_File_write(f, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE);
        }
        fprintf(stderr, "%c", new_line);
        // MPI_File_write(f, &new_line, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    fprintf(stderr, "\n\n");
    // if (MPI_File_close(&f) != MPI_SUCCESS)
    // {
    //     printf("Failure in closing the file.\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // }
    // debug(0, "I am here 3-1\n");
}

void saveGridMarkFile(grid G, char *filename)
{
    // MPI_File f;
    // if (MPI_File_open(MPI_COMM_SELF, "test_mark.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f) != MPI_SUCCESS)
    // {
    //     printf("Failure in opening the file.\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // };
    char c;
    char new_line = '\n';
    fprintf(stderr, "#GRID_MARK\n\n");
    for (int y = 0; y < G.Y; y++)
    {
        for (int x = 0; x < G.X; x++)
        {
            int m = G.mark[x][y];
            switch (m)
            {
            case M_NULL:
                c = ' ';
                break;
            case M_PATH:
                c = 'p';
                break;
            case M_USED:
                c = 'u';
                break;
            case M_FRONT:
                c = 'f';
                break;
            default:
                c = ' ';
            }
            fprintf(stderr, "%c", c);
            // MPI_File_write(f, &c, 1, MPI_CHAR, MPI_STATUS_IGNORE);
        }
        fprintf(stderr, "%c", new_line);
        // MPI_File_write(f, &new_line, 1, MPI_CHAR, MPI_STATUS_IGNORE);
    }
    fprintf(stderr, "\n\n");
    // if (MPI_File_close(&f) != MPI_SUCCESS)
    // {
    //     printf("Failure in closing the file.\n");
    //     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    // }
    // debug(0, "I am here 3-2\n");
}

void debug(int rank, char *format, ...)
{
    va_list args; // Variable argument list

    // Initialize variable argument list;
    // format is last argument before va_start (args, format)
    va_start(args, format);

    printf("%d:\t", rank);
    vprintf(format, args);
    //fflush(stdin);
    va_end(args);
}