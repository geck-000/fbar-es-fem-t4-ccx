/*     F-barES-FEM-T4: a pthread mutex for the F-barES-FEM-T4 node->element map.
 *
 *     ccx parallelises with pthreads and compiles the Fortran WITHOUT
 *     -fopenmp (see FFLAGS), so every !$omp directive in a .f file is an
 *     inert comment.  The one-time node->element map that U2 and U3 build
 *     over the base U4 tets must be guarded by a real mutex, or every
 *     assembly thread fills it at once.  A pthread mutex is what this
 *     threading model actually needs, and it carries the memory barrier
 *     that makes the finished map visible to every other thread.
 */

#include <pthread.h>

static pthread_mutex_t fbarmutex = PTHREAD_MUTEX_INITIALIZER;

void fbarlock_(void)   { pthread_mutex_lock(&fbarmutex); }
void fbarunlock_(void) { pthread_mutex_unlock(&fbarmutex); }
