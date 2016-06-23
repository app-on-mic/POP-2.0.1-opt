#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>

/*int backtrace (void **buffer, int size);*/
/*char ** backtrace_symbols (void *const *buffer, int size);*/
/*void backtrace_symbols_fd (void *const *buffer, int size, int fd);*/

/* Obtain a backtrace and print it to stdout. */
void
print_trace (void)
{
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace (array, 10);
  strings = backtrace_symbols (array, size);

  printf ("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
     printf ("%s\n", strings[i]);

  free (strings);
}

/* A dummy function to make the backtrace more interesting. */
void
dummy_function (void)
{
  print_trace ();
}

/*int
main (void)
{
  dummy_function ();
  return 0;
}*/
