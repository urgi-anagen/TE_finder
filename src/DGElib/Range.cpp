#include <Range.h>

void *Range::clone(void) const
{
  return (void*) new Range(this);
};



