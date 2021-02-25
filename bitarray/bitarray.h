/*
   Copyright (c) 2008 - 2021, Ilan Schnell; All Rights Reserved
   bitarray is published under the PSF license.

   Author: Ilan Schnell
*/
#define BITARRAY_VERSION  "1.7.0"

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#if PY_MAJOR_VERSION == 2
/* the Py_MIN macro was introduced in Python 3.3 */
#define Py_MIN(x, y)  (((x) > (y)) ? (y) : (x))
#define PySlice_GetIndicesEx(slice, len, start, stop, step, slicelength) \
    PySlice_GetIndicesEx(((PySliceObject *) slice),                      \
                         (len), (start), (stop), (step), (slicelength))
#define PyLong_FromLong  PyInt_FromLong
#endif

#include <vector>
#include "_allocator.hpp"
#include "bit.hpp"

typedef long long int idx_t;
typedef unsigned long long int word_type;

/* ob_size is the byte count of the buffer, not the number of elements.
   The number of elements (bits) is nbits. */
typedef struct {
    PyObject_VAR_HEAD
    std::vector<word_type, PyAlloc<word_type> >* ob_item;
    bit::bit_iterator<decltype(std::begin(*words))> bits;
    Py_ssize_t allocated;       /* how many bytes allocated */
    Py_ssize_t nbits;           /* length of bitarray, i.e. elements */
    int endian;                 /* bit endianness of bitarray */
    int ob_exports;             /* how many buffer exports */
    PyObject *weakreflist;      /* list of weak references */
} bitarrayobject;

/* --- bit endianness --- */
#define ENDIAN_LITTLE  0
#define ENDIAN_BIG     1

/* returns the endianness string from bitarrayobject */
#define ENDIAN_STR(o)  ((o)->endian == ENDIAN_LITTLE ? "little" : "big")

/* number of bits that can be stored in given bytes */
#define BITS(bytes)  ((bytes) << 3)

/* number of bytes necessary to store given bits */
#define BYTES(bits)  ((bits) == 0 ? 0 : (((bits) - 1) / 8 + 1))

#define BITMASK(endian, i)  \
    (((char) 1) << ((endian) == ENDIAN_LITTLE ? ((i) % 8) : (7 - (i) % 8)))

/* ------------ low level access to bits in bitarrayobject ------------- */

/* sets unused padding bits (within last byte of buffer) to 0,
   and return the number of padding bits -- self->nbits is unchanged */
static inline int
setunused(bitarrayobject* self)
{
    const idx_t n = BITS(Py_SIZE(self));
    idx_t i;
    int res = 0;

    for (i = self->nbits; i < n; i++) {
        self->bits[i] = bit::bit0;
        res++;
    }
    assert(res < 8);
    return res;
}

static unsigned char bitcount_lookup[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
};
