/*
   Copyright (c) 2008 - 2019, Ilan Schnell
   bitarray is published under the PSF license.

   This file is the C part of the bitarray package.
   All functionality is implemented here.

   Author: Ilan Schnell
*/

#include "reverse.hpp"
#include <memory>
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

#ifdef IS_PY3K
#define Py_TPFLAGS_HAVE_WEAKREFS  0
#endif

#if PY_MAJOR_VERSION == 3 || (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION == 7)
/* (new) buffer protocol */
#define WITH_BUFFER
#endif

#ifdef STDC_HEADERS
#include <stddef.h>
#else  /* !STDC_HEADERS */
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>      /* For size_t */
#endif /* HAVE_SYS_TYPES_H */
#endif /* !STDC_HEADERS */

#include <vector>
#include "bit.hpp"

typedef long long int idx_t;

/* throughout:  0 = little endian   1 = big endian */
#define DEFAULT_ENDIAN  1

/* Unlike the normal convention, ob_size is the byte count, not the number
   of elements.  The reason for doing this is that we can use our own
   special idx_t for the number of bits (which can exceed 2^32 on a 32 bit
   machine.  */
struct bitarrayobject {
    PyObject_VAR_HEAD
    std::vector<unsigned char>* words;
    bit::bit_iterator<decltype(std::begin(*words))> bits;
    long long nbits;   /* length of bitarray, i.e. elements */
    int endian;                 /* bit endianness of bitarray */
    int ob_exports;             /* how many buffer exports */
    PyObject *weakreflist;      /* list of weak references */
};

extern PyTypeObject Bitarraytype;
#define bitarray_Check(obj)  PyObject_TypeCheck(obj, &Bitarraytype)

/* number of bytes necessary to store given bits */
#define BYTES(bits)  (((bits) == 0) ? 0 : (((bits) - 1) / 8 + 1))
#define BITS(bytes)  ((idx_t) (bytes) << 3)


/* ------------ low level access to bits in bitarrayobject ------------- */
static int
check_overflow(idx_t nbits)
{
    assert(nbits >= 0);
    if (sizeof(void *) == 4) {  /* 32bit system */
        const idx_t max_bits = ((idx_t) 1) << 34;  /* 2^34 = 16 Gbits*/
        if (nbits > max_bits) {
            char buff[256];
            sprintf(buff, "cannot create bitarray of size %lld, "
                          "max size is %lld", nbits, max_bits);
            PyErr_SetString(PyExc_OverflowError, buff);
            return -1;
        }
    }
    return 0;
}

static int
resize(bitarrayobject* self, idx_t nbits)
{
    long long old_nbits = self->nbits;
    self->words->resize(BYTES(nbits));
    if (nbits < old_nbits) {
        self->words->shrink_to_fit();
    }
    self->nbits = nbits;
    Py_SIZE(self) = self->words->size();
    self->bits = bit::bit_iterator(std::begin(*(self->words)));
    return 0;
}

/* create new bitarray object without initialization of buffer */
static PyObject *
newbitarrayobject(PyTypeObject* type, idx_t nbits, int endian)
{
    bitarrayobject* obj;
    Py_ssize_t nbytes;

    if (check_overflow(nbits) < 0)
        return NULL;

    obj = (bitarrayobject *) type->tp_alloc(type, 0);
    if (obj == NULL)
        return NULL;

    nbytes = (Py_ssize_t) BYTES(nbits);
    Py_SIZE(obj) = nbytes;
    obj->nbits = nbits;
    obj->endian = endian;
    if (nbytes == 0) {
        obj->words = new std::vector<unsigned char>(0);
    }
    else {
        obj->words = new std::vector<unsigned char>(nbytes, 0);
        if (obj->words == NULL) {
            PyObject_Del(obj);
            PyErr_NoMemory();
            return NULL;
        }
    }
    obj->weakreflist = NULL;
    obj->bits = bit::bit_iterator(std::begin(*(obj->words)));
    return (PyObject *) obj;
}

static void
bitarray_dealloc(bitarrayobject* self)
{
    if (self->weakreflist != NULL)
        PyObject_ClearWeakRefs((PyObject *) self);

    if (self->words != NULL)
        delete self->words;

    Py_TYPE(self)->tp_free((PyObject *) self);
}

enum op_type {
    OP_and,
    OP_or,
    OP_xor,
};


// TODO handle endianess
/* copy n bits from other (starting at b) to self (starting at a) */
static void
copy_n(bitarrayobject* self, idx_t a,
       bitarrayobject* other, idx_t b, idx_t n)
{
    assert(0 <= n && n <= self->nbits && n <= other->nbits);
    assert(0 <= a && a <= self->nbits - n);
    assert(0 <= b && b <= other->nbits - n);
    if (n == 0)
        return;

    // TODO has to be a better way to do this...
    if (&*self->words <= &*other->words && &(*other->words)[0] < &(*((*self->words).end())))
        bit::copy_backward(
            bit::next(self->bits, b), 
            bit::next(self->bits, b + n), 
            bit::next(other->bits, a + n)
        );
    else {
        bit::copy(
            bit::next(self->bits, b), 
            bit::next(self->bits, b + n), 
            bit::next(other->bits, a)
        );
    }
    return;
}

// TODO replace with true bit vector
/* starting at start, delete n bits from self */
static int
delete_n(bitarrayobject* self, idx_t start, idx_t n)
{
    assert(0 <= start && start <= self->nbits);
    assert(0 <= n && n <= self->nbits - start);
    if (n == 0)
        return 0;

    copy_n(self, start + n, self, start, self->nbits - start - n);
    return resize(self, self->nbits - n);
}

/* starting at start, insert n (uninitialized) bits into self */
static int
insert_n(bitarrayobject* self, idx_t start, idx_t n)
{
    assert(0 <= start && start <= self->nbits);
    assert(n >= 0);
    if (n == 0)
        return 0;

    self->nbits += 1;
    if (resize(self, self->nbits) < 0) {
        return -1;
    }
    copy_n(self, start + n, self, start, self->nbits - start - n);
    return 0;
}

/* repeat self n times */
static int
repeat(bitarrayobject* self, idx_t n)
{
    idx_t nbits, i;

    if (n < 0) {
        return -1;
    }
    if (n > 1) {
        nbits = self->nbits;
        resize(self, n * self->nbits);
        for (i = 1; i < n; i++)
            copy_n(self, i * nbits, self, 0, nbits);
    }
    return 0;
}

/* sets ususet bits to 0, i.e. the ones in the last byte (if any),
   and return the number of bits set -- self->nbits is unchanged */
// TODO this could be optimized
static int
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

/* perform bitwise in-place operation */
static int
bitwise(bitarrayobject *self, PyObject *arg, enum op_type oper)
{
    bitarrayobject *other;
    Py_ssize_t i;

    if (!bitarray_Check(arg)) {
        PyErr_SetString(PyExc_TypeError,
                        "bitarray object expected for bitwise operation");
        return -1;
    }
    other = (bitarrayobject *) arg;
    if (self->nbits != other->nbits) {
        PyErr_SetString(PyExc_ValueError,
               "bitarrays of equal length expected for bitwise operation");
        return -1;
    }
    //setunused(self);
    //setunused(other);
    switch (oper) {
    case OP_and:
        for (i = 0; i < Py_SIZE(self); i++)
            (*self->words)[i] &= (*other->words)[i];
        break;
    case OP_or:
        for (i = 0; i < Py_SIZE(self); i++)
            (*self->words)[i] |= (*other->words)[i];
        break;
    case OP_xor:
        for (i = 0; i < Py_SIZE(self); i++)
            (*self->words)[i] ^= (*other->words)[i];
        break;
    default:  /* should never happen */
        return -1;
    }
    return 0;
}

/* set the bits from start to stop (excluding) in self to val */
static void
setrange(bitarrayobject* self, idx_t start, idx_t stop, int val)
{
    assert(0 <= start && start <= self->nbits);
    assert(0 <= stop && stop <= self->nbits);

    if (self->nbits == 0 || start >= stop)
        return;

    bit::fill(self->bits, self->bits + self->nbits, val ? bit::bit1 : bit::bit0);
    return;
}

/* returns number of 1 bits */
static idx_t
count(bitarrayobject* self, int val, idx_t start, idx_t stop)
{
    assert(0 <= start && start <= self->nbits);
    assert(0 <= stop && stop <= self->nbits);
    assert(0 <= val && val <= 1);
    assert(BYTES(stop) <= Py_SIZE(self));

    return bit::count(
        bit::next(self->bits, start), 
        bit::next(self->bits, stop), 
        val ? bit::bit1 : bit::bit0
    );
}

/* return index of first occurrence of vi, -1 when x is not in found. */
static idx_t
findfirst(bitarrayobject* self, int val, idx_t start, idx_t stop)
{
    assert(0 <= start && start <= self->nbits);
    assert(0 <= stop && stop <= self->nbits);
    assert(0 <= val && val <= 1);
    assert(BYTES(stop) <= Py_SIZE(self));

    if (self->nbits == 0 || start >= stop)
        return -1;

    auto ret_iter = bit::find(
        bit::next(self->bits, start),
        bit::next(self->bits, stop),
        val ? bit::bit1 : bit::bit0
    );

    return ret_iter != bit::next(self->bits, self->nbits) ? bit::distance(self->bits, ret_iter) : -1;
}

/* search for the first occurrence of bit vector pattern (in self), starting at start,
   and return its index (or -1 when not found)
*/
static idx_t
search(bitarrayobject* self, bitarrayobject* pattern, idx_t start)
{
    auto ret_iter = bit::search(
        bit::next(self->bits, start), bit::next(self->bits, self->nbits),
        pattern->bits, bit::next(pattern->bits, pattern->nbits)
    );
    return ret_iter != bit::next(self->bits, self->nbits) ? bit::distance(self->bits, ret_iter) : -1;
}

// TODO replace with bit:: function
static void
invert(bitarrayobject* self)
{
    Py_ssize_t i;
    for (i = 0; i < Py_SIZE(self); i++)
        (*self->words)[i] = ~(*self->words)[i];
}

static int
set_item(bitarrayobject* self, idx_t i, PyObject *v)
{
    long vi;

    assert(0 <= i && i < self->nbits);
    vi = PyObject_IsTrue(v);
    if (vi < 0)
        return -1;
    self->bits[i] = vi ? bit::bit1 : bit::bit0;
    return 0;
}

// TODO replace with true bit vector
// TODO catch error
static int
append_item(bitarrayobject* self, PyObject *item)
{
    if (self->words->size() == self->nbits) {
        self->words->push_back(static_cast<char>(0));
    }
    self->nbits += 1;
    return set_item(self, self->nbits - 1, item);
}

static PyObject *
unpack(bitarrayobject* self, char zero, char one)
{
    PyObject *result;
    Py_ssize_t i;
    char *str;

    if (self->nbits > PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "bitarray too large to unpack");
        return NULL;
    }
    str = (char *) PyMem_Malloc((size_t) self->nbits);
    if (str == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    for (i = 0; i < self->nbits; i++) {
        *(str + i) = self->bits[i] ? one : zero;
    }
    result = PyBytes_FromStringAndSize(str, (Py_ssize_t) self->nbits);
    PyMem_Free((void *) str);
    return result;
}

static int
extend_bitarray(bitarrayobject* self, bitarrayobject* other)
{
    idx_t n_sum;
    idx_t n_other_bits;

    if (other->nbits == 0)
        return 0;

    /* Note that other may be self.  Thus we take the size before we resize,
       ensuring we only copy the right parts of the array. */
    n_other_bits = other->nbits;
    n_sum = self->nbits + other->nbits;
    if (resize(self, n_sum) < 0)
        return -1;
    self->nbits = n_sum;

    copy_n(self, n_sum - n_other_bits, other, 0, n_other_bits);
    return 0;
}

// TODO need to understand more about this, can probably improve
static int
extend_iter(bitarrayobject* self, PyObject *iter)
{
    PyObject *item;

    assert(PyIter_Check(iter));
    while ((item = PyIter_Next(iter)) != NULL) {
        if (append_item(self, item) < 0) {
            Py_DECREF(item);
            return -1;
        }
        Py_DECREF(item);
    }
    if (PyErr_Occurred())
        return -1;

    return 0;
}

static int
extend_list(bitarrayobject* self, PyObject *list)
{
    PyObject *item;
    Py_ssize_t n, i;

    assert(PyList_Check(list));
    n = PyList_Size(list);
    if (n == 0)
        return 0;

    resize(self, self->nbits + n);

    for (i = 0; i < n; i++) {
        item = PyList_GetItem(list, i);
        if (item == NULL)
            return -1;
        if (set_item(self, self->nbits - n + i, item) < 0)
            return -1;
    }
    return 0;
}

static int
extend_tuple(bitarrayobject* self, PyObject *tuple)
{
    PyObject *item;
    Py_ssize_t n, i;

    assert(PyTuple_Check(tuple));
    n = PyTuple_Size(tuple);
    if (n == 0)
        return 0;

    resize(self, self->nbits + n);

    for (i = 0; i < n; i++) {
        item = PyTuple_GetItem(tuple, i);
        if (item == NULL)
            return -1;
        if (set_item(self, self->nbits - n + i, item) < 0)
            return -1;
    }
    return 0;
}

/* extend_bytes(): extend the bitarray from a PyBytes object, where each
   whole character is converted to a single bit */
enum conv_t {
    BYTES_01,    /*  '0' -> 0    '1'  -> 1   no other characters allowed */
    BYTES_RAW,   /*  0x00 -> 0   other -> 1                              */
};

static int
extend_bytes(bitarrayobject* self, PyObject *bytes, enum conv_t conv)
{
    Py_ssize_t strlen, i;
    char c, *str;
    auto bv = bit::bit0;

    assert(PyBytes_Check(bytes));
    strlen = PyBytes_Size(bytes);
    if (strlen == 0)
        return 0;

    if (resize(self, self->nbits + strlen) < 0)
        return -1;

    str = PyBytes_AsString(bytes);

    for (i = 0; i < strlen; i++) {
        c = *(str + i);
        /* depending on conv, map c to bit */
        switch (conv) {
        case BYTES_01:
            switch (c) {
            case '0': bv = bit::bit0; break;
            case '1': bv = bit::bit1; break;
            default:
                PyErr_Format(PyExc_ValueError,
                             "character must be '0' or '1', found '%c'", c);
                return -1;
            }
            break;
        case BYTES_RAW:
            bv = c ? bit::bit1 : bit::bit0;
            break;
        default:  /* should never happen */
            return -1;
        }
        self->bits[self->nbits - strlen + i] =  bv;
    }
    return 0;
}

static int
extend_rawbytes(bitarrayobject* self, PyObject *bytes)
{
    Py_ssize_t strlen;
    char *str;

    assert(PyBytes_Check(bytes) && self->nbits % 8 == 0);
    strlen = PyBytes_Size(bytes);
    if (strlen == 0)
        return 0;

    if (resize(self, self->nbits + BITS(strlen)) < 0)
        return -1;

    str = PyBytes_AsString(bytes);
    auto str_bit_it = bit::bit_iterator(reinterpret_cast<unsigned char*>(str));
    bit::copy(str_bit_it, bit::next(str_bit_it, BITS(strlen)), self->bits);
    return 0;
}

static int
extend_dispatch(bitarrayobject* self, PyObject *obj)
{
    PyObject *iter;
    int ret;

    /* dispatch on type */
    if (bitarray_Check(obj))                              /* bitarray */
        return extend_bitarray(self, (bitarrayobject* ) obj);

    if (PyList_Check(obj))                                    /* list */
        return extend_list(self, obj);

    if (PyTuple_Check(obj))                                  /* tuple */
        return extend_tuple(self, obj);

    if (PyBytes_Check(obj))                               /* bytes 01 */
        return extend_bytes(self, obj, BYTES_01);

    if (PyUnicode_Check(obj)) {                /* (unicode) string 01 */
        PyObject *bytes;
        bytes = PyUnicode_AsEncodedString(obj, NULL, NULL);
        ret = extend_bytes(self, bytes, BYTES_01);
        Py_DECREF(bytes);
        return ret;
    }

    if (PyIter_Check(obj))                                    /* iter */
        return extend_iter(self, obj);

    /* finally, try to get the iterator of the object */
    iter = PyObject_GetIter(obj);
    if (iter == NULL) {
        PyErr_SetString(PyExc_TypeError, "could not extend bitarray");
        return -1;
    }
    ret = extend_iter(self, iter);
    Py_DECREF(iter);
    return ret;
}

/* --------- helper functions NOT involving bitvector objects ------------ */

#define ENDIAN_STR(ba)  (((ba)->endian) ? "big" : "little")

#ifdef IS_PY3K
#define IS_INDEX(x)  (PyLong_Check(x) || PyIndex_Check(x))
#define IS_INT_OR_BOOL(x)  (PyBool_Check(x) || PyLong_Check(x))
#else  /* Py 2 */
#define IS_INDEX(x)  (PyInt_Check(x) || PyLong_Check(x) || PyIndex_Check(x))
#define IS_INT_OR_BOOL(x)  (PyBool_Check(x) || PyInt_Check(x) || \
                                               PyLong_Check(x))
#endif

/* given an PyLong (which must be 0 or 1), or a PyBool, return 0 or 1,
   or -1 on error */
static int
IntBool_AsInt(PyObject *v)
{
    long x;

    if (PyBool_Check(v))
        return PyObject_IsTrue(v);

#ifndef IS_PY3K
    if (PyInt_Check(v)) {
        x = PyInt_AsLong(v);
    }
    else
#endif
    if (PyLong_Check(v)) {
        x = PyLong_AsLong(v);
    }
    else {
        PyErr_SetString(PyExc_TypeError, "integer or bool expected");
        return -1;
    }

    if (x < 0 || x > 1) {
        PyErr_SetString(PyExc_ValueError,
                        "integer value between 0 and 1 expected");
        return -1;
    }
    return (int) x;
}

/* Normalize index (which may be negative), such that 0 <= i <= n */
static void
normalize_index(idx_t n, idx_t *i)
{
    if (*i < 0) {
        *i += n;
        if (*i < 0)
            *i = 0;
    }
    if (*i > n)
        *i = n;
}

/* Extract a slice index from a PyInt or PyLong or an object with the
   nb_index slot defined, and store in *i.
   However, this function returns -1 on error and 0 on success.
   This is almost _PyEval_SliceIndex() with Py_ssize_t replaced by idx_t
*/
static int
getIndex(PyObject *v, idx_t *i)
{
    idx_t x;

#ifndef IS_PY3K
    if (PyInt_Check(v)) {
        x = PyInt_AS_LONG(v);
    }
    else
#endif
    if (PyLong_Check(v)) {
        x = PyLong_AsLongLong(v);
    }
    else if (PyIndex_Check(v)) {
        x = PyNumber_AsSsize_t(v, NULL);
        if (x == -1 && PyErr_Occurred())
            return -1;
    }
    else {
        PyErr_SetString(PyExc_TypeError, "slice indices must be integers or "
                                         "None or have an __index__ method");
        return -1;
    }
    *i = x;
    return 0;
}

/* this is PySlice_GetIndicesEx() with Py_ssize_t replaced by idx_t */
static int
slice_GetIndicesEx(PySliceObject *r, idx_t length,
                   idx_t *start, idx_t *stop, idx_t *step, idx_t *slicelength)
{
    idx_t defstart, defstop;

    if (r->step == Py_None) {
        *step = 1;
    }
    else {
        if (getIndex(r->step, step) < 0)
            return -1;
        if (*step == 0) {
            PyErr_SetString(PyExc_ValueError, "slice step cannot be zero");
            return -1;
        }
    }
    defstart = *step < 0 ? length - 1 : 0;
    defstop = *step < 0 ? -1 : length;

    if (r->start == Py_None) {
        *start = defstart;
    }
    else {
        if (getIndex(r->start, start) < 0)
            return -1;
        if (*start < 0) *start += length;
        if (*start < 0) *start = (*step < 0) ? -1 : 0;
        if (*start >= length) *start = (*step < 0) ? length - 1 : length;
    }

    if (r->stop == Py_None) {
        *stop = defstop;
    }
    else {
        if (getIndex(r->stop, stop) < 0)
            return -1;
        if (*stop < 0) *stop += length;
        if (*stop < 0) *stop = -1;
        if (*stop > length) *stop = length;
    }

    if ((*step < 0 && *stop >= *start) || (*step > 0 && *start >= *stop)) {
        *slicelength = 0;
    }
    else if (*step < 0) {
        *slicelength = (*stop - *start + 1) / (*step) + 1;
    }
    else {
        *slicelength = (*stop - *start - 1) / (*step) + 1;
    }

    return 0;
}

/**************************************************************************
                         Implementation of API methods
 **************************************************************************/

static PyObject *
bitarray_length(bitarrayobject* self)
{
    return PyLong_FromLongLong(self->nbits);
}

PyDoc_STRVAR(length_doc,
"length() -> int\n\
\n\
Return the length, i.e. number of bits stored in the bitarray.\n\
This method is preferred over `__len__` (used when typing `len(a)`),\n\
since `__len__` will fail for a bitarray object with 2^31 or more elements\n\
on a 32bit machine, whereas this method will return the correct value,\n\
on 32bit and 64bit machines.");

PyDoc_STRVAR(len_doc,
"__len__() -> int\n\
\n\
Return the length, i.e. number of bits stored in the bitarray.\n\
This method will fail for a bitarray object with 2^31 or more elements\n\
on a 32bit machine.  Use bitarray.length() instead.");


static PyObject *
bitarray_copy(bitarrayobject* self)
{
    PyObject *res;

    res = newbitarrayobject(Py_TYPE(self), self->nbits, self->endian);
    if (res == NULL)
        return NULL;

    *((bitarrayobject*) res)->words = *self->words;
    return res;
}

PyDoc_STRVAR(copy_doc,
"copy() -> bitarray\n\
\n\
Return a copy of the bitarray.");


static PyObject *
bitarray_count(bitarrayobject *self, PyObject *args)
{
    PyObject *x = Py_True;
    idx_t start = 0, stop = self->nbits;
    long vi;

    if (!PyArg_ParseTuple(args, "|OLL:count", &x, &start, &stop))
        return NULL;

    vi = PyObject_IsTrue(x);
    if (vi < 0)
        return NULL;

    normalize_index(self->nbits, &start);
    normalize_index(self->nbits, &stop);

    return PyLong_FromLongLong(count(self, vi, start, stop));
}

PyDoc_STRVAR(count_doc,
"count(value=True, start=0, stop=<end of array>, /) -> int\n\
\n\
Count the number of occurrences of bool(value) in the bitarray.");


static PyObject *
bitarray_index(bitarrayobject *self, PyObject *args)
{
    PyObject *x;
    idx_t i, start = 0, stop = self->nbits;
    long vi;

    if (!PyArg_ParseTuple(args, "O|LL:index", &x, &start, &stop))
        return NULL;

    vi = PyObject_IsTrue(x);
    if (vi < 0)
        return NULL;

    normalize_index(self->nbits, &start);
    normalize_index(self->nbits, &stop);

    i = findfirst(self, vi, start, stop);
    if (i < 0) {
        PyErr_SetString(PyExc_ValueError, "index(x): x not in bitarray");
        return NULL;
    }
    return PyLong_FromLongLong(i);
}

PyDoc_STRVAR(index_doc,
"index(value, start=0, stop=<end of array>, /) -> int\n\
\n\
Return index of the first occurrence of `bool(value)` in the bitarray.\n\
Raises `ValueError` if the value is not present.");


static PyObject *
bitarray_extend(bitarrayobject *self, PyObject *obj)
{
    if (extend_dispatch(self, obj) < 0)
        return NULL;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(extend_doc,
"extend(iterable, /)\n\
\n\
Append bits to the end of the bitarray.  The objects which can be passed\n\
to this method are the same iterable objects which can given to a bitarray\n\
object upon initialization.");


static PyObject *
bitarray_contains(bitarrayobject *self, PyObject *x)
{
    long res;

    if (IS_INT_OR_BOOL(x)) {
        int vi;

        vi = IntBool_AsInt(x);
        if (vi < 0)
            return NULL;
        res = findfirst(self, vi, 0, self->nbits) >= 0;
    }
    else if (bitarray_Check(x)) {
        res = search(self, (bitarrayobject *) x, 0) >= 0;
    }
    else {
        PyErr_SetString(PyExc_TypeError, "bitarray or bool expected");
        return NULL;
    }
    return PyBool_FromLong(res);
}

PyDoc_STRVAR(contains_doc,
"__contains__(value, /) -> bool\n\
\n\
Return True if bitarray contains value, False otherwise.\n\
The value may be a boolean (or integer between 0 and 1), or a bitarray.");


static PyObject *
bitarray_search(bitarrayobject *self, PyObject *args)
{
    PyObject *list = NULL;   /* list of matching positions to be returned */
    PyObject *x, *item = NULL;
    Py_ssize_t limit = -1;
    bitarrayobject *xa;
    idx_t p;

    if (!PyArg_ParseTuple(args, "O|n:_search", &x, &limit))
        return NULL;

    if (!bitarray_Check(x)) {
        PyErr_SetString(PyExc_TypeError, "bitarray expected for search");
        return NULL;
    }
    xa = (bitarrayobject *) x;
    if (xa->nbits == 0) {
        PyErr_SetString(PyExc_ValueError, "can't search for empty bitarray");
        return NULL;
    }
    list = PyList_New(0);
    if (list == NULL)
        return NULL;
    if (xa->nbits > self->nbits || limit == 0)
        return list;

    p = 0;
    while (1) {
        p = search(self, xa, p);
        if (p < 0)
            break;
        item = PyLong_FromLongLong(p);
        p++;
        if (item == NULL || PyList_Append(list, item) < 0) {
            Py_XDECREF(item);
            Py_XDECREF(list);
            return NULL;
        }
        Py_DECREF(item);
        if (limit > 0 && PyList_Size(list) >= limit)
            break;
    }
    return list;
}

PyDoc_STRVAR(search_doc,
"search(bitarray, limit=<none>, /) -> list\n\
\n\
Searches for the given bitarray in self, and return the list of start\n\
positions.\n\
The optional argument limits the number of search results to the integer\n\
specified.  By default, all search results are returned.");

static PyObject *
bitarray_buffer_info(bitarrayobject *self)
{
    PyObject *res, *ptr;

    ptr = PyLong_FromVoidPtr(&(*self->words)[0]),
    res = Py_BuildValue("OLsiL",
                        ptr,
                        (idx_t) Py_SIZE(self),
                        ENDIAN_STR(self),
                        (int) (8*BYTES(self->nbits) - self->nbits),
                        (idx_t) self->words->capacity());
    Py_DECREF(ptr);
    return res;
}

PyDoc_STRVAR(buffer_info_doc,
"buffer_info() -> tuple\n\
\n\
Return a tuple (address, size, endianness, unused, allocated) giving the\n\
current memory address, the size (in bytes) used to hold the bitarray's\n\
contents, the bit endianness as a string, the number of unused bits\n\
(e.g. a bitarray of length 11 will have a buffer size of 2 bytes and\n\
5 unused bits), and the size (in bytes) of the allocated memory.");


static PyObject *
bitarray_endian(bitarrayobject *self)
{
#ifdef IS_PY3K
    return PyUnicode_FromString(ENDIAN_STR(self));
#else
    return PyString_FromString(ENDIAN_STR(self));
#endif
}

PyDoc_STRVAR(endian_doc,
"endian() -> str\n\
\n\
Return the bit endianness as a string (either `little` or `big`).");


static PyObject *
bitarray_append(bitarrayobject *self, PyObject *v)
{
    if (append_item(self, v) < 0)
        return NULL;

    Py_RETURN_NONE;
}

PyDoc_STRVAR(append_doc,
"append(item, /)\n\
\n\
Append the value `bool(item)` to the end of the bitarray.");


static PyObject *
bitarray_all(bitarrayobject *self)
{
    if (findfirst(self, 0, 0, self->nbits) >= 0)
        Py_RETURN_FALSE;
    else
        Py_RETURN_TRUE;
}

PyDoc_STRVAR(all_doc,
"all() -> bool\n\
\n\
Returns True when all bits in the array are True.");


static PyObject *
bitarray_any(bitarrayobject *self)
{
    if (findfirst(self, 1, 0, self->nbits) >= 0)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(any_doc,
"any() -> bool\n\
\n\
Returns True when any bit in the array is True.");


static PyObject *
bitarray_reduce(bitarrayobject *self)
{
    PyObject *dict, *repr = NULL, *result = NULL;
    char *str;

    dict = PyObject_GetAttrString((PyObject *) self, "__dict__");
    if (dict == NULL) {
        PyErr_Clear();
        dict = Py_None;
        Py_INCREF(dict);
    }
    /* the first byte indicates the number of unused bits at the end, and
       the rest of the bytes consist of the raw binary data */
    str = (char *) PyMem_Malloc(Py_SIZE(self) + 1);
    if (str == NULL) {
        PyErr_NoMemory();
        goto error;
    }
    str[0] = (char) setunused(self);
    memcpy(str + 1, &(*self->words)[0], self->words->size());
    repr = PyBytes_FromStringAndSize(str, self->words->size() + 1);
    if (repr == NULL)
        goto error;
    PyMem_Free((void *) str);
    result = Py_BuildValue("O(Os)O", Py_TYPE(self),
                           repr, ENDIAN_STR(self), dict);
error:
    Py_DECREF(dict);
    Py_XDECREF(repr);
    return result;
}

PyDoc_STRVAR(reduce_doc, "state information for pickling");

static PyObject *
bitarray_reverse(bitarrayobject *self)
{

    if (self->nbits < 2)
        Py_RETURN_NONE;

    bit::reverse(self->bits, bit::next(self->bits, self->nbits));
    Py_RETURN_NONE;
}

PyDoc_STRVAR(reverse_doc,
"reverse()\n\
\n\
Reverse the order of bits in the array (in-place).");


static PyObject *
bitarray_fill(bitarrayobject *self)
{
    long p;

    p = setunused(self);
    self->nbits += p;
#ifdef IS_PY3K
    return PyLong_FromLong(p);
#else
    return PyInt_FromLong(p);
#endif
}

PyDoc_STRVAR(fill_doc,
"fill() -> int\n\
\n\
Adds zeros to the end of the bitarray, such that the length of the bitarray\n\
will be a multiple of 8.  Returns the number of bits added (0..7).");

static PyObject *
bitarray_invert(bitarrayobject *self)
{
    invert(self);
    Py_RETURN_NONE;
}

PyDoc_STRVAR(invert_doc,
"invert()\n\
\n\
Invert all bits in the array (in-place),\n\
i.e. convert each 1-bit into a 0-bit and vice versa.");


static PyObject *
bitarray_bytereverse(bitarrayobject *self)
{
    for (auto it = self->words->begin(); it != self->words->end(); ++it) {
        *it = bit::_bitswap(*it);
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(bytereverse_doc,
"bytereverse()\n\
\n\
For all bytes representing the bitarray, reverse the bit order (in-place).\n\
Note: This method changes the actual machine values representing the\n\
bitarray; it does not change the endianness of the bitarray object.");

static PyObject *
bitarray_setall(bitarrayobject *self, PyObject *v)
{
    long vi;

    vi = PyObject_IsTrue(v);
    if (vi < 0)
        return NULL;
    setrange(self, 0, self->nbits, vi);

    Py_RETURN_NONE;
}

PyDoc_STRVAR(setall_doc,
"setall(value, /)\n\
\n\
Set all bits in the bitarray to `bool(value)`.");


// TODO will replace with bit::sort once complete (currently no custom comparator support)
static PyObject *
bitarray_sort(bitarrayobject *self, PyObject *args, PyObject *kwds)
{
    idx_t n, n0, n1;
    int reverse = 0;
    static char *kwlist[] = {(char*) "reverse", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i:sort", kwlist, &reverse))
        return NULL;

    n = self->nbits;
    n1 = count(self, 1, 0, n);

    if (reverse) {
        setrange(self, 0, n1, 1);
        setrange(self, n1, n, 0);
    }
    else {
        n0 = n - n1;
        setrange(self, 0, n0, 0);
        setrange(self, n0, n, 1);
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(sort_doc,
"sort(reverse=False)\n\
\n\
Sort the bits in the array (in-place).");

/* since too many details differ between the Python 2 and 3 implementation
   of this function, we choose to have two separate function implementation,
   even though this means some of the code is duplicated in the two versions
*/
#ifdef IS_PY3K
static PyObject *
bitarray_fromfile(bitarrayobject *self, PyObject *args)
{
    PyObject *f;
    Py_ssize_t newsize, nbytes = -1;
    PyObject *reader, *rargs, *result;
    size_t nread;
    idx_t t, p;

    if (!PyArg_ParseTuple(args, "O|n:fromfile", &f, &nbytes))
        return NULL;

    if (nbytes == 0)
        Py_RETURN_NONE;

    reader = PyObject_GetAttrString(f, "read");
    if (reader == NULL)
    {
        PyErr_SetString(PyExc_TypeError,
                        "first argument must be an open file");
        return NULL;
    }
    rargs = Py_BuildValue("(n)", nbytes);
    if (rargs == NULL) {
        Py_DECREF(reader);
        return NULL;
    }
    result = PyEval_CallObject(reader, rargs);
    if (result != NULL) {
        if (!PyBytes_Check(result)) {
            PyErr_SetString(PyExc_TypeError,
                            "first argument must be an open file");
            Py_DECREF(result);
            Py_DECREF(rargs);
            Py_DECREF(reader);
            return NULL;
        }
        nread = PyBytes_Size(result);

        t = self->nbits;
        p = setunused(self);
        self->nbits += p;

        newsize = Py_SIZE(self) + nread;

        if (resize(self, BITS(newsize)) < 0) {
            Py_DECREF(result);
            Py_DECREF(rargs);
            Py_DECREF(reader);
            return NULL;
        }
        memcpy(&(*self->words)[0] + (Py_SIZE(self) - nread),
               PyBytes_AS_STRING(result), nread);

        if (nbytes > 0 && nread < (size_t) nbytes) {
            PyErr_SetString(PyExc_EOFError, "not enough items read");
            return NULL;
        }
        if (delete_n(self, t, p) < 0)
            return NULL;
        Py_DECREF(result);
    }

    Py_DECREF(rargs);
    Py_DECREF(reader);
    Py_RETURN_NONE;
}
#else  /* Python 2 */
static PyObject *
bitarray_fromfile(bitarrayobject *self, PyObject *args)
{
    PyObject *f;
    FILE *fp;
    Py_ssize_t newsize, nbytes = -1;
    size_t nread;
    idx_t t, p;
    long cur;

    if (!PyArg_ParseTuple(args, "O|n:fromfile", &f, &nbytes))
        return NULL;

    fp = PyFile_AsFile(f);
    if (fp == NULL) {
        PyErr_SetString(PyExc_TypeError,
                        "first argument must be an open file");
        return NULL;
    }

    /* find number of bytes till EOF */
    if (nbytes < 0) {
        if ((cur = ftell(fp)) < 0)
            goto EOFerror;

        if (fseek(fp, 0L, SEEK_END) || (nbytes = ftell(fp)) < 0)
            goto EOFerror;

        nbytes -= cur;
        if (fseek(fp, cur, SEEK_SET)) {
        EOFerror:
            PyErr_SetString(PyExc_EOFError, "could not find EOF");
            return NULL;
        }
    }
    if (nbytes == 0)
        Py_RETURN_NONE;

    /* file exists and there are more than zero bytes to read */
    t = self->nbits;
    p = setunused(self);
    self->nbits += p;

    newsize = Py_SIZE(self) + nbytes;
    if (resize(self, BITS(newsize)) < 0)
        return NULL;

    nread = fread(&(*self->words)[0] + (Py_SIZE(self) - nbytes), 1, nbytes, fp);
    if (nread < (size_t) nbytes) {
        newsize -= nbytes - nread;
        if (resize(self, BITS(newsize)) < 0)
            return NULL;
        PyErr_SetString(PyExc_EOFError, "not enough items in file");
        return NULL;
    }

    if (delete_n(self, t, p) < 0)
        return NULL;
    Py_RETURN_NONE;
}
#endif

PyDoc_STRVAR(fromfile_doc,
"fromfile(f, n=<till EOF>, /)\n\
\n\
Read n bytes from the file object f and append them to the bitarray\n\
interpreted as machine values.  When n is omitted, as many bytes are\n\
read until EOF is reached.");

/* since too many details differ between the Python 2 and 3 implementation
   of this function, we choose to have two separate function implementation
*/
#ifdef IS_PY3K
static PyObject *
bitarray_tofile(bitarrayobject *self, PyObject *f)
{
    PyObject *writer, *value, *args, *result;

    if (f == NULL) {
        PyErr_SetString(PyExc_TypeError, "writeobject with NULL file");
        return NULL;
    }
    writer = PyObject_GetAttrString(f, "write");
    if (writer == NULL)
        return NULL;
    setunused(self);
    value = PyBytes_FromStringAndSize(reinterpret_cast<char*>(&(*self->words)[0]), Py_SIZE(self));
    if (value == NULL) {
        Py_DECREF(writer);
        return NULL;
    }
    args = PyTuple_Pack(1, value);
    if (args == NULL) {
        Py_DECREF(value);
        Py_DECREF(writer);
        return NULL;
    }
    result = PyEval_CallObject(writer, args);
    Py_DECREF(args);
    Py_DECREF(value);
    Py_DECREF(writer);
    if (result == NULL)
    {
        PyErr_SetString(PyExc_TypeError, "open file expected");
        return NULL;
    }
    Py_DECREF(result);
    Py_RETURN_NONE;
}
#else  /* Python 2 */
static PyObject *
bitarray_tofile(bitarrayobject *self, PyObject *f)
{
    FILE *fp;

    fp = PyFile_AsFile(f);
    if (fp == NULL) {
        PyErr_SetString(PyExc_TypeError, "open file expected");
        return NULL;
    }
    if (Py_SIZE(self) == 0)
        Py_RETURN_NONE;

    setunused(self);
    if (fwrite(&(*self->words)[0], 1, Py_SIZE(self), fp) !=
        (size_t) Py_SIZE(self))
    {
        PyErr_SetFromErrno(PyExc_IOError);
        clearerr(fp);
        return NULL;
    }
    Py_RETURN_NONE;
}
#endif

PyDoc_STRVAR(tofile_doc,
"tofile(f, /)\n\
\n\
Write all bits (as machine values) to the file object f.\n\
When the length of the bitarray is not a multiple of 8,\n\
the remaining bits (1..7) are set to 0.");

static PyObject *
bitarray_tolist(bitarrayobject *self)
{
    PyObject *list;
    idx_t i;

    list = PyList_New((Py_ssize_t) self->nbits);
    if (list == NULL)
        return NULL;

    for (i = 0; i < self->nbits; i++) {
        if (self->bits[i] == bit::bit1) {
            if (PyList_SetItem(list, (Py_ssize_t) i, 
                        self->bits[i] ? Py_True : Py_False) < 0) {
                return NULL;
            }
        }
    }
    return list;
}

PyDoc_STRVAR(tolist_doc,
"tolist() -> list\n\
\n\
Return an ordinary list with the items in the bitarray.\n\
Note that the list object being created will require 32 or 64 times more\n\
memory than the bitarray object, which may cause a memory error if the\n\
bitarray is very large.\n\
Also note that to extend a bitarray with elements from a list,\n\
use the extend method.");


static PyObject *
bitarray_frombytes(bitarrayobject *self, PyObject *bytes)
{
    idx_t t, p;

    if (!PyBytes_Check(bytes)) {
        PyErr_SetString(PyExc_TypeError, "bytes expected");
        return NULL;
    }

    /* Before we extend the raw bytes with the new data, we need to store
       the current size and pad the last byte, as our bitarray size might
       not be a multiple of 8.  After extending, we remove the padding
       bits again.  The same is done in bitarray_fromfile().
    */
    t = self->nbits;
    p = setunused(self);
    self->nbits += p;

    if (extend_rawbytes(self, bytes) < 0)
        return NULL;
    if (delete_n(self, t, p) < 0)
        return NULL;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(frombytes_doc,
"frombytes(bytes, /)\n\
\n\
Append from a byte string, interpreted as machine values.");


static PyObject *
bitarray_tobytes(bitarrayobject *self)
{
    setunused(self);
    return PyBytes_FromStringAndSize(reinterpret_cast<char*>(&(*self->words)[0]), Py_SIZE(self));
}

PyDoc_STRVAR(tobytes_doc,
"tobytes() -> bytes\n\
\n\
Return the byte representation of the bitarray.\n\
When the length of the bitarray is not a multiple of 8, the few remaining\n\
bits (1..7) are considered to be 0.");


static PyObject *
bitarray_to01(bitarrayobject *self)
{
#ifdef IS_PY3K
    PyObject *string;
    PyObject *unpacked;

    unpacked = unpack(self, '0', '1');
    if (unpacked == NULL)
        return NULL;
    string = PyUnicode_FromEncodedObject(unpacked, NULL, NULL);
    Py_DECREF(unpacked);
    return string;
#else
    return unpack(self, '0', '1');
#endif
}

PyDoc_STRVAR(to01_doc,
"to01() -> str\n\
\n\
Return a string containing '0's and '1's, representing the bits in the\n\
bitarray object.\n\
Note: To extend a bitarray from a string containing '0's and '1's,\n\
use the extend method.");


static PyObject *
bitarray_unpack(bitarrayobject *self, PyObject *args, PyObject *kwds)
{
    char zero = 0x00, one = 0xff;
    static char *kwlist[] = {(char*) "zero", (char*) "one", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|cc:unpack", kwlist,
                                     &zero, &one))
        return NULL;

    return unpack(self, zero, one);
}

PyDoc_STRVAR(unpack_doc,
"unpack(zero=b'\\x00', one=b'\\xff') -> bytes\n\
\n\
Return bytes containing one character for each bit in the bitarray,\n\
using the specified mapping.");


static PyObject *
bitarray_pack(bitarrayobject *self, PyObject *bytes)
{
    if (!PyBytes_Check(bytes)) {
        PyErr_SetString(PyExc_TypeError, "bytes expected");
        return NULL;
    }
    if (extend_bytes(self, bytes, BYTES_RAW) < 0)
        return NULL;

    Py_RETURN_NONE;
}
PyDoc_STRVAR(pack_doc,
"pack(bytes, /)\n\
\n\
Extend the bitarray from bytes, where each byte corresponds to a single\n\
bit.  The byte `b'\\x00'` maps to bit 0 and all other characters map to\n\
bit 1.\n\
This method, as well as the unpack method, are meant for efficient\n\
transfer of data between bitarray objects to other python objects\n\
(for example NumPy's ndarray object) which have a different memory view.");


static PyObject *
bitarray_repr(bitarrayobject *self)
{
    PyObject *bytes;
    PyObject *unpacked;
#ifdef IS_PY3K
    PyObject *decoded;
#endif

    if (self->nbits == 0) {
        bytes = PyBytes_FromString("bitarray()");
    }
    else {
        bytes = PyBytes_FromString("bitarray(\'");
        unpacked = unpack(self, '0', '1');
        if (unpacked == NULL)
            return NULL;
        PyBytes_ConcatAndDel(&bytes, unpacked);
        PyBytes_ConcatAndDel(&bytes, PyBytes_FromString("\')"));
    }
#ifdef IS_PY3K
    decoded = PyUnicode_FromEncodedObject(bytes, NULL, NULL);
    Py_DECREF(bytes);
    return decoded;
#else
    return bytes;  /* really a string in Python 2 */
#endif
}

static PyObject *
bitarray_insert(bitarrayobject *self, PyObject *args)
{
    idx_t i;
    PyObject *v;

    if (!PyArg_ParseTuple(args, "LO:insert", &i, &v))
        return NULL;

    normalize_index(self->nbits, &i);

    if (insert_n(self, i, 1) < 0)
        return NULL;
    if (set_item(self, i, v) < 0)
        return NULL;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(insert_doc,
"insert(index, value, /)\n\
\n\
Insert `bool(value)` into the bitarray before index.");


static PyObject *
bitarray_pop(bitarrayobject *self, PyObject *args)
{
    idx_t i = -1;

    if (!PyArg_ParseTuple(args, "|L:pop", &i))
        return NULL;

    if (self->nbits == 0) {
        /* special case -- most common failure cause */
        PyErr_SetString(PyExc_IndexError, "pop from empty bitarray");
        return NULL;
    }
    if (i < 0)
        i += self->nbits;

    if (i < 0 || i >= self->nbits) {
        PyErr_SetString(PyExc_IndexError, "pop index out of range");
        return NULL;
    }
    auto vi = self->bits[i];
    if (delete_n(self, i, 1) < 0)
        return NULL;
    return vi ? Py_True : Py_False;
}

PyDoc_STRVAR(pop_doc,
"pop(index=-1, /) -> item\n\
\n\
Return the i-th (default last) element and delete it from the bitarray.\n\
Raises `IndexError` if bitarray is empty or index is out of range.");


static PyObject *
bitarray_remove(bitarrayobject *self, PyObject *v)
{
    idx_t i;
    long vi;

    vi = PyObject_IsTrue(v);
    if (vi < 0)
        return NULL;

    i = findfirst(self, vi, 0, self->nbits);
    if (i < 0) {
        PyErr_SetString(PyExc_ValueError, "remove(x): x not in bitarray");
        return NULL;
    }
    if (delete_n(self, i, 1) < 0)
        return NULL;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(remove_doc,
"remove(value, /)\n\
\n\
Remove the first occurrence of `bool(value)` in the bitarray.\n\
Raises `ValueError` if item is not present.");


/* --------- special methods ----------- */

// TODO use pdep for step sizes. 
static PyObject *
bitarray_getitem(bitarrayobject *self, PyObject *a)
{
    PyObject *res;
    idx_t start, stop, step, slicelength, j, i = 0;

    if (IS_INDEX(a)) {
        if (getIndex(a, &i) < 0)
            return NULL;
        if (i < 0)
            i += self->nbits;
        if (i < 0 || i >= self->nbits) {
            PyErr_SetString(PyExc_IndexError, "bitarray index out of range");
            return NULL;
        }
        return self->bits[i] ? Py_True : Py_False;
    }
    if (PySlice_Check(a)) {
        if (slice_GetIndicesEx((PySliceObject *) a, self->nbits,
                               &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }
        res = newbitarrayobject(Py_TYPE(self), slicelength, self->endian);
        if (res == NULL)
            return NULL;

        for (i = 0, j = start; i < slicelength; i++, j += step)
            ((bitarrayobject *) res)->bits[i] = self->bits[i];

        return res;
    }
    PyErr_SetString(PyExc_TypeError, "index or slice expected");
    return NULL;
}

/* Sets the elements, specified by slice, in self to the value(s) given by v
   which is either a bitarray or a boolean.
*/
static int
setslice(bitarrayobject *self, PySliceObject *slice, PyObject *v)
{
    idx_t start, stop, step, slicelength, j, i = 0;

    if (slice_GetIndicesEx(slice, self->nbits,
                           &start, &stop, &step, &slicelength) < 0)
        return -1;

    if (bitarray_Check(v)) {
#define vv  ((bitarrayobject *) v)
        if (vv->nbits == slicelength) {
            for (i = 0, j = start; i < slicelength; i++, j += step)
                self->bits[j] = vv->bits[i];
            return 0;
        }
        if (step != 1) {
            char buff[256];
            sprintf(buff, "attempt to assign sequence of size %lld "
                          "to extended slice of size %lld",
                    vv->nbits, (idx_t) slicelength);
            PyErr_SetString(PyExc_ValueError, buff);
            return -1;
        }
        /* make self bigger or smaller */
        if (vv->nbits > slicelength) {
            if (insert_n(self, start, vv->nbits - slicelength) < 0)
                return -1;
        }
        else {
            if (delete_n(self, start, slicelength - vv->nbits) < 0)
                return -1;
        }
        /* copy the new values into self */
        copy_n(self, start, vv, 0, vv->nbits);
#undef vv
        return 0;
    }
    if (IS_INT_OR_BOOL(v)) {
        int vi;

        vi = IntBool_AsInt(v);
        if (vi < 0)
            return -1;
        for (i = 0, j = start; i < slicelength; i++, j += step)
            self->bits[j] = vi ? bit::bit1 : bit::bit0;
        return 0;
    }
    PyErr_SetString(PyExc_IndexError,
                    "bitarray or bool expected for slice assignment");
    return -1;
}

static PyObject *
bitarray_setitem(bitarrayobject *self, PyObject *args)
{
    PyObject *a, *v;
    idx_t i = 0;

    if (!PyArg_ParseTuple(args, "OO:__setitem__", &a, &v))
        return NULL;

    if (IS_INDEX(a)) {
        if (getIndex(a, &i) < 0)
            return NULL;
        if (i < 0)
            i += self->nbits;
        if (i < 0 || i >= self->nbits) {
            PyErr_SetString(PyExc_IndexError, "bitarray index out of range");
            return NULL;
        }
        if (set_item(self, i, v) < 0)
            return NULL;
        Py_RETURN_NONE;
    }
    if (PySlice_Check(a)) {
        if (setslice(self, (PySliceObject *) a, v) < 0)
            return NULL;
        Py_RETURN_NONE;
    }
    PyErr_SetString(PyExc_TypeError, "index or slice expected");
    return NULL;
}

static PyObject *
bitarray_delitem(bitarrayobject *self, PyObject *a)
{
    idx_t start, stop, step, slicelength, j, i = 0;

    if (IS_INDEX(a)) {
        if (getIndex(a, &i) < 0)
            return NULL;
        if (i < 0)
            i += self->nbits;
        if (i < 0 || i >= self->nbits) {
            PyErr_SetString(PyExc_IndexError, "bitarray index out of range");
            return NULL;
        }
        if (delete_n(self, i, 1) < 0)
            return NULL;
        Py_RETURN_NONE;
    }
    if (PySlice_Check(a)) {
        if (slice_GetIndicesEx((PySliceObject *) a, self->nbits,
                               &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }
        if (slicelength == 0)
            Py_RETURN_NONE;

        if (step < 0) {
            stop = start + 1;
            start = stop + step * (slicelength - 1) - 1;
            step = -step;
        }
        if (step == 1) {
            assert(stop - start == slicelength);
            if (delete_n(self, start, slicelength) < 0)
                return NULL;
            Py_RETURN_NONE;
        }
        /* this is the only complicated part when step > 1 */
        for (i = j = start; i < self->nbits; i++)
            if ((i - start) % step != 0 || i >= stop) {
                self->bits[j] = self->bits[i];
                j++;
            }
        if (resize(self, self->nbits - slicelength) < 0)
            return NULL;
        Py_RETURN_NONE;
    }
    PyErr_SetString(PyExc_TypeError, "index or slice expected");
    return NULL;
}

/* ---------- number methods ---------- */

static PyObject *
bitarray_add(bitarrayobject *self, PyObject *other)
{
    PyObject *res;

    res = bitarray_copy(self);
    if (extend_dispatch((bitarrayobject *) res, other) < 0) {
        Py_DECREF(res);
        return NULL;
    }
    return res;
}

static PyObject *
bitarray_iadd(bitarrayobject *self, PyObject *other)
{
    if (extend_dispatch(self, other) < 0)
        return NULL;
    Py_INCREF(self);
    return (PyObject *) self;
}

static PyObject *
bitarray_mul(bitarrayobject *self, PyObject *v)
{
    PyObject *res;
    idx_t vi = 0;

    if (!IS_INDEX(v)) {
        PyErr_SetString(PyExc_TypeError,
                        "integer value expected for bitarray repetition");
        return NULL;
    }
    if (getIndex(v, &vi) < 0)
        return NULL;
    res = bitarray_copy(self);
    if (repeat((bitarrayobject *) res, vi) < 0) {
        Py_DECREF(res);
        return NULL;
    }
    return res;
}

static PyObject *
bitarray_imul(bitarrayobject *self, PyObject *v)
{
    idx_t vi = 0;

    if (!IS_INDEX(v)) {
        PyErr_SetString(PyExc_TypeError,
            "integer value expected for in-place bitarray repetition");
        return NULL;
    }
    if (getIndex(v, &vi) < 0)
        return NULL;
    if (repeat(self, vi) < 0)
        return NULL;
    Py_INCREF(self);
    return (PyObject *) self;
}

static PyObject *
bitarray_cpinvert(bitarrayobject *self)
{
    PyObject *res;

    res = bitarray_copy(self);
    invert((bitarrayobject *) res);
    return res;
}

#define BITWISE_FUNC(oper)  \
static PyObject *                                                   \
bitarray_ ## oper (bitarrayobject *self, PyObject *other)           \
{                                                                   \
    PyObject *res;                                                  \
                                                                    \
    res = bitarray_copy(self);                                      \
    if (bitwise((bitarrayobject *) res, other, OP_ ## oper) < 0) {  \
        Py_DECREF(res);                                             \
        return NULL;                                                \
    }                                                               \
    return res;                                                     \
}

BITWISE_FUNC(and)
BITWISE_FUNC(or)
BITWISE_FUNC(xor)


#define BITWISE_IFUNC(oper)  \
static PyObject *                                            \
bitarray_i ## oper (bitarrayobject *self, PyObject *other)   \
{                                                            \
    if (bitwise(self, other, OP_ ## oper) < 0)               \
        return NULL;                                         \
    Py_INCREF(self);                                         \
    return (PyObject *) self;                                \
}

BITWISE_IFUNC(and)
BITWISE_IFUNC(or)
BITWISE_IFUNC(xor)

/******************* variable length encoding and decoding ***************/

static int
check_codedict(PyObject *codedict)
{
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    if (!PyDict_Check(codedict)) {
        PyErr_SetString(PyExc_TypeError, "dict expected");
        return -1;
    }
    if (PyDict_Size(codedict) == 0) {
        PyErr_SetString(PyExc_ValueError, "prefix code dict empty");
        return -1;
    }
    while (PyDict_Next(codedict, &pos, &key, &value)) {
        if (!bitarray_Check(value)) {
            PyErr_SetString(PyExc_TypeError,
                            "bitarray expected for dict value");
            return -1;
        }
        if (((bitarrayobject *) value)->nbits == 0) {
            PyErr_SetString(PyExc_ValueError, "non-empty bitarray expected");
            return -1;
        }
    }
    return 0;
}

static PyObject *
bitarray_encode(bitarrayobject *self, PyObject *args)
{
    PyObject *codedict, *iterable, *iter, *symbol, *bits;

    if (!PyArg_ParseTuple(args, "OO:encode", &codedict, &iterable))
        return NULL;

    if (check_codedict(codedict) < 0)
        return NULL;

    iter = PyObject_GetIter(iterable);
    if (iter == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return NULL;
    }
    /* extend self with the bitarrays from codedict */
    while ((symbol = PyIter_Next(iter)) != NULL) {
        bits = PyDict_GetItem(codedict, symbol);
        Py_DECREF(symbol);
        if (bits == NULL) {
            PyErr_SetString(PyExc_ValueError,
                            "symbol not defined in prefix code");
            goto error;
        }
        if (extend_bitarray(self, (bitarrayobject *) bits) < 0)
            goto error;
    }
    Py_DECREF(iter);
    if (PyErr_Occurred())
        return NULL;
    Py_RETURN_NONE;
error:
    Py_DECREF(iter);
    return NULL;
}

PyDoc_STRVAR(encode_doc,
"encode(code, iterable, /)\n\
\n\
Given a prefix code (a dict mapping symbols to bitarrays),\n\
iterate over the iterable object with symbols, and extend the bitarray\n\
with the corresponding bitarray for each symbols.");

/* Binary tree definition */
typedef struct _bin_node
{
    struct _bin_node *child[2];
    PyObject *symbol;
} binode;


static binode *
new_binode(void)
{
    binode *nd;

    nd = (binode *) PyMem_Malloc(sizeof(binode));
    if (nd == NULL) {
        PyErr_NoMemory();
        return NULL;
    }
    nd->child[0] = NULL;
    nd->child[1] = NULL;
    nd->symbol = NULL;
    return nd;
}

static void
delete_binode_tree(binode *tree)
{
    if (tree == NULL)
        return;

    delete_binode_tree(tree->child[0]);
    delete_binode_tree(tree->child[1]);
    PyMem_Free(tree);
}

static int
insert_symbol(binode *tree, bitarrayobject *ba, PyObject *symbol)
{
    binode *nd = tree, *prev;
    Py_ssize_t i;
    int k;

    for (i = 0; i < ba->nbits; i++) {
        k = ba->bits[i] ? 1 : 0;
        prev = nd;
        nd = nd->child[k];

        /* we cannot have already a symbol when branching to the new leaf */
        if (nd && nd->symbol)
            goto ambiguity;

        if (!nd) {
            nd = new_binode();
            if (nd == NULL)
                return -1;
            prev->child[k] = nd;
        }
    }
    /* the new leaf node cannot already have a symbol or children */
    if (nd->symbol || nd->child[0] || nd->child[1])
        goto ambiguity;

    nd->symbol = symbol;
    return 0;

 ambiguity:
    PyErr_SetString(PyExc_ValueError, "prefix code ambiguous");
    return -1;
}

static binode *
make_tree(PyObject *codedict)
{
    binode *tree;
    PyObject *symbol, *array;
    Py_ssize_t pos = 0;

    tree = new_binode();
    if (tree == NULL)
        return NULL;

    while (PyDict_Next(codedict, &pos, &symbol, &array)) {
        if (insert_symbol(tree, (bitarrayobject *) array, symbol) < 0) {
            delete_binode_tree(tree);
            return NULL;
        }
    }
    return tree;
}

/*
  Traverse tree using the branches corresponding to the bitarray `ba`,
  starting at *indexp.  Return the symbol at the leaf node, or NULL
  when the end of the bitarray has been reached, or on error (in which
  case the appropriate PyErr_SetString is set.
*/
static PyObject *
traverse_tree(binode *tree, bitarrayobject *ba, idx_t *indexp)
{
    binode *nd = tree;
    int k;

    while (*indexp < ba->nbits) {
        k = ba->bits[*indexp] ? 1 : 0;
        (*indexp)++;
        nd = nd->child[k];
        if (nd == NULL) {
            PyErr_SetString(PyExc_ValueError,
                            "prefix code does not match data in bitarray");
            return NULL;
        }
        if (nd->symbol)  /* leaf */
            return nd->symbol;
    }
    if (nd != tree)
        PyErr_SetString(PyExc_ValueError, "decoding not terminated");

    return NULL;
}

static PyObject *
bitarray_decode(bitarrayobject *self, PyObject *codedict)
{
    binode *tree, *nd;
    PyObject *list;
    Py_ssize_t i;
    int k;

    if (check_codedict(codedict) < 0)
        return NULL;

    tree = make_tree(codedict);
    if (tree == NULL || PyErr_Occurred())
        return NULL;

    nd = tree;
    list = PyList_New(0);
    if (list == NULL) {
        delete_binode_tree(tree);
        return NULL;
    }
    /* traverse tree (just like above) */
    for (i = 0; i < self->nbits; i++) {
        k = self->bits[i] ? 1 : 0;
        nd = nd->child[k];
        if (nd == NULL) {
            PyErr_SetString(PyExc_ValueError,
                            "prefix code does not match data in bitarray");
            goto error;
        }
        if (nd->symbol) {  /* leaf */
            if (PyList_Append(list, nd->symbol) < 0)
                goto error;
            nd = tree;
        }
    }
    if (nd != tree) {
        PyErr_SetString(PyExc_ValueError, "decoding not terminated");
        goto error;
    }
    delete_binode_tree(tree);
    return list;

error:
    delete_binode_tree(tree);
    Py_DECREF(list);
    return NULL;
}

PyDoc_STRVAR(decode_doc,
"decode(code, /) -> list\n\
\n\
Given a prefix code (a dict mapping symbols to bitarrays),\n\
decode the content of the bitarray and return it as a list of symbols.");


/*********************** (Bitarray) Decode Iterator *********************/


typedef struct {
    PyObject_HEAD
    bitarrayobject *bao;        /* bitarray we're searching in */
    binode *tree;               /* prefix tree containing symbols */
    idx_t index;                /* current index in bitarray */
} decodeiterobject;

static PyObject *
decodeiter_next(decodeiterobject *it);

static void
decodeiter_dealloc(decodeiterobject *it);

static int
decodeiter_traverse(decodeiterobject *it, visitproc visit, void *arg);

static PyTypeObject DecodeIter_Type = {
#ifdef IS_PY3K
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                                        /* ob_size */
#endif
    "bitarraydecodeiterator",                 /* tp_name */
    sizeof(decodeiterobject),                 /* tp_basicsize */
    0,                                        /* tp_itemsize */
    /* methods */
    (destructor) decodeiter_dealloc,          /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    PyObject_GenericGetAttr,                  /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
    0,                                        /* tp_doc */
    (traverseproc) decodeiter_traverse,       /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    PyObject_SelfIter,                        /* tp_iter */
    (iternextfunc) decodeiter_next,           /* tp_iternext */
    0,                                        /* tp_methods */
};

#define DecodeIter_Check(op)  PyObject_TypeCheck(op, &DecodeIter_Type)

/* create a new initialized bitarray search iterator object */
static PyObject *
bitarray_iterdecode(bitarrayobject *self, PyObject *codedict)
{
    decodeiterobject *it;  /* iterator to be returned */
    binode *tree;

    if (check_codedict(codedict) < 0)
        return NULL;

    tree = make_tree(codedict);
    if (tree == NULL || PyErr_Occurred())
        return NULL;

    it = PyObject_GC_New(decodeiterobject, &DecodeIter_Type);
    if (it == NULL)
        return NULL;

    it->tree = tree;

    Py_INCREF(self);
    it->bao = self;
    it->index = 0;
    PyObject_GC_Track(it);
    return (PyObject *) it;
}

PyDoc_STRVAR(iterdecode_doc,
"iterdecode(code, /) -> iterator\n\
\n\
Given a prefix code (a dict mapping symbols to bitarrays),\n\
decode the content of the bitarray and return an iterator over\n\
the symbols.");

static PyObject *
decodeiter_next(decodeiterobject *it)
{
    PyObject *symbol;

    assert(DecodeIter_Check(it));
    symbol = traverse_tree(it->tree, it->bao, &(it->index));
    if (symbol == NULL)  /* stop iteration OR error occured */
        return NULL;
    Py_INCREF(symbol);
    return symbol;
}

static void
decodeiter_dealloc(decodeiterobject *it)
{
    delete_binode_tree(it->tree);
    PyObject_GC_UnTrack(it);
    Py_XDECREF(it->bao);
    PyObject_GC_Del(it);
}

static int
decodeiter_traverse(decodeiterobject *it, visitproc visit, void *arg)
{
    Py_VISIT(it->bao);
    return 0;
}


/*********************** (Bitarray) Search Iterator *********************/

typedef struct {
    PyObject_HEAD
    bitarrayobject *bao;        /* bitarray we're searching in */
    bitarrayobject *xa;         /* bitarray being searched for */
    idx_t p;                    /* current search position */
} searchiterobject;


static PyObject *
searchiter_next(searchiterobject *it);

static void
searchiter_dealloc(searchiterobject *it);

static int
searchiter_traverse(searchiterobject *it, visitproc visit, void *arg);

static PyTypeObject SearchIter_Type = {
#ifdef IS_PY3K
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                                        /* ob_size */
#endif
    "bitarraysearchiterator",                 /* tp_name */
    sizeof(searchiterobject),                 /* tp_basicsize */
    0,                                        /* tp_itemsize */
    /* methods */
    (destructor) searchiter_dealloc,          /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    PyObject_GenericGetAttr,                  /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
    0,                                        /* tp_doc */
    (traverseproc) searchiter_traverse,       /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    PyObject_SelfIter,                        /* tp_iter */
    (iternextfunc) searchiter_next,           /* tp_iternext */
    0,                                        /* tp_methods */
};

#define SearchIter_Check(op)  PyObject_TypeCheck(op, &SearchIter_Type)

/* create a new initialized bitarray search iterator object */
static PyObject *
bitarray_itersearch(bitarrayobject *self, PyObject *x)
{
    searchiterobject *it;  /* iterator to be returned */
    bitarrayobject *xa;

    if (!bitarray_Check(x)) {
        PyErr_SetString(PyExc_TypeError, "bitarray expected for itersearch");
        return NULL;
    }
    xa = (bitarrayobject *) x;
    if (xa->nbits == 0) {
        PyErr_SetString(PyExc_ValueError, "can't search for empty bitarray");
        return NULL;
    }

    it = PyObject_GC_New(searchiterobject, &SearchIter_Type);
    if (it == NULL)
        return NULL;

    Py_INCREF(self);
    it->bao = self;
    Py_INCREF(xa);
    it->xa = xa;
    it->p = 0;  /* start search at position 0 */
    PyObject_GC_Track(it);
    return (PyObject *) it;
}

PyDoc_STRVAR(itersearch_doc,
"itersearch(bitarray, /) -> iterator\n\
\n\
Searches for the given a bitarray in self, and return an iterator over\n\
the start positions where bitarray matches self.");

static PyObject *
searchiter_next(searchiterobject *it)
{
    idx_t p;

    assert(SearchIter_Check(it));
    p = search(it->bao, it->xa, it->p);
    if (p < 0)  /* no more positions -- stop iteration */
        return NULL;
    it->p = p + 1;  /* next search position */
    return PyLong_FromLongLong(p);
}

static void
searchiter_dealloc(searchiterobject *it)
{
    PyObject_GC_UnTrack(it);
    Py_XDECREF(it->bao);
    Py_XDECREF(it->xa);
    PyObject_GC_Del(it);
}

static int
searchiter_traverse(searchiterobject *it, visitproc visit, void *arg)
{
    Py_VISIT(it->bao);
    return 0;
}

/*************************** Method definitions *************************/

static PyMethodDef
bitarray_methods[] = {
    {"all",          (PyCFunction) bitarray_all,         METH_NOARGS,
     all_doc},
    {"any",          (PyCFunction) bitarray_any,         METH_NOARGS,
     any_doc},
    {"append",       (PyCFunction) bitarray_append,      METH_O,
     append_doc},
    {"buffer_info",  (PyCFunction) bitarray_buffer_info, METH_NOARGS,
     buffer_info_doc},
    {"bytereverse",  (PyCFunction) bitarray_bytereverse, METH_NOARGS,
     bytereverse_doc},
    {"copy",         (PyCFunction) bitarray_copy,        METH_NOARGS,
     copy_doc},
    {"count",        (PyCFunction) bitarray_count,       METH_VARARGS,
     count_doc},
    {"decode",       (PyCFunction) bitarray_decode,      METH_O,
     decode_doc},
    {"iterdecode",   (PyCFunction) bitarray_iterdecode,  METH_O,
     iterdecode_doc},
    {"encode",       (PyCFunction) bitarray_encode,      METH_VARARGS,
     encode_doc},
    {"endian",       (PyCFunction) bitarray_endian,      METH_NOARGS,
     endian_doc},
    {"extend",       (PyCFunction) bitarray_extend,      METH_O,
     extend_doc},
    {"fill",         (PyCFunction) bitarray_fill,        METH_NOARGS,
     fill_doc},
    {"fromfile",     (PyCFunction) bitarray_fromfile,    METH_VARARGS,
     fromfile_doc},
    {"frombytes",    (PyCFunction) bitarray_frombytes,   METH_O,
     frombytes_doc},
    {"index",        (PyCFunction) bitarray_index,       METH_VARARGS,
     index_doc},
    {"insert",       (PyCFunction) bitarray_insert,      METH_VARARGS,
     insert_doc},
    {"invert",       (PyCFunction) bitarray_invert,      METH_NOARGS,
     invert_doc},
    {"length",       (PyCFunction) bitarray_length,      METH_NOARGS,
     length_doc},
    {"pack",         (PyCFunction) bitarray_pack,        METH_O,
     pack_doc},
    {"pop",          (PyCFunction) bitarray_pop,         METH_VARARGS,
     pop_doc},
    {"remove",       (PyCFunction) bitarray_remove,      METH_O,
     remove_doc},
    {"reverse",      (PyCFunction) bitarray_reverse,     METH_NOARGS,
     reverse_doc},
    {"setall",       (PyCFunction) bitarray_setall,      METH_O,
     setall_doc},
    {"search",       (PyCFunction) bitarray_search,      METH_VARARGS,
     search_doc},
    {"itersearch",   (PyCFunction) bitarray_itersearch,  METH_O,
     itersearch_doc},
    {"sort",         (PyCFunction) bitarray_sort,        METH_VARARGS |
                                                         METH_KEYWORDS,
     sort_doc},
    {"tofile",       (PyCFunction) bitarray_tofile,      METH_O,
     tofile_doc},
    {"tolist",       (PyCFunction) bitarray_tolist,      METH_NOARGS,
     tolist_doc},
    {"tobytes",      (PyCFunction) bitarray_tobytes,     METH_NOARGS,
     tobytes_doc},
    {"to01",         (PyCFunction) bitarray_to01,        METH_NOARGS,
     to01_doc},
    {"unpack",       (PyCFunction) bitarray_unpack,      METH_VARARGS |
                                                         METH_KEYWORDS,
     unpack_doc},

    /* special methods */
    {"__copy__",     (PyCFunction) bitarray_copy,        METH_NOARGS,
     copy_doc},
    {"__deepcopy__", (PyCFunction) bitarray_copy,        METH_O,
     copy_doc},
    {"__len__",      (PyCFunction) bitarray_length,      METH_NOARGS,
     len_doc},
    {"__contains__", (PyCFunction) bitarray_contains,    METH_O,
     contains_doc},
    {"__reduce__",   (PyCFunction) bitarray_reduce,      METH_NOARGS,
     reduce_doc},

    /* slice methods */
    {"__delitem__",  (PyCFunction) bitarray_delitem,     METH_O,       0},
    {"__getitem__",  (PyCFunction) bitarray_getitem,     METH_O,       0},
    {"__setitem__",  (PyCFunction) bitarray_setitem,     METH_VARARGS, 0},

    /* number methods */
    {"__add__",      (PyCFunction) bitarray_add,         METH_O,       0},
    {"__iadd__",     (PyCFunction) bitarray_iadd,        METH_O,       0},
    {"__mul__",      (PyCFunction) bitarray_mul,         METH_O,       0},
    {"__rmul__",     (PyCFunction) bitarray_mul,         METH_O,       0},
    {"__imul__",     (PyCFunction) bitarray_imul,        METH_O,       0},
    {"__and__",      (PyCFunction) bitarray_and,         METH_O,       0},
    {"__or__",       (PyCFunction) bitarray_or,          METH_O,       0},
    {"__xor__",      (PyCFunction) bitarray_xor,         METH_O,       0},
    {"__iand__",     (PyCFunction) bitarray_iand,        METH_O,       0},
    {"__ior__",      (PyCFunction) bitarray_ior,         METH_O,       0},
    {"__ixor__",     (PyCFunction) bitarray_ixor,        METH_O,       0},
    {"__invert__",   (PyCFunction) bitarray_cpinvert,    METH_NOARGS,  0},

    {NULL,           NULL}  /* sentinel */
};


static PyObject *
bitarray_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *a;  /* to be returned in some cases */
    PyObject *initial = NULL;
    char *endian_str = NULL;
    int endian;
    static char *kwlist[] = {(char*) "initial", (char*) "endian", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds,
                        "|Os:bitarray", kwlist, &initial, &endian_str))
        return NULL;

    if (endian_str == NULL) {
        endian = DEFAULT_ENDIAN;  /* use default value */
    }
    else if (strcmp(endian_str, "little") == 0) {
        endian = 0;
    }
    else if (strcmp(endian_str, "big") == 0) {
        endian = 1;
    }
    else {
        PyErr_SetString(PyExc_ValueError,
                        "endian must be 'little' or 'big'");
        return NULL;
    }

    /* no arg or None */
    if (initial == NULL || initial == Py_None)
        return newbitarrayobject(type, 0, endian);

    /* int, long */
    if (IS_INDEX(initial)) {
        idx_t nbits = 0;

        if (getIndex(initial, &nbits) < 0)
            return NULL;
        if (nbits < 0) {
            PyErr_SetString(PyExc_ValueError,
                            "cannot create bitarray with negative length");
            return NULL;
        }
        return newbitarrayobject(type, nbits, endian);
    }

    /* from bitarray itself */
    if (bitarray_Check(initial)) {
#define np  ((bitarrayobject *) initial)
        a = newbitarrayobject(type, np->nbits,
                              endian_str == NULL ? np->endian : endian);
        if (a == NULL)
            return NULL;
        memcpy(&(*((bitarrayobject *) a)->words)[0], &(*np->words)[0], Py_SIZE(np));
#undef np
        return a;
    }

    /* bytes */
    if (PyBytes_Check(initial)) {
        Py_ssize_t strlen;
        char *str;

        strlen = PyBytes_Size(initial);
        if (strlen == 0)        /* empty string */
            return newbitarrayobject(type, 0, endian);

        str = PyBytes_AsString(initial);
        if (0 <= str[0] && str[0] < 8) {
            /* when the first character is smaller than 8, it indicates the
               number of unused bits at the end, and rest of the bytes
               consist of the raw binary data, this is used for pickling */
            if (strlen == 1 && str[0] > 0) {
                PyErr_Format(PyExc_ValueError,
                             "did not expect 0x0%d", (int) str[0]);
                return NULL;
            }
            a = newbitarrayobject(type, BITS(strlen - 1) - ((idx_t) str[0]),
                                  endian);
            if (a == NULL)
                return NULL;
            memcpy(&(*((bitarrayobject *) a)->words)[0], str + 1, strlen - 1);
            return a;
        }
    }

#define CHECK_TYPE(type)  \
    if (Py ## type ## _Check(initial)) {                                  \
        PyErr_SetString(PyExc_TypeError,                                  \
                        "cannot create bitarray from " #type " object");  \
        return NULL;                                                      \
    }
CHECK_TYPE(Float)
CHECK_TYPE(Complex)
#undef CHECK_TYPE

    /* leave remaining type dispatch to the extend method */
    a = newbitarrayobject(type, 0, endian);
    if (a == NULL)
        return NULL;
    if (extend_dispatch((bitarrayobject *) a, initial) < 0) {
        Py_DECREF(a);
        return NULL;
    }
    return a;
}


static PyObject *
richcompare(PyObject *v, PyObject *w, int op)
{
    int cmp, vi, wi;
    idx_t i, vs, ws;

    if (!bitarray_Check(v) || !bitarray_Check(w)) {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }
#define va  ((bitarrayobject *) v)
#define wa  ((bitarrayobject *) w)
    vs = va->nbits;
    ws = wa->nbits;
    if (vs != ws) {
        /* shortcut for EQ/NE: if sizes differ, the bitarrays differ */
        if (op == Py_EQ)
            Py_RETURN_FALSE;
        if (op == Py_NE)
            Py_RETURN_TRUE;
    }

    /* to avoid uninitialized warning for some compilers */
    vi = wi = 0;
    /* search for the first index where items are different */
    for (i = 0; i < vs && i < ws; i++) {
        // TODO is > overloaded for bitref?
        vi = va->bits[i] ? 1 : 0;
        wi = wa->bits[i] ? 1 : 0;
        if (vi != wi) {
            /* we have an item that differs -- first, shortcut for EQ/NE */
            if (op == Py_EQ)
                Py_RETURN_FALSE;
            if (op == Py_NE)
                Py_RETURN_TRUE;
            /* compare the final item using the proper operator */
            switch (op) {
            case Py_LT: cmp = vi <  wi; break;
            case Py_LE: cmp = vi <= wi; break;
            case Py_EQ: cmp = vi == wi; break;
            case Py_NE: cmp = vi != wi; break;
            case Py_GT: cmp = vi >  wi; break;
            case Py_GE: cmp = vi >= wi; break;
            default: return NULL;  /* cannot happen */
            }
            return PyBool_FromLong((long) cmp);
        }
    }
#undef va
#undef wa

    /* no more items to compare -- compare sizes */
    switch (op) {
    case Py_LT: cmp = vs <  ws; break;
    case Py_LE: cmp = vs <= ws; break;
    case Py_EQ: cmp = vs == ws; break;
    case Py_NE: cmp = vs != ws; break;
    case Py_GT: cmp = vs >  ws; break;
    case Py_GE: cmp = vs >= ws; break;
    default: return NULL;  /* cannot happen */
    }
    return PyBool_FromLong((long) cmp);
}
/************************** Bitarray Iterator **************************/

typedef struct {
    PyObject_HEAD
    bitarrayobject *bao;        /* bitarray we're iterating over */
    idx_t index;                /* current index in bitarray */
} bitarrayiterobject;


static PyObject *
bitarrayiter_next(bitarrayiterobject *it);

static void
bitarrayiter_dealloc(bitarrayiterobject *it);

static int
bitarrayiter_traverse(bitarrayiterobject *it, visitproc visit, void *arg);

static PyTypeObject BitarrayIter_Type = {
#ifdef IS_PY3K
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                                        /* ob_size */
#endif
    "bitarrayiterator",                       /* tp_name */
    sizeof(bitarrayiterobject),               /* tp_basicsize */
    0,                                        /* tp_itemsize */
    /* methods */
    (destructor) bitarrayiter_dealloc,        /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    0,                                        /* tp_repr */
    0,                                        /* tp_as_number */
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    0,                                        /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    PyObject_GenericGetAttr,                  /* tp_getattro */
    0,                                        /* tp_setattro */
    0,                                        /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
    0,                                        /* tp_doc */
    (traverseproc) bitarrayiter_traverse,     /* tp_traverse */
    0,                                        /* tp_clear */
    0,                                        /* tp_richcompare */
    0,                                        /* tp_weaklistoffset */
    PyObject_SelfIter,                        /* tp_iter */
    (iternextfunc) bitarrayiter_next,         /* tp_iternext */
    0,                                        /* tp_methods */
};

#define BitarrayIter_Check(op)  PyObject_TypeCheck(op, &BitarrayIter_Type)

/* create a new initialized bitarray iterator object, this object is
   returned when calling item(a) */
static PyObject *
bitarray_iter(bitarrayobject *self)
{
    bitarrayiterobject *it;

    assert(bitarray_Check(self));
    it = PyObject_GC_New(bitarrayiterobject, &BitarrayIter_Type);
    if (it == NULL)
        return NULL;

    Py_INCREF(self);
    it->bao = self;
    it->index = 0;
    PyObject_GC_Track(it);
    return (PyObject *) it;
}

static PyObject *
bitarrayiter_next(bitarrayiterobject *it)
{
    long vi;

    assert(BitarrayIter_Check(it));
    if (it->index < it->bao->nbits) {
        vi = it->bao->bits[it->index] ? 1 : 0;
        it->index++;
        return PyBool_FromLong(vi);
    }
    return NULL;  /* stop iteration */
}

static void
bitarrayiter_dealloc(bitarrayiterobject *it)
{
    PyObject_GC_UnTrack(it);
    Py_XDECREF(it->bao);
    PyObject_GC_Del(it);
}

static int
bitarrayiter_traverse(bitarrayiterobject *it, visitproc visit, void *arg)
{
    Py_VISIT(it->bao);
    return 0;
}


/********************* Bitarray Buffer Interface ************************/
#ifdef WITH_BUFFER

#if PY_MAJOR_VERSION == 2       /* old buffer protocol */
static Py_ssize_t
bitarray_buffer_getreadbuf(bitarrayobject *self,
                           Py_ssize_t index, const void **ptr)
{
    if (index != 0) {
        PyErr_SetString(PyExc_SystemError, "accessing non-existent segment");
        return -1;
    }
    *ptr = (void *) self->ob_item;
    return Py_SIZE(self);
}

static Py_ssize_t
bitarray_buffer_getwritebuf(bitarrayobject *self,
                            Py_ssize_t index, const void **ptr)
{
    if (index != 0) {
        PyErr_SetString(PyExc_SystemError, "accessing non-existent segment");
        return -1;
    }
    *ptr = (void *) self->ob_item;
    return Py_SIZE(self);
}

static Py_ssize_t
bitarray_buffer_getsegcount(bitarrayobject *self, Py_ssize_t *lenp)
{
    if (lenp)
        *lenp = Py_SIZE(self);
    return 1;
}

static Py_ssize_t
bitarray_buffer_getcharbuf(bitarrayobject *self,
                           Py_ssize_t index, const char **ptr)
{
    if (index != 0) {
        PyErr_SetString(PyExc_SystemError, "accessing non-existent segment");
        return -1;
    }
    *ptr = self->ob_item;
    return Py_SIZE(self);
}

#endif

static int
bitarray_getbuffer(bitarrayobject *self, Py_buffer *view, int flags)
{
    int ret;
    void *ptr;

    if (view == NULL) {
        self->ob_exports++;
        return 0;
    }
    ptr = (void *) &(*self->words)[0];
    ret = PyBuffer_FillInfo(view, (PyObject *) self, ptr,
                            Py_SIZE(self), 0, flags);
    if (ret >= 0) {
        self->ob_exports++;
    }
    return ret;
}

static void
bitarray_releasebuffer(bitarrayobject *self, Py_buffer *view)
{
    self->ob_exports--;
}

static PyBufferProcs bitarray_as_buffer = {
#if PY_MAJOR_VERSION == 2   /* old buffer protocol */
    (readbufferproc) bitarray_buffer_getreadbuf,
    (writebufferproc) bitarray_buffer_getwritebuf,
    (segcountproc) bitarray_buffer_getsegcount,
    (charbufferproc) bitarray_buffer_getcharbuf,
#endif
    (getbufferproc) bitarray_getbuffer,
    (releasebufferproc) bitarray_releasebuffer,
};

#endif  /* WITH_BUFFER */

/************************** Bitarray Type *******************************/

PyTypeObject Bitarraytype = {
#ifdef IS_PY3K
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                                        /* ob_size */
#endif
    "bitarray._bitarray",                     /* tp_name */
    sizeof(bitarrayobject),                   /* tp_basicsize */
    0,                                        /* tp_itemsize */
    /* methods */
    (destructor) bitarray_dealloc,            /* tp_dealloc */
    0,                                        /* tp_print */
    0,                                        /* tp_getattr */
    0,                                        /* tp_setattr */
    0,                                        /* tp_compare */
    (reprfunc) bitarray_repr,                 /* tp_repr */
    0,                                        /* tp_as_number*/
    0,                                        /* tp_as_sequence */
    0,                                        /* tp_as_mapping */
    PyObject_HashNotImplemented,              /* tp_hash */
    0,                                        /* tp_call */
    0,                                        /* tp_str */
    PyObject_GenericGetAttr,                  /* tp_getattro */
    0,                                        /* tp_setattro */
#ifdef WITH_BUFFER
    &bitarray_as_buffer,                      /* tp_as_buffer */
#else
    0,                                        /* tp_as_buffer */
#endif
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_WEAKREFS
#if defined(WITH_BUFFER) && PY_MAJOR_VERSION == 2
    | Py_TPFLAGS_HAVE_NEWBUFFER
#endif
    ,                                         /* tp_flags */
    0,                                        /* tp_doc */
    0,                                        /* tp_traverse */
    0,                                        /* tp_clear */
    richcompare,                              /* tp_richcompare */
    offsetof(bitarrayobject, weakreflist),    /* tp_weaklistoffset */
    (getiterfunc) bitarray_iter,              /* tp_iter */
    0,                                        /* tp_iternext */
    bitarray_methods,                         /* tp_methods */
    0,                                        /* tp_members */
    0,                                        /* tp_getset */
    0,                                        /* tp_base */
    0,                                        /* tp_dict */
    0,                                        /* tp_descr_get */
    0,                                        /* tp_descr_set */
    0,                                        /* tp_dictoffset */
    0,                                        /* tp_init */
    PyType_GenericAlloc,                      /* tp_alloc */
    bitarray_new,                             /* tp_new */
    PyObject_Del,                             /* tp_free */
};

/*************************** Module functions **********************/

static PyObject *
bitdiff(PyObject *self, PyObject *args)
{
    PyObject *a, *b;
    Py_ssize_t i;
    idx_t res = 0;
    unsigned char c;

    if (!PyArg_ParseTuple(args, "OO:bitdiff", &a, &b))
        return NULL;
    if (!(bitarray_Check(a) && bitarray_Check(b))) {
        PyErr_SetString(PyExc_TypeError, "bitarray object expected");
        return NULL;
    }

#define aa  ((bitarrayobject *) a)
#define bb  ((bitarrayobject *) b)
    if (aa->nbits != bb->nbits) {
        PyErr_SetString(PyExc_ValueError,
                        "bitarrays of equal length expected");
        return NULL;
    }
    setunused(aa);
    setunused(bb);
    // TODO This can probably be done with a combination of accumulate and transform
    for (i = 0; i < Py_SIZE(aa); i++) {
        c = (*aa->words)[i] ^ (*bb->words)[i];
        res += bit::_popcnt(c);
    }
#undef aa
#undef bb
    return PyLong_FromLongLong(res);
}

PyDoc_STRVAR(bitdiff_doc,
"bitdiff(a, b, /) -> int\n\
\n\
Return the difference between two bitarrays a and b.\n\
This is function does the same as (a ^ b).count(), but is more memory\n\
efficient, as no intermediate bitarray object gets created.\n\
Deprecated since version 1.2.0, use `bitarray.util.count_xor()` instead.");


static PyObject *
bits2bytes(PyObject *self, PyObject *v)
{
    idx_t n = 0;

    if (!IS_INDEX(v)) {
        PyErr_SetString(PyExc_TypeError, "integer expected");
        return NULL;
    }
    if (getIndex(v, &n) < 0)
        return NULL;
    if (n < 0) {
        PyErr_SetString(PyExc_ValueError, "non-negative integer expected");
        return NULL;
    }
    return PyLong_FromLongLong(BYTES(n));
}

PyDoc_STRVAR(bits2bytes_doc,
"bits2bytes(n, /) -> int\n\
\n\
Return the number of bytes necessary to store n bits.");


static PyObject *
sysinfo(void)
{
    return Py_BuildValue("iiiiL",
                         (int) sizeof(void *),
                         (int) sizeof(size_t),
                         (int) sizeof(Py_ssize_t),
                         (int) sizeof(idx_t),
                         (idx_t) PY_SSIZE_T_MAX);
}

PyDoc_STRVAR(sysinfo_doc,
"_sysinfo() -> tuple\n\
\n\
tuple(sizeof(void *),\n\
      sizeof(size_t),\n\
      sizeof(Py_ssize_t),\n\
      sizeof(idx_t),\n\
      PY_SSIZE_T_MAX)");

/*
   In retrospect, I wish I had never added any modules functions here.
   These, and possibly many others should be part of a separate utility
   module.  Anyway, at this point (2019) it is too late to remove them,
   so I will just leave them here, but not any new ones.
*/
static PyMethodDef module_functions[] = {
    {"bitdiff",    (PyCFunction) bitdiff,    METH_VARARGS, bitdiff_doc   },
    {"bits2bytes", (PyCFunction) bits2bytes, METH_O,       bits2bytes_doc},
    {"_sysinfo",   (PyCFunction) sysinfo,    METH_NOARGS,  sysinfo_doc   },
    {NULL,         NULL}  /* sentinel */
};

/*********************** Install Module **************************/

#ifdef IS_PY3K
static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "_bitarray", 0, -1, module_functions,
};
PyMODINIT_FUNC
PyInit__bitarray(void)
#else
PyMODINIT_FUNC
init_bitarray(void)
#endif
{
    PyObject *m;

    Py_TYPE(&Bitarraytype) = &PyType_Type;
    Py_TYPE(&SearchIter_Type) = &PyType_Type;
    Py_TYPE(&DecodeIter_Type) = &PyType_Type;
    Py_TYPE(&BitarrayIter_Type) = &PyType_Type;
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
    if (m == NULL)
        return NULL;
#else
    m = Py_InitModule3("_bitarray", module_functions, 0);
    if (m == NULL)
        return;
#endif

    Py_INCREF((PyObject *) &Bitarraytype);
    PyModule_AddObject(m, "_bitarray", (PyObject *) &Bitarraytype);
#ifdef IS_PY3K
    return m;
#endif
}
