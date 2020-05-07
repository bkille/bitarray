#include "Python.h"
#include <limits>
#include <iostream>

template <class T>
class PyAlloc {
    public:
    // type definitions
    typedef T        value_type;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef T&       reference;
    typedef const T& const_reference;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;

    // rebind allocator to type U
    template <class U>
    struct rebind {
       typedef PyAlloc<U> other;
    };

    // return address of values
    pointer address (reference value) const {
       return &value;
    }
    const_pointer address (const_reference value) const {
       return &value;
    }

    /* constructors and destructor
    * - nothing to do because the allocator has no state
    */
    PyAlloc() throw() {
    }
    PyAlloc(const PyAlloc&) throw() {
    }
    template <class U>
     PyAlloc (const PyAlloc<U>&) throw() {
    }
    ~PyAlloc() throw() {
    }

    // return maximum number of elements that can be allocated
    size_type max_size () const throw() {
       return std::numeric_limits<std::size_t>::max() / sizeof(T);
    }

    // initialize elements of allocated storage p with value value
    // Allocate memory
    pointer allocate(size_type num, const_pointer /* hint */ = 0)
    {
       //std::cerr << "allocate " << num << " element(s)"
                 //<< " of size " << sizeof(T) << std::endl;
       if(num > max_size()) {
           throw std::bad_alloc();
       }
       return static_cast<pointer>(PyMem_New(value_type, num));
    }

    // destroy elements of initialized storage p
    void destroy (pointer p) {
       // destroy objects by calling their destructor
       p->~T();
    }

    // deallocate storage p of deleted elements
    void deallocate (pointer p, size_type num) {
       // print message and deallocate memory with global delete
       //std::cerr << "deallocate " << num << " element(s)"
                 //<< " of size " << sizeof(T)
                 //<< " at: " << (void*)p << std::endl;
       PyMem_Del((void*) p);
    }
};

// return that all specializations of this allocator are interchangeable
template <class T1, class T2>
bool operator== (const PyAlloc<T1>&,
                const PyAlloc<T2>&) throw() {
   return true;
}
template <class T1, class T2>
bool operator!= (const PyAlloc<T1>&,
                const PyAlloc<T2>&) throw() {
   return false;
}
