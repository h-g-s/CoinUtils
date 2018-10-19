#include <cstdlib>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <limits>
#include <cstddef>
#include <algorithm>
#include "vint_set.h"

#define VINT_SET_MIN_CAP 512
#define VINT_SET_FLUSH 1000000

struct _IntSet {
    std::vector<int> elements;
    size_t notUpdated;
};

IntSet* vint_set_create() {
    IntSet *iset = new IntSet;
    iset->elements.reserve(VINT_SET_MIN_CAP);
    iset->notUpdated = 0;
    return iset;
}

void vint_set_free(IntSet **_iset) {
    IntSet *iset = *_iset;
    delete iset;
    _iset = nullptr;
}

IntSet* vint_set_clone(const IntSet *iset) {
    IntSet *clone = new IntSet;
    clone->elements = iset->elements;
    clone->notUpdated = iset->notUpdated;
    return clone;
}

void vint_set_remove_duplicates(IntSet *iset) {
    std::sort(iset->elements.begin(), iset->elements.end());
    auto it = std::unique(iset->elements.begin(), iset->elements.end());
    iset->elements.resize(static_cast<unsigned long>(std::distance(iset->elements.begin(), it)));
    iset->notUpdated = 0;
}

void vint_set_add_using_original_indexes(IntSet *iset, const int elements[], int size, const int orig[]) {
    for(int i = 0; i < size; i++) {
        iset->elements.push_back(orig[elements[i]]);
    }
    iset->notUpdated += size;

    if(iset->notUpdated >= VINT_SET_FLUSH) {
        vint_set_remove_duplicates(iset);
    }
}

void vint_set_add(IntSet *iset, const int elements[], int size) {
    iset->elements.insert(iset->elements.end(), elements, elements + size);
    iset->notUpdated += size;

    if(iset->notUpdated >= VINT_SET_FLUSH) {
        vint_set_remove_duplicates(iset);
    }
}

const std::vector<int>& vint_set_get_elements( IntSet *iset ) {
    if(iset->notUpdated > 0) {
        vint_set_remove_duplicates(iset);
    }
#ifdef DEBUG
    if(!std::is_sorted(iset->elements.begin(), iset->elements.end())) {
        fprintf(stderr, "ERROR: vector is unsorted\n");
        exit(EXIT_FAILURE);
    }
#endif
    return iset->elements;
}

int vint_set_size(IntSet *iset) {
    if(iset->notUpdated > 0) {
        vint_set_remove_duplicates(iset);
    }
    return (int)iset->elements.size();
}

const int vint_set_find(IntSet *iset, int key) {
    if (iset->elements.empty())
        return -1;

    if(iset->notUpdated > 0) {
        vint_set_remove_duplicates(iset);
    }

    return bsearch_int(iset->elements, key);
}

const int vint_set_intersection(int intersection[], const int size, const int elements[], IntSet *is2) {
    if(is2->notUpdated > 0) {
        vint_set_remove_duplicates(is2);
    }
    int result = 0;
    for (int i = 0; i < size; i++) {
        if (vint_set_find(is2, elements[i]))
            intersection[result++] = elements[i];
    }

    return result;
}

int vint_set_equals(IntSet *is1, IntSet *is2) {
    if(is1->notUpdated > 0) {
        vint_set_remove_duplicates(is1);
    }
    if(is2->notUpdated > 0) {
        vint_set_remove_duplicates(is2);
    }
    if (is1->elements.size() != is2->elements.size())
        return 0;

    return (is1->elements != is2->elements);
}

inline int bsearch_int(const std::vector<int> &v, const int key) {
    register int l = 0;
    register int r = ((int)v.size()) - 1;
    register int m;
    while (l <= r) {
        m = (l + r) / 2;
        if (v[m] == key)
            return m;
        else {
            if (key < v[m])
                r = m - 1;
            else
                l = m + 1;
        }
    }

    return -1;
}

/* inserts given element in an appropriated position in a
 * vector keeping everything sorted */
inline void vint_insert_sort(const int key, std::vector<int> &v) {
#ifdef DEBUG
    for (size_t i = 1; i < v.size(); i++) {
        if (v[i - 1] > v[i]) {
            fprintf(stderr, "ERROR: passing an unsorted vector to function vint_insert_sort\n");
            exit(EXIT_FAILURE);
        }
    }
#endif
    if(!v.empty()) {
        /* doing a binary search */
        register int l = 0;
        register int r = ((int)v.size()) - 1;
        register int m;
        register int ip = -1;  /* insertion pos */
        while (l <= r) {
            m = (l + r) / 2;
            if (v[m] == key) {
                ip = m;
                goto FOUND_POS;
            } else {
                if (key < v[m])
                    r = m - 1;
                else
                    l = m + 1;
            }
        }
        FOUND_POS:
        if (ip == -1)
            ip = l;

        v.insert(v.begin() + ip, key);
    } else {
        v.push_back(key);
    }
}
