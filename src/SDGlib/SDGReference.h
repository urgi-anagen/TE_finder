/*
 * Copyright (c) Laboratoire de Dynamique du Genome et Evolution - Institut 
 * Jacques Monod, 1997. All rights reserved.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * this copyright and notice appears in all copies.
 *
 * This file was written by Hadi Quesneville at Laboratoire de Dynamique 
 * du Genome et Evolution, Institut Jacques Monod, 2 place Jussieu, 
 * 75251 Paris Cedex 05, France.
 * 
 * Contact: Hadi Quesneville
 * Laboratoire de Dynamique du Genome et Evolution,
 * Institut Jacques Monod,
 * 2, place Jussieu, 75251 Paris Cedex 05, France
 * e-mail: hq@ccr.jussieu.fr
 *
 *
 * Laboratoire de Dynamique du Genome et Evolution disclaims all warranties
 * with regard to this software.
 */
//----------------------------------------------------------------------------
// Reference.h 
//----------------------------------------------------------------------------

#ifndef SDGREFERENCE_H
#define SDGREFERENCE_H

#include <stddef.h>
#include <iostream>
#include <SDGError.h>

template<class T>
class SDGReference {

    struct Ref {
        T *ptr;
        unsigned count;

        Ref(const T *p) : count(1) {
            ptr = (T *) (p->clone());
        };
    };

protected:

    Ref *ref;

public:


    SDGReference(void) {
        ref = NULL;
    };

    SDGReference(const T *p) {
        ref = new Ref(p);
    };

    SDGReference(const T &obj) {
        ref = new Ref(&obj);
    };

    SDGReference(Ref *r) {
        if (r)
            r->count++;
        ref = r;
    };

    SDGReference(const SDGReference<T> &robj) {
        if (robj.ref)
            robj.ref->count++;
        if(ref==robj.ref){
            throw "SDGReference: Circular reference error !";
        }
        ref = robj.ref;
    };

    virtual ~SDGReference() {
        clear();
    };

    void clear() {
        if (ref && --(ref->count) == 0) {
            if (ref->ptr) delete ref->ptr;
            delete ref;
        }
        ref = NULL;
    }

    SDGReference<T> &operator=(const SDGReference<T> &robj) {
        if(ref==robj.ref){
            throw "SDGReference: Circular reference error !";
        }
        clear();
        if (robj.isAllocated()) {
            robj.ref->count++;
            ref = robj.ref;
        }
        return *this;
    };

    bool isAllocated() const {
        if (!ref) return false;
        if (ref->ptr) return true;
        else return false;
    };

    const T *getPointer() const {
        if (!ref)
            throw SDGException(this,
                               "SDGReference<T>::getPointer():  empty reference !!!");
        return ref->ptr;
    };

    T *getMutablePointer() {
        if (!ref)
            throw SDGException(this,
                               "SDGReference<T>::getPointer():  empty reference !!!");
        makeExclusif();
        return ref->ptr;
    };

    const T &getObject() const {
        if (!ref)
            throw SDGException(this,
                               "const SDGReference<T>::getObject():  empty reference !!!");
        return *(ref->ptr);
    };

    T &getMutableObject() {
        if (!ref)
            throw SDGException(this,
                               "const SDGReference<T>::getObject():  empty reference !!!");
        makeExclusif();
        return *(ref->ptr);
    };

    void makeExclusif() {
        if (ref) {
            Ref *new_ref = new Ref(ref->ptr);
            clear();
            ref = new_ref;
        } else
            ref = NULL;
    };

    friend int operator==(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) == *(r2.ref->ptr)) return 1;
        return 0;
    };

    friend int operator!=(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) != *(r2.ref->ptr)) return 1;
        return 0;
    };

    friend int operator>(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) > *(r2.ref->ptr)) return 1;
        return 0;
    };

    friend int operator>=(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) >= *(r2.ref->ptr)) return 1;
        return 0;
    };

    friend int operator<(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) < *(r2.ref->ptr)) return 1;
        return 0;
    };

    friend int operator<=(const SDGReference<T> &r1, const SDGReference<T> &r2) {
        if (*(r1.ref->ptr) <= *(r2.ref->ptr)) return 1;
        return 0;
    };


    friend std::ostream &operator<<(std::ostream &os, const SDGReference<T> &robj) {
        if (robj.ref) {
            os << "[count=" << robj.ref->count << "]";
            //	 return operator<<(os,*(robj.ref->ptr));
            return os;
        } else
            os << "(ref NULL!!)";
        return os;
    };

    friend std::istream &operator>>(std::istream &is, const SDGReference<T> &robj) {
        if (!robj.ref)
            throw SDGException(&robj,
                               "SDGReference<T>::operator>>():  empty reference !!!");

        return operator>>(is, *(robj.ref->ptr));
    };
};
#endif







