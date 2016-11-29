/**
 *
 * SkipList.h
 *
 */

#ifndef SKIPLIST_H
#define SKIPLIST_H

#include <limits>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

#include "Alea.h"

template <class T>
class SkipList
{
 public:

  //  typedef SkipList<T>::Node value_type;
  //  typedef SkipList<T>::Node* iterator;
  //  typedef const SkipList<T>::Node* const_iterator;


  struct Node
  {
    typedef Node* iterator;
    typedef const Node* const_iterator;
    typedef Node value_type;
    T value;
    long long key;
    std::vector< Node* > forward;
    
    Node(T val, long long k=std::numeric_limits<long long>::max(), 
	 unsigned l=0)
      :  value(val),key(k),forward(l+1){};
    Node(long long k=std::numeric_limits<long long>::max(), unsigned l=0) 
      : key(k),forward(l+1){};
  };

 private:

  Node nil;
  Node head;
  unsigned max_level;
  unsigned level;

 public:
  

  SkipList(unsigned maxl=32) : max_level(maxl),level(1)
    {
      head.forward.resize(maxl+1);
      head.key=INT_MIN;
      nil.forward[1]=&nil;
      for(typename std::vector< Node* >::iterator i=head.forward.begin();
	  i!=head.forward.end();i++)
	*i=&nil;
    };
  ~SkipList(void){clear();};


  void clear(void)
    {
      Node* x=head.forward[1];
      while(x!=&nil)
	{
	  Node* next=x->forward[1];
	  delete x;
	  x=next;
	};
      for(typename std::vector< Node* >::iterator i=head.forward.begin();
	  i!=head.forward.end();i++)
	*i=&nil;
    };

  Node& search(long long searchKey)
    {
      Node* x=&head;

      for(unsigned i=level;i>0; i--)
	  while(x->forward[i]->key<searchKey) x=x->forward[i];

      x=x->forward[1];
      if(x->key==searchKey) return *x; else return nil;
    };

  Node& searchAfter(long long searchKey)
    {
      Node* x=&head;

      for(unsigned i=level;i>0; i--)
	  while(x->forward[i]->key<searchKey) x=x->forward[i];

      x=x->forward[1];
      if(x->key==searchKey) return *(x->forward[1]); else return *x;
    };

  Node& searchBefore(long long searchKey)
    {
      Node* prev=&head,*x=&head;

      for(unsigned i=level;i>0; i--)
	  while(x->forward[i]->key<searchKey)
	    {
	      prev=x;
	      x=x->forward[i];
	    }
      
      if(x->key==searchKey) x=prev;
      if(x==&head) x=&nil;
      return *x;
    };

  void insert(long long searchKey, T newValue)
    {
      std::vector<Node*> update(max_level+1);

      Node* x=&head;
      for(unsigned i=level;i>0; i--)
	{
	  while(x->forward[i]->key<searchKey) x=x->forward[i];
	  update[i]=x;
	}
      x=x->forward[1];
      if(x->key==searchKey) x->value=newValue; 
      else
	{
	  unsigned newLevel=randomLevel();
	  if(newLevel>level)
	    {
	      for(unsigned i=level+1; i<=newLevel;i++)
		update[i]=&head;
	      level=newLevel;
	    }
	  x=new Node(newValue,searchKey,newLevel);
	  for(unsigned i=1;i<=newLevel;i++)
	    {
	      x->forward[i]=update[i]->forward[i];
	      update[i]->forward[i]=x;
	    }
	}
    };

  void remove(long long searchKey)
    {
      std::vector<Node*> update(max_level+1);

      Node* x=&head;
      for(unsigned i=level;i>0; i--)
	{
	  while(x->forward[i]->key<searchKey) x=x->forward[i];
	  update[i]=x;
	}
      x=x->forward[1];
      if(x->key==searchKey) 	
	{
	  for(unsigned i=1;i<=level;i++)
	    {
	      if(update[i]->forward[i]!=x) break;
	      update[i]->forward[i]=x->forward[i];
	    }
	  delete x;
	  while(level>1 && head.forward[level]==&nil)
	    level--;
	}
    };

  unsigned randomLevel()
    {
      unsigned lvl=1;
      while(alea.flip(0.5) && lvl<max_level) lvl++;
      return lvl;
    };
  

  bool isNil(const Node& n){return (n.key==nil.key)?true:false ;};
  Node getNil(void){return nil;};
  
  Node& first(void){return *head.forward[1];};
  Node& last(void){return searchBefore(nil.key);};
  Node& next(const Node& n)
    {
      if(isNil(n)) return nil;
      return *(n.forward[1]);
    };
  friend std::ostream& operator<<(std::ostream& os, SkipList& sl)
    {
      Node n=sl.first();
      while(!sl.isNil(n))
	{
	  os<<n.key<<"="<<n.value<<std::endl;
	  n=sl.next(n);
	}
      return os;
    };
  
  Node* begin(void){return head.forward[1];};
  const Node* begin(void) const {return head.forward[1];};

  Node* end(void){return &searchBefore(nil.key);};
  const Node* end(void) const {&searchBefore(nil.key);};
  /*
  iterator operator++(const iterator& it) {return &(next(*it));};
  const_iterator operator++(const const_iterator& it) const 
    {return &(next(*it));};
  */
};
#endif




