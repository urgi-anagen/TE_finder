/**
 *
 * Graph.h
 *
 */

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>
#include <SDGString.h>
#include <SDGError.h>

template<typename T>
class Graph
{
  struct Node
    {
      unsigned num_node;
      std::list<Node *> list_next;
      T value;
      Node(unsigned num, T val):num_node(num),value(val){};
      void view(void)
      {
	std::cout<<"node:"<<num_node<<" ->"<<value<<" next:";
	for(typename std::list< Node* >::iterator i=list_next.begin();i!=list_next.end();i++)
	  std::cout<<" "<<(*i)->num_node<<"("<<(*i)->value<<")";
	std::cout<<std::endl;
      };

      void toDot(std::ostream& out)
      {
	out<<"\t"<<value<<std::endl;
	for(typename std::list< Node* >::iterator i=list_next.begin();i!=list_next.end();i++)
	  out<<"\t"<<value<<" -- "<<(*i)->value<<std::endl;
      };
      // TODO code duplication with view
      void toStream(std::ostream& out)
      {
	out<<"node:"<<num_node<<" ->"<<value<<" next:";
	for(typename std::list< Node* >::iterator i=list_next.begin();i!=list_next.end();i++)
		out<<" "<<(*i)->num_node<<"("<<(*i)->value<<")";
	out<<std::endl;
      };
  };
  std::vector< Node* > vec_node;
  std::map< T,Node* > map_node;

  std::vector<T> bfs(unsigned num_node, std::vector<bool>& mark_all)
    {
      if(num_node>=vec_node.size())
	throw SDGException(this,
		  "Graph::bfs(unsigned, std::vector<bool>&): unknown node!");
      
      
      std::vector<bool> mark(vec_node.size(),false);
      std::vector<int> level(vec_node.size(),-1);
      std::vector<int> next(vec_node.size(),-1);
      std::vector<int> prev(vec_node.size(),-1);
      unsigned d=0;
      unsigned f=0;
      next[d]=num_node;
      mark[num_node]=true;
      level[num_node]=0;
      while(d<=f)
	{
	  unsigned x=next[d];
	  Node node_x=*(vec_node[x]);
	  if(node_x.list_next.empty()) d++;
	  for(typename std::list<Node*>::iterator i=node_x.list_next.begin();
	      i!=node_x.list_next.end();i++)
	    {
	      unsigned y=(*i)->num_node;
	      if(!mark[y])
		{
		  prev[y]=x;
		  mark[y]=true;
		  level[y]=level[x]+1;
		  next[++f]=y;
		}
	    }
	  d++;
	}
      std::vector<T> list_val;
      for(typename std::vector<Node*>::iterator i=vec_node.begin();
	  i!=vec_node.end();i++)
	{
	  if(mark[(*i)->num_node])
	    {
	      mark_all[(*i)->num_node]=true;
	      list_val.push_back((*i)->value);
	    }
	}
      return list_val;
    };


 public:

  Graph(void){};
  virtual ~Graph(void){clear();};

  void clear(void)
    {
      for(typename std::vector<Node*>::iterator i=vec_node.begin();
	  i!=vec_node.end();i++)
	delete *i;
    };

  unsigned size(void){return vec_node.size();};

  Node* add_node(const T& val)
    {
      Node *pnode;
      typename std::map<T,Node*>::iterator i=map_node.find(val);
      if(i==map_node.end())
	{
	  pnode=new Node(vec_node.size(),val);
	  vec_node.push_back(pnode);
	  map_node[val]=pnode;
	}
      else
	pnode=i->second;
      
      return pnode;
    };

  void add_edge(const T& val1 ,const T& val2)
    {
      Node *pnode1,*pnode2;
      pnode1=add_node(val1);
      
      pnode2=add_node(val2);
      
      typename std::list<Node*>::iterator j;
      for(j=pnode1->list_next.begin();
	  j!=pnode1->list_next.end();j++)
	{
	  if(*j==pnode2) break;
	}
      if(j==pnode1->list_next.end())
	pnode1->list_next.push_back(pnode2);
      
      for(j=pnode2->list_next.begin();
	  j!=pnode2->list_next.end();j++)
	{
	  if(*j==pnode1) break;
	}
      if(j==pnode2->list_next.end())
	pnode2->list_next.push_back(pnode1);
    };


  std::vector<T> bfs(const T& val)
    {
      typename std::map<T,Node*>::iterator i=map_node.find(val);
      if(i==map_node.end())
	throw SDGException(this,"Graph::bfs(const T&): unknown node!");
      
      unsigned num_node=i->second->num_node;
      
      std::vector<bool> mark(vec_node.size(),false);
      return bfs(num_node,mark);
    };

  void connexComp(typename std::vector< typename std::vector<T> >& grp)
    {
      std::vector<bool> mark(vec_node.size(),false);
      
      for(typename std::vector<Node*>::iterator i=vec_node.begin();
	  i!=vec_node.end();i++)
	{
	  if(!mark[(*i)->num_node])
	    grp.push_back(bfs((*i)->num_node,mark));
	}
    };

  void view(void)
    {
      for(typename std::vector<Node*>::iterator i=vec_node.begin();i!=vec_node.end();i++)
	{
	  (*i)->view();
	}
    };

  void toStream(std::ostream& out)
    {
      for(typename std::vector<Node*>::iterator i=vec_node.begin();i!=vec_node.end();i++)
	{
	  (*i)->toStream(out);
	}
    };

  void toDot(std::ostream& out)
    {
      std::vector< std::vector<T> > gr;
      connexComp(gr);

      unsigned count=0;
      for(typename std::vector< std::vector<T> >::iterator l=gr.begin();l!=gr.end();l++)
	{
	  out<<"graph  cluster"<<++count<<"{"<<std::endl;
	  out<<"\tsize=\"8.5,12\";"<<std::endl;
	  out<<"\tconcentrate=true;"<<std::endl;
	  out<<"\torientation=landscape;"<<std::endl;
	  out<<"\tcenter=true;"<<std::endl;
	  for(typename std::vector<T>::iterator i=l->begin();i!=l->end();i++)
	    {
	      map_node[*i]->toDot(out);
	    }
	  out<<"}"<<std::endl;
	}

    }

  void toDot(SDGString file)
    {
      std::ofstream graphfile(file);

      std::vector< std::vector<T> > gr;
      connexComp(gr);

      unsigned count=0;
      for(typename std::vector< std::vector<T> >::iterator l=gr.begin();l!=gr.end();l++)
	{
	  graphfile<<"graph  cluster"<<++count<<"{"<<std::endl;
	  graphfile<<"\tsize=\"8.5,12\";"<<std::endl;
	  graphfile<<"\tconcentrate=true;"<<std::endl;
	  graphfile<<"\torientation=landscape;"<<std::endl;
	  graphfile<<"\tcenter=true;"<<std::endl;
	  for(typename std::vector<T>::iterator i=l->begin();i!=l->end();i++)
	    {
	      map_node[*i]->toDot(graphfile);
	    }
	  graphfile<<"}"<<std::endl;
	}

      graphfile.close();
    };

};

#endif






