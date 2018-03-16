#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "split.hpp"

namespace strom
    {

    class Tree;
    class TreeManip;
    class Likelihood;
    class Updater;  //POLNEW

    class Node
        {
            friend class Tree;
            friend class TreeManip;
            friend class Likelihood;
            friend class Updater;   //POLNEW

        public:
                                        Node();
                                        ~Node();

                    Node *              getParent()     {return _parent;}
                    Node *              getLeftChild()  {return _left_child;}
                    Node *              getRightSib()   {return _right_sib;}
                    int                 getNumber()     {return _number;}
                    std::string         getName()       {return _name;}
                    double              getEdgeLength() {return _edge_length;}
                    Split               getSplit()      {return _split;}

                    void                setEdgeLength(double v) {_edge_length = v;}

            static const double          _smallest_edge_length; //POLNEW
            typedef std::vector<Node>    Vector;
            typedef std::vector<Node *>  PtrVector;

        private:

            void                clear();

            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            std::string         _name;
            double              _edge_length;
            //double              _x;
            //double              _y;
            //unsigned            _n;
            Split               _split;
        };

    inline Node::Node()
        {
        //std::cout << "Creating Node object" << std::endl;
        clear();
        }

    inline Node::~Node()
        {
        //std::cout << "Destroying Node object" << std::endl;
        }

    inline void Node::clear()
        {
        _left_child = 0;
        _right_sib = 0;
        _parent = 0;
        _number = 0;
        _name = "";
        _edge_length = 0.0;
        //_x = 0.0;
        //_y = 0.0;
        //_n = 0;
        }

    }
