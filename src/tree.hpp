#pragma once

#include <memory>
#include <iostream>
#include "node.hpp"

namespace strom
    {

    class TreeManip;
    //class Likelihood;

    class Tree
        {

        friend class TreeManip;
        //friend class Likelihood;

        public:

                                        Tree();
                                        ~Tree();

        private:

            void                        clear();

            bool                        _is_rooted;
            Node *                      _root;
            unsigned                    _nleaves;
            Node::PtrVector             _preorder;
            Node::Vector                _nodes;

        public:

            typedef std::shared_ptr< Tree > SharedPtr;
        };

    inline Tree::Tree()
        {
        std::cout << "Constructing a Tree" << std::endl;
        clear();
        }

    inline Tree::~Tree()
        {
        std::cout << "Destroying a Tree" << std::endl;
        }

    inline void Tree::clear()
        {
        _is_rooted = false;
        _root = 0;
        _nodes.clear();
        _preorder.clear();
        }

    }
