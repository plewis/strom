#pragma once

#include "updater.hpp"

namespace strom
    {

    class Chain;

    class TreeUpdater : public Updater
        {
        friend class Chain;

        public:

                                                TreeUpdater();
            virtual                             ~TreeUpdater(); //POLWARN

            virtual double                      calcLogPrior() const;
            double                              calcLogTopologyPrior() const;
            //POLNEW deleted calcLogEdgeLengthPrior

        private:

            virtual void                        revert();
            virtual void                        proposeNewState();
            virtual void                        pullCurrentStateFromModel();
            virtual void                        pushCurrentStateToModel() const;

            virtual void                        reset();

            double                              _orig_edgelen_top;
            double                              _orig_edgelen_middle;
            double                              _orig_edgelen_bottom;

            double                              _new_edgelen_top;
            double                              _new_edgelen_middle;
            double                              _new_edgelen_bottom;

            unsigned                            _case;
            bool                                _topology_changed;
            Node *                              _x;
            Node *                              _y;
            Node *                              _a;
            Node *                              _b;
            Node *                              _c;
            Node *                              _d;

        public:
            typedef std::shared_ptr< TreeUpdater > SharedPtr;
        };

inline TreeUpdater::TreeUpdater()
    {
    // std::cout << "Creating a TreeUpdater" << std::endl;
    Updater::clear();
    reset();    //POLWARN
    }

inline TreeUpdater::~TreeUpdater()
    {
    // std::cout << "Destroying a TreeUpdater" << std::endl;
    }

inline void TreeUpdater::reset()
    {
    _topology_changed       = false;
    _orig_edgelen_top       = 0.0;
    _orig_edgelen_middle    = 0.0;
    _orig_edgelen_bottom    = 0.0;
    _new_edgelen_top        = 0.0;
    _new_edgelen_middle     = 0.0;
    _new_edgelen_bottom     = 0.0;
    _log_hastings_ratio     = 0.0;
    _case                   = 0;
    _x                      = 0;
    _y                      = 0;
    _a                      = 0;
    _b                      = 0;
    _c                      = 0;
    _d                      = 0;
    }

inline double TreeUpdater::calcLogTopologyPrior() const
    {
    typename Tree::SharedPtr tree = _tree_manipulator->getTree();
    assert(tree);
    double n = tree->numLeaves();
    if (tree->isRooted())
        n += 1.0;
    double log_num_topologies = lgamma(2.0*n - 5.0 + 1.0) - (n - 3.0)*log(2.0) - lgamma(n - 3.0 + 1.0);
    return -log_num_topologies;
    }

inline double TreeUpdater::calcLogPrior() const
    {
    double log_topology_prior = calcLogTopologyPrior();
    double log_edge_length_prior = Updater::calcEdgeLengthPrior(); //POLNEW
    return log_topology_prior + log_edge_length_prior;
    }

inline void TreeUpdater::pullCurrentStateFromModel()
    {
    }

inline void TreeUpdater::pushCurrentStateToModel() const
    {
    }

inline void TreeUpdater::proposeNewState()
    {
    _case = 0;
    _topology_changed = false;

    // Choose random internal node x and let a and b equal its left and right children, d its sibling, y its parent, and
    // c y's parent. Thus, x and y are the vertices at the end of the chosen internal edge, a and b are attached to x,
    // and c and d are attached to y:
    //
    //  a     b
    //   \   /
    //    \ /
    //     x     d
    //      \   /
    //       \ /
    //        y
    //        |
    //        |
    //        c
    //
    _x = _tree_manipulator->randomInternalEdge(_lot->uniform());

    _a = _x->getLeftChild();
    _b = _a->getRightSib();
    _y = _x->getParent();
    _c = _y->getParent();
    _d = 0;
    if (_x == _y->getLeftChild())
        _d = _x->getRightSib();
    else
        _d = _y->getLeftChild();

    // Choose focal 3-edge segment to shrink or grow
    bool a_on_path = (_lot->uniform() < 0.5);
    if (a_on_path)
        _orig_edgelen_top = _a->getEdgeLength();
    else
        _orig_edgelen_top = _b->getEdgeLength();

    _orig_edgelen_middle = _x->getEdgeLength();

    bool c_on_path = (_lot->uniform() < 0.5);
    if (c_on_path)
        _orig_edgelen_bottom = _y->getEdgeLength();
    else
        _orig_edgelen_bottom = _d->getEdgeLength();

    double m = exp(_lambda*(_lot->uniform() - 0.5));
    _log_hastings_ratio = 3.0*log(m);

    _new_edgelen_top    = m*_orig_edgelen_top;
    _new_edgelen_middle = m*_orig_edgelen_middle;
    _new_edgelen_bottom = m*_orig_edgelen_bottom;

    // Decide where along focal path (starting from top) to place moved node
    double new_focal_path_length = _new_edgelen_top + _new_edgelen_middle + _new_edgelen_bottom;
    double u = _lot->uniform();
    double new_attachment_point = u*new_focal_path_length;
    if (new_attachment_point <= Node::_smallest_edge_length)
        new_attachment_point = Node::_smallest_edge_length;
    else if (new_focal_path_length - new_attachment_point <= Node::_smallest_edge_length)
        new_attachment_point = new_focal_path_length - Node::_smallest_edge_length;

    // Decide which node to move, and whether the move involves a topology change
    u = _lot->uniform();
    if (a_on_path && c_on_path)
        {
        if (u < 0.5)
            {
            _case = 1;

            // (a)    b*
            //   \   /
            //    \ /
            //    (x)    d
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //       (c)

            if (new_attachment_point > _new_edgelen_top + _new_edgelen_middle)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_b, _d);
                _a->setEdgeLength(_new_edgelen_top + _new_edgelen_middle);
                _x->setEdgeLength(new_attachment_point - _a->getEdgeLength());
                _y->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            else
                {
                _topology_changed = false;
                _a->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top + _new_edgelen_middle - new_attachment_point);
                _y->setEdgeLength(_new_edgelen_bottom);
                }
            }
        else
            {
            _case = 2;

            // (a)    b
            //   \   /
            //    \ /
            //    (x)    d*
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //       (c)

            if (new_attachment_point < _new_edgelen_top)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_b, _d);
                _a->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top - new_attachment_point);
                _y->setEdgeLength(_new_edgelen_middle + _new_edgelen_bottom);
                }
            else
                {
                _topology_changed = false;
                _a->setEdgeLength(_new_edgelen_top);
                _x->setEdgeLength(new_attachment_point - _new_edgelen_top);
                _y->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            }
        }
    else if (!a_on_path && c_on_path)
        {
        if (u < 0.5)
            {
            _case = 3;

            //  a*   (b)
            //   \   /
            //    \ /
            //    (x)    d
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //       (c)

            if (new_attachment_point > _new_edgelen_top + _new_edgelen_middle)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_a, _d);
                _b->setEdgeLength(_new_edgelen_top + _new_edgelen_middle);
                _x->setEdgeLength(new_attachment_point - _b->getEdgeLength());
                _y->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            else
                {
                _topology_changed = false;
                _b->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top + _new_edgelen_middle - new_attachment_point);
                _y->setEdgeLength(_new_edgelen_bottom);
                }
            }
        else
            {
            _case = 4;

            //  a    (b)
            //   \   /
            //    \ /
            //    (x)    d*
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //       (c)

            if (new_attachment_point < _new_edgelen_top)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_a, _d);
                _b->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top - new_attachment_point);
                _y->setEdgeLength(_new_edgelen_middle + _new_edgelen_bottom);
                }
            else
                {
                _topology_changed = false;
                _b->setEdgeLength(_new_edgelen_top);
                _x->setEdgeLength(new_attachment_point - _new_edgelen_top);
                _y->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            }
        }
    else if (a_on_path && !c_on_path)
        {
        if (u < 0.5)
            {
            _case = 5;

            // (a)    b*
            //   \   /
            //    \ /
            //    (x)   (d)
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //        c

            if (new_attachment_point > _new_edgelen_top + _new_edgelen_middle)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_a, _d);
                _a->setEdgeLength(_new_edgelen_top + _new_edgelen_middle);
                _x->setEdgeLength(new_attachment_point - _a->getEdgeLength());
                _d->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            else
                {
                _topology_changed = false;
                _a->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top + _new_edgelen_middle - new_attachment_point);
                _d->setEdgeLength(_new_edgelen_bottom);
                }
            }
        else
            {
            _case = 6;

            // (a)    b
            //   \   /
            //    \ /
            //    (x)   (d)
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //        c*

            if (new_attachment_point < _new_edgelen_top)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_a, _d);
                _d->setEdgeLength(_new_edgelen_bottom + _new_edgelen_middle);
                _x->setEdgeLength(_new_edgelen_top - new_attachment_point);
                _a->setEdgeLength(new_attachment_point);
                }
            else
                {
                _topology_changed = false;
                _a->setEdgeLength(_new_edgelen_top);
                _x->setEdgeLength(new_attachment_point - _new_edgelen_top);
                _d->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            }
        }
    else
        {
        if (u < 0.5)
            {
            _case = 7;

            //  a*   (b)
            //   \   /
            //    \ /
            //    (x)   (d)
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //        c

            if (new_attachment_point > _new_edgelen_top + _new_edgelen_middle)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_b, _d);
                _d->setEdgeLength(new_focal_path_length - new_attachment_point);
                _x->setEdgeLength(new_attachment_point - _new_edgelen_top - _new_edgelen_middle);
                _b->setEdgeLength(_new_edgelen_top + _new_edgelen_middle);
                }
            else
                {
                _topology_changed = false;
                _b->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top + _new_edgelen_middle - new_attachment_point);
                _d->setEdgeLength(_new_edgelen_bottom);
                }
            }
        else
            {
            _case = 8;

            //  a    (b)
            //   \   /
            //    \ /
            //    (x)   (d)
            //      \   /
            //       \ /
            //       (y)
            //        |
            //        |
            //        c*

            if (new_attachment_point < _new_edgelen_top)
                {
                _topology_changed = true;
                _tree_manipulator->nniNodeSwap(_b, _d);
                _b->setEdgeLength(new_attachment_point);
                _x->setEdgeLength(_new_edgelen_top - new_attachment_point);
                _d->setEdgeLength(_new_edgelen_middle + _new_edgelen_bottom);
                }
            else
                {
                _topology_changed = false;
                _b->setEdgeLength(_new_edgelen_top);
                _x->setEdgeLength(new_attachment_point - _new_edgelen_top);
                _d->setEdgeLength(new_focal_path_length - new_attachment_point);
                }
            }
        }

    assert(_a->getEdgeLength() >= Node::_smallest_edge_length);
    assert(_b->getEdgeLength() >= Node::_smallest_edge_length);
    assert(_d->getEdgeLength() >= Node::_smallest_edge_length);
    assert(_x->getEdgeLength() >= Node::_smallest_edge_length);
    assert(_y->getEdgeLength() >= Node::_smallest_edge_length);
    }

inline void TreeUpdater::revert()
    {
    assert(_case > 0 && _case < 9);
    if (_case == 1 || _case == 2)
        {
        if (_topology_changed)
            {
            _tree_manipulator->nniNodeSwap(_d, _b);
            _a->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _y->setEdgeLength(_orig_edgelen_bottom);
            }
        else
            {
            _a->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _y->setEdgeLength(_orig_edgelen_bottom);
            }
        }
    else if (_case == 3 || _case == 4)
        {
        if (_topology_changed)
            {
            _tree_manipulator->nniNodeSwap(_d, _a);
            _b->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _y->setEdgeLength(_orig_edgelen_bottom);
            }
        else
            {
            _b->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _y->setEdgeLength(_orig_edgelen_bottom);
            }
        }
    else if (_case == 5 || _case == 6)
        {
        if (_topology_changed)
            {
            _tree_manipulator->nniNodeSwap(_d, _a);
            _a->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _d->setEdgeLength(_orig_edgelen_bottom);
            }
        else
            {
            _a->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _d->setEdgeLength(_orig_edgelen_bottom);
            }
        }
    else if (_case == 7 || _case == 8)
        {
        if (_topology_changed)
            {
            _tree_manipulator->nniNodeSwap(_d, _b);
            _b->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _d->setEdgeLength(_orig_edgelen_bottom);
            }
        else
            {
            _b->setEdgeLength(_orig_edgelen_top);
            _x->setEdgeLength(_orig_edgelen_middle);
            _d->setEdgeLength(_orig_edgelen_bottom);
            }
        }
    }

}
