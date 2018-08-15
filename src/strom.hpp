#pragma once

#include <iostream>
#include <memory>
#include "tree_summary.hpp"
#include "data.hpp"
//#include "likelihood.hpp"
#include <boost/program_options.hpp>

namespace strom {

class Strom
    {
    public:
                            Strom();
                            ~Strom();

        void                clear();
        void                processCommandLineOptions(int argc, const char * argv[]);
        void                run();

    private:
    
        void                        createOutgroup();


        std::string                 _output_file_name;
        std::vector<std::string>    _tree_file_names;
        std::vector<std::string>    _outgroup_names;
        //std::string                 _outgroup_name;
        unsigned                    _treefile_to_plot;
        unsigned                    _trees_to_skip;

        Data::SharedPtr             _data;
        //Likelihood::SharedPtr       _likelihood;
        TreeSummary::SharedPtr      _tree_summary;

        static std::string          _program_name;
        static unsigned             _major_version;
        static unsigned             _minor_version;

    };

inline Strom::Strom()
    {
    //std::cout << "Constructing a Strom" << std::endl;
    clear();
    }

inline Strom::~Strom()
    {
    //std::cout << "Destroying a Strom" << std::endl;
    }

inline void Strom::clear()
    {
    _output_file_name = "";
    //_tree_file_names = "";
    //_outgroup_name = "";
    _treefile_to_plot = 0;
    _trees_to_skip = 0;
    _data           = nullptr;
    //_likelihood     = nullptr;
    _tree_summary   = nullptr;
    }

inline void Strom::processCommandLineOptions(int argc, const char * argv[])
    {
    boost::program_options::variables_map       vm;
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("outfile,o",  boost::program_options::value(&_output_file_name), "name of file to which output is saved (if not specified, defaults to console output)")
        ("treefile,t",  boost::program_options::value(&_tree_file_names)->multitoken(), "name of tree file in NEXUS format")
        ("outgroup",  boost::program_options::value(&_outgroup_names)->multitoken(), "names of taxa to use as outgroup  (if only one name is specified, the tree will be rooted at that leaf)")
        //("outgroup",  boost::program_options::value(&_outgroup_name), "name of taxon to use as outgroup for tree drawing")
        ("plot",  boost::program_options::value(&_treefile_to_plot), "index of tree file to plot (first is index 0)")
        ("skip",  boost::program_options::value(&_trees_to_skip), "number of trees to skip at the beginning of each tree file")
        ;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    try
        {
        const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("strom.conf", desc, false);
        boost::program_options::store(parsed, vm);
        }
    catch(boost::program_options::reading_file & x)
        {
        std::cout << "Note: configuration file (strom.conf) not found" << std::endl;
        }
    boost::program_options::notify(vm);

    // If user specified --help on command line, output usage summary and quit
    if (vm.count("help") > 0)
        {
        std::cout << desc << "\n";
        std::exit(1);
        }

    // If user specified --version on command line, output version and quit
    if (vm.count("version") > 0)
        {
        std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
        std::exit(1);
        }

    // If user specified --plot on command line, check to make sure it is valid
    if (vm.count("plot") > 0)
        {
        if (_treefile_to_plot == 0)
            {
            std::cerr << "plot must be greater than zero" << std::endl;
            std::exit(1);
            }
        }
        
    if (vm.count("treefile") > 0)
        {
        std::cerr << _tree_file_names.size() << " tree file names were specified:\n";
        for (auto nm : _tree_file_names) {
            std::cerr << "  " << nm << std::endl;
            }
        }
        
    if (vm.count("outfile") > 0)
        {
        if (_output_file_name.size() == 0) {
            std::cerr << "output file name specified has zero length" << std::endl;
            std::exit(1);
            }
        //TODO check to make sure file _output_file_name does not exist
        std::cout << "output will be saved to a file named " << _output_file_name << std::endl;
        }
        
    if (vm.count("outgroup") > 0)
        {
        if (_outgroup_names.size() == 1)
            {
            std::cout << "outgroup comprises just one taxon:" << std::endl;
            std::cout << "  " << _outgroup_names[0] << std::endl;
            std::cout << "tree will be rooted at this leaf" << std::endl;
            }
        else
            {
            std::cout << "outgroup comprises these " << _outgroup_names.size() << " taxa:" << std::endl;
            for (auto nm : _outgroup_names)
                {
                std::cout << "  " << nm << std::endl;
                }
            std::cout << "tree will be rooted at the midpoint of the edge below the MRCA of the outgroup taxa" << std::endl;
            }
        }
    else
        std::cout << "No outgroup was specified" << std::endl;
    }
    
inline void Strom::run()
    {
    std::streambuf * buf;
    std::ofstream of;

    if (_output_file_name.size() > 0) {
        of.open(_output_file_name.c_str());
        buf = of.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    std::ostream outf(buf);

    outf << "Starting..." << std::endl;

    try
        {
        // Read in trees
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefiles(outf, _tree_file_names, _trees_to_skip, _treefile_to_plot, _outgroup_names);
        _tree_summary->showSummary(outf);
        }
    catch (XStrom & x)
        {
        std::cerr << "Strom encountered a problem:\n  " << x.what() << std::endl;
        }

    outf << "\nFinished!" << std::endl;
    }

} // namespace strom
